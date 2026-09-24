"""Joint fleet MILP: choose the CARRIERS and assign ALL 298 targets in one optimisation (stage 10).

Every previous fleet layer was two-stage and lost: either routes were grown independently and then set-covered
(they overlap -- 82 carrier routes needed 20 craft for 293 targets), or a partition was fixed first and routes grown
into it (they collapse -- Lloyd bases 17-22 -> 8-14 flybys).  Neither optimises  J = sum_c J_i(tank_c) + misses.

With the stage-10 insertion surrogate (tools/insert_model.py, 1225 real settles) the fleet cost is explicit, so the
whole thing is one MILP:

    y_c in {0,1}          open carrier c (library results/s7/libopt.pkl, bases settled by `bases`)
    x_{t,c,b} in {0,1}    target t is the (band b)-th insertion on carrier c
    z_t in {0,1}          target t is missed
    F_c = bf_c y_c + sum kg_{tcb} x_{tcb}                       propellant above the 600 kg dry mass
    J_c >= alpha_m F_c + beta_m y_c    for tangents m           epigraph of the CONVEX J_i(600+F) = 1+x+x^2
    sum_{c: t in base_c} y_c + sum_{c,b} x_{tcb} + z_t >= 1     coverage (a carrier's DP base flies for free)
    sum_t x_{tcb} <= C_b y_c                                    band capacity
    min sum_c J_c + sum_t z_t

The band index b is the SLOT the insertion occupies: kg_{tcb} = predict_kg(d_tc, depth = n0_c + m_b, ...), so the
surrogate's depth^0.91 law is carried exactly (the solver fills cheap low slots first and, by the rearrangement
inequality, puts the FARTHEST target in the shallowest slot -- which is also where predict_ok is highest).
N is an OUTPUT.  The LP relaxation is a genuine lower bound on what this carrier library can do.

Commands
  bases  lib.pkl outdir            settle every carrier's DP base, record base fuel + candidate table
  milp   outdir                    build and solve the MILP (LP bound first), write plan.json
  build  outdir                    realise the plan: real insertions with run_ialns.w_insert, save an IFleet

MEASURED, 2026-09-21 (results/s10/milp/calibration.json) -- READ THIS BEFORE TRUSTING A PLAN
------------------------------------------------------------------------------------------
The N=10 plan was built for real (134 of 177 planned insertions settled, 76 %):

    predicted sum J_i 15.53 covering 287   ->   SETTLED sum J_i 18.15 covering 244   (J 28.5 -> 74.2)

and the per-insertion calibration is the point:

  * the surrogate underprices an ASSIGNED insertion by x1.94 in aggregate (median ratio 1.97);
    shrinking the base-derived margin by n0/depth (--mshrink, now on by default) recovers a
    quarter of it (x1.57 residual), the rest is the freedom premium itself;
  * on these 134 assigned insertions predict_kg scores R2(log kg) 0.010 and median |err| 21.6 kg,
    against 0.651 / 18.4 kg on its own free/stratified sample, and against 26.5 kg for the constant
    "every insertion costs 37 kg".  (Part of that collapse is range restriction -- the MILP only
    uses d 0.03-0.12 AU while the fit spans 0.004-0.35 -- but the aggregate bias is not.)
  * end to end, same library and same settler: free choice (tools/carrier_farm.py, results/s7/farmO)
    gives 96 routes at a MEDIAN 0.0395 J/flyby, depth 34, tank 968 kg; MILP assignment gives
    0.0617 J/flyby, depth 29, tank 1347 kg.  x1.56 in J per flyby, everything else held fixed.

So the fleet MILP is well posed and its bounds are real, but its COST MATRIX is not: closest
approach is a good price for a target a route chooses and a poor one for a target it is given.

Usage: fleet_milp.py bases results/s7/libopt.pkl results/s10/milp
       fleet_milp.py milp results/s10/milp --pmin 0.3 --dmax 0.10
       fleet_milp.py greedy results/s10/milp --fix-n 10 --plan gplan_n10.json
       fleet_milp.py build results/s10/milp --plan gplan_n10.json
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pickle, pathlib, argparse, collections, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI, carrier_dp as CD
import insert_model as IM
from ctoc14.search import UNREACHABLE
from ctoc14.constants import DAY, T_MISSION

ALL = [t for t in range(1, 301) if t not in UNREACHABLE]
RIDGE_KG = 14.8                      # tank grows this much per flyby along the real depth/tank ridge
_A = None


# ------------------------------------------------------------------ stage 1: settle the carrier bases
def w_base(job):
    """Settle carrier c's DP base (leave-one-out retries), then tabulate every target's closest approach to the
    SETTLED base trajectory (the libopt 'neigh' was measured on the un-settled seed)."""
    k, c, dmax = job
    from ctoc14.globalopt import insert_block
    full = list(c['base'])
    if len(full) < 4:
        return k, None
    ip = None
    for drop in [None] + list(range(len(full))):
        base = [b for i, b in enumerate(full) if i != drop]
        try:
            q = CD.seed_ip(RI.eph(), c['tL'], np.array(c['vinf']),
                           np.array([b[1] for b in base]), np.array([b[2] for b in base]))
            if insert_block(q, [(b[0], b[1]) for b in base], log=RI.QUIET).max() > 1e4:
                continue
            if RI.settle(q, 100) <= 150:
                ip = q; break
        except Exception:
            continue
    if ip is None:
        return k, None
    st = RI.ist(ip); tank = float(ip.tank())
    if not np.isfinite(tank) or tank < 600 or tank > 1400:
        return k, None
    have = set(int(x) for x in st['asts'])
    left = [t for t in ALL if t not in have]
    cands = RI.w_cands(('h', st, left, dmax))
    by = {}
    for x in sorted(cands, key=lambda x: x['dist']):
        by.setdefault(x['ast'], []).append((float(x['dist']), float(x['t'])))
    by = {a: v[:3] for a, v in by.items()}
    return k, dict(k=k, tank=tank, base=sorted(have), tf=[float(x) for x in st['tf']],
                   asts=[int(x) for x in st['asts']], cands=by, st=st)


def cmd_bases(a):
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    (out / 'bases').mkdir(exist_ok=True)
    say = RI.logger(out / 'log.txt')
    lib = pickle.load(open(a.lib, 'rb'))
    say(f'settling {len(lib)} carrier bases (dmax {a.dmax} AU)')
    jobs = [(i, c, a.dmax) for i, c in enumerate(lib)]
    recs = {}; t0 = time.time(); nfail = 0
    with mp.get_context('fork').Pool(a.nproc) as pl:
        for k, r in pl.imap_unordered(w_base, jobs, chunksize=1):
            if r is None:
                nfail += 1; continue
            np.savez(out / 'bases' / f'base_{k:04d}.npz', **r['st'])
            recs[k] = {q: v for q, v in r.items() if q != 'st'}
            if len(recs) % 20 == 0:
                say(f'  {len(recs)} settled, {nfail} failed ({time.time()-t0:.0f} s)')
    say(f'{len(recs)} bases settled, {nfail} failed in {time.time()-t0:.0f} s')
    tk = np.array([r['tank'] for r in recs.values()]); nb = np.array([len(r['base']) for r in recs.values()])
    say(f'  base depth {nb.min()}-{nb.max()} (median {np.median(nb):.0f}); '
        f'base tank {tk.min():.0f}-{tk.max():.0f} kg (median {np.median(tk):.0f}); '
        f'union of bases covers {len(set(t for r in recs.values() for t in r["base"]))} targets')
    json.dump({str(k): v for k, v in recs.items()}, open(out / 'bases.json', 'w'))
    say(f'-> {out}/bases.json')


# ------------------------------------------------------------------ stage 2: the MILP
def _margin(tf, t):
    return max(3.0, min(400.0, float(np.abs(np.asarray(tf) - t).min()) / DAY))


def build_pairs(recs, a):
    """(carrier, target, band) columns with their surrogate kg cost.  Band b covers added-slots
    [b*W+1 .. (b+1)*W]; the representative added index is its midpoint, so depth = n0 + m_rep - 1."""
    W = a.band
    nband = int(np.ceil(a.max_added / W))
    cols = []                                          # (ci, targ, band, kg, d, t, p)
    for ci, r in enumerate(recs):
        n0 = len(r['base']); tf = np.array(r['tf'])
        for t, lst in r['cands'].items():
            d, te = lst[0]
            if d > a.dmax:
                continue
            mg0 = _margin(tf, te)
            for b in range(nband):
                m_rep = b * W + (W + 1) / 2.0
                dep = n0 + m_rep - 1.0
                tank = r['tank'] + RIDGE_KG * (m_rep - 1.0)
                # the base is SPARSE (17 flybys / 10 yr -> 215 d gaps); as the route fills, the gap the
                # insertion lands in shrinks like n0/depth, and kg ~ margin^-0.744.  Pricing at the base
                # margin underprices a full route by up to 1.7x, so shrink it with the depth.
                mg = max(3.0, mg0 * (n0 / dep) ** a.mshrink)
                p = IM.predict_ok(d, dep, tank, mg)
                if p < a.pmin:
                    break                              # p falls with depth: no deeper band will pass either
                kg = a.kcal * max(a.kg_floor, IM.predict_kg(d, dep, tank, mg))
                if a.risk:
                    kg = kg / max(p, 0.05)             # risk-loaded: price the expected number of attempts
                cols.append((ci, int(t), b, float(kg), float(d), float(te), float(p)))
    return cols, nband


def tangents(a):
    """Tangent lines of J_i(600+F) = 1 + F/1400 + (F/1400)^2 at F = 0, step, ... : J >= alpha F + beta."""
    out = []
    for F0 in np.arange(0.0, a.max_fuel + 1e-9, a.tan_step):
        x0 = F0 / 1400.0
        al = (1.0 + 2.0 * x0) / 1400.0
        be = (1.0 + x0 + x0 * x0) - al * F0
        out.append((float(al), float(be)))
    return out


def solve_milp(recs, cols, nband, a, say, relax=False, fix_n=None):
    from scipy.optimize import milp, linprog, LinearConstraint, Bounds
    from scipy.sparse import coo_matrix
    nc = len(recs); nt = len(ALL); ti = {t: i for i, t in enumerate(ALL)}
    W = a.band
    ny, nx, nz, nj = nc, len(cols), nt, nc
    iy = 0; ix = iy + ny; iz = ix + nx; ij = iz + nz; nvar = ij + nj
    TG = tangents(a)
    R = []; C = []; V = []; lo = []; hi = []; row = 0

    # 1. coverage  sum_{c: t in base} y + sum x + z >= 1
    base_of = [set(r['base']) for r in recs]
    for t in ALL:
        for c in range(nc):
            if t in base_of[c]: R.append(row); C.append(iy + c); V.append(1.0)
        R.append(row); C.append(iz + ti[t]); V.append(1.0)
        lo.append(1.0); hi.append(np.inf); row += 1
    cov_row = {t: ti[t] for t in ALL}
    for j, (ci, t, b, kg, d, te, p) in enumerate(cols):
        R.append(cov_row[t]); C.append(ix + j); V.append(1.0)

    # 2. band capacity  sum_t x_{c,b} - W y_c <= 0
    brow = {}
    for c in range(nc):
        for b in range(nband):
            brow[(c, b)] = row
            R.append(row); C.append(iy + c); V.append(-float(W))
            lo.append(-np.inf); hi.append(0.0); row += 1
    for j, (ci, t, b, kg, d, te, p) in enumerate(cols):
        R.append(brow[(ci, b)]); C.append(ix + j); V.append(1.0)

    # 3. J epigraph   J_c - alpha F_c - beta y_c >= 0,  F_c = bf_c y_c + sum kg x
    bycar = collections.defaultdict(list)
    for j, (ci, t, b, kg, d, te, p) in enumerate(cols):
        bycar[ci].append((j, kg))
    for c in range(nc):
        bf = recs[c]['tank'] - 600.0
        for (al, be) in TG:
            R.append(row); C.append(ij + c); V.append(1.0)
            R.append(row); C.append(iy + c); V.append(-(al * bf + be))
            for j, kg in bycar[c]:
                R.append(row); C.append(ix + j); V.append(-al * kg)
            lo.append(0.0); hi.append(np.inf); row += 1

    # 4. tank cap  bf y + sum kg x <= max_fuel
    for c in range(nc):
        bf = recs[c]['tank'] - 600.0
        R.append(row); C.append(iy + c); V.append(bf)
        for j, kg in bycar[c]:
            R.append(row); C.append(ix + j); V.append(kg)
        lo.append(-np.inf); hi.append(a.max_fuel); row += 1

    # 5. optional craft-count equality
    if fix_n is not None:
        for c in range(nc): R.append(row); C.append(iy + c); V.append(1.0)
        lo.append(float(fix_n)); hi.append(float(fix_n)); row += 1

    A = coo_matrix((V, (R, C)), shape=(row, nvar)).tocsc()
    obj = np.zeros(nvar); obj[iz:iz + nz] = 1.0; obj[ij:ij + nj] = 1.0
    bounds = Bounds(np.zeros(nvar), np.concatenate([np.ones(ny + nx + nz), np.full(nj, np.inf)]))
    say(f'  MILP {nvar} vars ({ny} y, {nx} x, {nz} z, {nj} J), {row} rows, {len(V)} nnz')
    if relax:                                       # integrality 0 -> plain LP through the same HiGHS path
        res = milp(obj, constraints=LinearConstraint(A, lo, hi), integrality=np.zeros(nvar),
                   bounds=bounds, options=dict(time_limit=a.tlim, presolve=True))
    else:
        integ = np.concatenate([np.ones(ny + nx + nz), np.zeros(nj)])
        res = milp(obj, constraints=LinearConstraint(A, lo, hi), integrality=integ, bounds=bounds,
                   options=dict(time_limit=a.tlim, mip_rel_gap=a.gap, presolve=True))
    if res.x is None:
        say(f'  FAILED: {res.message}'); return None
    x = res.x
    out = dict(obj=float(res.fun), msg=str(res.message)[:70],
               bound=float(getattr(res, 'mip_dual_bound', res.fun) or res.fun))
    out['open'] = [c for c in range(nc) if x[iy + c] > 0.5]
    out['miss'] = [t for t in ALL if x[iz + ti[t]] > 0.5]
    asg = collections.defaultdict(list)
    for j, (ci, t, b, kg, d, te, p) in enumerate(cols):
        if x[ix + j] > 0.5:
            asg[ci].append(dict(ast=t, band=b, kg=kg, d=d, t=te, p=p))
    out['assign'] = {str(c): sorted(v, key=lambda q: q['band']) for c, v in asg.items()}
    out['J_c'] = {str(c): float(x[ij + c]) for c in out['open']}
    fuel = {}
    for c in out['open']:
        fuel[str(c)] = recs[c]['tank'] - 600.0 + sum(q['kg'] for q in asg.get(c, []))
    out['fuel'] = fuel
    return out


def cmd_milp(a):
    out = pathlib.Path(a.out); say = RI.logger(out / f'log{a.tag}.txt')
    raw = json.load(open(out / 'bases.json'))
    recs = []
    for k, r in sorted(raw.items(), key=lambda kv: int(kv[0])):
        r['cands'] = {int(t): v for t, v in r['cands'].items()}
        r['lib'] = int(k); recs.append(r)
    if a.ncar and a.ncar < len(recs):
        # diverse subset: greedily maximise NEW base coverage (the library is a local-search cloud, so the
        # top-value carriers are near-duplicates and a value cut alone would give the MILP nothing to combine)
        pick = []; seen = set(); pool = list(recs)
        while pool and len(pick) < a.ncar:
            j = max(range(len(pool)), key=lambda i: (len(set(pool[i]['base']) - seen), -pool[i]['tank']))
            pick.append(pool.pop(j)); seen |= set(pick[-1]['base'])
        recs = pick
    say(f'=== MILP over {len(recs)} carriers, dmax {a.dmax} AU, pmin {a.pmin}, band {a.band}, '
        f'max_added {a.max_added}, risk {bool(a.risk)}')
    cols, nband = build_pairs(recs, a)
    say(f'  {len(cols)} (carrier,target,band) columns, {nband} bands; '
        f'{len(set((c[0], c[1]) for c in cols))} distinct (carrier,target) pairs; '
        f'reachable targets {len(set(c[1] for c in cols) | set(t for r in recs for t in r["base"]))}/{len(ALL)}')
    t0 = time.time()
    lp = solve_milp(recs, cols, nband, a, say, relax=True)
    if lp: say(f'  LP bound  sum J_i {lp["obj"]:.4f}  -> J {lp["obj"]+2:.4f}   ({time.time()-t0:.0f} s)')
    t0 = time.time()
    r = solve_milp(recs, cols, nband, a, say, relax=False, fix_n=a.fix_n)
    if r is None: return
    r['lp_bound'] = lp['obj'] if lp else None
    nfb = {c: len(recs[c]['base']) + len(r['assign'].get(str(c), [])) for c in r['open']}
    sumJ = sum(IM.cost_J(600 + r['fuel'][str(c)]) for c in r['open'])
    say(f'  MILP {r["msg"]}: {len(r["open"])} craft, {len(ALL)-len(r["miss"])} covered, '
        f'predicted sum J_i {sumJ:.4f}, J {sumJ + len(r["miss"]) + 2:.4f}  '
        f'(HiGHS obj {r["obj"]:.4f}, dual bound {r["bound"]:.4f})   ({time.time()-t0:.0f} s)')
    for c in sorted(r['open'], key=lambda c: -nfb[c]):
        A = r['assign'].get(str(c), [])
        say(f'   carrier {recs[c]["lib"]:4d}: base {len(recs[c]["base"]):2d} @ {recs[c]["tank"]:.0f} kg '
            f'+ {len(A):2d} assigned = {nfb[c]:2d} fb, predicted tank {600+r["fuel"][str(c)]:.0f} kg, '
            f'J_i {IM.cost_J(600+r["fuel"][str(c)]):.4f} ({IM.cost_J(600+r["fuel"][str(c)])/max(nfb[c],1):.4f} J/fb); '
            f'kg ' + ' '.join(f'{q["kg"]:.0f}' for q in A))
    if r['miss']: say(f'   missed: {r["miss"]}')
    plan = dict(carriers=[dict(lib=recs[c]['lib'], tank=recs[c]['tank'], base=recs[c]['base'],
                               assign=r['assign'].get(str(c), []), fuel=r['fuel'][str(c)]) for c in r['open']],
                miss=r['miss'], pred_sumJ=sumJ, lp_bound=r['lp_bound'], args=vars(a))
    json.dump(plan, open(out / a.plan, 'w'), indent=1)
    say(f'-> {out}/{a.plan}')


# ------------------------------------------------------------------ stage 2b: greedy on the SAME objective
# HiGHS needs hours to close this MILP, so the plan that actually gets built comes from a vectorised
# regret-insertion heuristic minimising exactly the MILP objective.  The MILP's job is then the BOUND.
def _mats(recs, a):
    """K[c,t] = B d^N (margin/60)^-S  and D[c,t] = d, MG[c,t] = margin; inf where the pair does not exist."""
    m = IM.model()['cost']; nt = len(ALL); ti = {t: i for i, t in enumerate(ALL)}
    nc = len(recs)
    K = np.full((nc, nt), np.inf); D = np.full((nc, nt), np.inf); MG = np.full((nc, nt), 60.0)
    BASE = np.zeros((nc, nt), bool)
    for c, r in enumerate(recs):
        for t in r['base']:
            if t in ti: BASE[c, ti[t]] = True
        tf = np.array(r['tf'])
        for t, lst in r['cands'].items():
            if t not in ti: continue
            d, te = lst[0]
            if d > a.dmax: continue
            mg = _margin(tf, te)
            K[c, ti[t]] = max(a.kg_floor, m['B'] * d ** m['N'])
            D[c, ti[t]] = d; MG[c, ti[t]] = mg
    return K, D, MG, BASE


def assign(sel, K, D, MG, BASE, recs, a):
    """Regret insertion of every uncovered target into the chosen carriers, with the surrogate's depth law.
    Returns (sum J_i + misses, tank[], depth[], plan per carrier, missed targets)."""
    m = IM.model()['cost']; w = IM.model()['success']['w']
    nt = K.shape[1]; S = np.asarray(sel)
    tank = np.array([recs[c]['tank'] for c in S]); dep = np.array([float(len(recs[c]['base'])) for c in S])
    n0 = np.array([float(len(recs[c]['base'])) for c in S])
    need = ~BASE[S].any(0)
    Ks = K[S]; Ds = D[S]; Ms0 = MG[S]
    out = {c: [] for c in S}; miss = []
    while need.any():
        Ms = np.maximum(3.0, Ms0 * ((n0 / np.maximum(dep, 1.0)) ** a.mshrink)[:, None])
        f = ((tank / 900.0) ** m['P_tank'] * (np.maximum(dep, 1.0) / 30.0) ** m['Q_depth'])
        kg = a.kcal * Ks * f[:, None] * (Ms / 60.0) ** (-m['S_margin'])
        z = (w[0] + w[1] * Ds + w[2] * dep[:, None] / 10.0 + w[3] * tank[:, None] / 1000.0
             + w[4] * np.log10(np.maximum(Ms, 1.0)))
        p = 1.0 / (1.0 + np.exp(-np.clip(z, -40, 40)))
        dJ = IM.cost_J(tank[:, None] + kg) - IM.cost_J(tank[:, None])
        dJ = np.where((p >= a.pmin) & np.isfinite(kg) & (tank[:, None] + kg <= 600 + a.max_fuel), dJ, np.inf)
        dJ[:, ~need] = np.inf
        best = dJ.min(0)
        if not np.isfinite(best).any():
            miss += [ALL[j] for j in np.nonzero(need)[0]]; break
        srt = np.sort(dJ, axis=0)
        second = srt[1] if dJ.shape[0] > 1 else np.full(nt, 1.0)
        regret = np.where(np.isfinite(best), np.minimum(second, 1.0) - best, -np.inf)
        j = int(np.argmax(np.where(best < 1.0, regret, -np.inf)))
        if not np.isfinite(best[j]) or best[j] >= 1.0:          # nothing worth placing: miss the rest
            for jj in np.nonzero(need)[0]:
                if not np.isfinite(best[jj]) or best[jj] >= 1.0: miss.append(ALL[jj]); need[jj] = False
            continue
        i = int(np.argmin(dJ[:, j]))
        out[S[i]].append(dict(ast=ALL[j], kg=float(kg[i, j]), d=float(Ds[i, j]),
                              t=float(recs[S[i]]['cands'][ALL[j]][0][1]), band=len(out[S[i]]), p=float(p[i, j])))
        tank[i] += kg[i, j]; dep[i] += 1.0; need[j] = False
    J = float(sum(IM.cost_J(x) for x in tank) + len(miss))
    return J, tank, dep, out, miss


def cmd_greedy(a):
    out = pathlib.Path(a.out); say = RI.logger(out / f'log{a.tag}.txt')
    raw = json.load(open(out / 'bases.json'))
    recs = []
    for k, r in sorted(raw.items(), key=lambda kv: int(kv[0])):
        r['cands'] = {int(t): v for t, v in r['cands'].items()}; r['lib'] = int(k); recs.append(r)
    K, D, MG, BASE = _mats(recs, a)
    say(f'=== GREEDY over {len(recs)} carriers, dmax {a.dmax}, pmin {a.pmin}; '
        f'{int(np.isfinite(K).sum())} pairs')
    rng = np.random.default_rng(a.seed)
    best = None
    for N in ([a.fix_n] if a.fix_n else range(a.nmin, a.nmax + 1)):
        # greedy construction: add the carrier that most reduces the objective
        sel = []
        while len(sel) < N:
            cand = [c for c in range(len(recs)) if c not in sel]
            sc = [(assign(sel + [c], K, D, MG, BASE, recs, a)[0], c) for c in cand]
            sc.sort(); sel.append(sc[0][1])
        J = assign(sel, K, D, MG, BASE, recs, a)[0]
        # local search: swap one carrier
        for it in range(a.ls_rounds):
            moved = False
            for pos in range(N):
                base = [c for i, c in enumerate(sel) if i != pos]
                sc = [(assign(base + [c], K, D, MG, BASE, recs, a)[0], c)
                      for c in range(len(recs)) if c not in base]
                sc.sort()
                if sc[0][0] < J - 1e-6:
                    J = sc[0][0]; sel = base + [sc[0][1]]; moved = True
            if not moved: break
        Jf, tank, dep, plan, miss = assign(sel, K, D, MG, BASE, recs, a)
        say(f'  N {N}: sum J_i + miss {Jf:.4f} -> J {Jf+2:.4f}; covered {len(ALL)-len(miss)}, '
            f'tanks ' + ' '.join(f'{x:.0f}' for x in tank) + f'; depths ' + ' '.join(f'{x:.0f}' for x in dep))
        if best is None or Jf < best[0]:
            best = (Jf, N, sel, tank, dep, plan, miss)
    Jf, N, sel, tank, dep, plan, miss = best
    sumJ = float(sum(IM.cost_J(x) for x in tank))
    say(f'BEST N {N}: predicted sum J_i {sumJ:.4f}, misses {len(miss)}, J {sumJ+len(miss)+2:.4f}')
    for i, c in enumerate(sel):
        A = plan[c]
        say(f'   carrier {recs[c]["lib"]:4d}: base {len(recs[c]["base"]):2d} @ {recs[c]["tank"]:.0f} kg + '
            f'{len(A):2d} assigned = {int(dep[i])} fb, tank {tank[i]:.0f} kg, J_i {IM.cost_J(tank[i]):.4f} '
            f'({IM.cost_J(tank[i])/max(dep[i],1):.4f} J/fb); d ' + ' '.join(f'{q["d"]:.3f}' for q in A))
    if miss: say(f'   missed {len(miss)}: {miss}')
    pl = dict(carriers=[dict(lib=recs[c]['lib'], tank=recs[c]['tank'], base=recs[c]['base'],
                             assign=plan[c], fuel=float(tank[i] - 600.0)) for i, c in enumerate(sel)],
              miss=miss, pred_sumJ=sumJ, lp_bound=None, args={k: v for k, v in vars(a).items()})
    json.dump(pl, open(out / a.plan, 'w'), indent=1)
    say(f'-> {out}/{a.plan}')


# ------------------------------------------------------------------ stage 3: build it for real
def w_build(job):
    ci, lib, st, items, a = job
    st = {k: (float(v) if k == 'tL' else v) for k, v in st.items()}
    tank = float(RI.ipr(st).tank())
    log = []; done = []; failed = []
    t0 = time.time()
    for it in items:
        if tank > a.max_tank or time.time() - t0 > a.tmax:
            failed.append(it['ast']); continue
        X = int(it['ast']); got = False
        # candidate epochs: the planned one first, then the local minima on the CURRENT trajectory
        eps = [it['t']]
        try:
            cc = sorted(RI.w_cands(('h', st, [X], a.dmax)), key=lambda q: q['dist'])[:a.ntry]
            eps += [q['t'] for q in cc if abs(q['t'] - it['t']) > 5 * DAY]
        except Exception:
            pass
        for te in eps[:a.ntry]:
            r = RI.w_insert(('h', st, X, te))
            if r.get('ok') and r['tank'] <= a.max_tank:
                log.append((X, round(r['tank'] - tank, 1), round(it['kg'], 1), round(it['d'], 4)))
                st, tank = r['st'], r['tank']; got = True; break
        if got: done.append(X)
        else: failed.append(X)
    return ci, dict(lib=lib, st=st, tank=tank, done=done, failed=failed, log=log, sec=time.time() - t0)


def cmd_build(a):
    out = pathlib.Path(a.out); say = RI.logger(out / 'log.txt')
    plan = json.load(open(out / a.plan))
    say(f'=== BUILD {len(plan["carriers"])} routes, predicted sum J_i {plan["pred_sumJ"]:.4f}')
    jobs = []
    for ci, c in enumerate(plan['carriers']):
        z = np.load(out / 'bases' / f'base_{c["lib"]:04d}.npz'); st = {k: z[k] for k in z.files}
        st['tL'] = float(st['tL'])
        jobs.append((ci, c['lib'], st, c['assign'], a))
    res = {}
    with mp.get_context('fork').Pool(min(a.nproc, len(jobs))) as pl:
        for ci, r in pl.imap_unordered(w_build, jobs):
            res[ci] = r
            pc = plan['carriers'][ci]
            nfb = len(r['st']['asts'])
            say(f'  carrier {r["lib"]:4d}: {len(r["done"])}/{len(pc["assign"])} inserted -> {nfb} fb @ '
                f'{r["tank"]:.0f} kg (predicted {600+pc["fuel"]:.0f}), J_i {RI.cost(r["tank"]):.4f} '
                f'({RI.cost(r["tank"])/nfb:.4f} J/fb), failed {r["failed"]}  ({r["sec"]:.0f} s)')
            for X, real, pred, d in r['log']:
                say(f'      {X:4d} d {d:.3f} AU: predicted {pred:6.1f} kg, real {real:6.1f} kg')
    fl = RI.IFleet()
    mdir = out / 'routes'; mdir.mkdir(exist_ok=True)
    for f in mdir.glob('route_*.npz'): f.unlink()
    for ci in sorted(res):
        r = res[ci]
        fl.routes[f'{ci:02d}'] = dict(st=r['st'], tank=r['tank'])
        np.savez(mdir / f'route_{ci:02d}.npz', **r['st'])
    fl.save(out / 'fleet', note='joint fleet MILP (stage 10)')
    cov = fl.coverage()
    sumJ = sum(RI.cost(x['tank']) for x in fl.routes.values())
    say(f'BUILT: {len(fl.routes)} craft, covered {len(cov)}/{len(ALL)}, settled sum J_i {sumJ:.4f}, '
        f'J {fl.J():.4f}  (MILP predicted sum J_i {plan["pred_sumJ"]:.4f}, gap {sumJ-plan["pred_sumJ"]:+.4f})')
    miss = sorted(set(ALL) - set(cov))
    say(f'  uncovered {len(miss)}: {miss}')
    say(f'  -> {out}/fleet, {mdir}/route_*.npz')
    json.dump(dict(sumJ=sumJ, J=fl.J(), covered=len(cov), n=len(fl.routes), pred_sumJ=plan['pred_sumJ'],
                   miss=miss, per_route={n: dict(tank=x['tank'], fb=len(x['st']['asts']))
                                         for n, x in fl.routes.items()}),
              open(out / 'built.json', 'w'), indent=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['bases', 'milp', 'greedy', 'build'])
    ap.add_argument('lib', nargs='?', default='results/s7/libopt.pkl')
    ap.add_argument('out', nargs='?', default='results/s10/milp')
    ap.add_argument('--nproc', type=int, default=10)
    ap.add_argument('--dmax', type=float, default=0.10, help='max closest approach admitted as a column [AU]')
    ap.add_argument('--pmin', type=float, default=0.30, help='predict_ok floor for a column')
    ap.add_argument('--band', type=int, default=4, help='slots per depth band')
    ap.add_argument('--max-added', type=int, default=32)
    ap.add_argument('--max-fuel', type=float, default=700.0, help='kg above dry mass (tank cap 600+this)')
    ap.add_argument('--tan-step', type=float, default=60.0)
    ap.add_argument('--kg-floor', type=float, default=3.0)
    ap.add_argument('--risk', action='store_true', help='price kg/p instead of kg')
    ap.add_argument('--ncar', type=int, default=0, help='keep only the best N carriers')
    ap.add_argument('--fix-n', type=int, default=None)
    ap.add_argument('--tlim', type=float, default=600.0)
    ap.add_argument('--gap', type=float, default=0.005)
    ap.add_argument('--plan', default='plan.json'); ap.add_argument('--tag', default='')
    ap.add_argument('--max-tank', type=float, default=1400.0)
    ap.add_argument('--ntry', type=int, default=3)
    ap.add_argument('--tmax', type=float, default=1800.0)
    ap.add_argument('--nmin', type=int, default=8); ap.add_argument('--nmax', type=int, default=14)
    ap.add_argument('--ls-rounds', type=int, default=3); ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--kcal', type=float, default=1.0,
                    help='multiplier on predict_kg; the MEASURED value for assigned insertions is 1.57 '
                         'with --mshrink 1 (1.94 without) -- see the module docstring')
    ap.add_argument('--mshrink', type=float, default=1.0,
                    help='exponent on (n0/depth) applied to the base-derived margin as the route fills')
    a = ap.parse_args()
    if a.cmd == 'bases':
        a.dmax = max(a.dmax, 0.25)
        cmd_bases(a)
    elif a.cmd == 'milp':
        cmd_milp(a)
    elif a.cmd == 'greedy':
        cmd_greedy(a)
    else:
        cmd_build(a)


if __name__ == '__main__':
    main()
