"""Stage 14, track A, stage 1 -- THE CEILING TEST: a beam over a prescribed target pool whose state velocity is the
impulsive TWIN's velocity at the last flyby, not the planner's own Lambert/linear arrival velocity.

State   a SETTLED impulsive-twin prefix (ctoc14/impulsive.ImpulsiveProblem: launch v_inf + 20-d impulse nodes + Kepler
        arcs, flyby epochs free) and its end state (t, r, v) at the last flyby, read from ImpulsiveProblem.integrate.
Expand  ctoc14.search.expand (Lambert 15-400 d with tof refinement + LinLeg 20-600 d, production BP settings) from a
        search.State built on the twin end state; prize 1 on the pool, every other target excluded.
Settle  child twin = parent's settled twin + one flyby (globalopt.add_flyby semantics: the flyby list stays sorted) +
        regular nodes up to the new epoch + a leg seed (Lambert junction impulse 1 d after the flyby from the coasted
        state, or the LinLeg thrust profile, cheaper first, the other on failure), then run_ialns.settle (restore, L1 SCP
        over the WHOLE prefix, restore, thrust-cap homotopy).  Settled = max miss <= 150 km.
Score   search.State.score: t_last/T + w_fuel F/1400 - prizes, F = 1600 (1 - exp(-dv_twin/ve)) (the planner's m0-1600
        propellant convention applied to the twin's delta-v; w_fuel = CM.w_fuel(480, 1600) = 0.52 as in production).
Modes
  free    depth-synchronous beam.  Children of all beam states are ranked by the leg-model score, deduplicated on
          (visited set, 5-d bin), settled in rank order until --beam succeed or --tries attempts are spent; survivors are
          re-ranked by their TWIN score.  Roots = search.roots (launch grid, v_inf <= 4): exact Lambert arcs, zero thrust.
  guided  diagnostic along the SOURCE route's own order: at every step, from the myopic prefix twin, is the next source
          target among expand()'s children (production and relaxed limits)?  single-leg costs at the source epoch; then
          the flyby is added and settled (at the admissible child's epoch, else the source epoch) and the chain goes on.
Outputs <out>/log.txt, <out>/ckpt.pkl (restart point, rewritten after every level), <out>/result.json,
        <out>/best.npz (deepest settled twin, run_ialns ist format)
Usage   s14_twinbeam.py free   SRC_FLEET ROUTE OUT [--dv 1.2 --drmax 0.15 --beam 20 --tries 40]
        s14_twinbeam.py guided SRC_FLEET ROUTE OUT
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, signal, pathlib, argparse, pickle, warnings
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from greedy_cover import BP, ALL, CM, _timeout
from ctoc14.search import State, expand, roots, mask_of, UNREACHABLE
from ctoc14.impulsive import ImpulsiveProblem, LIN_MISS_MAX
from ctoc14.kepler import propagate_batch
from ctoc14.lambert import lambert, lambert_best
from ctoc14.linleg import LinLeg
from ctoc14.globalopt import irls_min_fuel
from ctoc14.constants import DAY, AU, VE, TMAX, T_MISSION, M_DRY, FUEL_MAX, VINF_MAX
warnings.filterwarnings('ignore')

M0P = 1600.0              # planner launch mass (production)
DTAU = 20 * DAY           # twin node spacing
TOL = 150.0               # km: settled
LAG = 1.0 * DAY           # junction impulse lag after the departure flyby (impulsive.from_tour default)


def params(dv, drmax, mc=60):
    """Production planner parameters (linchpin_replan / greedy_cover) with leg cap dv [km/s] and linear drmax [AU]."""
    return BP(beam=20, w_fuel=CM.w_fuel(480.0, M0P), m_margin=40.0, dv_max=dv,
              lin_tofs=np.arange(20, 600 + 1e-9, 10) * DAY, tofs=np.arange(15, 400 + 1e-9, 5) * DAY,
              lin_drmax=drmax * AU, vinf_cap=4.0, tof_refine=True, max_depth=60, max_children=mc)


# ------------------------------------------------------------------ twin prefix states
def centres(tL, t_hi):
    edges = np.arange(tL, t_hi + DTAU, DTAU)
    return 0.5 * (edges[:-1] + edges[1:])


def pack(ip, nreg, miss, how=''):
    """Twin state dict from a settled ImpulsiveProblem."""
    Yf, _ = ip.integrate()
    j = int(np.argmax(ip.tf))
    o = np.argsort(ip.tf)
    return dict(st=RI.ist(ip), nreg=int(nreg), tank=float(ip.tank()), dv=float(ip.dv()), miss=float(miss),
                t_end=float(ip.tf[j]), r_end=Yf[j, :3].copy(), v_end=Yf[j, 3:6].copy(),
                asts=[int(ip.asts[i]) for i in o], tfs=[float(ip.tf[i]) for i in o], how=how)


def root_state(tL, vinf, ast, t1):
    """Exact zero-thrust twin of a launch leg (|v_inf| <= 4: the Lambert arc itself)."""
    cen = centres(tL, t1)
    ip = ImpulsiveProblem(RI.eph(), tL, vinf, cen, np.zeros((len(cen), 3)), [t1], [int(ast)])
    Yf, _ = ip.integrate(); dm, _ = ip.misses(Yf)
    return pack(ip, len(cen), float(np.linalg.norm(dm, axis=1).max()), 'root')


def search_state(ts_, excl_mask, P):
    m = M0P * np.exp(-ts_['dv'] / VE)
    n = len(ts_['asts'])
    return State(t=ts_['t_end'], r=np.array(ts_['r_end']), v=np.array(ts_['v_end']), m=float(m),
                 visited=excl_mask | mask_of(ts_['asts']),
                 seq=tuple((a, t, 0.0, 0.0) for a, t in zip(ts_['asts'], ts_['tfs'])), fuel=float(M0P - m),
                 t_launch=float(ts_['st']['tL']), vinf=np.array(ts_['st']['vinf']), m0=M0P, rare=float(n))


def twin_score(ts_, P):
    F = M0P * (1.0 - np.exp(-ts_['dv'] / VE))
    return P.w_t * ts_['t_end'] / T_MISSION + P.w_fuel * F / FUEL_MAX - len(ts_['asts'])


def seeds(t0, r0, v0, R, tof):
    """Leg seeds (cost km/s, impulse epochs, impulses) for (t0, r0, v0) -> R at t0 + tof, cheapest first:
    Lambert (0-2 rev, best junction) from the state coasted LAG past the flyby; LinLeg distributed profile."""
    out = []
    lag = min(LAG, 0.1 * tof)
    rc, vc = propagate_batch(r0[None], v0[None], np.array([lag]))
    v1, _, _, _ = lambert_best(rc, R[None], np.array([tof - lag]), vc, nrev_max=2)
    if np.all(np.isfinite(v1[0])):
        d = v1[0] - vc[0]
        out.append(('lam', float(np.linalg.norm(d)), np.array([t0 + lag]), d[None].copy()))
    try:
        L = LinLeg(r0, v0, np.array([tof]), dtau=10 * DAY)
        c, _, _, U = L.solve(0, R[None], TMAX / 1000.0 * 1e-3, ub=0.9, return_U=True)
        if np.isfinite(c[0]) and L.check(0, U[0], R) < LIN_MISS_MAX:
            K = int(L.K[0]); keep = [k for k in range(K) if np.linalg.norm(U[0][k]) * L.dtau > 1e-6]
            if keep:
                out.append(('lin', float(c[0]), t0 + L.taus[keep], U[0][keep] * L.dtau))
    except Exception:
        pass
    out.sort(key=lambda x: x[1])
    return out


def child_problem(par, ast, t_new, seed):
    st = par['st']; tL = float(st['tL'])
    cen = centres(tL, t_new); nreg = par['nreg']
    add_t = [cen[nreg:]]; add_T = [np.zeros((max(0, len(cen) - nreg), 3))]
    if seed is not None:
        add_t.append(np.asarray(seed[2], float)); add_T.append(np.asarray(seed[3], float).reshape(-1, 3))
    ts = np.concatenate([np.asarray(st['ts'], float)] + add_t)
    Ts = np.concatenate([np.asarray(st['Ts'], float).reshape(-1, 3)] + add_T)
    o = np.argsort(ts, kind='stable'); ts = ts[o]; Ts = Ts[o]
    tf = np.append(np.asarray(st['tf'], float), t_new); asts = [int(a) for a in st['asts']] + [int(ast)]
    o2 = np.argsort(tf, kind='stable')
    ip = ImpulsiveProblem(RI.eph(), tL, np.asarray(st['vinf'], float), ts, Ts, tf[o2], [asts[i] for i in o2])
    return ip, max(nreg, len(cen))


def settle_child(job):
    """job = (parent twin state, ast, t_new, iters, timeout[, lin]). lin = linearised-twin step (T, dt, dvinf) on the
    zero-impulse extension (lin_price): used as the FIRST seed when given. Returns (twin state or None, info dict)."""
    par, ast, t_new, iters, timeout = job[:5]; lin = job[5] if len(job) > 5 else None
    E = RI.eph(); tic = time.time()
    R = E.ast_states_at(np.array([ast - 1]), np.array([t_new]))[0][0]
    sd = seeds(par['t_end'], np.asarray(par['r_end']), np.asarray(par['v_end']), R, t_new - par['t_end'])
    info = dict(ast=int(ast), t=float(t_new), seeds=[(s[0], round(s[1], 3)) for s in sd], tries=[])
    attempts = ([('lintwin', None)] if lin is not None else []) + [(s[0], s) for s in sd]
    for kind, s in attempts[:2]:
        if kind == 'lintwin':
            ip, nreg = child_problem(par, ast, t_new, None)
            ip.Ts = np.array(lin['T'], float); ip.tf = np.minimum(ip.tf + np.asarray(lin['dt']) * DAY, T_MISSION - 2 * DAY)
            vi = ip.vinf + np.asarray(lin['dvinf']); nv = np.linalg.norm(vi)
            ip.vinf = vi * min(1.0, (VINF_MAX - 1e-6) / nv) if nv > 0 else vi
        else:
            ip, nreg = child_problem(par, ast, t_new, s)
        try:
            signal.signal(signal.SIGALRM, _timeout); signal.setitimer(signal.ITIMER_REAL, timeout)
            miss = float(RI.settle(ip, iters))
        except Exception as e:
            info['tries'].append((kind, type(e).__name__)); continue
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
        info['tries'].append((kind, round(miss, 1)))
        if miss <= TOL:
            out = pack(ip, nreg, miss, kind)
            info.update(ok=True, tank=out['tank'], dv_marg=round(out['dv'] - par['dv'], 4), dt_s=round(time.time() - tic, 1))
            return out, info
    info.update(ok=False, dt_s=round(time.time() - tic, 1))
    return None, info


# ------------------------------------------------------------------ linearised-twin leg pricing (gate lintwin)
def epoch_candidates(par, remaining, K=3, sep=30 * DAY):
    """Candidate epochs per remaining target: the K lowest local minima (>= sep apart) of the single-leg cost curves
    from the twin end state -- single-rev Lambert junction over 15-400 d / 5 d and LinLeg L1 cost over 20-600 d / 10 d,
    WITHOUT any dv / eta / drmax gate (they only locate encounter windows). Returns [(ast, t, single-leg cost)]."""
    E = RI.eph(); t0 = par['t_end']; r0 = np.asarray(par['r_end']); v0 = np.asarray(par['v_end'])
    idx = np.array(sorted(remaining), int) - 1
    if len(idx) == 0:
        return []
    curves = {int(i) + 1: [] for i in idx}
    tl = np.arange(15, 400 + 1e-9, 5) * DAY; tl = tl[t0 + tl <= T_MISSION - 2 * DAY]
    tn = np.arange(20, 600 + 1e-9, 10) * DAY; tn = tn[t0 + tn <= T_MISSION - 2 * DAY]
    nX = len(idx)
    if len(tl):
        R = E.ast_states_at(idx[None, :], (t0 + tl)[:, None])[0]
        v1, _ = lambert(np.broadcast_to(r0, (len(tl) * nX, 3)), R.reshape(-1, 3), np.repeat(tl, nX))
        dv = np.linalg.norm(v1 - v0, axis=-1).reshape(len(tl), nX)
        dv = np.where(np.isfinite(dv), dv, np.inf)
        for q in range(nX): curves[int(idx[q]) + 1].append((tl, dv[:, q]))
    if len(tn):
        L = LinLeg(r0, v0, tn, dtau=10 * DAY); acc = 50.0 * TMAX / M0P * 1e-3       # loose bound: epochs only
        C = np.full((len(tn), nX), np.inf)
        for j in range(len(tn)):
            R = E.ast_states_at(idx, np.full(nX, t0 + tn[j]))[0]
            c, _, _ = L.solve(j, R, acc, ub=0.8); C[j] = c
        for q in range(nX): curves[int(idx[q]) + 1].append((tn, C[:, q]))
    out = []
    for X, cl in curves.items():
        mins = []
        for tt, cc in cl:
            for i in range(len(tt)):
                if not np.isfinite(cc[i]): continue
                lo = cc[i - 1] if i > 0 else np.inf; hi = cc[i + 1] if i + 1 < len(tt) else np.inf
                if cc[i] <= lo and cc[i] <= hi: mins.append((float(cc[i]), float(t0 + tt[i])))
        mins.sort(); keep = []
        for c, t in mins:
            if all(abs(t - u) >= sep for _, u in keep): keep.append((c, t))
            if len(keep) >= K: break
        out += [(X, t, c) for c, t in keep]
    return out


def lin_price(par, ast, t):
    """Marginal delta-v [km/s] of adding (ast, t) to the settled prefix by the LINEARISED TWIN: the zero-impulse
    extension is linearised (all earlier impulses, flyby epochs and v_inf free) and one globalopt.irls_min_fuel
    min-L1 solve is taken.  Returns (marginal dv, step dict for settle_child, coast miss AU)."""
    ip, _ = child_problem(par, ast, t, None)
    Yf, Ys = ip.integrate(); dm, vr = ip.misses(Yf)
    j = [i for i, a in enumerate(ip.asts) if a == ast][-1]; cm = float(np.linalg.norm(dm[j]) / AU)
    B, Lv = ip.linearise(Yf, Ys)
    # |v_inf| at the 4 km/s cap: the settle clips v_inf, so the linear price must not use it (frozen, conservative)
    wv = 1e6 if np.linalg.norm(ip.vinf) >= VINF_MAX - 0.05 else 1.0
    T, dt, dvi = irls_min_fuel(B, Lv, vr, dm, ip.Ts, ip.vinf, ip.tcap, rho=0.0, w_v=wv)
    # nonlinear check of the full linear step: the residual miss after ONE integration tells whether the linearised
    # price can be trusted (coast misses of 0.4-2.4 AU give garbage prices of 0.1-1.3 km/s that never settle)
    vi = ip.vinf + dvi; nv = np.linalg.norm(vi)
    if nv > VINF_MAX - 1e-6: vi = vi * (VINF_MAX - 1e-6) / nv
    tf2 = np.minimum(ip.tf + dt * DAY, T_MISSION - 2 * DAY)
    Yf2, _ = ip.integrate(Ts=T, tf=tf2, vinf=vi); dm2, _ = ip.misses(Yf2, tf2)
    res = float(np.linalg.norm(dm2, axis=1).max())
    return float(np.linalg.norm(T, axis=1).sum()) - par['dv'], dict(T=T, dt=dt, dvinf=dvi, res_km=res), cm


def _price_parent(job):
    """All lintwin candidates of one parent: [(score, ast, t, marginal dv, step, coast miss, single-leg cost)]."""
    par, pool, K, cap, wt, wf, res_max = job
    rem = [x for x in pool if x not in set(par['asts'])]
    out = []
    for X, t, c1 in epoch_candidates(par, rem, K):
        try:
            dvm, lin, cm = lin_price(par, X, t)
        except Exception:
            continue
        if not np.isfinite(dvm) or dvm > cap or not (lin['res_km'] <= res_max):
            continue
        F = M0P * (1.0 - np.exp(-(par['dv'] + max(dvm, 0.0)) / VE))
        sc = wt * t / T_MISSION + wf * F / FUEL_MAX - (len(par['asts']) + 1)
        out.append((sc, X, t, dvm, lin, cm, c1))
    return out


def _pool_map(fn, jobs, nproc):
    if nproc <= 1 or len(jobs) <= 1:
        return [fn(j) for j in jobs]
    import multiprocessing as mp
    with mp.get_context('fork').Pool(min(nproc, len(jobs))) as pool:
        return pool.map(fn, jobs, chunksize=1)


def save_ckpt(out, obj):
    tmp = out / 'ckpt.pkl.tmp'
    with open(tmp, 'wb') as f:
        pickle.dump(obj, f)
    os.replace(tmp, out / 'ckpt.pkl')


def best_of(a, b):
    """Deeper wins; equal depth: lighter tank."""
    if b is None: return a
    if a is None: return b
    ka, kb = (len(a['asts']), -a['tank']), (len(b['asts']), -b['tank'])
    return a if ka >= kb else b


# ------------------------------------------------------------------ free twin beam
def run_free(a, src, pool, say, out):
    E = RI.eph(); P = params(a.dv, a.drmax, a.mc)
    P.prize = np.zeros(300)
    for t in pool: P.prize[t - 1] = 1.0
    excl = sorted(set(ALL) - set(pool)); excl_mask = mask_of(set(excl) | UNREACHABLE)
    g = [float(x) for x in a.grid.split(',')]; grid = np.arange(g[0], g[1] + 1e-9, g[2]) * DAY
    ck = out / 'ckpt.pkl'
    if ck.exists():
        S = pickle.load(open(ck, 'rb'))
        say(f'RESUME from level {S["level"]}: beam {len(S["beam"])}, best {len(S["best"]["asts"])} @ {S["best"]["tank"]:.1f} kg')
    else:
        rts = roots(E, grid, M0P, excl, P)
        rts.sort(key=lambda s: s.score(P)); seen = set(); b2 = []
        for s in rts:
            key = (s.seq[0][0], int(s.t_launch // (20 * DAY)))
            if key in seen: continue
            seen.add(key); b2.append(s)
        b2 = b2[:a.roots]
        beam = [root_state(s.t_launch, s.vinf, s.seq[0][0], s.seq[0][1]) for s in b2]
        best = None
        for b in beam: best = best_of(best, b)
        S = dict(level=1, beam=beam, best=best, hist=[], t_start=time.time(), wall=0.0,
                 settles=0, settle_ok=0, settle_s=0.0, end=None)
        say(f'roots: {len(rts)} launch legs -> {len(beam)} kept (max root miss {max(b["miss"] for b in beam):.1f} km)')
        save_ckpt(out, S)
    while S['end'] is None:
        tic = time.time(); beam = S['beam']
        kids = []                                  # (score, parent, ast, t, visited mask, extra)
        if a.gate == 'leg':
            for i, ts_ in enumerate(beam):
                for c in expand(E, search_state(ts_, excl_mask, P), P):
                    kids.append((c.score(P), i, int(c.seq[-1][0]), float(c.seq[-1][1]), c.visited,
                                 dict(dv_leg=round(c.seq[-1][2], 3))))
        else:
            res = _pool_map(_price_parent, [(ts_, pool, a.k_epochs, a.lin_cap, P.w_t, P.w_fuel, a.res_max) for ts_ in beam],
                            a.nproc)
            for i, (ts_, lst) in enumerate(zip(beam, res)):
                vm = excl_mask | mask_of(ts_['asts'])
                for sc, X, t, dvm, lin, cm, c1 in lst:
                    kids.append((sc, i, int(X), float(t), vm | (1 << (int(X) - 1)),
                                 dict(dv_lin=round(dvm, 4), coast_au=round(cm, 4), dv_leg1=round(c1, 3),
                                      res_km=round(lin['res_km']), lin=lin)))
        t_gen = time.time() - tic
        kids.sort(key=lambda x: x[0])
        seen = set(); cand = []
        for kd in kids:
            key = (kd[4], int(kd[3] // (5 * DAY)))
            if key in seen: continue
            seen.add(key); cand.append(kd)
        new = []; infos = []; tried = 0; nk = len(kids); nu = len(cand)
        while cand and len(new) < a.beam and tried < a.tries:
            take = cand[:max(1, min(a.nproc, a.tries - tried))]; cand = cand[len(take):]
            jobs = [(beam[i], X, t, a.iters, a.timeout) + ((ex['lin'],) if 'lin' in ex else ()) for sc, i, X, t, vm, ex in take]
            res = _pool_map(settle_child, jobs, a.nproc)
            tried += len(take)
            for (sc, i, X, t, vm, ex), (tw, info) in zip(take, res):
                info['parent'] = i; info['score'] = round(sc, 5)
                info.update({k: v for k, v in ex.items() if k != 'lin'})
                infos.append(info)
                S['settles'] += 1; S['settle_s'] += info['dt_s']
                if tw is not None:
                    S['settle_ok'] += 1; new.append(tw)
        new.sort(key=lambda t: twin_score(t, P))
        seen = set(); nb = []
        for t in new:
            key = (mask_of(t['asts']), int(t['t_end'] // (5 * DAY)))
            if key in seen: continue
            seen.add(key); nb.append(t)
        nb = nb[:a.beam]
        depth = S['level'] + 1
        rec = dict(depth=depth, parents=len(beam), children=nk, unique=nu, tried=tried, settled=len(new), kept=len(nb),
                   wall_s=round(time.time() - tic, 1), gen_s=round(t_gen, 1),
                   calib=[(i.get('dv_lin', i.get('dv_leg')), i.get('dv_marg')) for i in infos if i.get('ok')][:40],
                   parent_t_end_d=[round(min(b['t_end'] for b in beam) / DAY, 1), round(max(b['t_end'] for b in beam) / DAY, 1)],
                   fails=[i for i in infos if not i.get('ok')][:12])
        if nb:
            rec.update(t_end_d=[round(min(b['t_end'] for b in nb) / DAY, 1), round(max(b['t_end'] for b in nb) / DAY, 1)],
                       tank=[round(min(b['tank'] for b in nb), 1), round(max(b['tank'] for b in nb), 1)])
            for t in nb: S['best'] = best_of(S['best'], t)
        S['hist'].append(rec)
        say(f'depth {depth:2d}: {len(beam)} parents -> {nk} children ({nu} unique) -> tried {tried}, settled {len(new)}, '
            f'kept {len(nb)}' + (f'; t_end {rec["t_end_d"][0]:.0f}-{rec["t_end_d"][1]:.0f} d, tank '
                                 f'{rec["tank"][0]:.0f}-{rec["tank"][1]:.0f}' if nb else '') +
            f'; best {len(S["best"]["asts"])} @ {S["best"]["tank"]:.1f} kg ({rec["wall_s"]:.0f} s)')
        if not nb:
            if nk == 0:
                tl = T_MISSION - max(b['t_end'] for b in beam)
                S['end'] = dict(kind='no_children', depth=S['level'], time_left_d=round(tl / DAY, 1))
            else:
                S['end'] = dict(kind='settle_wall', depth=S['level'], tried=tried)
        elif len(nb[0]['asts']) >= len(pool):
            S['end'] = dict(kind='pool_exhausted', depth=depth)
        S['beam'] = nb if nb else beam; S['level'] = depth if nb else S['level']
        S['wall'] = S.get('wall', 0.0) + (time.time() - tic)
        save_ckpt(out, S)
    return S


# ------------------------------------------------------------------ guided diagnostic
def leg_costs(ts_, ast, t2, acc):
    """legcheck-style single-leg costs from the twin end state to ast at epoch t2."""
    E = RI.eph(); r1 = np.asarray(ts_['r_end']); v1 = np.asarray(ts_['v_end']); tof = t2 - ts_['t_end']
    if tof <= 0:
        return dict(tof_d=round(tof / DAY, 1))
    R2 = E.ast_states_at(np.array([ast - 1]), np.array([t2]))[0][0]
    rc, _ = propagate_batch(r1[None], v1[None], np.array([tof])); miss = float(np.linalg.norm(rc[0] - R2) / AU)
    vl, _ = lambert(r1[None], R2[None], np.array([tof])); dv1 = float(np.linalg.norm(vl[0] - v1))
    vb, _, nrev, _ = lambert_best(r1[None], R2[None], np.array([tof]), v1[None], nrev_max=2)
    dvb = float(np.linalg.norm(vb[0] - v1))
    L = LinLeg(r1, v1, np.array([tof]), dtau=10 * DAY); c, _, _ = L.solve(0, R2[None], acc, ub=0.8)
    return dict(tof_d=round(tof / DAY, 1), coast_miss_au=round(miss, 4), lam1_dv=round(dv1, 3), lam_best_dv=round(dvb, 3),
                lin_dv=round(float(c[0]), 3) if np.isfinite(c[0]) else None)


def run_guided(a, src, say, out):
    E = RI.eph()
    Pp = params(1.2, 0.15, mc=1000); Pr = params(2.0, 0.25, mc=1000)
    st = src; o = np.argsort(st['tf']); asts = [int(st['asts'][i]) for i in o]; tfs = [float(st['tf'][i]) for i in o]
    pool = sorted(asts)
    for P in (Pp, Pr):
        P.prize = np.zeros(300)
        for t in pool: P.prize[t - 1] = 1.0
    excl_mask = mask_of((set(ALL) - set(pool)) | UNREACHABLE)
    lc = {}
    try:
        L = json.load(open(ROOT / 'results/s13/linchpin/legcheck_t10d.json'))['routes'][a.route]['legs']
        lc = {l['k']: l for l in L}
    except Exception:
        pass
    ck = out / 'ckpt.pkl'
    if ck.exists():
        S = pickle.load(open(ck, 'rb')); say(f'RESUME guided at step {S["k"]}')
    else:
        tL = float(st['tL']); rE, vE = E.earth_state(tL)
        R1 = E.ast_states_at(np.array([asts[0] - 1]), np.array([tfs[0]]))[0][0]
        v1, _ = lambert(rE[None], R1[None], np.array([tfs[0] - tL]))
        vinf = v1[0] - vE; nv = np.linalg.norm(vinf)
        root = root_state(tL, vinf * min(1.0, VINF_MAX / nv), asts[0], tfs[0])
        if nv > VINF_MAX:
            ip = RI.ipr(root['st']); miss = RI.settle(ip, a.iters); root = pack(ip, root['nreg'], miss, 'root_settled')
        S = dict(k=1, cur=root, steps=[dict(k=0, ast=asts[0], vinf=round(float(nv), 3), miss=round(root['miss'], 1))],
                 skipped=[], t_start=time.time(), wall=0.0)
        say(f'guided root: ast {asts[0]} at {tfs[0] / DAY:.1f} d, |v_inf| {nv:.3f} km/s')
        save_ckpt(out, S)
    while S['k'] < len(asts):
        tic = time.time(); k = S['k']; cur = S['cur']; tgt = asts[k]; t_src = tfs[k]
        rec = dict(k=k, ast=tgt, t_src_d=round(t_src / DAY, 1), t_prev_twin_d=round(cur['t_end'] / DAY, 1),
                   src_prev_d=round(tfs[k - 1] / DAY, 1))
        s = search_state(cur, excl_mask, Pp)
        acc = TMAX / s.m * 1e-3
        rec['myopic'] = leg_costs(cur, tgt, t_src, acc)
        opts = []
        for tag, P in (('prod', Pp), ('rel', Pr)):
            ch = [c for c in expand(E, s, P) if c.seq[-1][0] == tgt]
            if ch:
                cb = min(ch, key=lambda c: abs(c.seq[-1][1] - t_src))
                cc = min(ch, key=lambda c: c.seq[-1][2])
                rec[tag] = dict(ok=True, n=len(ch), dv_near=round(cb.seq[-1][2], 3), dt_near_d=round((cb.seq[-1][1] - t_src) / DAY, 1),
                                dv_min=round(cc.seq[-1][2], 3), dt_min_d=round((cc.seq[-1][1] - t_src) / DAY, 1))
                opts.append((tag, cb.seq[-1][1]))
            else:
                rec[tag] = dict(ok=False)
        # lintwin view of the true next leg: nearest single-leg encounter epoch, linearised-twin price there and at the
        # source epoch (all earlier impulses free)
        try:
            ec = epoch_candidates(cur, [tgt], K=5)
            if ec:
                X_, te, c1 = min(ec, key=lambda z: abs(z[1] - t_src))
                rec['lt_epoch_dt_d'] = round((te - t_src) / DAY, 1); rec['lt_epoch_leg1'] = round(c1, 3)
                rec['lt_price_epoch'] = round(lin_price(cur, tgt, te)[0], 4)
                rec['lt_epochs_d'] = [round((z[1] - t_src) / DAY, 1) for z in ec]
            rec['lt_price_src'] = round(lin_price(cur, tgt, t_src)[0], 4)
        except Exception as e:
            rec['lt_error'] = repr(e)[:80]
        l = lc.get(k)
        if l:
            rec['fullroute'] = dict(ok_planner=l['ok_planner'], ok_relaxed=l['ok_relaxed'], lam1_dv=l['lam1_dv'],
                                    lin_dv=l['lin_dv'], coast_miss_au=l['coast_miss_au'], dv_twin=l['dv_twin'])
        # extend the twin: admissible child's epoch first, then the source epoch
        tlist = []
        for tag, t in opts + [('src', t_src)]:
            if t > cur['t_end'] + 5 * DAY and all(abs(t - u) > 0.5 * DAY for _, u in tlist):
                tlist.append((tag, t))
        got = None; tries = []
        for tag, t in tlist[:2]:
            tw, info = settle_child((cur, tgt, t, a.iters, a.timeout))
            tries.append(dict(at=tag, t_d=round(t / DAY, 1), **{kk: info.get(kk) for kk in ('ok', 'tank', 'dt_s', 'tries', 'seeds')}))
            if tw is not None:
                got = tw; break
        rec['settle'] = tries
        if got is not None:
            rec.update(settled=True, tank=round(got['tank'], 2), dtank=round(got['tank'] - cur['tank'], 2),
                       dv_marg=round(got['dv'] - cur['dv'], 4),
                       t_twin_d=round(max(got['tfs']) / DAY, 1))
            S['cur'] = got
        else:
            rec['settled'] = False; S['skipped'].append(tgt)
        rec['wall_s'] = round(time.time() - tic, 1)
        S['steps'].append(rec); S['k'] = k + 1; S['wall'] += time.time() - tic
        say(f'leg {k:2d} -> {tgt:3d} (src {t_src / DAY:6.1f} d, tof_src {(t_src - tfs[k - 1]) / DAY:5.1f}): myopic lam1 '
            f'{rec["myopic"].get("lam1_dv")} lin {rec["myopic"].get("lin_dv")} miss {rec["myopic"].get("coast_miss_au")} | '
            f'prod {"Y" if rec["prod"]["ok"] else "n"} rel {"Y" if rec["rel"]["ok"] else "n"} | full-route '
            f'{"Y" if l and l["ok_planner"] else "n"}{"Y" if l and l["ok_relaxed"] else "n"} | twin '
            + (f'OK tank {rec["tank"]:.1f} (+{rec["dtank"]:.1f})' if rec['settled'] else 'FAIL') + f' ({rec["wall_s"]:.0f} s)')
        save_ckpt(out, S)
    return S


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['free', 'guided']); ap.add_argument('src'); ap.add_argument('route'); ap.add_argument('out')
    ap.add_argument('--dv', type=float, default=1.2); ap.add_argument('--drmax', type=float, default=0.15)
    ap.add_argument('--beam', type=int, default=20); ap.add_argument('--tries', type=int, default=40)
    ap.add_argument('--roots', type=int, default=60); ap.add_argument('--mc', type=int, default=60)
    ap.add_argument('--grid', default='0,800,20'); ap.add_argument('--iters', type=int, default=100)
    ap.add_argument('--timeout', type=float, default=240.0); ap.add_argument('--nproc', type=int, default=1)
    ap.add_argument('--gate', default='leg', choices=['leg', 'lintwin'],
                    help='leg: children from search.expand at --dv/--drmax (the spec); lintwin: every remaining pool '
                         'target at its single-leg encounter epochs, priced by the linearised whole-prefix twin')
    ap.add_argument('--k-epochs', type=int, default=3); ap.add_argument('--lin-cap', type=float, default=3.0)
    ap.add_argument('--res-max', type=float, default=3.0e6,
                    help='lintwin: max miss [km] left after one integration of the linear step (linearisation validity)')
    a = ap.parse_args()
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt')
    F = RI.IFleet(a.src); src = F.routes[a.route]['st']; tank0 = F.routes[a.route]['tank']
    pool = sorted(int(x) for x in src['asts'])
    say(f'=== s14_twinbeam {a.mode} route {a.route}: pool {len(pool)} (own targets), source {len(pool)} fb @ {tank0:.1f} kg; '
        + (f'gate {a.gate} dv {a.dv} drmax {a.drmax} beam {a.beam} tries {a.tries} roots {a.roots}'
           + (f' k_epochs {a.k_epochs} lin_cap {a.lin_cap} res_max {a.res_max:.3g}' if a.gate == 'lintwin' else '')
           if a.mode == 'free' else ''))
    t0 = time.time()
    if a.mode == 'free':
        S = run_free(a, src, pool, say, out)
        b = S['best']; got = sorted(set(b['asts']) & set(pool))
        o = np.argsort(src['tf']); src_t = [float(src['tf'][i]) for i in o]
        res = dict(mode='free', gate=a.gate, route=a.route, pool_size=len(pool), dv=a.dv, drmax=a.drmax, beam=a.beam, tries=a.tries,
                   reflown=len(got), frac=round(len(got) / len(pool), 4), tank=round(b['tank'], 2), src_tank=round(tank0, 2),
                   miss_km=round(b['miss'], 1), t_launch_d=round(float(b['st']['tL']) / DAY, 1),
                   t_end_d=round(max(b['tfs']) / DAY, 1), src_t_end_d=round(src_t[-1] / DAY, 1),
                   pace_d_per_fb=round((max(b['tfs']) - float(b['st']['tL'])) / DAY / len(b['asts']), 1),
                   src_pace_d_per_fb=round((src_t[-1] - float(src['tL'])) / DAY / len(pool), 1),
                   missing=sorted(set(pool) - set(got)), order=b['asts'], epochs_d=[round(t / DAY, 1) for t in b['tfs']],
                   src_epochs_d=[round(t / DAY, 1) for t in src_t], end=S['end'],
                   settles=S['settles'], settle_ok=S['settle_ok'], settle_cpu_s=round(S['settle_s'], 1),
                   wall_s=round(S['wall'], 1), hist=S['hist'])
        np.savez(out / 'best.npz', **b['st'])
        say(f'RESULT route {a.route}: re-flown {len(got)}/{len(pool)} at twin tank {b["tank"]:.1f} kg (source {tank0:.1f}); '
            f'end {S["end"]}; settles {S["settles"]} ({S["settle_ok"]} ok), wall {S["wall"]:.0f} s')
    else:
        S = run_guided(a, src, say, out)
        st = S['steps'][1:]
        n = len(st)
        cnt = lambda f: sum(1 for r in st if f(r))
        res = dict(mode='guided', route=a.route, legs=n, pool_size=len(pool), src_tank=round(tank0, 2),
                   adm_prod=cnt(lambda r: r['prod']['ok']), adm_rel=cnt(lambda r: r['rel']['ok']),
                   adm_full_prod=cnt(lambda r: r.get('fullroute', {}).get('ok_planner')),
                   adm_full_rel=cnt(lambda r: r.get('fullroute', {}).get('ok_relaxed')),
                   settled=cnt(lambda r: r['settled']), skipped=S['skipped'],
                   final_n=len(S['cur']['asts']), final_tank=round(S['cur']['tank'], 2), wall_s=round(S['wall'], 1),
                   steps=S['steps'])
        runs = []; c = 0
        for r in st:
            c = c + 1 if r['prod']['ok'] else 0; runs.append(c)
        res['longest_prod_run'] = max(runs) if runs else 0
        np.savez(out / 'final.npz', **S['cur']['st'])
        say(f'RESULT guided {a.route}: {n} legs; next source target admissible from the myopic twin: prod {res["adm_prod"]}, '
            f'relaxed {res["adm_rel"]} (full-route twin: {res["adm_full_prod"]}, {res["adm_full_rel"]}); twin settled '
            f'{res["settled"]}/{n}; final {res["final_n"]} fb @ {res["final_tank"]:.1f} kg (source {tank0:.1f}); '
            f'wall {S["wall"]:.0f} s')
    json.dump(res, open(out / 'result.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
