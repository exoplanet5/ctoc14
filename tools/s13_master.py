"""Stage 13 (PCC-8G) [M]: exact fleet selection over every settled route column on disk.
docs/stage13_solver_plan.md section 4.4.  A copy of tools/route_cover.py (`solve`) with its defects fixed:

  * miss-checked columns: load() integrates the stored impulsive state and drops it if the flyby miss > 150 km, the
    tank is outside [600, --max-tank], or it has no flybys (route_cover.w_load never checked the stored miss);
  * sum y <= N: the craft count is a CONSTRAINT (--N 8,9 solves each), and the LP value (same constraints) is reported;
  * --exact re-settles the picks and REMOVES a pick that fails before re-solving (route_cover kept failed picks in
    the column list and could re-pick them; it also saved the ORIGINAL state with the RE-SETTLED tank);
  * --strip: the fleet_dedupe bulk rule (RI.w_remove) when the picked fleet covers a target twice.

    min  sum_r cost(tank_r) y_r + sum_t z_t   s.t.  sum_{r: t in r} y_r + z_t >= 1 (t in the 298),  sum_r y_r <= N

usage: s13_master.py OUT PAT [PAT ...] [--N 8,9] [--exact] [--strip] [--tlim 300] [--gap 1e-3] [--nproc 4]
                     [--max-tank 1400] [--exclude GLOB,GLOB]
Outputs: OUT/N<N>/fleet (RI.IFleet dir), OUT/N<N>/picked.json, OUT/master.json, OUT/log.txt."""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, glob, time, fnmatch, pathlib, argparse, collections, multiprocessing as mp
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
_CWD0 = os.getcwd()
os.chdir(ROOT)
import numpy as np
import run_ialns as RI
from ctoc14.search import UNREACHABLE
ALL = [t for t in range(1, 301) if t not in UNREACHABLE]
MISS_OK = 150.0


def _state(f):
    z = np.load(f); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL'])
    return st


def load(job):
    """-> (column dict or None, reason).  Column: f, tank, miss, asts (reachable targets, sorted)."""
    f, max_tank = job
    try:
        st = _state(f)
        if 'asts' not in st or len(st['asts']) == 0:
            return None, 'empty'
        ip = RI.ipr(st); tank = float(ip.tank())
        if not np.isfinite(tank) or tank < 600.0 or tank > max_tank:
            return None, 'tank'
        Yf, _ = ip.integrate(); dm, _ = ip.misses(Yf)
        miss = float(np.linalg.norm(dm, axis=1).max()) if len(dm) else 0.0
        if not np.isfinite(miss) or miss > MISS_OK:
            return None, 'miss'
        asts = sorted(set(int(x) for x in st['asts']) - set(UNREACHABLE))
        if not asts:
            return None, 'empty'
        return dict(f=str(f), tank=tank, miss=miss, asts=asts), 'ok'
    except Exception as e:
        return None, f'error {type(e).__name__}'


def w_resettle(f):
    st = _state(f); ip = RI.ipr(st)
    try:
        miss = float(RI.settle(ip, 100))
    except Exception as e:
        return f, None, float('inf'), type(e).__name__
    return f, RI.ist(ip), miss, None


def _matrix(cols):
    from scipy.sparse import coo_matrix
    nr = len(cols); ti = {t: i for i, t in enumerate(ALL)}; nt = len(ALL)
    R = []; C = []
    for j, x in enumerate(cols):
        for t in x['asts']:
            R.append(ti[t]); C.append(j)
    R += list(range(nt)); C += list(range(nr, nr + nt))
    A = coo_matrix((np.ones(len(R)), (R, C)), shape=(nt, nr + nt)).tocsr()
    c = np.concatenate([[RI.cost(x['tank']) for x in cols], np.ones(nt)])
    card = np.concatenate([np.ones(nr), np.zeros(nt)])[None, :]
    return A, c, card, nr, nt


def lp_value(cols, N):
    from scipy.optimize import linprog
    from scipy.sparse import vstack, csr_matrix
    A, c, card, nr, nt = _matrix(cols)
    Aub = vstack([-A, csr_matrix(card)]).tocsr(); bub = np.concatenate([-np.ones(nt), [N]])
    res = linprog(c, A_ub=Aub, b_ub=bub, bounds=(0, 1), method='highs')
    return float(res.fun) + 2.0 if res.status == 0 else float('nan')


def solve(cols, N, say, tlim=300.0, gap=1e-3):
    """HiGHS MILP; returns dict(pick, covered, miss, sumJi, J, status, mip_gap) or None."""
    from scipy.optimize import milp, LinearConstraint, Bounds
    A, c, card, nr, nt = _matrix(cols)
    cons = [LinearConstraint(A, np.ones(nt), np.full(nt, np.inf))]
    if N is not None:
        cons.append(LinearConstraint(card, -np.inf, N))
    res = milp(c, constraints=cons, integrality=np.ones(nr + nt), bounds=Bounds(0, 1),
               options=dict(time_limit=tlim, mip_rel_gap=gap))
    if res.x is None:
        say(f'  MILP N<={N} failed: {res.message}'); return None
    pick = [j for j in range(nr) if res.x[j] > 0.5]
    cov = set(t for j in pick for t in cols[j]['asts'])
    sj = sum(RI.cost(cols[j]['tank']) for j in pick)
    out = dict(pick=pick, covered=len(cov), miss=nt - len(cov), sumJi=sj, J=sj + (nt - len(cov)) + 2.0,
               status=res.message[:60], mip_gap=float(getattr(res, 'mip_gap', float('nan')) or 0.0))
    say(f'  MILP N<={N}: {res.message[:40]}: {len(pick)} craft, {len(cov)} covered, sum J_i {sj:.4f}, '
        f'J {out["J"]:.4f} (incl. 2 permanent misses)')
    return out


def solve_budget(cols, N, B, say, tlim=300.0, gap=1e-3):
    """Capacity frontier: MAX coverage with sum y <= N and sum cost(tank) y <= B (ties: the cheaper fleet).
    min -sum_t c_t + 1e-3 sum_r cost_r y_r   s.t.  c_t <= sum_{r: t in r} y_r,  sum y <= N,  sum cost y <= B."""
    from scipy.optimize import milp, LinearConstraint, Bounds
    from scipy.sparse import hstack, eye, csr_matrix
    A, c, card, nr, nt = _matrix(cols)
    Ay = A[:, :nr]; cost = c[:nr]
    M = hstack([Ay, -eye(nt)]).tocsr()                                   # sum_r y_r - c_t >= 0
    obj = np.concatenate([1e-3 * cost, -np.ones(nt)])
    cons = [LinearConstraint(M, np.zeros(nt), np.full(nt, np.inf)),
            LinearConstraint(csr_matrix(np.concatenate([np.ones(nr), np.zeros(nt)])[None, :]), -np.inf, N),
            LinearConstraint(csr_matrix(np.concatenate([cost, np.zeros(nt)])[None, :]), -np.inf, B)]
    res = milp(obj, constraints=cons, integrality=np.ones(nr + nt), bounds=Bounds(0, 1),
               options=dict(time_limit=tlim, mip_rel_gap=gap))
    if res.x is None:
        say(f'  MILP N<={N}, sum J_i<={B}: failed: {res.message}'); return None
    pick = [j for j in range(nr) if res.x[j] > 0.5]
    cov = set(t for j in pick for t in cols[j]['asts'])
    sj = sum(RI.cost(cols[j]['tank']) for j in pick)
    say(f'  MILP N<={N}, sum J_i<={B}: {res.message[:30]}: {len(pick)} craft, {len(cov)} covered, sum J_i {sj:.4f}')
    return dict(pick=pick, covered=len(cov), miss=nt - len(cov), sumJi=sj, J=sj + (nt - len(cov)) + 2.0,
                status=res.message[:60])


def strip(fl, pmap, say):
    """fleet_dedupe --bulk rule: keep every duplicated target on the host whose removal saves least, strip it from
    the others in one RI.w_remove per route (only removals saving > 0.5 kg); a route whose bulk removal does not
    settle is kept as is.  Coverage never shrinks (one host always keeps each target)."""
    cov = fl.coverage(); dup = {t: sorted(h) for t, h in cov.items() if len(h) > 1}
    if not dup:
        return dict(dup=0)
    jobs = [(h, fl.routes[h]['st'], [t]) for t, H in dup.items() for h in H]
    res = pmap(RI.w_remove, jobs); g = {}
    for r in res:
        if r.get('ok'):
            g[(r['name'], int(r['asts'][0]))] = fl.routes[r['name']]['tank'] - r['tank']
    drop = collections.defaultdict(list)
    for t, H in dup.items():
        pr = sorted((g.get((h, t), -1e9), h) for h in H)
        for _, h in pr[1:]:
            if g.get((h, t), 0.0) > 0.5:
                drop[h].append(t)
    n0 = len(cov); rep = dict(dup=len(dup), trials=len(jobs), drops={h: v for h, v in drop.items()}, applied={})
    out2 = pmap(RI.w_remove, [(h, fl.routes[h]['st'], v) for h, v in sorted(drop.items())])
    for r in out2:
        h = r['name']
        if not r.get('ok') or r['tank'] >= fl.routes[h]['tank']:
            say(f'  strip {h}: bulk removal of {len(r["asts"])} did not settle/improve -- kept'); continue
        say(f'  strip {h}: -{len(r["asts"])} flybys, tank {fl.routes[h]["tank"]:.1f} -> {r["tank"]:.1f} kg')
        rep['applied'][h] = (len(r['asts']), round(fl.routes[h]['tank'], 1), round(r['tank'], 1))
        fl.routes[h] = dict(st=r['st'], tank=r['tank'])
    assert len(fl.coverage()) == n0, 'strip shrank coverage'
    return rep


def exact_loop(C, solve_fn, pmap, verified, dropped, say, rounds=6):
    """Re-settle the picks of solve_fn(C) (RI.settle, 100 iters); a pick that does not re-settle is REMOVED from C and
    the problem re-solved; a lighter re-settled state replaces the column's state and tank."""
    r = solve_fn(C)
    for it in range(rounds):
        if r is None:
            return None
        todo = [C[j]['f'] for j in r['pick'] if C[j]['f'] not in verified]
        if not todo:
            return r
        say(f'  exact: re-settling {len(todo)} picks')
        changed = False
        for f, st, miss, err in pmap(w_resettle, todo):
            j = next(i for i, x in enumerate(C) if x['f'] == f)
            if st is None or not (miss <= MISS_OK):
                say(f'  exact: {f} does NOT re-settle (miss {miss:.0f} km{", " + err if err else ""}) -> dropped')
                dropped.append(f); C.pop(j); changed = True; continue
            tank = float(RI.ipr(st).tank())
            verified[f] = dict(miss=miss, tank=tank)
            if tank < C[j]['tank'] - 1e-6:
                if tank < C[j]['tank'] - 1.0:
                    say(f'  exact: {f}: tank {C[j]["tank"]:.1f} -> {tank:.1f} kg'); changed = True
                C[j]['tank'] = tank; C[j]['st'] = st
        if not changed:
            return r
        r = solve_fn(C)
    return r


def emit(d, r, C, a, pmap, say, note):
    """Save the picked fleet (the state that matches the tank), optionally stripped; returns the report record."""
    fl = RI.IFleet(); picks = []
    for i, j in enumerate(sorted(r['pick'], key=lambda j: (-len(C[j]['asts']), C[j]['tank']))):
        st = C[j].get('st') or _state(C[j]['f'])
        fl.routes[f'{i + 1:02d}'] = dict(st=st, tank=float(C[j]['tank']))
        picks.append(dict(f=C[j]['f'], n=len(C[j]['asts']), tank=round(C[j]['tank'], 3),
                          resettled=bool(C[j].get('st') is not None)))
    srep = None
    dup = sorted(t for t, h in fl.coverage().items() if len(h) > 1)     # multiply covered, before any strip
    if a.strip and dup:
        say(f'  strip: the picked fleet covers {len(dup)} targets more than once')
        srep = strip(fl, pmap, say)
    fl.save(d / 'fleet', note=note)
    json.dump(picks, open(d / 'picked.json', 'w'), indent=1)
    sj = sum(RI.cost(x['tank']) for x in fl.routes.values()); cv = len(fl.coverage())
    say(f'  {note}: saved {len(fl.routes)} routes -> {d}/fleet; covered {cv}, sum J_i {sj:.4f}, '
        f'J {sj + 298 - cv + 2:.4f}; depths {[len(x["st"]["asts"]) for x in fl.routes.values()]}')
    return dict(craft=len(fl.routes), covered=cv, sumJi=round(sj, 4), J=round(sj + 298 - cv + 2, 4), status=r['status'],
                picks=picks, dup_targets=dup, strip=srep, fleet=str(d / 'fleet'))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('pats', nargs='+')
    ap.add_argument('--N', default='8,9'); ap.add_argument('--exact', action='store_true')
    ap.add_argument('--strip', action='store_true'); ap.add_argument('--tlim', type=float, default=300.0)
    ap.add_argument('--gap', type=float, default=1e-3); ap.add_argument('--nproc', type=int, default=4)
    ap.add_argument('--max-tank', type=float, default=1400.0); ap.add_argument('--exclude', default='*/relocated/*')
    ap.add_argument('--budgets', default='', help='also the capacity frontier: max coverage at sum J_i <= B, e.g. 10.55,10.6,10.75')
    a = ap.parse_args()
    if a.nproc > 4:
        raise SystemExit('CPU budget: --nproc must be <= 4')
    out = pathlib.Path(a.out) if os.path.isabs(a.out) else pathlib.Path(_CWD0) / a.out
    out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt'); tic = time.time()
    excl = [e for e in a.exclude.split(',') if e]
    files = sorted(set(f for p in a.pats for f in glob.glob(p, recursive=True)
                       if not any(fnmatch.fnmatch(f, e) for e in excl)))
    say(f's13_master: {len(files)} route files from {a.pats} (excluded {excl}); N {a.N}; budgets {a.budgets or "-"}')
    RI.eph()
    pool = mp.get_context('fork').Pool(a.nproc) if a.nproc > 1 else None
    pmap = (lambda f, js: pool.map(f, js, chunksize=4)) if pool else (lambda f, js: [f(j) for j in js])
    try:
        res = pmap(load, [(f, a.max_tank) for f in files])
        why = collections.Counter(r for _, r in res)
        cols = [c for c, r in res if c is not None]
        best = {}
        for x in cols:                               # dedupe: same target set -> lightest tank
            k = frozenset(x['asts'])
            if k not in best or x['tank'] < best[k]['tank']:
                best[k] = x
        cols = sorted(best.values(), key=lambda x: (-len(x['asts']), x['tank']))
        union = set(t for x in cols for t in x['asts'])
        say(f'{len(cols)} distinct columns (load: {dict(why)}); depth {len(cols[0]["asts"])}..{len(cols[-1]["asts"])}; '
            f'union covers {len(union)} of {len(ALL)}')
        report = dict(files=len(files), load=dict(why), columns=len(cols), union=len(union), N={}, frontier={})
        for N in [int(x) for x in a.N.split(',')]:
            say(f'--- N <= {N}')
            lp = lp_value(cols, N); say(f'  LP value (sum y <= {N}): J >= {lp:.4f}')
            C = [dict(x) for x in cols]; verified = {}; dropped = []
            r = exact_loop(C, lambda C: solve(C, N, say, a.tlim, a.gap), pmap, verified, dropped, say,
                           rounds=6 if a.exact else 0) if a.exact else solve(C, N, say, a.tlim, a.gap)
            if r is None:
                report['N'][N] = dict(lp=lp, failed=True); continue
            rec = emit(out / f'N{N}', r, C, a, pmap, say, f's13_master N<={N}')
            rec.update(lp=round(lp, 4), dropped=dropped); report['N'][N] = rec
            for B in [float(x) for x in a.budgets.split(',') if x]:
                rb = exact_loop(C, lambda C: solve_budget(C, N, B, say, a.tlim, a.gap), pmap, verified, dropped, say,
                                rounds=6) if a.exact else solve_budget(C, N, B, say, a.tlim, a.gap)
                if rb is None:
                    continue
                rec = emit(out / f'N{N}_B{B:g}', rb, C, a, pmap, say, f's13_master N<={N} sum J_i<={B:g}')
                report['frontier'].setdefault(str(N), {})[f'{B:g}'] = rec
    finally:
        if pool is not None:
            pool.close(); pool.join()
    report['wall_s'] = round(time.time() - tic, 1)
    json.dump(report, open(out / 'master.json', 'w'), indent=1)
    say(f'master done ({report["wall_s"]} s) -> {out}/master.json')


if __name__ == '__main__':
    main()
