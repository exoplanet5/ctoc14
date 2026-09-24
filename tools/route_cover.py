"""Set cover over EVERY impulsive route on disk (stage 7, section 9).

The six fleets built so far leave ~50-70 targets each but only 3 targets in common: the leftovers are not a hard
family, they are whatever each construction happened to drop.  Routes are independent spacecraft (J = sum_i J_i +
N_miss, no coupling), so any subset of routes ever built is a legal fleet and the fleet problem is a plain set cover:

    min  sum_r cost(tank_r) y_r + sum_t z_t      s.t.  sum_{r: t in r} y_r + z_t >= 1,  y, z binary

N is an OUTPUT, not a constraint.  --exact re-settles the chosen routes (tanks in old dirs were settled with older
code) and re-solves once with the corrected costs.

Usage: route_cover.py 'results/**/route_*.npz' [--out results/s7/cover] [--max-tank 1400]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, glob, time, pathlib, argparse, collections, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from ctoc14.search import UNREACHABLE
ALL = [t for t in range(1, 301) if t not in UNREACHABLE]


def w_load(f):
    try:
        z = np.load(f); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL'])
        if 'asts' not in st or len(st['asts']) == 0: return None
        ip = RI.ipr(st); tank = float(ip.tank())
        if not np.isfinite(tank) or tank < 600 or tank > 2000: return None
        return dict(f=str(f), tank=tank, asts=sorted(set(int(x) for x in st['asts']) - set(UNREACHABLE)))
    except Exception:
        return None


def w_settle(f):
    z = np.load(f); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL'])
    ip = RI.ipr(st)
    try:
        miss = RI.settle(ip, 100)
    except Exception:
        return str(f), None
    return str(f), (RI.ist(ip), float(ip.tank()), miss)


def lp_duals(cols, say):
    """LP relaxation of the set cover; duals pi_t = what covering target t is worth to the fleet [J units]."""
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    nr = len(cols); ti = {t: i for i, t in enumerate(ALL)}; nt = len(ALL)
    c = np.concatenate([[RI.cost(x['tank']) for x in cols], np.ones(nt)])
    R = []; C = []
    for j, x in enumerate(cols):
        for t in x['asts']: R.append(ti[t]); C.append(j)
    R += list(range(nt)); C += list(range(nr, nr + nt))
    A = coo_matrix((-np.ones(len(R)), (R, C)), shape=(nt, nr + nt)).tocsr()
    res = linprog(c, A_ub=A, b_ub=-np.ones(nt), bounds=(0, 1), method='highs')
    pi = {t: float(max(0.0, -res.ineqlin.marginals[ti[t]])) for t in ALL}
    say(f'LP bound {res.fun + 2:.3f}; duals: median {np.median(list(pi.values())):.4f}, '
        f'max {max(pi.values()):.3f}, >= 0.05: {sum(1 for v in pi.values() if v >= 0.05)}')
    return pi, float(res.fun) + 2


def solve(cols, say, tlim=300.0, gap=0.001):
    from scipy.optimize import milp, LinearConstraint, Bounds
    from scipy.sparse import coo_matrix
    nr = len(cols); ti = {t: i for i, t in enumerate(ALL)}; nt = len(ALL)
    c = np.concatenate([[RI.cost(x['tank']) for x in cols], np.ones(nt)])
    R = []; C = []
    for j, x in enumerate(cols):
        for t in x['asts']:
            R.append(ti[t]); C.append(j)
    R += list(range(nt)); C += list(range(nr, nr + nt))
    A = coo_matrix((np.ones(len(R)), (R, C)), shape=(nt, nr + nt)).tocsr()
    res = milp(c, constraints=LinearConstraint(A, np.ones(nt), np.full(nt, np.inf)),
               integrality=np.ones(nr + nt), bounds=Bounds(0, 1),
               options=dict(time_limit=tlim, mip_rel_gap=gap))
    if res.x is None: say(f'MILP failed: {res.message}'); return None
    pick = [j for j in range(nr) if res.x[j] > 0.5]
    miss = int(round(res.x[nr:].sum()))
    say(f'MILP {res.message[:40]}: {len(pick)} craft, {nt - miss} covered, '
        f'sum J_i {sum(RI.cost(cols[j]["tank"]) for j in pick):.3f}, J {res.fun + 2:.3f} (incl. 2 permanent misses)')
    return pick, miss, float(res.fun) + 2


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('pats', nargs='+'); ap.add_argument('--out', default='results/s7/cover')
    ap.add_argument('--max-tank', type=float, default=1400.0); ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--tlim', type=float, default=300.0); ap.add_argument('--exact', action='store_true')
    ap.add_argument('--duals', action='store_true', help='also write the LP duals for column generation')
    a = ap.parse_args()
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt')
    files = sorted(set(f for p in a.pats for f in glob.glob(p, recursive=True)))
    say(f'{len(files)} route files')
    with mp.get_context('fork').Pool(a.nproc) as pl:
        cols = [x for x in pl.map(w_load, files, chunksize=8) if x and x['tank'] <= a.max_tank and x['asts']]
        best = {}
        for x in cols:                                     # dedupe: same target set -> lightest tank
            k = frozenset(x['asts'])
            if k not in best or x['tank'] < best[k]['tank']: best[k] = x
        cols = sorted(best.values(), key=lambda x: -len(x['asts']))
        say(f'{len(cols)} distinct routes; depth {len(cols[0]["asts"])}..{len(cols[-1]["asts"])}; '
            f'union covers {len(set(t for x in cols for t in x["asts"]))} of {len(ALL)}')
        if a.duals:
            pi, lb = lp_duals(cols, say)
            json.dump({str(k): v for k, v in pi.items()}, open(out / 'duals.json', 'w'))
            say(f'-> {out}/duals.json')
        r = solve(cols, say, a.tlim)
        if r is None: return
        pick, miss, J = r
        if a.exact:
            say('re-settling the chosen routes with current code...')
            res = dict(pl.map(w_settle, [cols[j]['f'] for j in pick]))
            fixed = []
            for j in pick:
                v = res[cols[j]['f']]
                if v is None or v[2] > 150: say(f'  {cols[j]["f"]}: does NOT settle -> dropped'); continue
                st, tank, ms = v
                fixed.append(dict(f=cols[j]['f'], tank=tank, asts=sorted(set(int(x) for x in st['asts']) - set(UNREACHABLE)), st=st))
                if abs(tank - cols[j]['tank']) > 1: say(f'  {cols[j]["f"]}: tank {cols[j]["tank"]:.0f} -> {tank:.0f} kg')
            keep = [x for x in cols if x['f'] not in set(f['f'] for f in fixed)] + [{k: v for k, v in f.items() if k != 'st'} for f in fixed]
            say('re-solving with corrected tanks')
            r2 = solve(keep, say, a.tlim)
            if r2: pick, miss, J = r2[0], r2[1], r2[2]; cols = keep
    fl = RI.IFleet()
    for k, j in enumerate(pick):
        z = np.load(cols[j]['f']); st = {q: z[q] for q in z.files}; st['tL'] = float(st['tL'])
        fl.routes[f'{k:02d}'] = dict(st=st, tank=cols[j]['tank'])
    fl.save(out / 'fleet', note='set cover over all routes')
    json.dump([cols[j]['f'] for j in pick], open(out / 'picked.json', 'w'), indent=1)
    say(f'saved {len(pick)} routes -> {out}/fleet   (covered {len(fl.coverage())}, J {fl.J():.3f})')


if __name__ == '__main__':
    main()
