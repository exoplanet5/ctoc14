"""Column-library statistics for the stage 6 gates: depth distribution per pool fraction, the reach matrix
(which skeleton can reach which easy target, and in how many columns) and the scarce easy targets (reachable by
<= --scarce skeletons).  Also reports the LP relaxation bound of the one-column-per-skeleton set-cover, which is
the honest ceiling of the column library before integrality.
Usage: col_stats.py results/n8s6/B/cols_*.jsonl [--scarce 2] [--skel-dir results/n8/skel8b] [--lp]
"""
import sys, json, glob, argparse, pathlib, collections
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
from ctoc14.search import UNREACHABLE

ap = argparse.ArgumentParser(); ap.add_argument('cols', nargs='+'); ap.add_argument('--scarce', type=int, default=2)
ap.add_argument('--skel-dir', default='results/n8/skel8b'); ap.add_argument('--lp', action='store_true')
ap.add_argument('--tank-cal', type=float, default=1.05); ap.add_argument('--json', default='')
a = ap.parse_args()
ALL = [t for t in range(1, 301) if t not in UNREACHABLE]
cost = lambda tank: 1.0 + (tank - 600.0) / 1400.0 + ((tank - 600.0) / 1400.0) ** 2

cols = []
for pat in a.cols:
    for f in sorted(glob.glob(pat)):
        cols += [json.loads(l) for l in open(f)]
if not cols:
    print('no columns'); sys.exit(1)
import run_ialns as RI
skel = RI.IFleet(a.skel_dir)
hard = {}
for n, r in skel.routes.items():
    for t in r['st']['asts']: hard[int(t)] = n
easy = [t for t in ALL if t not in hard]
skels = sorted(set(c['skel'] for c in cols))
print(f'{len(cols)} columns over {len(skels)} skeletons; {len(easy)} easy targets, {len(hard)} hard')

# depth (easy targets per column) and tank per skeleton
print('\nper skeleton: columns, easy-depth min/med/max, tank med, distinct easy reached')
for s in skels:
    cs = [c for c in cols if c['skel'] == s]
    d = np.array([len([t for t in c['targets'] if t not in hard]) for c in cs])
    tk = np.array([c['tank'] * a.tank_cal for c in cs])
    reach = set(t for c in cs for t in c['targets'] if t not in hard)
    print(f'  {s}: {len(cs):5d} cols, easy {d.min():2d}/{int(np.median(d)):2d}/{d.max():2d}, '
          f'tank {np.median(tk):4.0f}, reaches {len(reach):3d} easy')

# reach matrix
R = {t: set() for t in easy}; NC = collections.Counter()
for c in cols:
    for t in c['targets']:
        if t not in hard: R[t].add(c['skel']); NC[t] += 1
nre = [t for t in easy if not R[t]]
scarce = sorted(t for t in easy if 0 < len(R[t]) <= a.scarce)
print(f'\nunreachable by any column ({len(nre)}): {nre}')
print(f'scarce (<= {a.scarce} skeletons) ({len(scarce)}): ' + ' '.join(f'{t}[{"".join(sorted(R[t]))}]' for t in scarce))
h = collections.Counter(len(R[t]) for t in easy)
print('skeletons able to reach a target: ' + ' '.join(f'{k}:{h[k]}' for k in sorted(h)))

# capacity check: best single column per skeleton, and the sum of the 8 deepest (upper bound on coverage)
best = {s: max((len([t for t in c['targets'] if t not in hard]) for c in cols if c['skel'] == s), default=0) for s in skels}
print(f'\ndeepest easy column per skeleton: ' + ' '.join(f'{s}:{best[s]}' for s in skels) +
      f'  sum {sum(best.values())} vs {len(easy)} needed')

if a.lp:
    from scipy.optimize import linprog
    from scipy.sparse import lil_matrix
    nc, nt = len(cols), len(easy); tid = {t: i for i, t in enumerate(easy)}
    c_obj = np.concatenate([[cost(c['tank'] * a.tank_cal) for c in cols], np.full(nt, -1.0)])
    A = lil_matrix((len(skels) + nt, nc + nt)); lb = np.zeros(len(skels) + nt); ub = np.zeros(len(skels) + nt)
    for i, s in enumerate(skels):
        for j, c in enumerate(cols):
            if c['skel'] == s: A[i, j] = 1.0
        lb[i] = ub[i] = 1.0
    for t in easy:
        r = len(skels) + tid[t]; A[r, nc + tid[t]] = 1.0
        for j, c in enumerate(cols):
            if t in c['targets']: A[r, j] = -1.0
        lb[r] = -np.inf; ub[r] = 0.0
    from scipy.optimize import LinearConstraint
    res = linprog(c_obj, A_ub=None, bounds=[(0, 1)] * (nc + nt),
                  A_eq=None, method='highs', constraints=None) if False else None
    import scipy.optimize as so
    res = so.linprog(c_obj, A_ub=A.tocsr()[len(skels):], b_ub=ub[len(skels):],
                     A_eq=A.tocsr()[:len(skels)], b_eq=ub[:len(skels)], bounds=[(0, 1)] * (nc + nt), method='highs')
    if res.status == 0:
        y = res.x[nc:]
        print(f'\nLP relaxation: coverage {y.sum():.1f}/{nt} easy (+{len(hard)} hard), objective {res.fun:.3f}')
if a.json:
    json.dump(dict(reach={str(t): sorted(R[t]) for t in easy}, scarce=scarce, unreachable=nre,
                   best_depth=best), open(a.json, 'w'), indent=1)
