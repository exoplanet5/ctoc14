"""Pick one column per skeleton from tools/skel_fill.py --columns output: minimise sum J_i(tank) + 1 J per uncovered
target (the true objective), one column per skeleton (MIP, scipy milp).  Writes the chosen columns to out.json.
Usage: skel_select.py columns.jsonl out.json [--exclude skel:col:n,...] [--miss 1.0] [--time 120]"""
import sys, json, argparse, pathlib
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
from ctoc14.search import UNREACHABLE
M_DRY = 600.0
cost = lambda tank: 1.0 + (tank - M_DRY) / 1400.0 + ((tank - M_DRY) / 1400.0) ** 2

ap = argparse.ArgumentParser(); ap.add_argument('cols'); ap.add_argument('out')
ap.add_argument('--exclude', default=''); ap.add_argument('--miss', type=float, default=1.0)
ap.add_argument('--time', type=float, default=120.0); ap.add_argument('--tank-cal', type=float, default=1.0,
                help='multiply planner tank by this before costing (twin/planner ratio)')
a = ap.parse_args()
excl = set(a.exclude.split(',')) if a.exclude else set()
cols = [json.loads(l) for l in open(a.cols)]
cols = [c for c in cols if f'{c["skel"]}:{c["col"]}:{c["n"]}' not in excl]
skels = sorted(set(c['skel'] for c in cols)); ALL = [t for t in range(1, 301) if t not in UNREACHABLE]
tid = {t: i for i, t in enumerate(ALL)}
nc, nt = len(cols), len(ALL)
c_obj = np.concatenate([[cost(c['tank'] * a.tank_cal) for c in cols], np.full(nt, -a.miss)])   # min sum cost - miss*y
A = lil_matrix((len(skels) + nt, nc + nt))
lb = np.zeros(len(skels) + nt); ub = np.zeros(len(skels) + nt)
for i, s in enumerate(skels):
    for j, c in enumerate(cols):
        if c['skel'] == s: A[i, j] = 1.0
    lb[i] = ub[i] = 1.0
for t in ALL:
    r = len(skels) + tid[t]
    A[r, nc + tid[t]] = 1.0
    for j, c in enumerate(cols):
        if t in c['targets']: A[r, j] = -1.0
    lb[r] = -np.inf; ub[r] = 0.0                   # y_t - sum x_c(t) <= 0
res = milp(c_obj, constraints=LinearConstraint(A.tocsr(), lb, ub), integrality=np.ones(nc + nt),
           bounds=Bounds(0, 1), options=dict(time_limit=a.time, disp=False))
if res.x is None:
    print('no solution', res.message); sys.exit(1)
x = np.round(res.x[:nc]).astype(int); chosen = [c for j, c in enumerate(cols) if x[j]]
cov = set(t for c in chosen for t in c['targets']); left = [t for t in ALL if t not in cov]
sumJ = sum(cost(c['tank'] * a.tank_cal) for c in chosen)
print(f'{nc} columns, {len(skels)} skeletons -> chosen ' + ' '.join(f'{c["skel"]}:c{c["col"]}:{c["n"]}@{c["tank"]:.0f}' for c in chosen))
print(f'covered {len(cov)}/{nt}, leftovers {len(left)}: {left}')
print(f'sum J_i {sumJ:.3f} -> J {sumJ + 2 + len(left):.3f} with misses, {sumJ + 2:.3f} if all inserted free  (status {res.status}, {res.message})')
json.dump(dict(chosen=chosen, covered=sorted(cov), leftovers=left, sum_Ji=sumJ), open(a.out, 'w'))
