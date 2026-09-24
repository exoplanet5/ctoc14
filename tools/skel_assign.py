"""Assign every easy target to a skeleton from the MEASURED reach matrix of a column library (stage 6 C).

Stage 5's geometric pre-assignment (`skel_fill --assign`, closest approach of the asteroid to the skeleton
trajectory) failed because distance is a poor proxy for "this route can actually fly there in its schedule".
A column library answers the same question by measurement: w[s][t] = how many of skeleton s's columns contain
target t (0 = s never reaches t in any tour).  The assignment is a transportation problem

    max  sum_{s,t} w_hat[s][t] x[s][t]      s.t.  sum_s x[s][t] <= 1,  sum_t x[s][t] <= cap
    (w_hat = w / max_s w, so a target contributes its RELATIVE preference, and scarce targets -- reachable by
     one or two skeletons -- dominate because nobody else can claim them)

solved as a MIP; unassigned targets stay in the shared pad that `skel_fill --pools` offers every route at
--pad-prize.  Writes {route: [ids]} for `skel_fill --pools`.

Usage: skel_assign.py 'results/n8s6/B/cols_*.jsonl' pools.json [--cap 32] [--skel-dir results/n8/skel8b]
"""
import sys, json, glob, argparse, pathlib, collections
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
from ctoc14.search import UNREACHABLE

ap = argparse.ArgumentParser(); ap.add_argument('cols', nargs='+'); ap.add_argument('out')
ap.add_argument('--cap', type=int, default=32, help='max targets assigned to one skeleton (its measured capacity)')
ap.add_argument('--skel-dir', default='results/n8/skel8b'); ap.add_argument('--time', type=float, default=120.0)
ap.add_argument('--depth-weight', type=float, default=0.0, help='also weight a column by its depth')
a = ap.parse_args()

cols = []
for pat in a.cols:
    for f in sorted(glob.glob(pat)): cols += [json.loads(l) for l in open(f)]
if not cols: print('no columns'); sys.exit(1)
import run_ialns as RI
skel = RI.IFleet(a.skel_dir)
hard = {int(t): n for n, r in skel.routes.items() for t in r['st']['asts']}
ALL = [t for t in range(1, 301) if t not in UNREACHABLE]
easy = [t for t in ALL if t not in hard]
skels = sorted(set(c['skel'] for c in cols))
W = {s: collections.Counter() for s in skels}
for c in cols:
    w = 1.0 + a.depth_weight * len([t for t in c['targets'] if t not in hard])
    for t in c['targets']:
        if t not in hard: W[c['skel']][t] += w
tot = {s: sum(W[s].values()) or 1.0 for s in skels}
print(f'{len(cols)} columns, {len(skels)} skeletons, {len(easy)} easy targets')

pairs = [(s, t) for s in skels for t in easy if W[s][t] > 0]
idx = {p: i for i, p in enumerate(pairs)}
if not pairs: print('empty reach matrix'); sys.exit(1)
# relative preference: a target reachable only by one skeleton scores 1 for it, 0 elsewhere
wmax = {t: max(W[s][t] for s in skels) for t in easy}
obj = np.array([-(W[s][t] / wmax[t]) for s, t in pairs])
A = lil_matrix((len(easy) + len(skels), len(pairs)))
lb = np.zeros(len(easy) + len(skels)); ub = np.zeros(len(easy) + len(skels))
for i, t in enumerate(easy):
    for s in skels:
        if (s, t) in idx: A[i, idx[(s, t)]] = 1.0
    lb[i] = -np.inf; ub[i] = 1.0
for j, s in enumerate(skels):
    r = len(easy) + j
    for t in easy:
        if (s, t) in idx: A[r, idx[(s, t)]] = 1.0
    lb[r] = -np.inf; ub[r] = a.cap
res = milp(obj, constraints=LinearConstraint(A.tocsr(), lb, ub), integrality=np.ones(len(pairs)),
           bounds=Bounds(0, 1), options=dict(time_limit=a.time, disp=False))
if res.x is None: print('no solution', res.message); sys.exit(1)
x = np.round(res.x).astype(int)
pools = {s: [] for s in skels}
for (s, t), v in zip(pairs, x):
    if v: pools[s].append(t)
done = set(t for v in pools.values() for t in v)
left = [t for t in easy if t not in done]
nre = [t for t in easy if wmax.get(t, 0) == 0]
print('assigned: ' + ' '.join(f'{s}:{len(pools[s])}' for s in skels) + f'  total {len(done)}/{len(easy)}')
print(f'unassigned ({len(left)}): {left}')
print(f'reachable by no column ({len(nre)}): {nre}')
h = collections.Counter(sum(1 for s in skels if W[s][t] > 0) for t in easy)
print('skeletons reaching a target: ' + ' '.join(f'{k}:{h[k]}' for k in sorted(h)))
json.dump({s: sorted(v) for s, v in pools.items()}, open(a.out, 'w'), indent=1)
print(f'-> {a.out}')
