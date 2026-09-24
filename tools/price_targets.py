"""LP-relaxation pricing of the set cover over a route pool: dual value of each target's coverage constraint = how much J
the pool would save per unit of extra coverage of that target (0 = cheaply covered, ~1 = effectively missed).
Writes {id: weight} json to be used as --rare-weights in run_pool_campaign.py.
Usage: price_targets.py out.json pool.jsonl [pool2.jsonl ...] [--max-craft 12] [--fuel-margin 0.3] [--reserve 20]"""
import sys, json, os, glob, pathlib, argparse, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from scipy.optimize import linprog
from scipy.sparse import lil_matrix, csr_matrix
from ctoc14.constants import cost_sc, M_DRY, M0_MAX
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('src', nargs='+'); ap.add_argument('--max-craft', type=int, default=12)
ap.add_argument('--fuel-margin', type=float, default=0.3); ap.add_argument('--reserve', type=float, default=20.0); ap.add_argument('--min-targets', type=int, default=3)
a = ap.parse_args()
best = {}
def cost(t): return float(cost_sc(min(M0_MAX, M_DRY + t['fuel_est'] * (1 + a.fuel_margin) + a.reserve)))
for s in a.src:
    files = sorted(glob.glob(f'{s}/tour_*.json')) if os.path.isdir(s) else [s]
    for f in files:
        lines = open(f) if f.endswith('.jsonl') else [open(f).read()]
        for line in lines:
            t = json.loads(line); S = frozenset(l['ast'] for l in t['legs'])
            if len(S) < a.min_targets: continue
            c = cost(t)
            if S not in best or c < best[S]: best[S] = c
items = list(best.items()); n = len(items); print(f'{n} distinct route sets')
targets = [t for t in range(1, 301) if t not in (131, 144)]; tidx = {t: i for i, t in enumerate(targets)}; nT = len(targets)
# min sum c_r y_r - sum z_a  s.t. z_a - sum_{r∋a} y_r <= 0, sum y_r <= max_craft, 0<=y,z<=1
c = np.concatenate([np.array([c for _, c in items]), -np.ones(nT)])
A = lil_matrix((nT + 1, n + nT))
for j, (S, _) in enumerate(items):
    for t in S:
        if t in tidx: A[tidx[t], j] = -1.0
for i in range(nT): A[i, n + i] = 1.0
A[nT, :n] = 1.0
b = np.concatenate([np.zeros(nT), [a.max_craft]])
res = linprog(c, A_ub=csr_matrix(A), b_ub=b, bounds=[(0, 1)] * (n + nT), method='highs')
duals = -res.ineqlin.marginals[:nT]          # >= 0: value of relaxing coverage of target a
z = res.x[n:]; y = res.x[:n]
print(f'LP value {res.fun + nT + 2:.3f} (J incl. 131/144), fractional craft {y.sum():.2f}; targets with z<0.99: {(z < 0.99).sum()}')
w = {int(t): float(duals[i]) for i, t in enumerate(targets)}
order = sorted(w, key=lambda t: -w[t]); print('highest-priced targets:', [(t, round(w[t], 2)) for t in order[:30]])
json.dump(w, open(a.out, 'w'), indent=0)
