"""Set-packing selection of disjoint tours from a library (planned tours: tour JSONs).
min sum_j J_j y_j + N_miss  s.t. each asteroid in at most one selected tour.  Usage: select_tours_milp.py glob [glob ...]"""
import sys, json, glob, pathlib, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix
from ctoc14.constants import cost_sc
files = sorted({f for g in sys.argv[1:] for f in glob.glob(g)})
tours = []
for f in files:
    try: t = json.load(open(f))
    except Exception: continue
    if isinstance(t, dict) and 'legs' in t and t['legs']:
        tours.append((f, t['m0'], sorted({l['ast'] for l in t['legs']})))
n = len(tours); print(f'{n} tours in library')
c = np.array([cost_sc(m0) - len(S) for _, m0, S in tours])          # J_j - |S_j|  (misses are 300 - sum|S_j|)
A = lil_matrix((300, n))
for j, (_, _, S) in enumerate(tours):
    for a in S: A[a - 1, j] = 1
res = milp(c, constraints=LinearConstraint(A.tocsr(), 0, 1), integrality=np.ones(n), bounds=Bounds(0, 1))
y = np.round(res.x).astype(int); sel = [j for j in range(n) if y[j]]
cov = set().union(*[set(tours[j][2]) for j in sel]) if sel else set()
J = sum(cost_sc(tours[j][1]) for j in sel) + 300 - len(cov)
print(f'selected {len(sel)} tours, covered {len(cov)}, planned J = {J:.3f}')
for j in sel: print(f'  {tours[j][0]}  m0={tours[j][1]:.0f} n={len(tours[j][2])}')
