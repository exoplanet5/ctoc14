"""Max coverage of a target subset S by K pool routes (honest pool + s16 best routes), HiGHS MILP, time-limited.
S = tail set (r5-r9 targets) optionally minus the 17 cheap deep-absorbable ones (q_tailabsorb.txt)."""
import json, sys, csv, numpy as np, time
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix, csr_matrix
T = {int(r['ast']): r for r in csv.DictReader(open('/Users/mickey/solarsystem/ctoc14/results/s17/redteam/targets.csv'))}
tail = sorted(k for k, v in T.items() if v['cls'] == 'tail')
absorb = [4, 300, 265, 142, 176, 14, 182, 130, 39, 184, 213, 18, 230, 56, 242, 26, 6]
K = int(sys.argv[1]); mode = sys.argv[2]; tl = float(sys.argv[3]) if len(sys.argv) > 3 else 120
S = [x for x in tail if (mode == 'all' or x not in absorb)]
Sidx = {x: i for i, x in enumerate(S)}
L = [json.loads(l) for l in open('/Users/mickey/solarsystem/ctoc14/results/s16/honest/index.jsonl')]
cols = []
for r in L:
    if r['status'] == 'rejected': continue
    z = np.load(f'/Users/mickey/solarsystem/ctoc14/results/s16/honest/route_{r["key"]}.npz')
    a = sorted(set(int(x) for x in z['asts']) & set(S))
    if len(a) >= 12: cols.append((r['key'], a, float(r['tank_out']), len(z['asts'])))
for f in ['r5', 'r6', 'r7', 'r8', 'r9']:
    z = np.load(f'/Users/mickey/solarsystem/ctoc14/results/s16/best/fleet/route_{f}.npz')
    cols.append((f, sorted(set(int(x) for x in z['asts']) & set(S)), 0.0, len(z['asts'])))
# dedupe by target set, keep lightest
d = {}
for c in cols:
    k = tuple(c[1])
    if k not in d or c[2] < d[k][2]: d[k] = c
cols = list(d.values()); nC = len(cols); nS = len(S)
print(f'|S| {nS}, columns {nC}, K {K}', flush=True)
# vars: x_j (nC), y_t (nS); max sum y  - 1e-5 * tank
c = np.concatenate([np.array([1e-6 * cc[2] for cc in cols]), -np.ones(nS)])
A = lil_matrix((nS + 1, nC + nS))
for j, cc in enumerate(cols):
    for x in cc[1]: A[Sidx[x], j] = -1.0
for t in range(nS): A[t, nC + t] = 1.0
A[nS, :nC] = 1.0
lb = np.concatenate([np.full(nS, -np.inf), [0]]); ub = np.concatenate([np.zeros(nS), [K]])
tic = time.time()
res = milp(c, constraints=LinearConstraint(csr_matrix(A), lb, ub), integrality=np.concatenate([np.ones(nC), np.zeros(nS)]),
           bounds=Bounds(0, 1), options=dict(time_limit=tl, disp=False, mip_rel_gap=0))
x = res.x[:nC] > 0.5; cov = int(round(res.x[nC:].sum()))
pk = [cols[j] for j in np.where(x)[0]]
print(f'status {res.status} {res.message}; covered {cov}/{nS} with {len(pk)} routes; bound {-res.mip_dual_bound if hasattr(res, "mip_dual_bound") and res.mip_dual_bound is not None else "?"}; {time.time() - tic:.0f} s')
for p in pk: print('  ', p[0], 'n', p[3], 'in S', len(p[1]), 'tank', round(p[2]))
