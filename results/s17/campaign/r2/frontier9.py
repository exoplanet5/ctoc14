"""R2 CAMPAIGN quick check: 9-craft (and 8-craft) MIN-COST coverage frontier over parents vs parents+children.
Question: do stitched children cut the 9-craft cost enough to pay (-30 kg by 09-24 12:00, -77 kg by 09-25 12:00)?
min sum cost_j  s.t. sum_j x_j <= N, y_t <= sum_{j∋t} x_j, sum y_t >= K.  Pre-registered in PREREG.txt.
usage: frontier9.py parents|children   -> r2/frontier_<tag>.json"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib
sys.path.insert(0, '/Users/mickey/solarsystem/ctoc14/results/s17/campaign')
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix, csr_matrix, hstack, identity
import xover_scan as X
tag = sys.argv[1]; OUTF = X.OUT / 'r2' / f'frontier_{tag}.json'
cols = X.load_cols(tag == 'children')
nr = len(cols); print(tag, nr, 'cols', sum(c['kind'] == 'child' for c in cols), 'children', flush=True)
R = []; C = []
for j, c in enumerate(cols):
    for t in c['set']:
        if t in X.TI:
            R.append(X.TI[t]); C.append(j)
A = coo_matrix((np.ones(len(R)), (R, C)), shape=(X.NT, nr)).tocsr()
cost = np.array([c['cost'] for c in cols])
res = []
for N, K, tl in ((9, 298, 420), (9, 296, 300), (9, 294, 300), (8, 286, 300)):
    M = hstack([A, -identity(X.NT, format='csr')]).tocsr()
    cons = [LinearConstraint(M, 0, np.inf),
            LinearConstraint(csr_matrix(np.concatenate([np.ones(nr), np.zeros(X.NT)])[None, :]), -np.inf, N),
            LinearConstraint(csr_matrix(np.concatenate([np.zeros(nr), np.ones(X.NT)])[None, :]), K, np.inf)]
    obj = np.concatenate([cost, np.zeros(X.NT)])
    t0 = time.time()
    r = milp(obj, constraints=cons, integrality=np.concatenate([np.ones(nr), np.zeros(X.NT)]), bounds=Bounds(0, 1),
             options=dict(time_limit=tl, mip_rel_gap=1e-5))
    out = dict(N=N, K=K, msg=str(r.message)[:70], sec=round(time.time() - t0),
               bound=None if getattr(r, 'mip_dual_bound', None) is None else float(r.mip_dual_bound))
    if r.x is not None:
        pick = [j for j in range(nr) if r.x[j] > 0.5]
        cov = set().union(*[cols[j]['set'] for j in pick]) & set(X.ALL)
        out.update(sumJi=float(cost[pick].sum()), covered=len(cov), misses=sorted(set(X.ALL) - cov),
                   picks=[dict(kind=cols[j]['kind'], n=len(cols[j]['set']), tank=round(cols[j]['tank'], 1),
                               **({k: cols[j][k] for k in ('a', 'b', 'k', 'L', 'dvj')} if cols[j]['kind'] == 'child'
                                  else dict(f=cols[j]['f']))) for j in pick])
    print(json.dumps({k: v for k, v in out.items() if k != 'picks'}), flush=True)
    res.append(out); json.dump(res, open(OUTF, 'w'), indent=1)
print('done', flush=True)
