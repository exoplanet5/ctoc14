import os, sys, glob, pickle, pathlib
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
SP = pathlib.Path(__file__).resolve().parent
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
os.chdir(ROOT)
import numpy as np
from route_cover import w_load, ALL
import run_ialns as RI
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix
corpus = pickle.load(open(SP / 'corpus.pkl', 'rb'))
new = [x for x in (w_load(f) for f in sorted(glob.glob(str(SP / 'passes/*/route_*.npz')))) if x]
print('corpus', len(corpus), 'new settled pass routes', len(new))
ti = {t: i for i, t in enumerate(ALL)}; nt = len(ALL)
def solve(cols, K, label):
    nr = len(cols); R = []; C = []
    for j, x in enumerate(cols):
        for t in x['asts']: R.append(ti[t]); C.append(j)
    R += list(range(nt)); C += list(range(nr, nr + nt))
    A = coo_matrix((np.ones(len(R)), (R, C)), shape=(nt, nr + nt)).tocsr()
    card = np.concatenate([np.ones(nr), np.zeros(nt)])[None, :]
    Ji = np.array([RI.cost(x['tank']) for x in cols])
    c = np.concatenate([Ji, np.ones(nt)])
    res = milp(c, constraints=[LinearConstraint(A, np.ones(nt), np.full(nt, np.inf)), LinearConstraint(card, 0, K)],
               integrality=np.ones(nr + nt), bounds=Bounds(0, 1), options=dict(time_limit=300, mip_rel_gap=1e-4))
    pick = [j for j in range(nr) if res.x[j] > 0.5]
    cov = set(t for j in pick for t in cols[j]['asts']); sJ = sum(Ji[j] for j in pick)
    src = ['pass' if 'passes' in cols[j]['f'] else 'corpus' for j in pick]
    print(f'{label} K<={K}: {len(pick)} routes, covered {len(cov)}, slots {sum(len(cols[j]["asts"]) for j in pick)}, sumJi {sJ:.3f}, J {sJ + nt - len(cov) + 2:.3f}; '
          f'depths {[len(cols[j]["asts"]) for j in pick]} tanks {[round(cols[j]["tank"]) for j in pick]} src {src}', flush=True)
for K in (8, 9):
    solve(new, K, 'passes-only')
    solve(corpus + new, K, 'corpus+passes')
