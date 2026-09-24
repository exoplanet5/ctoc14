"""Leftovers of a pass fleet and their closest approach to every host (RI.w_cands, dmax 0.35 AU). Read-only."""
import os, sys, json, pathlib
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
os.chdir(ROOT)
import numpy as np, multiprocessing as mp
import run_ialns as RI
from greedy_cover import ALL
src = sys.argv[1]; out = sys.argv[2]
fl = RI.IFleet(src); cov = fl.coverage(); U = sorted(set(ALL) - set(cov))
with mp.get_context('fork').Pool(8) as pl:
    res = pl.map(RI.w_cands, [(n, r['st'], U, 0.35) for n, r in fl.routes.items()])
best = {u: [] for u in U}
for n, lst in zip(fl.routes, res):
    for c in lst:
        best[c['ast']].append((round(c['dist'], 4), n, round(c['t'] / 86400.0, 1)))
rows = {u: sorted(v)[:4] for u, v in best.items()}
print(fl.summary()); print(len(U), 'leftovers')
d1 = sorted((v[0][0] if v else 9.0) for v in rows.values())
for th in (0.03, 0.05, 0.08, 0.12, 0.2, 0.35):
    print(f'  nearest host within {th:.2f} AU: {sum(d <= th for d in d1)}')
from collections import Counter
print('nearest host counts:', Counter(v[0][1] for v in rows.values() if v))
json.dump(dict(src=src, leftovers=U, near={str(u): v for u, v in rows.items()}), open(out, 'w'), indent=1)
