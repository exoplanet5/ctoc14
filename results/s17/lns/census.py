"""s17 LNS census: forced (lin_price) insertion price of every r9 (and r8) target into every other E-fleet route,
approach minima <= 0.25 AU (s15b_close._price_host = w_cands + lin_price, 3 epoch offsets).  Read-only on tools."""
import os
for _v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, multiprocessing as mp
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT/'tools')):
    sys.path.insert(0, _p)
os.chdir(ROOT)
import numpy as np
import run_ialns as RI
import s15b_close as CL
from ctoc14.constants import DAY
F = RI.IFleet('results/s16/best/fleet')
vict = sys.argv[1] if len(sys.argv) > 1 else 'r9'
T = sorted(int(a) for a in F.routes[vict]['st']['asts'])
jobs = [(k, r['st'], T, 0.25) for k, r in sorted(F.routes.items()) if k != vict]
tic = time.time()
with mp.get_context('fork').Pool(2) as pool:
    R = pool.map(CL._price_host, jobs, chunksize=1)
out = [c for lst in R for c in lst]
json.dump(dict(victim=vict, targets=T, cands=out, wall=time.time()-tic), open(f'results/s17/lns/census_{vict}.json','w'), indent=1)
print(f'{len(out)} cands in {time.time()-tic:.0f} s')
