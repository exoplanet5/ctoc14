import os
os.environ['OMP_NUM_THREADS']='1'
import sys, time, numpy as np
sys.path.insert(0,'.'); sys.path.insert(0,'tools')
from run_ialns import IFleet, ipr, w_cands, settle, cost
from ctoc14.globalopt import insert_homotopy
from ctoc14.constants import DAY
fleet = IFleet('results/newgen/ifleet3'); X = int(sys.argv[1]); hosts = sys.argv[2].split(',')
for h in hosts:
    cands = sorted(w_cands((h, fleet.routes[h]['st'], [X], 0.6)), key=lambda c: c['dist'])[:3]
    for c in cands:
        ip = ipr(fleet.routes[h]['st']); J0 = cost(ip.tank()); tic = time.time()
        msgs = []
        d = insert_homotopy(ip, X, c['t'], stages=10, iters_stage=6, tol_km=100, verbose=True, log=msgs.append)
        hom = [m for m in msgs if 'homotopy' in m]
        print(f'{h} t={c["t"]/DAY:.0f} d dist {c["dist"]:.3f} AU: homotopy end miss {d.max():.0f} km; stages: ' + ' | '.join(m.split(':')[1].strip()[:28] for m in hom[::max(1,len(hom)//5)]), flush=True)
        if d.max() < 1e4:
            miss = settle(ip, 100); print(f'   settle miss {miss:.0f}, dJ {cost(ip.tank())-J0:+.4f} ({time.time()-tic:.0f} s)', flush=True)
