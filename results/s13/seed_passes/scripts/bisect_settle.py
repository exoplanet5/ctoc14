"""Measurement only: find the first prefix length at which a planner tour stops settling."""
import os, sys, json, pathlib, signal, time
for _v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'): os.environ.setdefault(_v,'1')
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/'tools'))
import numpy as np, run_ialns as RI
from ctoc14.impulsive import settle_tour, from_tour
H = np.asarray(json.load(open(ROOT/'results/n8/prizes_H.json')), float)
rows = [json.loads(l) for l in open(sys.argv[1])]
r = max((x for x in rows if len(x['targets']) == 36), key=lambda x: -x['fuel_est'])
legs = r['legs']
print('tour', len(legs), 'launch', r['t_launch']/86400, 'legs:', [(l['ast'], round(l['tof']/86400), round(l['dv_est'],2), 'H' if H[l['ast']-1]>=1.5 else '') for l in legs])
def ok(k):
    t = dict(r); t['legs'] = legs[:k]; t0 = time.time()
    ip0 = from_tour(RI.eph(), t)
    tank0 = float(ip0.tank()) if ip0 is not None else float('nan')
    ip, miss, lag = settle_tour(RI.eph(), t, lambda ip: RI.settle(ip, 100))
    print(f'  prefix {k}: seed tank {tank0:.0f}, settled {ip is not None}, miss {miss:.0f} km, {time.time()-t0:.0f} s', flush=True)
    return ip is not None
lo, hi = 1, len(legs)
while hi - lo > 1:
    mid = (lo + hi) // 2
    if ok(mid): lo = mid
    else: hi = mid
print('last settling prefix', lo, 'first failing', hi, 'failing leg', legs[hi-1])
