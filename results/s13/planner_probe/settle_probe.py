"""Stage-13 settle probe: take stored deep planner tours (jsonl columns) and settle each into the impulsive twin with
impulsive.LIN_FALLBACK = 2.5 (current default) and = None (pre-09-20 behaviour).  Records ok / miss / twin tank /
wall time.  Usage: settle_probe.py out.json file.jsonl:idx[,idx] [file2.jsonl:idx ...]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, signal, pathlib, warnings
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from greedy_cover import _timeout
import ctoc14.impulsive as IM
from ctoc14.impulsive import settle_tour, from_tour
from ctoc14.lambert import lambert_best
from ctoc14.constants import DAY
warnings.filterwarnings('ignore')

out = sys.argv[1]; jobs = []
for spec in sys.argv[2:]:
    f, idx = spec.rsplit(':', 1)
    L = [json.loads(l) for l in open(ROOT / f)]
    for i in idx.split(','):
        r = L[int(i)]; jobs.append((f, int(i), r.get('tour', r)))
res = []
for f, i, tour in jobs:
    n = len(tour['legs']); m0 = tour['m0']; fuel = tour['fuel_est']
    # count legs whose Lambert junction exceeds the fallback threshold (the legs LIN_FALLBACK re-seeds)
    for lf in (2.5, None):
        IM.LIN_FALLBACK = lf
        rec = dict(file=f, idx=i, n=n, planner_fuel=round(fuel, 1), lin_fallback=lf)
        ip0 = from_tour(RI.eph(), tour)
        rec['seed_tank'] = round(float(ip0.tank()), 1) if ip0 is not None else None
        t1 = time.time()
        try:
            signal.signal(signal.SIGALRM, _timeout); signal.setitimer(signal.ITIMER_REAL, 300.0)
            ip, miss, lag = settle_tour(RI.eph(), tour, lambda ip: RI.settle(ip, 100))
        except Exception as e:
            ip, miss, lag = None, float('nan'), None; rec['error'] = type(e).__name__
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
        rec['t_s'] = round(time.time() - t1, 1)
        if ip is not None and miss <= 150:
            tank = float(ip.tank())
            rec.update(ok=True, miss_km=round(miss, 1), lag_d=lag / DAY, tank_twin=round(tank, 1),
                       kg_per_fb=round((tank - 600) / n, 2), dv_per_fb=round(ip.dv() / n, 4))
        else:
            rec.update(ok=False, miss_km=(float(miss) if np.isfinite(miss) else None))
        print(rec, flush=True); res.append(rec)
        json.dump(res, open(out, 'w'), indent=1)
