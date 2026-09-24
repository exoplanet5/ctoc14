"""Settle saved joint-planner tours (<stem>_tours.json) with the stock settle_tour and with the wait-aware seed.
Usage: settle_joint_tours.py stem [stem ...]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, signal, pathlib, warnings
HERE = pathlib.Path(__file__).resolve().parent; ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools')); sys.path.insert(0, str(HERE))
import numpy as np
import run_ialns as RI
from greedy_cover import _timeout
from ctoc14.impulsive import settle_tour, from_tour
from twin_wait import settle_tour_wait, from_tour_wait
from ctoc14.constants import DAY
warnings.filterwarnings('ignore')
out = []
for stem in sys.argv[1:]:
    tours = json.load(open(HERE / f'{stem}_tours.json'))
    for i, tour in enumerate(tours):
        n = len(tour['legs'])
        waits = sum(1 for k in range(1, n) if tour['legs'][k]['t_flyby'] - tour['legs'][k - 1]['t_flyby'] > tour['legs'][k]['tof'] + DAY)
        for name, fn, seed in (('stock', settle_tour, from_tour), ('wait-aware', settle_tour_wait, from_tour_wait)):
            ip0 = seed(RI.eph(), tour); rec = dict(stem=stem, craft=i, n=n, waits=waits, seed=name,
                                                    seed_tank=round(float(ip0.tank()), 1) if ip0 is not None else None)
            t1 = time.time()
            try:
                signal.signal(signal.SIGALRM, _timeout); signal.setitimer(signal.ITIMER_REAL, 240.0)
                ip, miss, lag = fn(RI.eph(), tour, lambda ip: RI.settle(ip, 100))
            except Exception as e:
                ip, miss = None, float('nan'); rec['error'] = type(e).__name__
            finally:
                signal.setitimer(signal.ITIMER_REAL, 0)
            rec['t_s'] = round(time.time() - t1, 1)
            if ip is not None and miss <= 150:
                tank = float(ip.tank()); rec.update(ok=True, tank_twin=round(tank, 1), kg_per_fb=round((tank - 600) / n, 2),
                                                    dv_per_fb=round(ip.dv() / n, 4), planner_fuel=round(tour['fuel_est'], 1))
                np.savez(HERE / f'{stem}_c{i}_{name}.npz', **RI.ist(ip))
            else:
                rec.update(ok=False, miss_km=float(miss) if np.isfinite(miss) else None)
            print(rec, flush=True); out.append(rec)
json.dump(out, open(HERE / f'settle_joint_{"_".join(sys.argv[1:])[:60]}.json', 'w'), indent=1)
