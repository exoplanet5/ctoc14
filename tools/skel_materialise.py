"""Settle the columns chosen by skel_select.py into an IFleet (impulsive twins).  A column that fails to settle is
reported (re-run skel_select.py with --exclude skel:col:n).
Usage: skel_materialise.py chosen.json out_dir [--iters 100] [--timeout 240]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, signal, pathlib, argparse, warnings
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from greedy_cover import _timeout
from ctoc14.impulsive import settle_tour
warnings.filterwarnings('ignore')
ap = argparse.ArgumentParser(); ap.add_argument('chosen'); ap.add_argument('out')
ap.add_argument('--iters', type=int, default=100); ap.add_argument('--timeout', type=float, default=240.0)
a = ap.parse_args()
out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt'); RI.eph()
ch = json.load(open(a.chosen)); fleet = RI.IFleet(); failed = []
for c in ch['chosen']:
    tic = time.time()
    try:
        signal.signal(signal.SIGALRM, _timeout); signal.setitimer(signal.ITIMER_REAL, a.timeout)
        ip, miss, lag = settle_tour(RI.eph(), c['tour'], lambda ip: RI.settle(ip, a.iters))
    except Exception as e:
        ip = None; miss = float('inf'); say(f'{c["skel"]}:c{c["col"]}:{c["n"]} raised {type(e).__name__}')
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    if ip is None or miss > 150.0:
        say(f'{c["skel"]}:c{c["col"]}:{c["n"]} FAILED (miss {miss:.0f}, {time.time()-tic:.0f} s)'); failed.append(f'{c["skel"]}:{c["col"]}:{c["n"]}'); continue
    tank = float(ip.tank()); fleet.routes[c['skel']] = dict(st=RI.ist(ip), tank=tank)
    say(f'{c["skel"]}:c{c["col"]}:{c["n"]} planner {c["tank"]:.0f} -> twin {tank:.0f} kg, J_i {RI.cost(tank):.3f} ({time.time()-tic:.0f} s)')
fleet.save(out, note='materialised'); say(f'fleet: {fleet.summary()}')
if failed:
    say('FAILED columns (exclude and re-select): ' + ','.join(failed))
json.dump(dict(failed=failed), open(out / 'materialise.json', 'w'))
