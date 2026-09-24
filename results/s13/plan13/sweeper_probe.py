"""Measurement: how many of a pass's leftovers ONE extra (9th) route flies (strict pool = leftovers only).
usage: sweeper_probe.py LEFTOVER_JSON OUT_JSON DVMAX DRMAX NPT MC"""
import os, sys, json, time, signal, pathlib
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
os.chdir(ROOT)
import numpy as np
import run_ialns as RI
import ctoc14.impulsive as IM
from greedy_cover import BP, CM, DeepCollect, ALL, _timeout
from ctoc14.search import beam_search, tour_json
from ctoc14.constants import DAY, AU
IM.LIN_FALLBACK = 2.5
L = json.load(open(sys.argv[1]))['leftovers']; out = sys.argv[2]
dv, dr, npt, mc = float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
prize = np.zeros(300)
for t in L: prize[t - 1] = 1.0
P = BP(beam=100, w_fuel=CM.w_fuel(480, 1600), m_margin=40.0, dv_max=dv, tofs=np.arange(15, 401, 5) * DAY,
       lin_tofs=np.arange(20, 601, 10) * DAY, lin_drmax=dr * AU, vinf_cap=4.0, tof_refine=True, max_depth=60,
       n_per_target=npt, max_children=mc)
P.prize = prize; P.collect = DeepCollect(4)
tic = time.time()
best, beam = beam_search(RI.eph(), excluded=sorted(set(ALL) - set(L)), m0=1600.0, t_launch_grid=np.arange(0, 1501, 20) * DAY,
                         P=P, n_proc=1, verbose=False)
byn = {}
for s in list(P.collect) + list(beam):
    k = len(s.seq)
    if k not in byn or s.fuel < byn[k].fuel: byn[k] = s
front = {k: round(float(CM.tank(byn[k].fuel, 1600.0))) for k in sorted(byn)}
res = dict(pool=len(L), dv=dv, dr=dr, npt=npt, mc=mc, beam_s=round(time.time() - tic), front=front, settled={})
for k in sorted(byn, reverse=True)[:3]:
    try:
        signal.signal(signal.SIGALRM, _timeout); signal.setitimer(signal.ITIMER_REAL, 240)
        ip, miss, lag = IM.settle_tour(RI.eph(), tour_json(byn[k]), lambda ip: RI.settle(ip, 100))
    except Exception:
        ip = None; miss = float('inf')
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    res['settled'][k] = None if (ip is None or miss > 150) else round(float(ip.tank()), 1)
print(json.dumps(res)); json.dump(res, open(out, 'w'), indent=1)
