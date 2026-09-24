"""Read-only probe: does keeping VELOCITY diversity in the beam raise sparse-pool depth?
The stock beam keeps n_per_target=2 TOF variants per (parent, target) and dedups on (visited, 5-d time bin), so two
states with the same targets and time but different arrival velocities collapse to one.  Here: n_per_target NPT,
max_children MC, dedup key (visited, 5-d bin, velocity bin of VB km/s).  Pool = gc16 remainder after routes 1-6 (87)
or 1-7 (67), strict exclusion, production legs.  usage: probe_vdiv.py OUT MODE NPT MC VB [BEAM]"""
import os, sys, json, time, signal, pathlib
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
_CWD0 = os.getcwd(); os.chdir(ROOT)
import numpy as np
import run_ialns as RI
import ctoc14.impulsive as IM
import ctoc14.search as S
from greedy_cover import BP, CM, DeepCollect, ALL, _timeout
from ctoc14.constants import DAY, AU
IM.LIN_FALLBACK = 2.5
out = os.path.join(_CWD0, sys.argv[1]); mode = sys.argv[2]; npt = int(sys.argv[3]); mc = int(sys.argv[4]); vb = float(sys.argv[5])
beamw = int(sys.argv[6]) if len(sys.argv) > 6 else 100


def beam_loop_v(eph, beam, P):
    beam = sorted(beam, key=lambda s: s.score(P))[:P.beam * 3]
    best = max(beam, key=lambda s: (s.n(), -s.score(P)))
    for depth in range(1, P.max_depth):
        children = [c for s in beam for c in S.expand(eph, s, P)]
        if not children: break
        children.sort(key=lambda s: s.score(P))
        seen = set(); nb = []
        for c in children:
            key = (c.visited, int(c.t // (5 * DAY))) + ((tuple(np.round(c.v / vb).astype(int)),) if vb > 0 else ())
            if key in seen: continue
            seen.add(key); nb.append(c)
            if len(nb) >= P.beam: break
        beam = nb; P.collect.extend(beam)
        if beam[0].n() > best.n(): best = beam[0]
    return best, beam


gc = RI.IFleet('results/newgen/gc16'); cov6 = set(); cov7 = set()
for n, r in gc.routes.items():
    s = set(int(x) for x in r['st']['asts']); cov7 |= s
    if n != '07': cov6 |= s
covd = cov6 if mode == 'p87' else cov7
pool = sorted(set(ALL) - covd)
prize = np.zeros(300)
for t in pool: prize[t - 1] = 1.0
P = BP(beam=beamw, w_fuel=CM.w_fuel(480, 1600), m_margin=40.0, dv_max=1.2, tofs=np.arange(15, 401, 5) * DAY,
       lin_tofs=np.arange(20, 601, 10) * DAY, lin_drmax=0.15 * AU, vinf_cap=4.0, tof_refine=True, max_depth=60,
       n_per_target=npt, max_children=mc)
P.prize = prize; P.collect = DeepCollect(8)
tic = time.time()
roots = S.roots(RI.eph(), np.arange(0, 801, 20) * DAY, 1600.0, sorted(covd), P)
roots.sort(key=lambda s: s.score(P)); seen = set(); b2 = []
for s in roots:
    key = (s.seq[0][0], int(s.t_launch // (20 * DAY)))
    if key in seen: continue
    seen.add(key); b2.append(s)
best, beam = beam_loop_v(RI.eph(), b2, P)
states = list(P.collect) + list(beam)
deep = max(states, key=lambda s: (len(s.seq), -s.fuel))
res = dict(mode=mode, npt=npt, mc=mc, vb=vb, beam=beamw, beam_s=round(time.time() - tic), n=len(deep.seq),
           planner=round(float(CM.tank(deep.fuel, 1600.0))))
try:
    signal.signal(signal.SIGALRM, _timeout); signal.setitimer(signal.ITIMER_REAL, 240)
    ip, miss, lag = IM.settle_tour(RI.eph(), S.tour_json(deep), lambda ip: RI.settle(ip, 100))
except Exception:
    ip = None; miss = float('inf')
finally:
    signal.setitimer(signal.ITIMER_REAL, 0)
res['twin'] = round(float(ip.tank())) if (ip is not None and miss <= 150) else None
json.dump(res, open(out, 'w')); print(json.dumps(res), flush=True)
