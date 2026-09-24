"""One designed phase curve: how many distinct asteroids, at what Delta-v, for a few orbit planes?"""
import sys, time, numpy as np
sys.path.insert(0, '.')
from ctoc14.kepler import Ephemeris
from ctoc14.phasedp import PhaseDP, events
from ctoc14.search import UNREACHABLE
E = Ephemeris(); ids = [k for k in range(1, 301) if k not in UNREACHABLE]
tic = time.time(); ev = events(E, ids)
print(f'{len(ev["t"])} events over {len(set(ev["ast"]))} asteroids ({time.time()-tic:.0f} s); '
      f'|z| p50 {np.percentile(abs(ev["z"]),50):.3f} p90 {np.percentile(abs(ev["z"]),90):.3f} AU', flush=True)
np.savez('results/newgen/scratch/events2.npz', **ev)
dp = PhaseDP(ev)
prize = np.ones(len(ids))
for iv, ztol in [((0.0, 0.0), 0.010), ((0.0, 0.0), 0.02), ((0.03, 0.0), 0.02), ((0.0, 0.03), 0.02), ((0.02, 0.02), 0.02)]:
    n = dp.set_plane(iv, ztol)
    tic = time.time(); val, path = dp.solve(prize, lam=0.05)
    col = dp.collected(path); asts = {int(dp.ev['ast'][j]) for j in col}; dv = dp.dv_of(path)
    print(f'plane {iv} ztol {ztol}: {n} admissible events -> value {val:.2f}, {len(col)} events, '
          f'{len(asts)} asteroids, dv {dv:.1f} km/s ({time.time()-tic:.0f} s)', flush=True)
