"""How many asteroids can ONE phase curve cover, and at what Delta-v? (pricing DP, all prizes = 1)"""
import sys, time, numpy as np
sys.path.insert(0, '.'); sys.path.insert(0, 'tools')
from ctoc14.kepler import Ephemeris
from ctoc14.phasedp import catalogue, PhaseDP
from ctoc14.search import UNREACHABLE
from ctoc14.constants import DAY
E = Ephemeris(); ids = [k for k in range(1, 301) if k not in UNREACHABLE]
tic = time.time(); cat = catalogue(E, ids); print(f'catalogue: {len(cat["k"])} samples, {len(set(cat["ast"]))} asteroids ({time.time()-tic:.0f} s)', flush=True)
np.savez('results/newgen/scratch/cat.npz', **{k: v for k, v in cat.items() if k != 'dt'}, dt=cat['dt'])
tic = time.time(); dp = PhaseDP(cat, step=5); print(f'DP grid: {dp.nt} steps x {dp.nphi} phi x {dp.ns} slopes ({time.time()-tic:.0f} s)', flush=True)
prize = np.ones(len(ids))
for lam in (0.02, 0.05, 0.1, 0.2):
    tic = time.time(); val, path = dp.solve(prize, lam)
    cov = dp.covered(path); dv = dp.dv_of(path)
    print(f'lam {lam}: value {val:.2f}, covers {len(cov)} asteroids, dv {dv:.1f} km/s '
          f'({len(cov)/max(dv,1e-9):.2f} targets/(km/s), {time.time()-tic:.0f} s)', flush=True)
