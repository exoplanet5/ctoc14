"""Self-consistent regime: near-circular ecliptic craft meeting asteroids at their node crossings.
r within r_tol of the craft's a (no eccentricity needed for the radius), phase slack +-2 deg (e ~ 0.02),
|z| <= z_tol (node crossing). Under these conditions the phase DP is sound."""
import sys, numpy as np
sys.path.insert(0, '.')
from ctoc14.phasedp import PhaseDP
from ctoc14.search import UNREACHABLE
ev = dict(np.load('results/newgen/scratch/events2.npz'))
ids = [k for k in range(1, 301) if k not in UNREACHABLE]
prize = np.ones(len(ids))
for r_tol, half, z_tol in [(0.03, 2.0, 0.005), (0.05, 3.0, 0.010), (0.08, 4.0, 0.015)]:
    dp = PhaseDP(ev, half_max=half); dp.r_tol = r_tol
    n = dp.set_plane((0.0, 0.0), z_tol)
    out = []
    for lam in (1.0, 2.0, 4.0, 8.0):
        val, path = dp.solve(prize, lam=lam)
        col = dp.collected(path); asts = {int(ev['ast'][j]) for j in col}; dv = dp.dv_of(path)
        out.append(f'lam {lam:4.1f}: {len(asts):3d} @ {dv:5.1f} km/s')
    print(f'r_tol {r_tol} half {half} z_tol {z_tol} ({n} events admissible) | ' + ' | '.join(out), flush=True)
