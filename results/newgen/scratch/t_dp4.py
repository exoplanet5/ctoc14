"""Calibrate the phase-slack cap against reality: the real fleet covers 24-37 asteroids per route at 11-19 km/s."""
import sys, numpy as np
sys.path.insert(0, '.')
from ctoc14.phasedp import PhaseDP
from ctoc14.search import UNREACHABLE
ev = dict(np.load('results/newgen/scratch/events2.npz'))
ids = [k for k in range(1, 301) if k not in UNREACHABLE]
prize = np.ones(len(ids))
for half in (1.0, 2.0, 4.0, 8.0):
    dp = PhaseDP(ev, half_max=half); dp.set_plane((0.0, 0.0), 0.015)
    out = []
    for lam in (1.0, 2.0, 4.0, 8.0, 16.0):
        val, path = dp.solve(prize, lam=lam)
        col = dp.collected(path); asts = {int(ev['ast'][j]) for j in col}; dv = dp.dv_of(path)
        out.append(f'lam {lam:4.1f}: {len(asts):3d} @ {dv:5.1f} km/s')
    print(f'phase slack +-{half:3.1f} deg | ' + ' | '.join(out), flush=True)
