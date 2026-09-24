"""Price sweep: coverage of one designed phase curve at a realistic Delta-v budget."""
import sys, time, numpy as np
sys.path.insert(0, '.')
from ctoc14.phasedp import PhaseDP
from ctoc14.search import UNREACHABLE
ev = dict(np.load('results/newgen/scratch/events2.npz'))
ids = [k for k in range(1, 301) if k not in UNREACHABLE]
dp = PhaseDP(ev); prize = np.ones(len(ids))
print('plane (0,0), z_tol 0.015')
dp.set_plane((0.0, 0.0), 0.015)
for lam in (0.2, 0.5, 1.0, 2.0, 4.0, 8.0):
    val, path = dp.solve(prize, lam=lam)
    col = dp.collected(path); asts = {int(ev['ast'][j]) for j in col}; dv = dp.dv_of(path)
    tank = 601.5 * np.exp(dv / 39.2266); x = (tank - 600) / 1400
    print(f'  lam {lam:4.1f}: {len(asts):3d} asteroids, dv {dv:6.1f} km/s -> tank {tank:6.0f} kg, '
          f'J_i {1 + x + x * x:.3f}, {len(asts) / max(dv, 1e-9):.2f} per km/s', flush=True)
