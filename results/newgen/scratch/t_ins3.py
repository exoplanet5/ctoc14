import os
os.environ['OMP_NUM_THREADS']='1'; os.environ['VECLIB_MAXIMUM_THREADS']='1'
import sys, time, numpy as np
sys.path.insert(0,'.')
from ctoc14.kepler import Ephemeris
from ctoc14.globalopt import problem_from_state, load_state, add_flyby, extend_arc, optimise, restore, to_massless, to_thrust, tighten
from ctoc14.constants import DAY
eph = Ephemeris()
host, ast, t, iters = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
gp0 = problem_from_state(eph, load_state(f'results/newgen/fleet_dd/craft_{host}.npz'))
print('before: m0', gp0.m0, 'fuel', gp0.fuel()); tic = time.time()
gp = to_massless(gp0); extend_arc(gp, t*DAY + 5*DAY); add_flyby(gp, ast, t*DAY)
d = restore(gp, tol_km=300, iters=40, rho=1e3); print(f'restored: fuel {gp.fuel():.1f} miss {d.max():.0f} ({time.time()-tic:.0f} s)')
optimise(gp, iters=iters, rho0=100.0)
d = restore(gp, tol_km=150, iters=20); print(f'massless opt: fuel {gp.fuel():.1f} dv {gp.dv():.3f} miss {d.max():.0f} ({time.time()-tic:.0f} s)')
gt = to_thrust(gp, margin=1.5); d = restore(gt, tol_km=150, iters=20); print(f'thrust: m0 {gt.m0:.1f} fuel {gt.fuel():.1f} miss {d.max():.0f}')
d = tighten(gt, 1.5); print(f'tight: m0 {gt.m0:.1f} fuel {gt.fuel():.1f} miss {d.max():.0f} ({time.time()-tic:.0f} s)')
