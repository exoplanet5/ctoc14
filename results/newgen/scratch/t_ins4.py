import os
os.environ['OMP_NUM_THREADS']='1'; os.environ['VECLIB_MAXIMUM_THREADS']='1'
import sys, time, numpy as np
sys.path.insert(0,'.')
from ctoc14.kepler import Ephemeris
from ctoc14.globalopt import problem_from_state, load_state, add_flyby, extend_arc, optimise, restore, to_massless, to_thrust, tighten
from ctoc14.constants import DAY, M_DRY
def cost(m0): x=(m0-600)/1400; return 1+x+x*x
eph = Ephemeris()
fleet, host, ast, t, iters, rho_r = sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6])
gp0 = problem_from_state(eph, load_state(f'{fleet}/craft_{host}.npz'))
tic = time.time(); log = []
gp = to_massless(gp0); f0 = gp.fuel(); extend_arc(gp, t*DAY + 5*DAY); add_flyby(gp, ast, t*DAY)
d = restore(gp, tol_km=300, iters=40, rho=rho_r, log=log.append); r_it = len(log)
print(f'[{host} {ast} rho_r={rho_r}] restored in {r_it} its: fuel {gp.fuel():.1f} (was {f0:.1f}) miss {d.max():.0f} ({time.time()-tic:.0f} s)', flush=True)
if d.max() > 300: sys.exit()
optimise(gp, iters=iters, rho0=100.0, log=lambda s: None)
d = restore(gp, tol_km=150, iters=20, log=lambda s: None)
gt = to_thrust(gp, margin=1.5); restore(gt, tol_km=150, iters=20, log=lambda s: None); d = tighten(gt, 1.5, log=lambda s: None)
print(f'[{host} {ast} rho_r={rho_r}] final m0 {gt.m0:.1f} (was {gp0.m0:.1f}) dJ {cost(gt.m0)-cost(gp0.m0):+.4f} miss {d.max():.0f} ({time.time()-tic:.0f} s)', flush=True)
