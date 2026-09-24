import os
os.environ['OMP_NUM_THREADS']='1'; os.environ['VECLIB_MAXIMUM_THREADS']='1'
import sys, time, numpy as np
sys.path.insert(0,'.')
from ctoc14.kepler import Ephemeris
from ctoc14.globalopt import problem_from_state, load_state, optimise, restore, to_thrust, tighten, write_craft
from ctoc14.impulsive import from_thrust_problem, to_exact
from ctoc14.constants import DAY
def cost(m0): x=(m0-600)/1400; return 1+x+x*x
eph = Ephemeris(); host = sys.argv[1]
gp = problem_from_state(eph, load_state(f'results/newgen/fleet_dd2/craft_{host}.npz'))
ip = from_thrust_problem(gp); restore(ip, tol_km=100, iters=30, rho=1e3, log=lambda s: None); optimise(ip, iters=150, rho0=100.0, log=lambda s: None); restore(ip, tol_km=100, iters=30, log=lambda s: None)
print(f'impulsive tank {ip.tank():.1f} (exact {gp.m0:.1f})'); tic = time.time()
gm = to_exact(ip); Yf,_ = gm.integrate(dense_samples=True); print('initial continuous miss', np.linalg.norm(gm.misses(Yf)[0],axis=1).max())
d = restore(gm, tol_km=300, iters=40, rho=1e3, log=lambda s: None); print(f'restored massless miss {d.max():.0f} fuel@{gm.m0:.0f} {gm.fuel():.1f} ({time.time()-tic:.0f} s)')
optimise(gm, iters=40, rho0=100.0, log=lambda s: None); restore(gm, tol_km=150, iters=20, log=lambda s: None)
gt = to_thrust(gm, 1.5); restore(gt, tol_km=150, iters=20, log=lambda s: None); d = tighten(gt, 1.5, log=lambda s: None)
print(f'exact: m0 {gt.m0:.1f} miss {d.max():.0f} ({time.time()-tic:.0f} s)')
rep, fuel, pk = write_craft(gt, f'results/newgen/scratch/imp_exact_{host}.txt', eph=eph)
print('valid', rep.ok, rep.errors[:3], 'flybys', len(rep.flybys), 'peak', pk, f'({time.time()-tic:.0f} s)')
