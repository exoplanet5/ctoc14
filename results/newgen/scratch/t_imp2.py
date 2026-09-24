import os
os.environ['OMP_NUM_THREADS']='1'; os.environ['VECLIB_MAXIMUM_THREADS']='1'
import sys, time, numpy as np
sys.path.insert(0,'.')
from ctoc14.kepler import Ephemeris
from ctoc14.globalopt import problem_from_state, load_state, optimise, restore, insert_homotopy
from ctoc14.impulsive import from_thrust_problem
from ctoc14.constants import DAY, VE
def cost(m0): x=(m0-600)/1400; return 1+x+x*x
eph = Ephemeris(); host, ast, t = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
gp = problem_from_state(eph, load_state(f'results/newgen/fleet_dd2/craft_{host}.npz'))
ip = from_thrust_problem(gp, dtau=20*DAY)
tic = time.time(); restore(ip, tol_km=100, iters=30, rho=1e3, log=lambda s: None); optimise(ip, iters=150, rho0=100.0, log=lambda s: None); restore(ip, tol_km=100, iters=30, log=lambda s: None)
J0 = cost(ip.tank()); print(f'{host} base: tank {ip.tank():.1f} (exact {gp.m0:.1f}) ({time.time()-tic:.1f} s)')
tic = time.time()
d = insert_homotopy(ip, ast, t*DAY, stages=10, iters_stage=6, tol_km=100)
print(f'  inserted {ast}: miss {d.max():.0f}, tank {ip.tank():.1f} ({time.time()-tic:.1f} s)')
optimise(ip, iters=150, rho0=100.0, log=lambda s: None); d = restore(ip, tol_km=100, iters=30, log=lambda s: None)
print(f'  optimised: tank {ip.tank():.1f}, dJ {cost(ip.tank())-J0:+.4f}, miss {d.max():.0f} ({time.time()-tic:.1f} s)')
