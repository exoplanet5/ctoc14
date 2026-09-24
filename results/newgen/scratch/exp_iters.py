import os
os.environ['OMP_NUM_THREADS']='1'
import sys, time, numpy as np
sys.path.insert(0,'.'); sys.path.insert(0,'tools')
from run_ialns import IFleet, ipr, eph, QUIET
from ctoc14.impulsive import to_exact
from ctoc14.globalopt import restore, optimise, to_thrust, tighten, problem_from_state, state_of
from ctoc14.constants import DAY
fleet = IFleet(sys.argv[1]); n = sys.argv[2]; iters = int(sys.argv[3])
ip = ipr(fleet.routes[n]['st']); gm = to_exact(ip, eph()); tic = time.time()
restore(gm, tol_km=300, iters=40, rho=1e3, log=QUIET)
for k in range(iters // 40):
    optimise(gm, iters=40, rho0=100.0, log=QUIET); restore(gm, tol_km=150, iters=20, log=QUIET)
    gt = to_thrust(gm, 1.5)
    print(f'{n}: after {40*(k+1)} its massless dv-equivalent tank {gt.m0:.1f} (impulsive {ip.tank():.1f}) ({time.time()-tic:.0f} s)', flush=True)
