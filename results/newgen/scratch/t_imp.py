import os
os.environ['OMP_NUM_THREADS']='1'; os.environ['VECLIB_MAXIMUM_THREADS']='1'
import sys, time, numpy as np
sys.path.insert(0,'.')
from ctoc14.kepler import Ephemeris
from ctoc14.globalopt import problem_from_state, load_state, optimise, restore
from ctoc14.impulsive import from_thrust_problem
from ctoc14.constants import DAY, VE
eph = Ephemeris(); host = sys.argv[1]; dtau = float(sys.argv[2])
gp = problem_from_state(eph, load_state(f'results/newgen/fleet_dd2/craft_{host}.npz'))
dv_exact = VE * np.log(gp.m0 / (gp.m0 - gp.fuel()))
ip = from_thrust_problem(gp, dtau=dtau*DAY)
tic = time.time(); Yf, Ys = ip.integrate(); dm, vr = ip.misses(Yf); print(f'{host}: exact dv {dv_exact:.3f} km/s (m0 {gp.m0:.1f}); binned impulses dv {ip.dv():.3f}, misses max {np.linalg.norm(dm,axis=1).max():.0f} km ({time.time()-tic:.3f} s)')
tic = time.time(); B, Lv = ip.linearise(Yf, Ys); print('linearise', time.time()-tic)
# FD check
m = len(ip.ts)//3; T2 = ip.Ts.copy(); T2[m,1] += 1e-3; Y2,_ = ip.integrate(Ts=T2)
print('lin err', np.linalg.norm((Y2[:,:3]-Yf[:,:3]) - B[:,m,:,1]*1e-3, axis=1).max(), 'of', np.linalg.norm(Y2[:,:3]-Yf[:,:3],axis=1).max())
tic = time.time()
d = restore(ip, tol_km=100, iters=30, rho=1e3, log=lambda s: None); print(f'restored: dv {ip.dv():.3f} miss {d.max():.0f} ({time.time()-tic:.1f} s)')
optimise(ip, iters=150, rho0=100.0, log=lambda s: None); d = restore(ip, tol_km=100, iters=30, log=lambda s: None)
print(f'optimised: dv {ip.dv():.3f} km/s, tank {ip.tank():.1f} kg (exact {gp.m0:.1f}), miss {d.max():.0f} ({time.time()-tic:.1f} s)')
