"""Are the routes converged? Longer SCP on two routes; report tank and the plane part."""
import os
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'): os.environ.setdefault(v,'1')
import sys, time, numpy as np
sys.path.insert(0, '.'); sys.path.insert(0, 'tools')
from run_ialns import IFleet, ipr, cost, QUIET
from ctoc14.globalopt import optimise, restore
from ctoc14.phasemodel import elements, K_PLANE, K_DRIFT, K_TOTAL, drift
F = IFleet('results/newgen/ifleet10')
for n in ('12', '09'):
    ip = ipr(F.routes[n]['st']); t0 = ip.tank()
    def parts(ip):
        Yf, Ys = ip.integrate(); el = elements(Ys[:, :3], Ys[:, 3:6])
        return (K_TOTAL * K_DRIFT * np.abs(np.diff(drift(el['a']))).sum(),
                K_TOTAL * K_PLANE * np.linalg.norm(np.diff(el['i_vec'], axis=0), axis=1).sum())
    d0, p0 = parts(ip)
    tic = time.time(); optimise(ip, iters=400, rho0=100.0, log=QUIET); restore(ip, tol_km=100, iters=30, log=QUIET)
    d1, p1 = parts(ip)
    print(f'route {n}: tank {t0:.1f} -> {ip.tank():.1f} kg (dJ {cost(ip.tank())-cost(t0):+.4f}), '
          f'drift part {d0:.2f} -> {d1:.2f}, plane part {p0:.2f} -> {p1:.2f} ({time.time()-tic:.0f} s)', flush=True)
