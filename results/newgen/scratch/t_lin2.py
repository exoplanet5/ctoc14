import sys, time, numpy as np
sys.path.insert(0,'.')
from ctoc14.kepler import Ephemeris
from ctoc14.globalopt import load_craft, GlobalProblem, irls_min_fuel
from ctoc14.constants import DAY, VE
eph = Ephemeris()
sc = int(sys.argv[1])
c = load_craft('results/CTOC14_Result_TEAM.txt', sc)
gp = GlobalProblem(eph, c, hstep=0.25*DAY)
z = np.load('results/newgen/scratch/lin_sc%d.npz'%sc)
B, Lv, dm, vr = z['B'], z['Lv'], z['dm'], z['vr']
for rho in [0.0, 1.0, 10.0]:
    tic=time.time()
    T, dt, dv = irls_min_fuel(B, Lv, vr, dm, gp.Ts, gp.vinf, 0.45, rho=rho)
    print(f'rho {rho}: linear-model fuel {gp.fuel(T):.1f} kg (ref {gp.fuel():.1f}); max|dt| {np.abs(dt).max():.2f} d; |dv| {np.linalg.norm(dv)*1e3:.1f} m/s; '
          f'max|dT| {np.linalg.norm(T-gp.Ts,axis=1).max():.3f} N; {time.time()-tic:.1f} s')
    Yf2,_ = gp.integrate(Ts=T, tf=gp.tf+dt*DAY, vinf=gp.vinf+dv)
    m2,_ = gp.misses(Yf2, gp.tf+dt*DAY)
    print('    nonlinear misses km: median %.0f max %.0f' % (np.median(np.linalg.norm(m2,axis=1)), np.linalg.norm(m2,axis=1).max()))
