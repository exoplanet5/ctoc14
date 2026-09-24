import sys, time, numpy as np
sys.path.insert(0,'.')
from ctoc14.kepler import Ephemeris
from ctoc14.globalopt import load_craft, GlobalProblem, irls_min_fuel
from ctoc14.constants import DAY
eph = Ephemeris()
c = load_craft('results/CTOC14_Result_TEAM.txt', 12)
gp = GlobalProblem(eph, c, hstep=0.25*DAY)
Yf, Ys = gp.integrate(dense_samples=True)
dm, vr = gp.misses(Yf)
tic=time.time(); B, Lv = gp.linearise(Yf, Ys); print('linearise', time.time()-tic)
for s in [100, 1500, 3000]:
    for ax in [0,1]:
        T2 = gp.Ts.copy(); T2[s, ax] += 0.002
        Yf2,_ = gp.integrate(Ts=T2, dense_samples=True)
        d_true = Yf2[:,0:3]-Yf[:,0:3]; d_lin = B[:, s, :, ax]*0.002
        print(f's={s} ax={ax}: max true {np.linalg.norm(d_true,axis=1).max():.1f} km, lin err {np.linalg.norm(d_true-d_lin,axis=1).max():.3f} km')
dvv = np.array([0, 1e-4, 0]); Yf3,_ = gp.integrate(vinf=gp.vinf+dvv, dense_samples=True)
print('vinf check', np.linalg.norm(Yf3[:,0:3]-Yf[:,0:3],axis=1).max(), np.linalg.norm((Yf3[:,0:3]-Yf[:,0:3]) - Lv@dvv,axis=1).max())
for rho in [1e2, 1e4]:
    T, dt, dv = irls_min_fuel(B, Lv, vr, dm, gp.Ts, gp.vinf, 0.45, rho=rho, w_t=1.0, w_v=1e6)
    Yf2,_ = gp.integrate(Ts=T, tf=gp.tf+dt*DAY, vinf=gp.vinf+dv, dense_samples=True)
    m2,_ = gp.misses(Yf2, gp.tf+dt*DAY)
    print(f'rho {rho:g}: fuel {gp.fuel(T):.2f}, max|dT| {np.linalg.norm(T-gp.Ts,axis=1).max():.2e}, true miss max {np.linalg.norm(m2,axis=1).max():.0f}')
