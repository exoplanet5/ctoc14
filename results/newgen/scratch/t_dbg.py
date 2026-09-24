import sys, time, numpy as np
sys.path.insert(0,'.')
from ctoc14.kepler import Ephemeris
from ctoc14.globalopt import load_craft, GlobalProblem, irls_min_fuel
from ctoc14.constants import DAY
eph = Ephemeris()
c = load_craft('results/CTOC14_Result_TEAM.txt', 12)
gp = GlobalProblem(eph, c, hstep=0.25*DAY)
z = np.load('results/newgen/scratch/lin_sc12.npz')
B, Lv, dm, vr = z['B'], z['Lv'], z['dm'], z['vr']
for rho in [1e2, 1e4, 1e6]:
  for wt, wv in [(1e-6,1e-2),(1.0,1e6)]:
    T, dt, dv = irls_min_fuel(B, Lv, vr, dm, gp.Ts, gp.vinf, 0.45, rho=rho, w_t=wt, w_v=wv)
    Yf2,_ = gp.integrate(Ts=T, tf=gp.tf+dt*DAY, vinf=gp.vinf+dv)
    m2,_ = gp.misses(Yf2, gp.tf+dt*DAY)
    # predicted
    pred = dm + np.einsum('jsac,sc->ja', B, T-gp.Ts) + Lv@dv + vr*dt[:,None]*DAY
    print(f'rho {rho:g} wt {wt} wv {wv}: fuel {gp.fuel(T):.2f}, max|dT| {np.linalg.norm(T-gp.Ts,axis=1).max():.2e}, max|dt| {np.abs(dt).max():.3f} d, |dv| {np.linalg.norm(dv)*1e3:.2f} m/s, pred miss {np.linalg.norm(pred,axis=1).max():.1f}, true miss max {np.linalg.norm(m2,axis=1).max():.0f} per-j {np.round(np.linalg.norm(m2,axis=1)[:8])}')
