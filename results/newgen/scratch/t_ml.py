import os
os.environ['OMP_NUM_THREADS']='1'; os.environ['VECLIB_MAXIMUM_THREADS']='1'
import sys, time, numpy as np
sys.path.insert(0,'.')
from ctoc14.kepler import Ephemeris
from ctoc14.globalopt import problem_from_state, load_state, add_flyby, extend_arc, irls_min_fuel, to_massless, to_thrust
from ctoc14.constants import DAY
eph = Ephemeris()
st = load_state('results/newgen/fleet_dd/craft_04.npz'); gp0 = problem_from_state(eph, st)
Yf0, _ = gp0.integrate(dense_samples=True); print('thrust-mode max miss', np.linalg.norm(gp0.misses(Yf0)[0],axis=1).max(), 'fuel', gp0.fuel())
gp = to_massless(gp0)
Yf, Ys = gp.integrate(dense_samples=True); print('massless max miss', np.linalg.norm(gp.misses(Yf)[0],axis=1).max(), 'fuel', gp.fuel())
back = to_thrust(gp, margin=1.5); Yb,_ = back.integrate(dense_samples=True); print('back to thrust: m0', back.m0, 'fuel', back.fuel(), 'max miss', np.linalg.norm(back.misses(Yb)[0],axis=1).max())
extend_arc(gp, 3173*DAY); add_flyby(gp, 11, 3168*DAY)
Yf, Ys = gp.integrate(dense_samples=True); dm, vr = gp.misses(Yf); d = np.linalg.norm(dm, axis=1); j = int(np.argmax(d))
B, Lv = gp.linearise(Yf, Ys)
for rho in [1e2, 1e5]:
    T, dt, dv = irls_min_fuel(B, Lv, vr, dm, gp.Ts, gp.vinf, gp.tcap, rho=rho)
    print(f'rho {rho:g}: max|dT| {np.linalg.norm(T-gp.Ts,axis=1).max():.3e}, fuel {gp.fuel(T):.1f}')
    for a in [0.1, 0.3, 1.0]:
        Ta = gp.Ts + a*(T-gp.Ts); tfa = gp.tf + a*dt*DAY; va = gp.vinf + a*dv
        Y2,_ = gp.integrate(Ts=Ta, tf=tfa, vinf=va, dense_samples=True); m2,_ = gp.misses(Y2, tfa); d2 = np.linalg.norm(m2,axis=1)
        print(f'   a={a}: new miss {d2[j]:.3e}; others max {np.delete(d2,j).max():.0f}')
