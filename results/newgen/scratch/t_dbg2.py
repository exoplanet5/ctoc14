import os
os.environ['OMP_NUM_THREADS']='1'; os.environ['VECLIB_MAXIMUM_THREADS']='1'
import sys, numpy as np
sys.path.insert(0,'.')
from ctoc14.kepler import Ephemeris
from ctoc14.globalopt import problem_from_state, load_state, add_flyby, extend_arc, irls_min_fuel
from ctoc14.constants import DAY
eph = Ephemeris()
st = load_state('results/newgen/fleet_dd/craft_04.npz'); gp = problem_from_state(eph, st)
extend_arc(gp, 3173*DAY); add_flyby(gp, 11, 3168*DAY)
Yf, Ys = gp.integrate(dense_samples=True); dm, vr = gp.misses(Yf); d = np.linalg.norm(dm, axis=1)
j = int(np.argmax(d)); print('new flyby index', j, 'miss', d[j], 'vrel', np.linalg.norm(vr[j]), 'tf', gp.tf[j]/DAY, 'ts range', gp.ts[0]/DAY, gp.ts[-1]/DAY)
print('miss dir . vrel', np.dot(dm[j], vr[j])/np.linalg.norm(dm[j])/np.linalg.norm(vr[j]))
B, Lv = gp.linearise(Yf, Ys)
for rho, wt, wv in [(1e2, 0, 1.0), (1e5, 0, 1.0), (1e7, 0, 1.0)]:
    T, dt, dv = irls_min_fuel(B, Lv, vr, dm, gp.Ts, gp.vinf, gp.tcap, rho=rho, w_v=wv)
    print(f'rho {rho:g} wt {wt} wv {wv}: max|dT| {np.linalg.norm(T-gp.Ts,axis=1).max():.3e} N, dt[j] {dt[j]:.2f} d, max|dt| other {np.abs(np.delete(dt,j)).max():.2f}, |dv| {np.linalg.norm(dv):.4f} km/s, fuel {gp.fuel(T):.1f}')
    for a in [0.01, 0.1, 0.3]:
        Ta = gp.Ts + a*(T-gp.Ts); tfa = gp.tf + a*dt*DAY; va = gp.vinf + a*dv
        Y2,_ = gp.integrate(Ts=Ta, tf=tfa, vinf=va, dense_samples=True); m2,_ = gp.misses(Y2, tfa); d2 = np.linalg.norm(m2,axis=1)
        pred = dm + a*(np.einsum('jsac,sc->ja', B, T-gp.Ts) + Lv@dv + vr*dt[:,None]*DAY)
        print(f'   a={a}: new miss {d2[j]:.3e} (pred {np.linalg.norm(pred[j]):.3e}); others max {np.delete(d2,j).max():.0f} (pred {np.linalg.norm(np.delete(pred,j,0),axis=1).max():.0f}); sum excess {np.maximum(0,d2-150).sum():.3e} vs {np.maximum(0,d-150).sum():.3e}')
