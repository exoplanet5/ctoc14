import os
os.environ['OMP_NUM_THREADS']='1'; os.environ['VECLIB_MAXIMUM_THREADS']='1'
import sys, numpy as np
sys.path.insert(0,'.')
from ctoc14.kepler import Ephemeris
from ctoc14.globalopt import problem_from_state, load_state, add_flyby, extend_arc, irls_min_fuel, to_massless
from ctoc14.constants import DAY
eph = Ephemeris()
fleet = sys.argv[1]
gp = to_massless(problem_from_state(eph, load_state(f'{fleet}/craft_04.npz')))
extend_arc(gp, 3173*DAY); add_flyby(gp, 11, 3168*DAY)
Yf, Ys = gp.integrate(dense_samples=True); dm, vr = gp.misses(Yf); d = np.linalg.norm(dm, axis=1)
j = int(np.argmax(d)); print('miss', d[j], 'others max', np.delete(d,j).max(), 'vrel', [round(float(np.linalg.norm(v)),2) for v in vr])
B, Lv = gp.linearise(Yf, Ys)
for rho in [1e3]:
    T, dt, dv = irls_min_fuel(B, Lv, vr, dm, gp.Ts, gp.vinf, gp.tcap, rho=rho, w_v=1.0+1e-4*rho)
    print(f'rho {rho:g}: max|dT| {np.linalg.norm(T-gp.Ts,axis=1).max():.3e}, dt {np.round(dt,2)}, |dv| {np.linalg.norm(dv):.4f}')
    for a in [0.01, 0.1, 0.5]:
        Ta = gp.Ts + a*(T-gp.Ts); tfa = gp.tf + a*dt*DAY; va = gp.vinf + a*dv
        Y2,_ = gp.integrate(Ts=Ta, tf=tfa, vinf=va, dense_samples=True); m2,_ = gp.misses(Y2, tfa); d2 = np.linalg.norm(m2,axis=1)
        pred = dm + a*(np.einsum('jsac,sc->ja', B, T-gp.Ts) + Lv@dv + vr*dt[:,None]*DAY)
        print(f'  a={a}: new {d2[j]:.3e} (pred {np.linalg.norm(pred[j]):.3e}); others {np.round(np.delete(d2,j)).astype(int).tolist()}')
