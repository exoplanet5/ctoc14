import sys, pathlib, numpy as np, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ctoc14.lambert import lambert
from ctoc14.kepler import Ephemeris, propagate_twobody
eph = Ephemeris()
# PDF example: Earth at t=0 -> asteroid 174 at t2, ballistic
t2 = 1.0008479152e+07
rE, vE = eph.earth_state(0.0); ra, va = eph.ast_state(173, t2)
v1, v2 = lambert(rE, ra, t2)
print('Lambert v1 =', v1, ' PDF v1 = [-29.389477847 -4.3166912970 -3.9417153392]  v_inf =', np.linalg.norm(v1 - vE))
r_chk, v_chk = propagate_twobody(rE, v1, t2); print('propagated end error km:', np.linalg.norm(r_chk - ra), ' v2 err m/s:', 1e3 * np.linalg.norm(v_chk - v2))
# random consistency tests incl. long way and near-180 deg
rng = np.random.default_rng(1)
N = 20000
t0 = rng.uniform(0, 4e8, N); tof = rng.uniform(20, 900, N) * 86400
k1 = rng.integers(0, 300, N); k2 = rng.integers(0, 300, N)
R1 = np.array([eph.ast_state(k, t)[0] for k, t in zip(k1, t0)]); R2 = np.array([eph.ast_state(k, t + d)[0] for k, t, d in zip(k2, t0, tof)])
tic = time.time(); V1, V2 = lambert(R1, R2, tof); dt = time.time() - tic
ok = np.isfinite(V1[:, 0])
print(f'{N} Lambert solves in {dt:.2f} s ({dt / N * 1e6:.1f} us each); solved {ok.sum()}')
err = np.array([np.linalg.norm(propagate_twobody(R1[i], V1[i], tof[i])[0] - R2[i]) for i in np.where(ok)[0][:2000]])
print('max endpoint error over 2000 checked solutions [km]:', err.max(), ' median', np.median(err))
