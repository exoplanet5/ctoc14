"""Check ephemeris + two-body propagation against the worked example in the CTOC14 problem statement (Sec. 6.2)."""
import numpy as np, sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ctoc14.kepler import Ephemeris, propagate_twobody
from ctoc14.constants import VE, DAY

eph = Ephemeris()
# line 1: launch, t=0
r1 = np.array([-1.9500328647e+07, 1.4581848347e+08, -7.6372048832e+03]); v1 = np.array([-2.9389477847e+01, -4.3166912970e+00, -3.9417153392e+00])
rE, vE = eph.earth_state(0.0)
print('Earth r(t0) error [km]  :', np.linalg.norm(rE - r1), '  Earth r =', rE)
print('launch v_inf [km/s]     :', np.linalg.norm(v1 - vE), ' vE =', vE)
# line 2: flyby of asteroid 174 at t = 1.0008479152e7 s (ballistic coast from launch)
t2 = 1.0008479152e+07
r2 = np.array([-1.1398611209e+08, -8.7213961221e+07, -1.6523817305e+07]); v2 = np.array([1.8269315350e+01, -2.4357060498e+01, 1.9760173753e+00])
rp, vp = propagate_twobody(r1, v1, t2)
print('coast to line 2: pos err [km] =', np.linalg.norm(rp - r2), ' vel err [m/s] =', 1e3 * np.linalg.norm(vp - v2))
ra, va = eph.ast_state(174 - 1, t2)
print('asteroid 174 at t2: distance to SC [km] =', np.linalg.norm(ra - r2), ' (must be <= 1000 km)')
# line 3: 10000 s later, thrust arc begins with T=0
t3 = 1.0008489152e+07
r3 = np.array([-1.1398592940e+08, -8.7214204791e+07, -1.6523797545e+07]); v3 = np.array([1.8269365516e+01, -2.4357022114e+01, 1.9760246476e+00])
rp3, vp3 = propagate_twobody(r2, v2, t3 - t2)
print('coast to line 3: pos err [km] =', np.linalg.norm(rp3 - r3), ' vel err [m/s] =', 1e3 * np.linalg.norm(vp3 - v3))
# lines 3-5: mass drop consistency for linear-ish ramp 0 -> |T4| over 86400 s
T4 = np.array([2.4534511893e-02, 2.8330454030e-02, 1.1593124092e-02]); T5 = np.array([3.9479106603e-02, 6.3443225699e-02, 2.3114772719e-02])
m3, m4, m5 = 2000.0, 1.9999567651e+03, 1.9998273580e+03
print('|T4| =', np.linalg.norm(T4), ' |T5| =', np.linalg.norm(T5), ' dm(3->4) =', m3 - m4, ' dm(4->5) =', m4 - m5)
print('trapezoid estimate dm(3->4) =', 0.5 * np.linalg.norm(T4) * DAY / (VE * 1e3), ' dm(4->5) =', 0.5 * (np.linalg.norm(T4) + np.linalg.norm(T5)) * DAY / (VE * 1e3))
# ballistic-flyby reachability check: is 174 reachable ballistically? show v_inf & tof
print('TOF to 174 [d] =', t2 / DAY)
