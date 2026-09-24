"""Validator interpolation + integration vs the PDF Sec 6.2 example lines: segment-wise from the PRINTED state of each line
(exactly what the checker does). Mass to 1e-7 kg (printed resolution), position to 10 m (printed resolution)."""
import sys, pathlib, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from ctoc14.validator import integrate_segment, validate
from ctoc14.kepler import propagate_twobody, Ephemeris
from ctoc14.thrust import thrust_interp
from ctoc14.constants import DAY, VE

L = {1: (0, 0.0000000000e+00, [-1.9500328647e+07, 1.4581848347e+08, -7.6372048832e+03], [-2.9389477847e+01, -4.3166912970e+00, -3.9417153392e+00], 2.0000000000e+03, [0, 0, 0], 0),
     2: (3, 1.0008479152e+07, [-1.1398611209e+08, -8.7213961221e+07, -1.6523817305e+07], [1.8269315350e+01, -2.4357060498e+01, 1.9760173753e+00], 2.0000000000e+03, [0, 0, 0], 174),
     3: (1, 1.0008489152e+07, [-1.1398592940e+08, -8.7214204791e+07, -1.6523797545e+07], [1.8269365516e+01, -2.4357022114e+01, 1.9760246476e+00], 2.0000000000e+03, [0, 0, 0], 0),
     4: (1, 1.0094889152e+07, [-1.1238880278e+08, -8.9304194192e+07, -1.6350356961e+07], [1.8700314237e+01, -2.4020840117e+01, 2.0387760044e+00], 1.9999567651e+03, [2.4534511893e-02, 2.8330454030e-02, 1.1593124092e-02], 0),
     5: (1, 1.0181289152e+07, [-1.1075466879e+08, -9.1364740090e+07, -1.6171502174e+07], [1.9125924514e+01, -2.3675430238e+01, 2.1013494039e+00], 1.9998273580e+03, [3.9479106603e-02, 6.3443225699e-02, 2.3114772719e-02], 0)}
# reconstructed 41-sample arc (docs/thrust_interpolation_rule.md)
t3 = L[3][1]; ts = t3 + np.arange(41) * DAY
mag = 0.5 * np.sin(np.pi * np.arange(41) / 40.0)
az = np.radians(40.107046 + 9.0 * np.arange(41)); el = np.radians(17.188734)
Ts = mag[:, None] * np.stack([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el) * np.ones(41)], axis=1)
print('reconstruction error of printed samples k=1,2 [N]:', np.abs(Ts[1] - L[4][5]).max(), np.abs(Ts[2] - L[5][5]).max())

def seg(a, b, arc):
    ev, t0, r0, v0, m0, _, _ = L[a]; _, t1, r1, v1, m1, _, _ = L[b]
    rr, vv, mm = integrate_segment(np.array(r0), np.array(v0), m0, t0, t1, ts if arc else None, Ts if arc else None)
    dr = np.linalg.norm(rr - np.array(r1)); dv = np.linalg.norm(vv - np.array(v1)); dm = mm - m1
    print(f'line {a}->{b} ({"arc" if arc else "coast"}): |dr| = {dr*1e3:8.2f} m  |dv| = {dv*1e6:8.3f} mm/s  dm = {dm:+.3e} kg   (printed resolution: 10 m, 1e-6 km/s, 1e-7 kg)')
    return dr, dv, dm

d12 = seg(1, 2, False); d23 = seg(2, 3, False); d34 = seg(3, 4, True); d45 = seg(4, 5, True)
# also: mass with the two candidate window readings, computed independently by quadrature of |T(t)|
from scipy.integrate import quad
def T_B(t):  # reading B: forward window {j, j+1, j+2, j+3} clipped
    n = len(ts)
    if t < ts[0] or t > ts[-1]: return np.zeros(3)
    j = min(max(int(np.searchsorted(ts, t, side='right') - 1), 0), n - 2); l = min(j, n - 4); idx = range(l, l + 4)
    out = np.zeros(3)
    for i in idx:
        w = 1.0
        for k in idx:
            if k != i: w *= (t - ts[k]) / (ts[i] - ts[k])
        out += w * Ts[i]
    return out
for name, f in [('A (code)', lambda t: np.linalg.norm(thrust_interp(ts, Ts, t))), ('B (forward window)', lambda t: np.linalg.norm(T_B(t)))]:
    dm34 = quad(f, ts[0], ts[1], epsabs=1e-12, limit=200)[0] / (VE * 1e3); dm45 = quad(f, ts[1], ts[2], epsabs=1e-12, limit=200)[0] / (VE * 1e3)
    print(f'reading {name:20s}: dm(3->4) = {dm34:.9f} (PDF {2000 - 1.9999567651e+03:.7f}, diff {dm34 - (2000 - 1.9999567651e+03):+.2e}); dm(4->5) = {dm45:.9f} (PDF {1.9999567651e+03 - 1.9998273580e+03:.7f}, diff {dm45 - (1.9999567651e+03 - 1.9998273580e+03):+.2e})')
ok = d34[2] < 1e-7 and d45[2] < 1e-7 and d34[0] < 0.01 and d45[0] < 0.01 and d12[0] < 0.01 and d23[0] < 0.01
print('PDF example reproduced to printed precision:', ok)
# asteroid 174 distance at the flyby line
eph = Ephemeris(); ra, _ = eph.ast_state(173, L[2][1]); print('asteroid 174 distance at line 2 [km]:', np.linalg.norm(np.array(L[2][2]) - ra))
