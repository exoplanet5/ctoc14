"""Reconstruct the PDF example thrust arc (|T| = 0.5 sin(pi t/40 d), direction rotating 9 deg/day about z, elevation ~17.2 deg,
41 samples at 86400 s) and check that the mass drop across the first two intervals matches the printed values to << 0.01 kg."""
import sys, pathlib, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ctoc14.thrust import thrust_interp, arc_max_thrust
from ctoc14.constants import VE, DAY
from scipy.integrate import quad
T4 = np.array([2.4534511893e-02, 2.8330454030e-02, 1.1593124092e-02]); T5 = np.array([3.9479106603e-02, 6.3443225699e-02, 2.3114772719e-02])
u4 = T4 / np.linalg.norm(T4); u5 = T5 / np.linalg.norm(T5)
az4 = np.degrees(np.arctan2(u4[1], u4[0])); az5 = np.degrees(np.arctan2(u5[1], u5[0])); el = np.degrees(np.arcsin(u4[2]))
print(f'|T4|={np.linalg.norm(T4):.6f} (0.5 sin 4.5deg = {0.5*np.sin(np.radians(4.5)):.6f}); az4={az4:.3f} az5={az5:.3f} (rate {az5-az4:.3f} deg/d) el={el:.3f}')
ts = np.arange(41) * DAY
mag = 0.5 * np.sin(np.pi * ts / (40 * DAY))
az = np.radians(az4 + (az5 - az4) * (np.arange(41) - 1)); elr = np.radians(el)
Ts = mag[:, None] * np.stack([np.cos(elr) * np.cos(az), np.cos(elr) * np.sin(az), np.sin(elr) * np.ones(41)], axis=1)
print('reconstructed T4, T5 error:', np.abs(Ts[1] - T4).max(), np.abs(Ts[2] - T5).max())
def dm(t0, t1):
    return quad(lambda t: np.linalg.norm(thrust_interp(ts, Ts, t)), t0, t1, epsabs=1e-10, limit=200)[0] / (VE * 1e3)
print(f'dm(3->4) model = {dm(0, DAY):.7f}  PDF = {2000.0 - 1.9999567651e+03:.7f}')
print(f'dm(4->5) model = {dm(DAY, 2*DAY):.7f}  PDF = {1.9999567651e+03 - 1.9998273580e+03:.7f}')
print('max |T| over arc =', arc_max_thrust(ts, Ts))
