"""Direct test: is a route's tangential Delta-v = k * total variation of its phase drift rate, and its normal Delta-v
= v * total variation of its orbit-normal direction? Measured on the craft's own trajectory (20 d sampling)."""
import sys, numpy as np
sys.path.insert(0, '.'); sys.path.insert(0, 'tools')
from run_ialns import IFleet, ipr, eph
from ctoc14.constants import AU, DAY, MU
YR = 365.25 * DAY; E = eph(); F = IFleet('results/newgen/ifleet10')
print('route | dv | dvT | k*TV(drift) | dvN | v*TV(normal) | drift range deg/yr')
K = 0.028; tot = []
for n, r in sorted(F.routes.items()):
    ip = ipr(r['st']); Yf, Ys = ip.integrate()
    rr = Ys[:, :3]; vv = Ys[:, 3:6]
    R = rr / np.linalg.norm(rr, axis=1)[:, None]; N = np.cross(rr, vv); Nn = N / np.linalg.norm(N, axis=1)[:, None]
    T = np.cross(Nn, R)
    dvT = np.abs((ip.Ts * T).sum(1)).sum(); dvN = np.abs((ip.Ts * Nn).sum(1)).sum(); dvR = np.abs((ip.Ts * R).sum(1)).sum()
    # semi-major axis along the route -> drift rate relative to Earth [deg/yr]
    v2 = (vv ** 2).sum(1); rn = np.linalg.norm(rr, axis=1)
    a = 1 / (2 / rn - v2 / MU)
    nn = np.degrees(np.sqrt(MU / a ** 3)) * YR; nE = np.degrees(np.sqrt(MU / (1.00091750 * AU) ** 3)) * YR
    drift = nn - nE
    tv_drift = np.abs(np.diff(drift)).sum()
    # orbit normal direction change [rad]
    cosang = np.clip((Nn[:-1] * Nn[1:]).sum(1), -1, 1)
    tv_norm = np.arccos(cosang).sum()
    v = 29.8
    tot.append((ip.dv(), dvT, K * tv_drift, dvN, v * tv_norm))
    print(f'{n} | {ip.dv():6.2f} | {dvT:6.2f} | {K * tv_drift:6.2f} | {dvN:5.2f} | {v * tv_norm:6.2f} | '
          f'{drift.min():7.1f} .. {drift.max():6.1f}')
tot = np.array(tot)
print(f'tangential: corr {np.corrcoef(tot[:,1], tot[:,2])[0,1]:.3f}, mean ratio pred/actual {np.mean(tot[:,2]/tot[:,1]):.2f}')
print(f'normal:     corr {np.corrcoef(tot[:,3], tot[:,4])[0,1]:.3f}, mean ratio pred/actual {np.mean(tot[:,4]/tot[:,3]):.2f}')
