"""How far can a flyby sit from the craft's mean-longitude line for free? (eccentricity slack 2e sin(L - peri))
and: is the piecewise-linear mean-longitude curve through the flybys a good model of the route's drift history?"""
import sys, numpy as np
sys.path.insert(0, '.'); sys.path.insert(0, 'tools')
from run_ialns import IFleet, ipr, eph
from ctoc14.phasemodel import elements, drift, K_DRIFT, N_E, YR
from ctoc14.constants import AU, DAY, MU
E = eph(); F = IFleet('results/newgen/ifleet10')
def wrap(x): return (x + 180) % 360 - 180
print('route | flybys | |true-mean longitude| p50 p90 max (deg) | e range | TV(drift) traj | TV via flyby points | ratio')
allslack = []
for n, r in sorted(F.routes.items()):
    ip = ipr(r['st']); Yf, Ys = ip.integrate()
    o = np.argsort(ip.tf)
    el_f = elements(Yf[o, :3], Yf[o, 3:6])
    th = np.degrees(np.arctan2(Yf[o, 1], Yf[o, 0]))
    slack = wrap(th - np.degrees(el_f['L']))
    allslack.append(slack)
    el_s = elements(Ys[:, :3], Ys[:, 3:6])
    tv_traj = np.abs(np.diff(drift(el_s['a']))).sum()
    # piecewise linear through the flyby points in (t, mean longitude relative to Earth)
    t = ip.tf[o] / YR
    Lrel = np.degrees(np.unwrap(el_f['L'])) - N_E * t
    sl = np.diff(Lrel) / np.diff(t)
    tv_pts = np.abs(np.diff(sl)).sum()
    print(f'{n} | {len(t):2d} | {np.percentile(np.abs(slack),50):5.1f} {np.percentile(np.abs(slack),90):5.1f} '
          f'{np.abs(slack).max():5.1f} | {el_f["e"].min():.3f}-{el_f["e"].max():.3f} | {tv_traj:7.1f} | {tv_pts:7.1f} | {tv_pts/tv_traj:.2f}')
s = np.abs(np.concatenate(allslack))
print(f'all flybys: |true - mean longitude| p50 {np.percentile(s,50):.1f} p90 {np.percentile(s,90):.1f} max {s.max():.1f} deg')
