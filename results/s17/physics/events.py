"""Node events of every target in the (t, lag) plane: minima of the distance to the 1-AU ecliptic circle (torus distance),
lag = heliocentric longitude of the target minus Earth's at that epoch. Also the lag of every s16a flyby."""
import os
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'): os.environ.setdefault(v,'1')
import sys, json, pathlib
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT/'tools'))
import numpy as np
from run_ialns import eph
from ctoc14.constants import AU, DAY, T_MISSION
OUT = ROOT/'results/s17/physics'
E = eph()
tt = np.arange(0, T_MISSION, 0.5*DAY)
rE, vE = E.earth_state(tt)
lamE = np.arctan2(rE[:,1], rE[:,0]); aE = 1.0009175
ev = []
for X in range(1, 301):
    r, v = E.ast_states_at(np.full(len(tt), X-1), tt)
    rho = np.hypot(r[:,0], r[:,1])/AU; z = r[:,2]/AU
    d = np.hypot(rho - aE, z)
    k = np.where((d[1:-1] <= d[:-2]) & (d[1:-1] < d[2:]) & (d[1:-1] < 0.15))[0] + 1
    for j in k:
        lam = np.arctan2(r[j,1], r[j,0])
        L = np.degrees((lam - lamE[j] + np.pi) % (2*np.pi) - np.pi)
        # relative speed vs local circular velocity (prograde, in ecliptic)
        vc = 29.78*np.array([-np.sin(lam), np.cos(lam), 0.0])
        ev.append(dict(ast=X, t=float(tt[j]/DAY), L=float(L), d=float(d[j]), rho=float(rho[j]), z=float(z[j]), lam=float(np.degrees(lam)),
                       vrel=float(np.linalg.norm(v[j]-vc))))
json.dump(ev, open(OUT/'events.json','w'))
# flyby lags of s16a
rows = json.load(open(OUT/'fleet_flybys.json'))
from run_ialns import ipr
import glob
fl = []
for f in sorted(glob.glob(str(ROOT/'results/s16/best/fleet/route_r*.npz'))):
    st = dict(np.load(f)); name = pathlib.Path(f).stem[6:]
    ip = ipr(st); Yf, _ = ip.integrate()
    for j, X in enumerate(ip.asts):
        t = ip.tf[j]; re, _ = E.earth_state(t); r = Yf[j,:3]
        L = np.degrees((np.arctan2(r[1],r[0]) - np.arctan2(re[1],re[0]) + np.pi) % (2*np.pi) - np.pi)
        fl.append(dict(route=name, ast=int(X), t=float(t/DAY), L=float(L)))
json.dump(fl, open(OUT/'fleet_lags.json','w'))
print(len(ev), 'events')
