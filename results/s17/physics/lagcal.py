"""Calibrate the (t, lag) phase-cost model on s16a: model dv = K * sum |slope changes| on the flyby (t, lag) sequence
(true-longitude lag, and mean-longitude lag from the craft's osculating elements) vs the twin dv."""
import os
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS'): os.environ.setdefault(v,'1')
import sys, json, glob, pathlib
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT/'tools'))
import numpy as np
from run_ialns import ipr, eph
from ctoc14.constants import MU, DAY
E = eph(); K = 0.0276  # km/s per deg/yr
def mean_lon(r, v):
    h = np.cross(r, v); rn = np.linalg.norm(r); a = 1/(2/rn - v@v/MU)
    ev = np.cross(v, h)/MU - r/rn; e = np.linalg.norm(ev)
    # true longitude (in ecliptic, small i) and equation of centre
    lt = np.arctan2(r[1], r[0]); pw = np.arctan2(ev[1], ev[0]); f = lt - pw
    Ean = 2*np.arctan2(np.sqrt(1-e)*np.sin(f/2), np.sqrt(1+e)*np.cos(f/2)); M = Ean - e*np.sin(Ean)
    return pw + M, lt
out = {}
for f in sorted(glob.glob(str(ROOT/'results/s16/best/fleet/route_r*.npz'))):
    nm = pathlib.Path(f).stem[6:]; ip = ipr(dict(np.load(f))); Yf, _ = ip.integrate(); o = np.argsort(ip.tf)
    rE0, vE0 = E.earth_state(ip.tL); lmE0, _ = mean_lon(rE0, vE0)
    rows = [(ip.tL/DAY, 0.0, 0.0)]
    for j in o:
        rE, vE = E.earth_state(ip.tf[j]); lmE, ltE = mean_lon(rE, vE)
        lm, lt = mean_lon(Yf[j,:3], Yf[j,3:6])
        rows.append((ip.tf[j]/DAY, np.degrees(lt-ltE), np.degrees(lm-lmE)))
    R = np.array(rows)
    res = {}
    for col, lab in ((1, 'true'), (2, 'mean')):
        L = np.degrees(np.unwrap(np.radians(R[:, col]))); t = R[:, 0]
        s = np.diff(L)/np.diff(t)*365.25
        res[lab] = K*np.abs(np.diff(s)).sum()
    out[nm] = dict(twin=ip.dv(), **res)
    print(f"{nm}: twin dv {ip.dv():6.2f} km/s | model on TRUE-lon lag {res['true']:7.2f} | on MEAN-lon lag {res['mean']:7.2f}")
json.dump(out, open(ROOT/'results/s17/physics/lagcal.json','w'), indent=1)
