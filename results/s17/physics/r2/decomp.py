"""R2 aux: decompose the closest approach of each survivor track to r9's hard-core targets into along-track / radial /
vertical offsets (local frame of the target at the approach epoch). Read-only on tools."""
import os
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'): os.environ.setdefault(v,'1')
import sys, json, pathlib
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT/'tools'))
os.chdir(ROOT)
import numpy as np
import run_ialns as RI
from ctoc14.constants import AU, DAY
E = RI.eph()
hard = [17,36,75,108,118,180,212,216,219,220,243,247,258,288]
out = []
for i in range(1, 9):
    st = dict(np.load(f'results/s16/best/fleet/route_r{i}.npz')); ip = RI.ipr(st)
    tt, r = RI.trajectory(ip, step=1*DAY)
    for X in hard:
        ra, va = E.ast_states_at(np.full(len(tt), X-1), tt)
        d = np.linalg.norm(r - ra, axis=1)/AU
        k = np.where((d[1:-1] <= d[:-2]) & (d[1:-1] < d[2:]) & (d[1:-1] < 0.25))[0] + 1
        for j in k:
            rr = ra[j]; uR = rr/np.linalg.norm(rr); h = np.cross(rr, va[j]); uN = h/np.linalg.norm(h); uT = np.cross(uN, uR)
            # frame of the ECLIPTIC circular motion at the target position
            ez = np.array([0,0,1.0]); uT2 = np.cross(ez, uR); uT2 /= np.linalg.norm(uT2); uR2 = np.cross(uT2, ez)
            dd = (r[j] - rr)/AU
            out.append(dict(host=f'r{i}', ast=X, t=float(tt[j]/DAY), dist=float(d[j]), d_rad=float(dd@uR2), d_along=float(dd@uT2), d_z=float(dd@ez),
                            r_craft=float(np.linalg.norm(r[j])/AU), r_ast=float(np.linalg.norm(rr)/AU)))
json.dump(out, open(ROOT/'results/s17/physics/r2/decomp.json','w'), indent=1)
best = {}
for o in out:
    if o['ast'] not in best or o['dist'] < best[o['ast']]['dist']: best[o['ast']] = o
print('best approach per hard target: dist | along radial vertical (AU)')
for X in hard:
    o = best.get(X)
    if o: print(f"{X:4d} {o['host']} t{o['t']:5.0f} {o['dist']:.3f} | {o['d_along']:+.3f} {o['d_rad']:+.3f} {o['d_z']:+.3f}")
    else: print(X, 'none < 0.25')
A = np.array([[abs(o['d_along']), abs(o['d_rad']), abs(o['d_z'])] for o in out if o['dist'] < 0.15])
print('approaches < 0.15 AU:', len(A), ' median |along| |rad| |z|', np.median(A, 0).round(3))
