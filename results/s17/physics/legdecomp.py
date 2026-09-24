"""(1) per-leg element-change decomposition of s16a twin dv: period (lag), in-plane shape (e-vector), plane (i-vector);
(2) track density: approach minima < 0.035/0.06 AU of each route's twin track to ALL targets, split own / other route."""
import os
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'): os.environ.setdefault(v,'1')
import sys, json, glob, pathlib
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT/'tools'))
import numpy as np
from run_ialns import ipr, eph, trajectory
from ctoc14.constants import AU, DAY, MU
OUT = ROOT/'results/s17/physics'
E = eph()
def elvec(r, v):
    h = np.cross(r, v); rn = np.linalg.norm(r)
    a = 1/(2/rn - v@v/MU)
    ev = np.cross(v, h)/MU - r/rn
    hn = h/np.linalg.norm(h)
    ivec = np.array([hn[0], hn[1]])  # ~ i*(sinO, -cosO): plane tilt vector (rad, small i)
    return a, ev[:2], ivec
own = {}; routes = {}
for f in sorted(glob.glob(str(ROOT/'results/s16/best/fleet/route_r*.npz'))):
    nm = pathlib.Path(f).stem[6:]; st = dict(np.load(f)); routes[nm] = st
    for X in st['asts']: own[int(X)] = nm
res = {}
print('route  sum_dv  | per-leg shares (km/s): period  e-vec  plane | sum of three / actual')
for nm, st in routes.items():
    ip = ipr(st); Yf, _ = ip.integrate(); o = np.argsort(ip.tf)
    # launch state
    X0 = (ip.rE, ip.vE + ip.vinf)
    states = [X0] + [(Yf[j,:3], Yf[j,3:6]) for j in o]
    dP = dE = dI = 0.0; tot = 0.0; legs = []
    tprev = ip.tL
    for k in range(1, len(states)):
        a0, e0, i0 = elvec(*states[k-1]); a1, e1, i1 = elvec(*states[k])
        v = 29.78
        per = v/2*abs(a1-a0)/((a0+a1)/2); ecc = v/2*np.linalg.norm(e1-e0); pl = v*np.linalg.norm(i1-i0)
        t = ip.tf[o[k-1]]; m = (ip.ts > tprev) & (ip.ts <= t); act = float(np.linalg.norm(ip.Ts[m],axis=1).sum()); tprev = t
        legs.append((per, ecc, pl, act)); dP += per; dE += ecc; dI += pl; tot += act
    legs = np.array(legs)
    res[nm] = dict(period=dP, ecc=dE, plane=dI, actual=tot, dv=ip.dv())
    print(f"{nm}: {ip.dv():6.2f} | {dP:6.2f} {dE:6.2f} {dI:6.2f} | {dP+dE+dI:6.2f} / {tot:6.2f}   per fb: {dP/len(o):.3f} {dE/len(o):.3f} {dI/len(o):.3f}")
# density
tgts = [X for X in range(1,301) if X not in (131,144)]
print('\nroute: approach minima of the twin track (whole mission) to targets: d<0.035 own/other/total | d<0.06 own/other/total | own targets never within 0.035 of other routes')
dens = {}
tracks = {}
for nm, st in routes.items():
    ip = ipr(st); tt, r = trajectory(ip); tracks[nm] = (tt, r)
    dmin = {}
    for X in tgts:
        ra, _ = E.ast_states_at(np.full(len(tt), X-1), tt)
        d = np.linalg.norm(r - ra, axis=1)/AU
        dmin[X] = float(d.min())
    dens[nm] = dmin
json.dump(dens, open(OUT/'track_dmin.json','w'))
for nm in routes:
    dm = dens[nm]
    c = lambda th, who: sum(1 for X in tgts if dm[X] < th and (own[X]==nm) == (who=='own'))
    print(f"{nm}: {c(0.035,'own'):3d} {c(0.035,'oth'):3d} {c(0.035,'own')+c(0.035,'oth'):3d} | {c(0.06,'own'):3d} {c(0.06,'oth'):3d} {c(0.06,'own')+c(0.06,'oth'):3d}")
json.dump(res, open(OUT/'leg_decomp.json','w'))
