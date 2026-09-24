"""Stage-17 PHYSICS diagnosis: per-flyby geometry of the s16a twin fleet + pool popularity of every target."""
import os
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'): os.environ.setdefault(v,'1')
import sys, json, glob, pathlib
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT/'tools'))
import numpy as np
from run_ialns import ipr, eph
from ctoc14.constants import AU, DAY, MU, VE, M_DRY
OUT = ROOT/'results/s17/physics'

def rv2el(r, v):
    r = np.asarray(r); v = np.asarray(v)
    h = np.cross(r, v); rn = np.linalg.norm(r); vn2 = v @ v
    a = 1.0/(2/rn - vn2/MU)
    evec = np.cross(v, h)/MU - r/rn; e = np.linalg.norm(evec)
    i = np.degrees(np.arccos(h[2]/np.linalg.norm(h)))
    return a/AU, e, i

E = eph()
el = E.elem  # a e i Om om M
A, Ee, I, Om, om = el[:,0], el[:,1], np.radians(el[:,2]), np.radians(el[:,3]), np.radians(el[:,4])
p = A*(1-Ee**2)
r_asc = p/(1+Ee*np.cos(om)); r_dsc = p/(1-Ee*np.cos(om))   # AU, node radii
q = A*(1-Ee); Q = A*(1+Ee)
P_yr = A**1.5
# closest node to 1 AU
dnode = np.minimum(abs(r_asc-1), abs(r_dsc-1))
tgt = dict(a=A, e=Ee, i=np.degrees(I), q=q, Q=Q, r_asc=r_asc, r_dsc=r_dsc, dnode=dnode, P=P_yr)
np.savez(OUT/'targets_geom.npz', **tgt)

fleet = sorted(glob.glob(str(ROOT/'results/s16/best/fleet/route_r*.npz')))
rows = []
for f in fleet:
    name = pathlib.Path(f).stem.replace('route_','')
    st = dict(np.load(f, allow_pickle=True)); ip = ipr(st)
    Yf, Ys = ip.integrate()
    ra, va = E.ast_states_at(np.array(ip.asts)-1, ip.tf)
    o = np.argsort(ip.tf)
    tprev = ip.tL
    for jj, j in enumerate(o):
        X = ip.asts[j]; t = ip.tf[j]
        # impulses in (tprev, t]
        m = (ip.ts > tprev) & (ip.ts <= t)
        dvleg = float(np.linalg.norm(ip.Ts[m], axis=1).sum())
        ac, ec, ic = rv2el(Yf[j,:3], Yf[j,3:6])
        r = Yf[j,:3]; vrel = np.linalg.norm(Yf[j,3:6]-va[j])
        rows.append(dict(route=name, k=jj, ast=int(X), t_d=float(t/DAY), leg_d=float((t-tprev)/DAY), dv_leg=dvleg,
                         r=float(np.linalg.norm(r)/AU), z=float(r[2]/AU), lat=float(np.degrees(np.arcsin(r[2]/np.linalg.norm(r)))),
                         vrel=float(vrel), ac=ac, ec=ec, ic=ic,
                         a=float(A[X-1]), e=float(Ee[X-1]), i=float(np.degrees(I[X-1])), q=float(q[X-1]), Q=float(Q[X-1]),
                         dnode=float(dnode[X-1])))
        tprev = t
    # trailing impulses after last flyby
    tank = ip.tank(); rows.append(dict(route=name, k=-1, tank=float(tank), n=len(ip.asts), tL=float(ip.tL/DAY), dv=float(ip.dv())))
json.dump(rows, open(OUT/'fleet_flybys.json','w'))

# pool popularity from the honest pool
cnt = np.zeros(301); cnt_good = np.zeros(301); cnt_deep = np.zeros(301)
nr = 0
for line in open(ROOT/'results/s16/honest/index.jsonl'):
    d = json.loads(line)
    if d.get('status') == 'rejected': continue
    fn = ROOT/'results/s16/honest'/f"route_{d['key']}.npz"
    if not fn.exists(): continue
    s = np.load(fn); asts = s['asts']; n = len(asts); tank = d.get('tank_out') or d.get('tank_in')
    kgfb = (tank-600)/n
    nr += 1
    cnt[asts] += 1
    if kgfb < 9.5: cnt_good[asts] += 1
    if n >= 36: cnt_deep[asts] += 1
np.savez(OUT/'pool_popularity.npz', cnt=cnt, cnt_good=cnt_good, cnt_deep=cnt_deep, nr=nr)
print('pool routes', nr)
