"""s17 redteam: target anatomy of the s16 best fleet + honest pool statistics.

Writes results/s17/redteam/targets.csv (one row per reachable target) and anatomy.json (summaries).
Per target: owner route in s16 best, pool frequency (all / deep n>=37 / tail-like n<=31), max pool depth containing it,
elements a e i Omega omega, node radii, |r_node - 1| min, Earth MOID (numeric), period, #node passages near 1 AU in the
window, flyby geometry in the fleet (craft r, ecliptic z, relative speed, epoch), and the fleet fuel law fit.
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, pathlib, collections
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)
import numpy as np
import run_ialns as RI
from ctoc14.constants import DAY, AU, T_MISSION, EARTH_ELEMENTS

OUT = ROOT / 'results/s17/redteam'
E = RI.eph()

# ---------------------------------------------------------------- fleet
F = {}
for f in sorted((ROOT / 'results/s16/best/fleet').glob('route_*.npz')):
    z = np.load(f); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL']); F[f.stem[6:]] = st
DEEP = ['r1', 'r2', 'r3', 'r4']; TAIL = ['r5', 'r6', 'r7', 'r8', 'r9']
owner = {}; tfly = {}
geo = {}
for k, st in F.items():
    ip = RI.ipr(st); Yf, _ = ip.integrate()
    for j, a in enumerate(ip.asts):
        owner[int(a)] = k; tfly[int(a)] = float(ip.tf[j]) / DAY
        ra, va = E.ast_states_at(np.array([int(a) - 1]), np.array([ip.tf[j]]))
        vr = float(np.linalg.norm(Yf[j, 3:6] - va[0]))
        geo[int(a)] = dict(r=float(np.linalg.norm(Yf[j, :3]) / AU), z=float(Yf[j, 2] / AU), vrel=vr)

# ---------------------------------------------------------------- honest pool
pool = [json.loads(l) for l in open(ROOT / 'results/s16/honest/index.jsonl')]
pool = [r for r in pool if r['status'] != 'rejected']
cnt_all = collections.Counter(); cnt_deep = collections.Counter(); cnt_short = collections.Counter(); maxdep = {}
prow = []
for r in pool:
    z = np.load(ROOT / f'results/s16/honest/route_{r["key"]}.npz')
    asts = [int(a) for a in z['asts']]; n = len(asts)
    tank = float(r['tank_out'])
    prow.append(dict(key=r['key'], n=n, tank=tank, asts=asts, tL=float(z['tL']) / DAY))
    for a in asts:
        cnt_all[a] += 1
        if n >= 37: cnt_deep[a] += 1
        if n <= 31: cnt_short[a] += 1
        maxdep[a] = max(maxdep.get(a, 0), n)
ndeep = sum(1 for p in prow if p['n'] >= 37)

# ---------------------------------------------------------------- elements and MOID
ids, el = E.ids, E.elem   # a e i Om om M
def rot(Om, i, om):
    cO, sO, ci, si, co, so = np.cos(Om), np.sin(Om), np.cos(i), np.sin(i), np.cos(om), np.sin(om)
    return np.array([[cO * co - sO * so * ci, -cO * so - sO * co * ci, sO * si],
                     [sO * co + cO * so * ci, -sO * so + cO * co * ci, -cO * si],
                     [so * si, co * si, ci]])
def orbit_pts(a, e, i, Om, om, n=720):
    nu = np.linspace(0, 2 * np.pi, n, endpoint=False)
    p = a * (1 - e * e); r = p / (1 + e * np.cos(nu))
    P = np.stack([r * np.cos(nu), r * np.sin(nu), 0 * nu], 1)
    return P @ rot(np.radians(Om), np.radians(i), np.radians(om)).T
Ea, Ee, Ei, EOm, Eom, _ = EARTH_ELEMENTS
EP = orbit_pts(Ea, Ee, Ei, EOm, Eom, 1440)
rows = []
for k in range(len(ids)):
    a, e, i, Om, om, M = el[k]
    X = int(ids[k])
    P = orbit_pts(a, e, i, Om, om, 1440)
    d = np.linalg.norm(P[:, None, :] - EP[None, ::4, :], axis=2)
    moid = float(d.min())
    p = a * (1 - e * e)
    rasc = p / (1 + e * np.cos(np.radians(om))); rdes = p / (1 - e * np.cos(np.radians(om)))
    dn = min(abs(rasc - 1.0), abs(rdes - 1.0))
    per = a ** 1.5
    syn = 1.0 / abs(1.0 - 1.0 / per) if abs(per - 1) > 1e-6 else np.inf
    # heliocentric longitude of the near-1-AU node
    lnode = Om if abs(rasc - 1) <= abs(rdes - 1) else (Om + 180.0) % 360
    rows.append(dict(ast=X, a=a, e=e, i=i, Om=Om, om=om, q=a * (1 - e), Q=a * (1 + e), rasc=rasc, rdes=rdes, dnode=dn,
                     moid=moid, P=per, syn=syn, lnode=lnode, nodepass=15.0 / per,
                     owner=owner.get(X, 'miss'), cls=('deep' if owner.get(X) in DEEP else 'tail' if owner.get(X) in TAIL else 'miss'),
                     t_fly=tfly.get(X, np.nan), pool_all=cnt_all[X], pool_deep=cnt_deep[X], pool_short=cnt_short[X],
                     maxdep=maxdep.get(X, 0), **{f'fly_{kk}': v for kk, v in geo.get(X, dict(r=np.nan, z=np.nan, vrel=np.nan)).items()}))

import csv
with open(OUT / 'targets.csv', 'w', newline='') as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

# ---------------------------------------------------------------- summaries
def med(xs):
    xs = [x for x in xs if np.isfinite(x)]
    return (round(float(np.median(xs)), 3), round(float(np.percentile(xs, 25)), 3), round(float(np.percentile(xs, 75)), 3)) if xs else None
S = {}
keys = ['a', 'e', 'i', 'q', 'Q', 'dnode', 'moid', 'P', 'syn', 'fly_r', 'fly_z', 'fly_vrel', 'pool_all', 'pool_deep', 'maxdep']
for c in ['deep', 'tail', 'miss']:
    rr = [r for r in rows if r['cls'] == c]
    S[c] = dict(n=len(rr), **{k: med([r[k] for r in rr]) for k in keys})
for k in TAIL + DEEP:
    rr = [r for r in rows if r['owner'] == k]
    S['route_' + k] = dict(n=len(rr), **{kk: med([r[kk] for r in rr]) for kk in ['i', 'e', 'dnode', 'fly_vrel', 'pool_deep', 'maxdep']})
# "unwanted": flown by tails and in NO deep pool route (or <= 1% of deep routes)
unw = [r for r in rows if r['cls'] == 'tail' and r['pool_deep'] <= max(1, 0.01 * ndeep)]
S['unwanted_def'] = f'tail-owned and in <= max(1, 1%) of the {ndeep} deep (n>=37) pool routes'
S['unwanted'] = dict(n=len(unw), asts=sorted(r['ast'] for r in unw), **{k: med([r[k] for r in unw]) for k in keys})
wanted_tail = [r for r in rows if r['cls'] == 'tail' and r['pool_deep'] > max(1, 0.01 * ndeep)]
S['tail_wanted'] = dict(n=len(wanted_tail), **{k: med([r[k] for r in wanted_tail]) for k in keys})
# fuel law over the honest pool: fuel = a + b n (twin tanks) by launch-day class
Pn = np.array([p['n'] for p in prow]); Pf = np.array([p['tank'] - 600 for p in prow])
m = Pn >= 10
A = np.vstack([np.ones(m.sum()), Pn[m]]).T; coef = np.linalg.lstsq(A, Pf[m], rcond=None)[0]
S['pool_fuel_fit_n>=10'] = dict(a_kg=round(float(coef[0]), 1), b_kg_per_fb=round(float(coef[1]), 2), npts=int(m.sum()))
# lower envelope: min fuel per depth bin
env = {}
for n in range(10, 49):
    s = Pf[Pn == n]
    if len(s): env[n] = dict(min=round(float(s.min()), 1), p10=round(float(np.percentile(s, 10)), 1), cnt=int(len(s)))
S['pool_fuel_envelope'] = env
# tail-union pool routes: routes whose targets are >= 80 % tail-owned
tailset = set(r['ast'] for r in rows if r['cls'] == 'tail')
tr = []
for p in prow:
    if p['n'] >= 15:
        fr = sum(1 for a in p['asts'] if a in tailset) / p['n']
        if fr >= 0.8:
            tr.append(dict(key=p['key'], n=p['n'], fuel=round(p['tank'] - 600, 1), kgfb=round((p['tank'] - 600) / p['n'], 2), frac=round(fr, 2)))
tr.sort(key=lambda x: x['kgfb'])
S['tail_union_routes'] = dict(count=len(tr), best=tr[:12])
# fleet anatomy
fa = {}
for k, st in F.items():
    tank = float(RI.ipr(st).tank()); n = len(st['asts'])
    fa[k] = dict(n=n, tank=round(tank, 1), fuel=round(tank - 600, 1), kgfb=round((tank - 600) / n, 2), tL=round(st['tL'] / DAY),
                 t_first=round(float(np.min(st['tf'])) / DAY), t_last=round(float(np.max(st['tf'])) / DAY))
S['fleet'] = fa
S['n_pool'] = len(prow); S['n_pool_deep'] = ndeep
json.dump(S, open(OUT / 'anatomy.json', 'w'), indent=1)
print(json.dumps({k: S[k] for k in ['deep', 'tail', 'unwanted', 'tail_wanted', 'pool_fuel_fit_n>=10', 'fleet', 'n_pool', 'n_pool_deep']}, indent=0))
print('tail-union routes', S['tail_union_routes'])
