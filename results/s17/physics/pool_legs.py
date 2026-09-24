"""Leg dv and leg time vs the RARITY (node events < 0.05 AU in 15 yr) of the arrival target, over the 2490 honest pool routes."""
import json, numpy as np, pathlib
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); OUT = ROOT/'results/s17/physics'
ev = json.load(open(OUT/'events.json'))
nev = np.zeros(301, int)
for e in ev:
    if e['d'] < 0.05: nev[e['ast']] += 1
DAY = 86400.0
recs = []
for line in open(ROOT/'results/s16/honest/index.jsonl'):
    d = json.loads(line)
    if d.get('status') == 'rejected': continue
    fn = ROOT/'results/s16/honest'/f"route_{d['key']}.npz"
    if not fn.exists(): continue
    s = np.load(fn); tf = s['tf']; asts = s['asts']; ts = s['ts']; Ts = s['Ts']; o = np.argsort(tf)
    n = len(asts); tprev = float(s['tL']); dvn = np.linalg.norm(Ts, axis=1)
    nrare = int((nev[asts] <= 4).sum())
    for j in o:
        m = (ts > tprev) & (ts <= tf[j])
        recs.append((n, nrare, int(asts[j]), nev[asts[j]], (tf[j]-tprev)/DAY, dvn[m].sum()))
        tprev = tf[j]
R = np.array(recs, float)
np.save(OUT/'pool_legs.npy', R)
print('legs', len(R))
for lab, sel in (('deep >=38', R[:,0] >= 38), ('mid 30-37', (R[:,0] >= 30) & (R[:,0] < 38)), ('shallow <30', R[:,0] < 30)):
    print(lab)
    for lo, hi in ((0,3),(3,5),(5,7),(7,10),(10,99)):
        s = sel & (R[:,3] >= lo) & (R[:,3] < hi)
        print(f"   n_ev [{lo:2d},{hi:2d}): legs {s.sum():6d}  med leg_d {np.median(R[s,4]):5.0f}  med dv {np.median(R[s,5]):.3f}  mean dv {R[s,5].mean():.3f}")
# per route: fraction of rare targets vs depth
routes = {}
for r in recs: routes.setdefault((r[0], r[1]), 0)
import collections
byn = collections.defaultdict(list)
for line in open(ROOT/'results/s16/honest/index.jsonl'):
    d = json.loads(line)
    if d.get('status') == 'rejected': continue
    fn = ROOT/'results/s16/honest'/f"route_{d['key']}.npz"
    if not fn.exists(): continue
    s = np.load(fn); a = s['asts']; n = len(a); tank = d.get('tank_out') or d.get('tank_in')
    byn[min(n//5*5, 45)].append(((nev[a] <= 4).sum(), (tank-600)/n))
print('route depth band: mean # rare(<=4 ev) targets, max # rare, med kg/fb')
for k in sorted(byn):
    v = np.array(byn[k]); print(f"  {k:2d}-{k+4}: routes {len(v):4d}  rare mean {v[:,0].mean():5.1f}  max {v[:,0].max():3.0f}  kg/fb {np.median(v[:,1]):5.1f}")
print('total targets with <=4 events:', int(sum(1 for X in range(1,301) if X not in (131,144) and nev[X] <= 4)))
