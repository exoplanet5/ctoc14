import json, numpy as np, pathlib
from scipy.stats import spearmanr
OUT = pathlib.Path('/Users/mickey/solarsystem/ctoc14/results/s17/physics')
ev = json.load(open(OUT/'events.json')); fl = json.load(open(OUT/'fleet_lags.json')); rows = json.load(open(OUT/'fleet_flybys.json'))
pop = np.load(OUT/'pool_popularity.npz'); cd = pop['cnt_deep']; cg = pop['cnt_good']
g = np.load(OUT/'targets_geom.npz')
route_of = {r['ast']: r['route'] for r in fl}
tg = [X for X in range(1,301) if X not in (131,144)]
nev = {X: sum(1 for e in ev if e['ast']==X and e['d']<0.05) for X in tg}
nev10 = {X: sum(1 for e in ev if e['ast']==X and e['d']<0.10) for X in tg}
x = np.array([nev[X] for X in tg]); y = np.array([cd[X] for X in tg])
print('spearman(n_events<0.05, deep-route count)', spearmanr(x, y))
print('spearman(a, deep count)', spearmanr(g['a'][np.array(tg)-1], y))
print('spearman(i, deep count)', spearmanr(g['i'][np.array(tg)-1], y))
print('spearman(dnode, deep count)', spearmanr(g['dnode'][np.array(tg)-1], y))
print('n_events<0.05 histogram:', np.bincount(np.minimum(x, 20)))
tail = {'r5','r6','r7','r8','r9'}
for lo, hi in ((0,3),(3,5),(5,7),(7,10),(10,15),(15,99)):
    S = [X for X in tg if lo <= nev[X] < hi]
    ft = np.mean([route_of[X] in tail for X in S]) if S else np.nan
    # leg duration / dv into these targets
    L = [r for r in rows if r['k']>=0 and lo <= nev.get(r['ast'],0) < hi]
    print(f"n_ev [{lo},{hi}): {len(S):3d} targets  frac in tails {ft:.2f}  med leg_d {np.median([r['leg_d'] for r in L]):5.0f}  med dv_leg {np.median([r['dv_leg'] for r in L]):.3f}  med a {np.median(g['a'][np.array(S)-1]):.2f}")
# second node
el = np.loadtxt('/Users/mickey/solarsystem/ctoc14/MEA.txt', comments='#')
ra, rd = g['r_asc'], g['r_dsc']
far = np.where(np.abs(ra-1) < np.abs(rd-1), rd, ra)   # the node NOT near 1 AU
sp = json.load(open(OUT/'popularity_split.json'))
for nm in ('low','high'):
    S = np.array(sp[nm]) - 1
    f = far[S]
    print(nm, 'other-node radius quartiles', np.percentile(f,[10,25,50,75,90]).round(2), ' frac other node in [0.75,1.30]:', np.mean((f>0.75)&(f<1.30)).round(2),
          ' i med', np.median(g['i'][S]).round(1))
