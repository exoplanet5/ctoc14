"""Fillers (>=10 node events) vs depth and cost in the honest pool; which classes carry the overlap between deep routes."""
import json, numpy as np, pathlib, itertools
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); OUT = ROOT/'results/s17/physics'
ev = json.load(open(OUT/'events.json'))
nev = np.zeros(301, int)
for e in ev:
    if e['d'] < 0.05: nev[e['ast']] += 1
cls = lambda X: 'rare' if nev[X] <= 4 else ('fill' if nev[X] >= 10 else 'med')
TG = [X for X in range(1,301) if X not in (131,144)]
print('class sizes', {c: sum(1 for X in TG if cls(X)==c) for c in ('rare','med','fill')})
R = []; seen = set()
for line in open(ROOT/'results/s16/honest/index.jsonl'):
    d = json.loads(line)
    if d.get('status') == 'rejected': continue
    fn = ROOT/'results/s16/honest'/f"route_{d['key']}.npz"
    if not fn.exists(): continue
    a = tuple(sorted(int(x) for x in np.load(fn)['asts']))
    if a in seen: continue
    seen.add(a); tank = d.get('tank_out') or d.get('tank_in')
    R.append((a, tank))
D = [(a, t) for a, t in R if len(a) >= 36]
print('deep distinct', len(D))
nf = np.array([sum(cls(X)=='fill' for X in a) for a, _ in D]); kg = np.array([(t-600)/len(a) for a, t in D]); dep = np.array([len(a) for a, _ in D])
for lo, hi in ((0,8),(8,11),(11,14),(14,17),(17,20),(20,40)):
    s = (nf >= lo) & (nf < hi)
    if s.sum(): print(f"  fillers {lo:2d}-{hi-1:2d}: routes {s.sum():4d}  med depth {np.median(dep[s]):4.0f}  med kg/fb {np.median(kg[s]):5.2f}  best {kg[s].min():5.2f}")
print('corr(fillers, kg/fb) = %.2f; corr(fillers, depth) = %.2f' % (np.corrcoef(nf, kg)[0,1], np.corrcoef(nf, dep)[0,1]))
# overlap composition over random pairs of deep routes
rng = np.random.default_rng(0); comp = {'rare':0,'med':0,'fill':0}; npair = 0
for _ in range(20000):
    i, j = rng.choice(len(D), 2, replace=False)
    sh = set(D[i][0]) & set(D[j][0])
    for X in sh: comp[cls(X)] += 1
    npair += 1
tot = sum(comp.values())
base = {c: sum(1 for X in TG if cls(X)==c)/len(TG) for c in comp}
print('shared targets per deep pair: %.1f; composition' % (tot/npair), {c: round(comp[c]/tot,2) for c in comp}, ' base', {c: round(base[c],2) for c in base})
# demand per class
dem = {c: [] for c in comp}
cntd = np.zeros(301)
for a, _ in D:
    for X in a: cntd[X] += 1
for X in TG: dem[cls(X)].append(cntd[X])
print('deep-route demand per target (mean # of deep routes carrying it):', {c: round(np.mean(v),1) for c, v in dem.items()})
