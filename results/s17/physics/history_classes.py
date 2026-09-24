"""Re-read older fleets in the pin / medium / filler classes (node events < 0.05 AU in 15 yr: <=4 / 5-9 / >=10)."""
import json, numpy as np, pathlib
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); OUT = ROOT/'results/s17/physics'
ev = json.load(open(OUT/'events.json')); nev = np.zeros(301, int)
for e in ev:
    if e['d'] < 0.05: nev[e['ast']] += 1
cls = lambda X: 'pin' if nev[X] <= 4 else ('fill' if nev[X] >= 10 else 'med')
C = lambda L: {c: sum(cls(x)==c for x in L) for c in ('pin','med','fill')}
fl = json.load(open(OUT/'fleet_lags.json'))
print('s16a routes (results/s16/best/fleet):')
for R in sorted(set(r['route'] for r in fl)):
    a = [r['ast'] for r in fl if r['route']==R]; print('  ', R, len(a), C(a))
print('stage-5 fill2 (results/n8/fill2), in fill order:')
for r in ['05','09','13','01','03','06','10','12']:
    a = [int(x) for x in np.load(ROOT/f'results/n8/fill2/route_{r}.npz')['asts']]; print('  ', r, len(a), C(a))
left5 = [2,5,12,17,31,32,40,41,53,57,61,64,69,75,93,94,100,104,107,110,120,124,129,130,133,134,138,140,150,152,156,159,165,167,171,172,175,177,182,198,204,209,215,223,225,234,250,253,263,266,274,275,277,283,286,287,288]
print('  fill2 leftovers', len(left5), C(left5))
left8 = [1,3,6,10,11,12,13,16,18,22,23,25,27,31,34,38,40,50,66,69,75,77,82,84,88,92,95,101,108,110,114,116,119,123,124,127,132,138,149,152,154,156,157,167,168,173,180,182,184,185,189,191,192,216,217,221,223,224,228,235,236,237,246,248,250,251,253,254,263,266,274,276,278,282,287,295]
print('stage-8 J9grow leftovers (results/s8/J9grow.out)', len(left8), C(left8))
iso = [64,119,137,138,157,158,168,216]
print('stage-14 isolated 8: n_events', {x: int(nev[x]) for x in iso})
TG = [X for X in range(1,301) if X not in (131,144)]; print('catalogue', C(TG))
