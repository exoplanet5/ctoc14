import json, numpy as np, pathlib
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); OUT = ROOT/'results/s17/physics'
ev = json.load(open(OUT/'events.json'))
nev = np.zeros(301, int)
for e in ev:
    if e['d'] < 0.05: nev[e['ast']] += 1
tg = [X for X in range(1,301) if X not in (131,144)]
rare = [X for X in tg if nev[X] <= 4]
pop = np.load(OUT/'pool_popularity.npz'); cd = pop['cnt_deep']
print('rare targets:', len(rare), ' cnt_deep (routes >=36 fb carrying it) quartiles', np.percentile(cd[rare],[0,10,25,50,75,90]),
      ' # with 0:', int((cd[rare]==0).sum()), ' <5:', int((cd[rare]<5).sum()))
R = []
for line in open(ROOT/'results/s16/honest/index.jsonl'):
    d = json.loads(line)
    if d.get('status') == 'rejected': continue
    fn = ROOT/'results/s16/honest'/f"route_{d['key']}.npz"
    if not fn.exists(): continue
    s = np.load(fn); a = s['asts']; n = len(a); tank = d.get('tank_out') or d.get('tank_in')
    tf = np.sort(s['tf']); span = (tf[-1] - float(s['tL']))/86400
    R.append((n, int((nev[a] <= 4).sum()), (tank-600)/n, tank, span/n))
R = np.array(R)
D = R[R[:,0] >= 36]
print('deep routes >=36:', len(D))
for k in range(0, 16, 2):
    s = (D[:,1] >= k) & (D[:,1] < k+2)
    if s.sum(): print(f"  rare {k:2d}-{k+1:2d}: routes {s.sum():4d}  med depth {np.median(D[s,0]):4.0f}  med kg/fb {np.median(D[s,2]):5.2f}  best kg/fb {D[s,2].min():5.2f}  med tank {np.median(D[s,3]):6.0f}  pace d/fb {np.median(D[s,4]):4.0f}")
c = np.polyfit(D[:,1], D[:,2], 1); print('kg/fb ~ %.3f * n_rare + %.2f' % tuple(c))
