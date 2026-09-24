import json, numpy as np, pathlib
OUT = pathlib.Path('/Users/mickey/solarsystem/ctoc14/results/s17/physics')
fl = json.load(open(OUT/'fleet_lags.json')); ev = json.load(open(OUT/'events.json'))
pop = np.load(OUT/'pool_popularity.npz'); cd = pop['cnt_deep']; cg = pop['cnt_good']; ca = pop['cnt']
route_of = {r['ast']: r['route'] for r in fl}
print('per-route lag curve (flyby lag deg vs t):')
for R in sorted(set(r['route'] for r in fl)):
    F = sorted([r for r in fl if r['route']==R], key=lambda r: r['t'])
    L = np.unwrap(np.radians([r['L'] for r in F])); L = np.degrees(L); t = np.array([r['t'] for r in F])
    dL = np.diff(L); dt = np.diff(t)
    print(f"{R}: n{len(F)} t {t[0]:5.0f}-{t[-1]:5.0f} L start {L[0]:7.1f} end {L[-1]:7.1f} range {L.min():7.1f}..{L.max():7.1f}  "
          f"total|dL| {np.abs(dL).sum():6.0f} deg  mean slope {(L[-1]-L[0])/(t[-1]-t[0])*365.25:6.1f} deg/yr")
    print('    L:', ' '.join(f"{x:.0f}" for x in L))
# per target event stats
E = {}
for e in ev: E.setdefault(e['ast'], []).append(e)
print()
def stats(X, dmax):
    es = [e for e in E.get(X, []) if e['d'] < dmax]
    if not es: return 0, np.nan
    L = np.radians([e['L'] for e in es]); R = np.abs(np.mean(np.exp(1j*L)))  # circular concentration
    return len(es), R
deep = {'r1','r2','r3','r4'}
for dmax in (0.05, 0.1):
    print(f'dmax {dmax}: group  n_events(med)  lag concentration R (med)   frac R>0.8')
    for grp, S in (('deep', deep), ('tail', {'r5','r6','r7','r8','r9'})):
        xs = [X for X in route_of if route_of[X] in S]
        st = np.array([stats(X, dmax) for X in xs])
        print(f'   {grp}: {len(xs)} tg  events med {np.nanmedian(st[:,0]):.1f} (IQR {np.nanpercentile(st[:,0],25):.0f}-{np.nanpercentile(st[:,0],75):.0f})  R med {np.nanmedian(st[:,1]):.2f}  R>0.8 {np.mean(st[:,1]>0.8):.2f}  zero-ev {np.mean(st[:,0]==0):.2f}')
# popularity quartiles by deep-route count
allX = [X for X in range(1,301) if X not in (131,144)]
cdv = np.array([cd[X] for X in allX]); order = np.argsort(cdv)
low = [allX[i] for i in order[:99]]; high = [allX[i] for i in order[-99:]]
el = np.loadtxt('/Users/mickey/solarsystem/ctoc14/MEA.txt', comments='#')
for nm, S in (('least-deep-wanted 99', low), ('most-deep-wanted 99', high)):
    st = np.array([stats(X, 0.1) for X in S]); st5 = np.array([stats(X, 0.05) for X in S])
    a = el[np.array(S)-1,1]; e = el[np.array(S)-1,2]; i = el[np.array(S)-1,3]
    tails = np.mean([route_of.get(X,'') in {'r5','r6','r7','r8','r9'} for X in S])
    vr = [np.median([x['vrel'] for x in E[X]]) if X in E else np.nan for X in S]
    print(f"{nm}: cnt_deep med {np.median([cd[X] for X in S]):.0f}; ev<0.1 med {np.nanmedian(st[:,0]):.0f}, ev<0.05 med {np.nanmedian(st5[:,0]):.0f}; R med {np.nanmedian(st[:,1]):.2f}; "
          f"a {np.median(a):.2f} e {np.median(e):.2f} i {np.median(i):.1f}; vrel {np.nanmedian(vr):.1f}; frac in s16a tails {tails:.2f}")
json.dump(dict(low=low, high=high), open(OUT/'popularity_split.json','w'))
