import json, numpy as np, pathlib
OUT = pathlib.Path('/Users/mickey/solarsystem/ctoc14/results/s17/physics')
rows = json.load(open(OUT/'fleet_flybys.json'))
fb = [r for r in rows if r['k'] >= 0]; hd = {r['route']: r for r in rows if r['k'] == -1}
pop = np.load(OUT/'pool_popularity.npz'); cnt = pop['cnt']; cg = pop['cnt_good']; cd = pop['cnt_deep']
deep = {'r1','r2','r3','r4'}
def med(x): return np.median(x) if len(x) else np.nan
print('route  n  tank  kg/fb  tL | med: leg_d dv_leg vrel  r   |lat| ic  ec  ac | tgt a e i dnode | pop all/good/deep')
for R in sorted(hd):
    F = [r for r in fb if r['route']==R]
    g = lambda k: np.array([r[k] for r in F])
    ast = g('ast')
    print(f"{R:3s} {hd[R]['n']:3d} {hd[R]['tank']:6.1f} {(hd[R]['tank']-600)/hd[R]['n']:5.1f} {hd[R]['tL']:6.0f} | "
          f"{med(g('leg_d')):5.0f} {med(g('dv_leg')):.3f} {med(g('vrel')):5.1f} {med(g('r')):.2f} {med(abs(g('lat'))):4.1f} "
          f"{med(g('ic')):4.1f} {med(g('ec')):.2f} {med(g('ac')):.2f} | {med(g('a')):.2f} {med(g('e')):.2f} {med(g('i')):4.1f} {med(g('dnode')):.3f} | "
          f"{med(cnt[ast]):5.0f} {med(cg[ast]):4.0f} {med(cd[ast]):4.0f}  t_end {g('t_d').max():.0f}")
for grp, S in (('DEEP', deep), ('TAIL', {'r5','r6','r7','r8','r9'})):
    F = [r for r in fb if r['route'] in S]
    g = lambda k: np.array([r[k] for r in F])
    print(grp, len(F), 'quartiles')
    for k in ('leg_d','dv_leg','vrel','r','lat','ic','ec','ac','a','e','i','q','Q','dnode'):
        x = g(k); x = abs(x) if k=='lat' else x
        print(f"   {k:6s} {np.percentile(x,25):8.3f} {np.percentile(x,50):8.3f} {np.percentile(x,75):8.3f}")
    ast = g('ast')
    print('   pop  ', np.percentile(cnt[ast],[25,50,75]), 'good', np.percentile(cg[ast],[25,50,75]), 'deep', np.percentile(cd[ast],[25,50,75]))
    print('   dv_leg sum km/s', g('dv_leg').sum(), ' per fb', g('dv_leg').mean())
