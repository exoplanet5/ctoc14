import os
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'): os.environ.setdefault(v,'1')
import sys, multiprocessing as mp; sys.path.insert(0,'.'); sys.path.insert(0,'tools')
import numpy as np, run_ialns as RI
from ctoc14.constants import DAY
fl=RI.IFleet('results/s14/probe_cover/N9/fleet')
ALL=[t for t in range(1,301) if t not in (131,144)]
cov=set(fl.coverage()); miss=[t for t in ALL if t not in cov]
sj=sum(RI.cost(r['tank']) for r in fl.routes.values())
print(f'N9 fleet: {len(fl.routes)} craft, {len(cov)} covered, sum J_i {sj:.4f}; missing {len(miss)}: {miss}')
with mp.get_context('fork').Pool(2) as pl:
    C=[c for lst in pl.map(RI.w_cands,[(n,r['st'],miss,0.30) for n,r in fl.routes.items()]) for c in lst]
best={}
for c in C:
    st=fl.routes[c['host']]['st']
    if np.abs(np.asarray(st['tf'])-c['t']).min()<=10*DAY: continue
    best.setdefault(c['ast'],[]).append(c)
jobs=[]
for t in miss:
    cs=sorted(best.get(t,[]),key=lambda c:c['dist'])[:3]
    print(f'  target {t:3d}: {len(best.get(t,[]))} approaches; nearest ' + ', '.join(f"{c['host']}@{c['dist']:.3f}AU" for c in cs))
    for c in cs: jobs.append((c['host'],fl.routes[c['host']]['st'],c['ast'],c['t']))
with mp.get_context('fork').Pool(2) as pl: R=pl.map(RI.w_insert,jobs,chunksize=1)
bt={}
for r in R:
    if not r.get('ok'): continue
    dJ=RI.cost(r['tank'])-RI.cost(fl.routes[r['name']]['tank'])
    if r['ast'] not in bt or dJ<bt[r['ast']][0]: bt[r['ast']]=(dJ,r['name'],r['tank']-fl.routes[r['name']]['tank'])
print('\ncheapest single insertion per missing target (each priced alone):')
tot=0; n=0
for t in miss:
    if t in bt: dJ,h,dk=bt[t]; print(f'  {t:3d}: +{dJ:.4f} J into {h} (+{dk:.0f} kg)'); tot+=dJ; n+=1
    else: print(f'  {t:3d}: NOT placeable (no settling insertion)')
un=len(miss)-n
print(f'\nplaced {n}/{len(miss)} for +{tot:.4f} J; {un} unplaceable (+{un}.0 J as misses)')
J=sj+tot+un+2
print(f'=> closed 9-craft J ~ {J:.3f}   (bars: 09-23 14.053 | 09-24 14.003 | 09-25 13.952 | 09-27 13.851)')
print('   caveat: prices are independent; several into one host escalate (measured), so this is optimistic')
