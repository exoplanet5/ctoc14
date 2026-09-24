import os
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'): os.environ.setdefault(v,'1')
import sys, json, multiprocessing as mp; sys.path.insert(0,'.'); sys.path.insert(0,'tools')
import numpy as np, run_ialns as RI
from ctoc14.constants import DAY
F=json.load(open('results/s15/select/r02f/frontier.json'))
def load(f):
    z=np.load(f); st={k:z[k] for k in z.files}; st['tL']=float(st['tL']); return st
out={}
for K in ('291','292','293','294','295','296'):
    fr=F[K]; sts={f'{i}':load(p) for i,p in enumerate(fr['picks'])}
    tanks={k:float(RI.ipr(s).tank()) for k,s in sts.items()}
    with mp.get_context('fork').Pool(2) as pl:
        C=[c for lst in pl.map(RI.w_cands,[(k,s,fr['misses'],0.30) for k,s in sts.items()]) for c in lst]
    print(f"\n{K} covered, sum J_i {fr['sumJi']}: depths {[len(s['asts']) for s in sts.values()]} tanks {[round(t) for t in tanks.values()]}")
    rows={}
    for t in fr['misses']:
        cs=[c for c in C if c['ast']==t and np.abs(np.asarray(sts[c['host']]['tf'])-c['t']).min()>10*DAY]
        cs.sort(key=lambda c:c['dist'])
        rows[t]=[(c['host'],round(c['dist'],3),len(sts[c['host']]['asts'])) for c in cs[:3]]
        print(f"   miss {t:3d}: " + ('; '.join(f"route {h} (depth {d}) {x:.3f} AU" for h,x,d in rows[t]) if rows[t] else 'NO approach within 0.30 AU'))
    out[K]=rows
json.dump(out,open('results/s15/miss_dist.json','w'),indent=1)
