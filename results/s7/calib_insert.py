import os
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'): os.environ.setdefault(v,'1')
import sys, json, multiprocessing as mp; sys.path.insert(0,'.'); sys.path.insert(0,'tools')
import numpy as np, run_ialns as RI
from ctoc14.constants import DAY
rows=[]
pool=mp.get_context('fork').Pool(2)
for rn in ('01','03'):
    fl=RI.IFleet('results/s7/base8'); st=fl.routes[rn]['st']; tank=fl.routes[rn]['tank']
    have=set(int(x) for x in st['asts']); T=[t for t in range(1,301) if t not in have and t not in (131,144)]
    c=RI.w_cands((rn,st,T,0.25)); best={}
    for x in sorted(c,key=lambda x:x['dist']):
        if np.abs(np.asarray(st['tf'])-x['t']).min()>12*DAY: best.setdefault(x['ast'],x)
    b=sorted(best.values(),key=lambda x:x['dist']); sel=b[::max(1,len(b)//24)][:24]
    res=pool.map(RI.w_insert,[(rn,st,x['ast'],x['t']) for x in sel])
    for x,r in zip(sel,res): rows.append((x['dist'], (r['tank']-tank) if r.get('ok') else np.nan))
rows.sort()
R=np.array(rows)
for lo,hi in ((0,.03),(.03,.06),(.06,.1),(.1,.15),(.15,.25)):
    m=(R[:,0]>=lo)&(R[:,0]<hi); v=R[m,1]
    print(f'dist {lo:.2f}-{hi:.2f} AU: n {m.sum():2d}, ok {np.isfinite(v).sum():2d}, dkg median {np.nanmedian(v) if np.isfinite(v).any() else float("nan"):.0f}, all {np.round(v[np.isfinite(v)]).astype(int).tolist()}')
