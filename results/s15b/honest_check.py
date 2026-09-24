import os
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'): os.environ.setdefault(v,'1')
import sys, multiprocessing as mp; sys.path.insert(0,'.'); sys.path.insert(0,'tools')
import numpy as np, run_ialns as RI
from ctoc14.constants import DAY
D=sys.argv[1]; fl=RI.IFleet(D)
def job(n):
    st=fl.routes[n]['st']; ip=RI.ipr(st); t0=float(ip.tank())
    ts=np.asarray(ip.ts); dt=np.diff(ts)/DAY; h=ip.h/DAY
    irr=int(np.sum(np.abs(dt-h)>0.01)); cap=float((np.linalg.norm(ip.Ts,axis=1)/ip.tcap).max())
    Yf,_=ip.integrate(); dm,_=ip.misses(Yf); m0=float(np.linalg.norm(dm,axis=1).max())
    # the windowed cap: no 20-d window may hold more than one cap's worth of impulse
    w=0.0
    for s in ts:
        sel=(ts>=s)&(ts<s+ip.h-1)
        w=max(w,float(np.linalg.norm(ip.Ts[sel],axis=1).sum()/ip.tcap))
    miss=RI.settle(ip,100); t1=float(ip.tank())
    return n,len(ip.asts),irr,cap,w,m0,t0,t1,miss
with mp.get_context('fork').Pool(9) as pl: R=pl.map(job,sorted(fl.routes))
cov=set().union(*[set(int(x) for x in fl.routes[n]['st']['asts']) for n in fl.routes]) - {131,144}
print(f"{'route':5} {'fb':>3} {'irreg':>5} {'max imp/cap':>11} {'max 20d window/cap':>18} {'miss km':>8} {'tank':>7} {'re-settled':>10}")
S0=S1=0
for n,fb,irr,cap,w,m0,t0,t1,miss in R:
    print(f"{n:5} {fb:3d} {irr:5d} {cap:11.3f} {w:18.3f} {m0:8.1f} {t0:7.1f} {t1:10.1f} (miss {miss:.0f})")
    S0+=RI.cost(t0); S1+=RI.cost(t1)
print(f"covered {len(cov)} distinct reachable; sum J_i stored {S0:.4f}, re-settled {S1:.4f}")
