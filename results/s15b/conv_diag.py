import os
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'): os.environ.setdefault(v,'1')
import sys, multiprocessing as mp; sys.path.insert(0,'.'); sys.path.insert(0,'tools')
import numpy as np, run_ialns as RI
from ctoc14.impulsive import to_exact
from ctoc14.constants import DAY
fl=RI.IFleet('results/s15b/best_closed/fleet')
def job(n):
    st=fl.routes[n]['st']; ip=RI.ipr(st)
    Yf,_=ip.integrate(); dm,_=ip.misses(Yf); twin_miss=float(np.linalg.norm(dm,axis=1).max())
    ts=np.asarray(ip.ts); dt=np.diff(ts)/DAY; h=ip.h/DAY
    irregular=int(np.sum(np.abs(dt-h)>0.5)); imp=np.linalg.norm(ip.Ts,axis=1)
    cap=float(ip.tcap) if np.isscalar(ip.tcap) or np.ndim(ip.tcap)==0 else float(np.max(ip.tcap))
    gm=to_exact(ip,RI.eph()); Yg,_=gm.integrate(dense_samples=True); dg,_=gm.misses(Yg)
    dg=np.linalg.norm(dg,axis=1)
    return n,len(ip.asts),len(ts),irregular,float(np.min(dt)),float(imp.max()),cap,twin_miss,float(dg.max()),float(np.median(dg))
with mp.get_context('fork').Pool(9) as pl: R=pl.map(job,sorted(fl.routes))
print(f"{'route':5} {'fb':>3} {'nodes':>5} {'irreg':>5} {'min dt d':>8} {'max imp':>8} {'tcap':>6} {'twin miss':>9} {'to_exact miss max / median km':>30}  export")
ok={'s2','s1','t1','h05'}
for n,fb,nn,irr,mdt,mi,cap,tm,gmax,gmed in R:
    print(f"{n:5} {fb:3d} {nn:5d} {irr:5d} {mdt:8.3f} {mi:8.3f} {cap:6.3f} {tm:9.1f} {gmax:15.3e} / {gmed:10.3e}   {'OK' if n in ok else 'FAILED'}")
