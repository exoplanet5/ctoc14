"""How deep can one beam route go as a function of the planner mass m0 (the fuel cap m0-640 limits the reachable tank)?"""
import os
for _v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'): os.environ.setdefault(_v,'1')
import sys, time, pathlib, json, argparse
ROOT = pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/'tools'))
import numpy as np
import run_ialns as RI
from ctoc14.search import Params, beam_search, tour_json, UNREACHABLE
from ctoc14.colgen import CostModel
from ctoc14.constants import DAY, AU
ap=argparse.ArgumentParser(); ap.add_argument('--m0', default='1000,1600,2000'); ap.add_argument('--beam',type=int,default=300)
ap.add_argument('--nproc',type=int,default=8); ap.add_argument('--wt',type=float,default=1.0); ap.add_argument('--out',default='')
a=ap.parse_args()
CM=CostModel(s=0.70,reserve=2.0)
ALL=[t for t in range(1,301) if t not in UNREACHABLE]
prize=np.zeros(300)
for t in ALL: prize[t-1]=1.0
eph=RI.eph()
rows=[]
for m0 in [float(x) for x in a.m0.split(',')]:
    wf=CM.w_fuel(0.5*(m0-640),m0)
    P=Params(beam=a.beam,w_fuel=wf,m_margin=40.0,dv_max=1.2,lin_tofs=np.arange(20,600+1e-9,10)*DAY,
             lin_drmax=0.15*AU,vinf_cap=2.0,tof_refine=True,prize=prize,w_t=a.wt)
    tic=time.time()
    best,beam=beam_search(eph,excluded=[],m0=m0,t_launch_grid=np.arange(0,800+1e-9,20)*DAY,P=P,n_proc=a.nproc,verbose=False)
    n=len(best.seq); tank=float(CM.tank(best.fuel,m0)); c=float(CM.cost(best.fuel,m0))
    print(f'm0 {m0:5.0f} (w_fuel {wf:.3f}): best {n:3d} flybys, fuel {best.fuel:6.0f} kg -> tank {tank:7.1f}, J_i {c:.3f}, '
          f'J/flyby {c/max(n,1):.4f}, launch {best.t_launch/DAY:.0f} d, end {best.t/DAY/365.25:.2f} yr  ({time.time()-tic:.0f} s)',flush=True)
    # depth profile of the final beam
    ns=sorted(len(s.seq) for s in beam)
    print(f'         final beam depth p10/p50/p90 {ns[len(ns)//10]}/{ns[len(ns)//2]}/{ns[9*len(ns)//10]}',flush=True)
    rows.append(dict(m0=m0,n=n,fuel=float(best.fuel),tank=tank,cost=c))
    if a.out:
        t=tour_json(best); t['tag']=f'm0test:{m0:.0f}'; t['targets']=sorted(x[0] for x in best.seq)
        open(a.out,'a').write(json.dumps(t)+'\n')
