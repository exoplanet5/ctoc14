"""Single full-set beam under different lever settings (launch v_inf cap, per-leg dv cap, beam width, time weight)."""
import os
for _v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'): os.environ.setdefault(_v,'1')
import sys, time, json, pathlib, argparse
ROOT = pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/'tools'))
import numpy as np
import run_ialns as RI
from ctoc14.search import Params, beam_search, tour_json, UNREACHABLE
from ctoc14.colgen import CostModel
from ctoc14.impulsive import settle_tour
from ctoc14.constants import DAY, AU
ap=argparse.ArgumentParser(); ap.add_argument('--m0',type=float,default=1600.0); ap.add_argument('--beam',type=int,default=300)
ap.add_argument('--nproc',type=int,default=8); ap.add_argument('--vinf',default='2.0,3.0,4.0'); ap.add_argument('--dvmax',default='1.2')
ap.add_argument('--tofmax',default='400'); ap.add_argument('--lintofmax',default='600'); ap.add_argument('--wt',type=float,default=1.0); ap.add_argument('--settle',action='store_true'); ap.add_argument('--out',default='')
a=ap.parse_args()
CM=CostModel(s=0.70,reserve=2.0); ALL=[t for t in range(1,301) if t not in UNREACHABLE]
prize=np.zeros(300)
for t in ALL: prize[t-1]=1.0
eph=RI.eph(); wf=CM.w_fuel(0.5*(a.m0-640),a.m0)
for vc in [float(x) for x in a.vinf.split(',')]:
  for tm in [float(x) for x in a.tofmax.split(',')]:
   for ltm in [float(x) for x in a.lintofmax.split(',')]:
    for dm in [float(x) for x in a.dvmax.split(',')]:
        P=Params(beam=a.beam,w_fuel=wf,m_margin=40.0,dv_max=dm,lin_tofs=np.arange(20,ltm+1e-9,10)*DAY,
                 tofs=np.arange(15,tm+1e-9,5)*DAY,
                 lin_drmax=0.15*AU,vinf_cap=vc,tof_refine=True,prize=prize,w_t=a.wt)
        tic=time.time()
        best,beam=beam_search(eph,excluded=[],m0=a.m0,t_launch_grid=np.arange(0,800+1e-9,20)*DAY,P=P,n_proc=a.nproc,verbose=False)
        n=len(best.seq); tank=float(CM.tank(best.fuel,a.m0)); c=float(CM.cost(best.fuel,a.m0))
        msg=(f'vinf {vc:.1f} dv_max {dm:.1f} tof<={tm:.0f} lintof<={ltm:.0f}: {n:3d} flybys, fuel {best.fuel:6.0f} -> planner tank {tank:7.1f} '
             f'(J_i {c:.3f}, {c/max(n,1):.4f}/flyby), launch {best.t_launch/DAY:.0f} d, end {best.t/DAY/365.25:.2f} yr, '
             f'|vinf| {np.linalg.norm(best.vinf):.2f} km/s  ({time.time()-tic:.0f} s)')
        if a.settle:
            ip,miss,lag=settle_tour(eph,tour_json(best),lambda ip: RI.settle(ip,100))
            msg += (f'  TWIN tank {ip.tank():.0f} (J_i {RI.cost(ip.tank()):.3f}, {RI.cost(ip.tank())/n:.4f}/flyby), miss {miss:.0f} km'
                    if ip is not None else f'  TWIN FAILED (miss {miss:.0f})')
        print(msg,flush=True)
        if a.out:
            t=tour_json(best); t['tag']=f'lever:{vc}:{dm}'; t['targets']=sorted(x[0] for x in best.seq)
            open(a.out,'a').write(json.dumps(t)+'\n')
