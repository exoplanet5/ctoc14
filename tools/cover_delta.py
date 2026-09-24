"""Does ADDING genuinely new routes move the exact set cover?  (stage 11)

Baseline pool (results/s11/cols_base.json) vs baseline + each new farm directory, exact MILP + LP bound each time.
This is the honest version of the pool-scaling test: random subsamples of a FIXED pool saturate by construction.
"""
import os
for _v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'): os.environ.setdefault(_v,'1')
import sys,json,glob,argparse,pathlib,time
import numpy as np
ROOT=pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/'tools'))
import run_ialns as RI, route_cover as RC, cover_curve as CV, pack_limit as PL
ap=argparse.ArgumentParser(); ap.add_argument('--base',default='results/s11/cols_base.json')
ap.add_argument('--add',nargs='*',default=[]); ap.add_argument('--nproc',type=int,default=4)
ap.add_argument('--tlim',type=float,default=200.0); ap.add_argument('--out',default='results/s11/delta.json')
a=ap.parse_args()
base=json.load(open(a.base))
rows=[]
def report(tag,cols):
    c=PL.dedupe(cols)
    r=CV.solve_quiet(c,a.tlim)
    eff=[x for x in c if RI.cost(x['tank'])/len(x['asts'])<=0.040]
    k,cov,sj=PL.maxpack(eff,12,90.0) if len(eff)>=4 else (0,0,0.0)
    row=dict(tag=tag,n=len(c),J=r['J'],lp=r['lp'],craft=r['n'],covered=r['covered'],sumJ=r['sumJ'],
             n_eff=len(eff),maxdisjoint=k,disjoint_cov=cov,pick=r['pick'])
    rows.append(row); json.dump(rows,open(a.out,'w'),indent=1)
    print('%-28s n=%4d  J %.3f  LP %.3f  craft %2d  cov %3d  sumJ %.3f | eff %3d, max disjoint %d cov %3d'%(
        tag,len(c),r['J'],r['lp'],r['n'],r['covered'],r['sumJ'],len(eff),k,cov),flush=True)
report('baseline',base)
acc=list(base)
for d in a.add:
    files=sorted(glob.glob(d+'/route_*.npz'))
    if not files: print('  (no routes in %s)'%d); continue
    import multiprocessing as mp
    with mp.get_context('fork').Pool(a.nproc) as pl:
        new=[x for x in pl.map(RC.w_load,files) if x and x['tank']<=1400 and x['asts']]
    print('  +%d routes from %s'%(len(new),d))
    acc=acc+new
    report('base+'+pathlib.Path(d).name,acc)
