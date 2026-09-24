"""Grow-overlapping -> prune -> re-settle, over the WHOLE efficient pool (stage 11).

The exact set cover is forced to pick routes that already partition; the pool's most efficient routes overlap and are
therefore unusable to it.  This picks K routes from the efficient pool to MAXIMISE coverage with overlap allowed,
then gives every duplicated target to one host only, re-settles, and reports the resulting sum J_i.
The measured refund of a removed flyby on a deep route is only ~13 kg, so this usually loses -- it is the test of
whether any overlapping combination refunds more than the partition penalty costs.

Usage: overlap_prune.py cols.json [--K 10] [--thr 0.040] [--nproc 6]
"""
import os
for _v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'): os.environ.setdefault(_v,'1')
import sys,json,argparse,pathlib,collections,multiprocessing as mp
import numpy as np
ROOT=pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/'tools'))
import run_ialns as RI, pack_limit as PL
from ctoc14.search import UNREACHABLE
from scipy.optimize import milp,LinearConstraint,Bounds
from scipy.sparse import coo_matrix
ALL=[t for t in range(1,301) if t not in UNREACHABLE]
ap=argparse.ArgumentParser(); ap.add_argument('cols',nargs='?',default='results/s11/cols_base.json')
ap.add_argument('--K',type=int,default=10); ap.add_argument('--thr',type=float,default=0.040)
ap.add_argument('--nproc',type=int,default=6); ap.add_argument('--lam',type=float,default=0.0,
                help='penalty on sum J_i in the selection objective (0 = pure max coverage)')
a=ap.parse_args()
cols=PL.dedupe(json.load(open(a.cols)))
sub=[x for x in cols if RI.cost(x['tank'])/len(x['asts'])<=a.thr]
nr=len(sub); ti={t:i for i,t in enumerate(ALL)}; nt=len(ALL)
R=[];C=[]
for j,x in enumerate(sub):
    for t in x['asts']: R.append(ti[t]); C.append(j)
A=coo_matrix((np.ones(len(R)),(R,C)),shape=(nt,nr)).tocsr()
from scipy.sparse import hstack, csr_matrix, eye
M=hstack([A,-eye(nt,format='csr')]).tocsr()
obj=np.concatenate([a.lam*np.array([RI.cost(x['tank']) for x in sub]),-np.ones(nt)])
res=milp(obj,constraints=[LinearConstraint(M,0,np.inf),
                          LinearConstraint(csr_matrix(np.concatenate([np.ones(nr),np.zeros(nt)])[None,:]),-np.inf,a.K)],
         integrality=np.ones(nr+nt),bounds=Bounds(0,1),options=dict(time_limit=300,mip_rel_gap=0.001))
pick=[j for j in range(nr) if res.x[j]>0.5]
sets={j:set(sub[j]['asts']) for j in pick}
cov=set().union(*sets.values()); tot=sum(len(s) for s in sets.values())
print(f'{len(pick)} efficient routes (thr {a.thr}): covered {len(cov)}/298, {tot} flybys, surplus {tot-len(cov)}, '
      f'sum J_i {sum(RI.cost(sub[j]["tank"]) for j in pick):.3f}')
own=collections.defaultdict(list)
for j,s in sets.items():
    for t in s: own[t].append(j)
drop=collections.defaultdict(list)
for t,v in own.items():
    if len(v)>1:
        keep=min(v,key=lambda j: sub[j]['tank'])
        for j in v:
            if j!=keep: drop[j].append(t)
jobs=[]
for j in pick:
    z=np.load(sub[j]['f']); st={k:z[k] for k in z.files}; st['tL']=float(st['tL'])
    if drop[j]: jobs.append((str(j),st,drop[j]))
print('dropping: '+' '.join(f'{n}:-{len(d)}' for n,_,d in jobs))
with mp.get_context('fork').Pool(a.nproc) as pl: out=pl.map(RI.w_remove,jobs)
res_by={r['name']:r for r in out}
tot=0.0; nfb=0; okall=True
for j in pick:
    r=res_by.get(str(j))
    if r is None: tot+=RI.cost(sub[j]['tank']); nfb+=len(sets[j]); continue
    if not r.get('ok'): print(f'  {j}: removal did not settle'); okall=False; tot+=RI.cost(sub[j]['tank']); nfb+=len(sets[j]); continue
    t0=sub[j]['tank']; t1=r['tank']; k=len(r['st']['asts'])
    print(f'  route {j}: {len(sets[j])} -> {k} fb, {t0:.0f} -> {t1:.0f} kg, J_i {RI.cost(t0):.3f} -> {RI.cost(t1):.3f}')
    tot+=RI.cost(t1); nfb+=k
print(f'pruned: {nfb} flybys, covered {len(cov)}, sum J_i {tot:.3f} -> J {tot+2+(298-len(cov)):.3f}  (banked 14.355)')
