"""Does the idealised carrier-capacity bound rise with library size?  (stage 11)"""
import os
for _v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'): os.environ.setdefault(_v,'1')
import sys,pickle,argparse,pathlib
import numpy as np
ROOT=pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/'tools'))
import carrier_capacity as CC
ap=argparse.ArgumentParser(); ap.add_argument('libs',nargs='+')
ap.add_argument('--K',type=int,default=8); ap.add_argument('--cap',type=int,default=40)
ap.add_argument('--dmax',type=float,default=0.10); ap.add_argument('--tlim',type=float,default=150.0)
ap.add_argument('--sizes',default='60,120,240,480,0'); ap.add_argument('--reps',type=int,default=2)
a=ap.parse_args()
lib=[]
for f in a.libs: lib+=pickle.load(open(f,'rb'))
C=CC.cand_sets(lib,a.dmax)
print(f'{len(lib)} carriers, dmax {a.dmax}, K={a.K}, cap={a.cap}')
rng=np.random.default_rng(0)
for s in [int(x) for x in a.sizes.split(',')]:
    if s==0 or s>=len(C): s=len(C)
    out=[]
    for rep in range(1 if s>=len(C) else a.reps):
        sub=C if s>=len(C) else [C[i] for i in rng.choice(len(C),s,replace=False)]
        r=CC.capacity(sub,a.K,a.cap,a.tlim)
        out.append(r[0] if r else -1)
    print(f'  n_carriers {s:4d}: max covered {out}',flush=True)
    if s==len(C): break
