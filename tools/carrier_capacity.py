"""IDEALISED upper bound on what any route farm can get out of a carrier library (stage 11).

A route grown from carrier c can only ever fly targets that carrier can reach: its DP base plus the targets whose
closest approach to its seed trajectory is under `dmax` (the measured insertion price is ~12 kg per 0.01 AU, so
d < 0.06 AU is "cheap" and d < 0.10 AU is "affordable").  Ignore timing, ignore the settle, ignore the tank, and ask
the purely combinatorial question:

    pick K carriers and assign every target to ONE of them, at most `cap` targets per carrier -- is 298 coverable?

Whatever a farm + set cover can ever achieve from this library is bounded by this number.  If the bound is already
short of 298, the library, not the farm, is the wall.

Usage: carrier_capacity.py lib.pkl [lib2.pkl ...] [--K 8,9,10] [--cap 40] [--dmax 0.06]
"""
import os
for _v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'): os.environ.setdefault(_v,'1')
import sys, pickle, argparse, pathlib
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/'tools'))
from ctoc14.search import UNREACHABLE
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix, csr_matrix, hstack, vstack
ALL=[t for t in range(1,301) if t not in UNREACHABLE]


def cand_sets(lib, dmax):
    return [set(t for t,_,_ in c['base']) | set(t for t,(d,_) in c['neigh'].items() if d < dmax) for c in lib]


def capacity(C, K, cap, tlim=300.0):
    """max targets coverable by K carriers, <= cap targets each. y binary, x continuous (flow -> integral)."""
    nc=len(C); ti={t:i for i,t in enumerate(ALL)}; nt=len(ALL)
    pairs=[(j,ti[t]) for j,S in enumerate(C) for t in S if t in ti]
    npair=len(pairs)
    # vars: x (npair) | y (nc) | z (nt)
    rows=[];colsA=[];vals=[]
    for p,(j,i) in enumerate(pairs):
        rows.append(i); colsA.append(p); vals.append(1.0)
    Az=coo_matrix((vals,(rows,colsA)),shape=(nt,npair)).tocsr()
    con_cov=hstack([Az, csr_matrix((nt,nc)), -1.0*csr_matrix(np.eye(nt))]).tocsr()   # sum_j x - z >= 0
    rows=[];colsB=[];vals=[]
    for p,(j,i) in enumerate(pairs):
        rows.append(j); colsB.append(p); vals.append(1.0)
    Ac=coo_matrix((vals,(rows,colsB)),shape=(nc,npair)).tocsr()
    con_cap=hstack([Ac, -cap*csr_matrix(np.eye(nc)), csr_matrix((nc,nt))]).tocsr()   # sum_t x - cap*y <= 0
    con_K=csr_matrix(np.concatenate([np.zeros(npair),np.ones(nc),np.zeros(nt)])[None,:])
    # x_p <= y_j
    rows=[];colsC=[];vals=[]
    for p,(j,i) in enumerate(pairs):
        rows+= [p,p]; colsC+=[p, npair+j]; vals+=[1.0,-1.0]
    con_xy=coo_matrix((vals,(rows,colsC)),shape=(npair,npair+nc+nt)).tocsr()
    obj=np.concatenate([np.zeros(npair),np.zeros(nc),-np.ones(nt)])
    integ=np.concatenate([np.zeros(npair),np.ones(nc),np.zeros(nt)])
    res=milp(obj,constraints=[LinearConstraint(con_cov,0,np.inf),LinearConstraint(con_cap,-np.inf,0),
                              LinearConstraint(con_K,-np.inf,K),LinearConstraint(con_xy,-np.inf,0)],
             integrality=integ,bounds=Bounds(0,1),options=dict(time_limit=tlim,mip_rel_gap=0.001))
    if res.x is None: return None
    return int(round(-res.fun)), [j for j in range(nc) if res.x[npair+j]>0.5]


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('libs',nargs='+')
    ap.add_argument('--K',default='8,9,10'); ap.add_argument('--cap',type=int,default=40)
    ap.add_argument('--dmax',default='0.06,0.10'); ap.add_argument('--tlim',type=float,default=300.0)
    a=ap.parse_args()
    lib=[]
    for f in a.libs: lib+=pickle.load(open(f,'rb'))
    print(f'{len(lib)} carriers')
    for dmax in [float(s) for s in a.dmax.split(',')]:
        C=cand_sets(lib,dmax)
        print(f'  dmax {dmax}: candidate set median {np.median([len(s) for s in C]):.0f}, union {len(set(t for s in C for t in s))}')
        for K in [int(s) for s in a.K.split(',')]:
            r=capacity(C,K,a.cap,a.tlim)
            print(f'    K={K:2d} cap={a.cap}: max covered {r[0] if r else "fail"}/298',flush=True)


if __name__=='__main__':
    main()
