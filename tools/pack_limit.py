"""Why the set cover stalls: the DISJOINTNESS structure of the route pool (stage 11).

Two measurements per efficiency threshold thr (J per flyby):
  * the disjointness graph of the routes with J/fb <= thr: edge = the two routes share no target.
    Its density p and its max clique k (exact MILP) -- k is the largest number of EFFICIENT routes that
    can fly together in one fleet without wasting a flyby.
  * the exact set cover restricted to that sub-pool -> the J such a fleet actually reaches.

Usage: pack_limit.py cols.json [--thr ...] [--sizes ...]
"""
import os
for _v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'): os.environ.setdefault(_v,'1')
import json, sys, itertools, argparse, pathlib
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/'tools'))
import run_ialns as RI
from ctoc14.search import UNREACHABLE
from scipy.optimize import milp, LinearConstraint, Bounds, linprog
from scipy.sparse import coo_matrix
ALL=[t for t in range(1,301) if t not in UNREACHABLE]


def dedupe(cols):
    best={}
    for x in cols:
        k=frozenset(x['asts'])
        if k not in best or x['tank']<best[k]['tank']: best[k]=x
    return list(best.values())


def maxpack(sub, nmax=12, tlim=120.0):
    """max targets covered by <= nmax pairwise-DISJOINT routes"""
    nr=len(sub); ti={t:i for i,t in enumerate(ALL)}; nt=len(ALL)
    R=[];C=[]
    for j,x in enumerate(sub):
        for t in x['asts']: R.append(ti[t]); C.append(j)
    A=coo_matrix((np.ones(len(R)),(R,C)),shape=(nt,nr)).tocsr()
    cnt=np.array([len(x['asts']) for x in sub],float)
    cons=[LinearConstraint(A,-np.inf,np.ones(nt)),LinearConstraint(np.ones((1,nr)),-np.inf,nmax)]
    res=milp(-cnt,constraints=cons,integrality=np.ones(nr),bounds=Bounds(0,1),options=dict(time_limit=tlim,mip_rel_gap=0.002))
    if res.x is None: return 0,0,0.0
    pick=[j for j in range(nr) if res.x[j]>0.5]
    return len(pick), int(round(-res.fun)), float(sum(RI.cost(sub[j]['tank']) for j in pick))


def cover(sub, tlim=120.0):
    nr=len(sub); ti={t:i for i,t in enumerate(ALL)}; nt=len(ALL)
    c=np.concatenate([[RI.cost(x['tank']) for x in sub],np.ones(nt)])
    R=[];C=[]
    for j,x in enumerate(sub):
        for t in x['asts']: R.append(ti[t]); C.append(j)
    R+=list(range(nt)); C+=list(range(nr,nr+nt))
    A=coo_matrix((np.ones(len(R)),(R,C)),shape=(nt,nr+nt)).tocsr()
    res=milp(c,constraints=LinearConstraint(A,np.ones(nt),np.full(nt,np.inf)),integrality=np.ones(nr+nt),
             bounds=Bounds(0,1),options=dict(time_limit=tlim,mip_rel_gap=0.0005))
    if res.x is None: return None
    pick=[j for j in range(nr) if res.x[j]>0.5]
    return float(res.fun)+2, len(pick), int(round(res.x[nr:].sum()))


def density(sub, cap=400, seed=0):
    S=[set(x['asts']) for x in sub]
    if len(S)>cap:
        rng=np.random.default_rng(seed); S=[S[i] for i in rng.choice(len(S),cap,replace=False)]
    dis=[1 if not (a&b) else 0 for a,b in itertools.combinations(S,2)]
    return float(np.mean(dis)) if dis else 0.0


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('cols',nargs='+')
    ap.add_argument('--thr',default='0.0340,0.0369,0.0380,0.0400,0.0420,0.0450,1.0')
    ap.add_argument('--tag',default='')
    a=ap.parse_args()
    cols=[]
    for f in a.cols: cols+=json.load(open(f))
    cols=dedupe(cols)
    print(f'{len(cols)} distinct routes {a.tag}')
    print('thr      n   density  maxpack(<=12)            cover J   craft  miss')
    for thr in [float(s) for s in a.thr.split(',')]:
        sub=[x for x in cols if RI.cost(x['tank'])/len(x['asts'])<=thr]
        if len(sub)<4: continue
        p=density(sub)
        k,cov,sj=maxpack(sub)
        cv=cover(sub)
        print('%.4f %4d  %.4f   %2d routes %3d tgt sumJ %6.3f  %s'%(thr,len(sub),p,k,cov,sj,
              ('J %.3f  N %2d  miss %d'%cv if cv else 'n/a')),flush=True)


if __name__=='__main__':
    main()
