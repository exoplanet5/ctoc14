"""Fast greedy+swap version of tools/carrier_capacity.py (stage 11): max targets coverable by K carriers,
<= cap targets each, a target assignable to a carrier only if its closest approach is under dmax.
Exact assignment for a fixed carrier set (max-flow / Hopcroft on the b-matching); greedy + 1-swap over carriers."""
import os
for _v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'): os.environ.setdefault(_v,'1')
import sys,pickle,argparse,pathlib
import numpy as np
ROOT=pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/'tools'))
from ctoc14.search import UNREACHABLE
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_flow
ALL=[t for t in range(1,301) if t not in UNREACHABLE]
TI={t:i for i,t in enumerate(ALL)}


def cand_sets(lib,dmax):
    return [set(t for t,_,_ in c['base'])|set(t for t,(d,_) in c['neigh'].items() if d<dmax) for c in lib]


def assign(sets,cap):
    """max b-matching: source -> carrier(cap) -> target(1) -> sink"""
    K=len(sets); nt=len(ALL); N=1+K+nt+1; snk=N-1
    rows=[];cols=[];vals=[]
    for j in range(K): rows.append(0); cols.append(1+j); vals.append(cap)
    for j,S in enumerate(sets):
        for t in S:
            if t in TI: rows.append(1+j); cols.append(1+K+TI[t]); vals.append(1)
    for i in range(nt): rows.append(1+K+i); cols.append(snk); vals.append(1)
    g=csr_matrix((np.array(vals,dtype=np.int32),(rows,cols)),shape=(N,N))
    return int(maximum_flow(g,0,snk).flow_value)


def solve(C,K,cap,iters=40,seed=0,verbose=False):
    rng=np.random.default_rng(seed)
    cur=[]
    for _ in range(K):                                    # greedy
        best=(-1,None)
        for j in range(len(C)):
            if j in cur: continue
            v=assign([C[x] for x in cur]+[C[j]],cap)
            if v>best[0]: best=(v,j)
        cur.append(best[1])
    val=assign([C[x] for x in cur],cap)
    for it in range(iters):                               # 1-swap
        improved=False
        order=rng.permutation(len(C))
        for p in range(K):
            for j in order[:200]:
                j=int(j)
                if j in cur: continue
                trial=cur[:p]+cur[p+1:]+[j]
                v=assign([C[x] for x in trial],cap)
                if v>val: cur,val=trial,v; improved=True; break
            if improved: break
        if not improved: break
    return val,cur


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('libs',nargs='+')
    ap.add_argument('--K',default='8,10'); ap.add_argument('--cap',type=int,default=40)
    ap.add_argument('--dmax',default='0.06,0.10,0.15')
    ap.add_argument('--sizes',default=''); ap.add_argument('--reps',type=int,default=2)
    a=ap.parse_args()
    lib=[]
    for f in a.libs: lib+=pickle.load(open(f,'rb'))
    print(f'{len(lib)} carriers')
    for dmax in [float(s) for s in a.dmax.split(',')]:
        C=cand_sets(lib,dmax)
        print(f' dmax {dmax}: cand median {np.median([len(s) for s in C]):.0f}, union {len(set(t for s in C for t in s))}')
        for K in [int(s) for s in a.K.split(',')]:
            if a.sizes:
                rng=np.random.default_rng(0)
                for s in [int(x) for x in a.sizes.split(',')]:
                    s=len(C) if (s==0 or s>len(C)) else s
                    vals=[]
                    for rep in range(1 if s>=len(C) else a.reps):
                        sub=C if s>=len(C) else [C[i] for i in rng.choice(len(C),s,replace=False)]
                        vals.append(solve(sub,K,a.cap)[0])
                    print(f'   K={K} cap={a.cap} n_carriers {s:4d}: covered {vals}',flush=True)
                    if s>=len(C): break
            else:
                v,pick=solve(C,K,a.cap)
                print(f'   K={K} cap={a.cap}: covered {v}/298',flush=True)


if __name__=='__main__':
    main()
