"""J vs the number of COMPLETE FLEETS in the pool (stage 11).

The exact cover can only use a route if the rest of the fleet covers its complement, so the useful unit of pool growth
is a whole near-disjoint fleet, not a route.  Group every route by the directory it was built in, subsample k of those
groups, and solve the cover on their union.  The slope of J(k) sizes the "mass-produce partitions" campaign; compare
with the J(n_routes) curve from cover_curve.py, which plateaus.
"""
import os
for _v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'): os.environ.setdefault(_v,'1')
import sys,json,pathlib,argparse,collections,re
import numpy as np
ROOT=pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/'tools'))
import run_ialns as RI, cover_curve as CV, pack_limit as PL
ap=argparse.ArgumentParser(); ap.add_argument('cols',nargs='?',default='results/s11/cols_base.json')
ap.add_argument('--reps',type=int,default=5); ap.add_argument('--tlim',type=float,default=60.0)
ap.add_argument('--min-size',type=int,default=5)
a=ap.parse_args()
cols=json.load(open(a.cols))
g=collections.defaultdict(list)
for x in cols: g[str(pathlib.Path(x['f']).parent)].append(x)
# keep only directories that look like a FLEET (a partition attempt): >= min-size routes, little internal overlap
fleets={}
for d,v in g.items():
    if len(v)<a.min_size: continue
    tot=sum(len(x['asts']) for x in v); cov=len(set(t for x in v for t in x['asts']))
    if tot-cov > 0.25*cov: continue            # a farm directory, not a fleet
    fleets[d]=v
print(f'{len(fleets)} fleet-like directories, {sum(len(v) for v in fleets.values())} routes '
      f'(of {len(g)} directories / {len(cols)} routes)')
keys=sorted(fleets)
rng=np.random.default_rng(0)
print('k_fleets  n_routes   mean J    min J   craft')
for k in [2,4,6,8,12,16,20,24,28,len(keys)]:
    if k>len(keys): continue
    Js=[];ns=[];cr=[]
    for rep in range(1 if k>=len(keys) else a.reps):
        sel=keys if k>=len(keys) else [keys[i] for i in rng.choice(len(keys),k,replace=False)]
        sub=PL.dedupe([x for d in sel for x in fleets[d]])
        if len(set(t for x in sub for t in x['asts']))<298: continue
        r=CV.solve_quiet(sub,a.tlim)
        if r: Js.append(r['J']); ns.append(len(sub)); cr.append(r['n'])
    if Js: print('%5d     %6.0f   %7.3f  %7.3f   %.1f'%(k,np.mean(ns),np.mean(Js),min(Js),np.mean(cr)),flush=True)
    if k>=len(keys): break
