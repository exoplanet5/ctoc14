"""Grow-overlapping -> assign -> PRUNE: the deepened fleet covers 298 with 67 surplus flybys.
Each surplus flyby is one a route need not fly. Assign every target to exactly one route (cheapest craft wins),
drop the rest, re-settle, and see whether sum J_i falls below the cover's 12.355."""
import os
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'): os.environ.setdefault(v,'1')
import sys, collections, multiprocessing as mp; sys.path.insert(0,'.'); sys.path.insert(0,'tools')
import numpy as np, run_ialns as RI
fl = RI.IFleet('results/s7/deep1/fleet')
A = {n: set(int(x) for x in r['st']['asts']) for n, r in fl.routes.items()}
own = collections.defaultdict(list)
for n, s in A.items():
    for t in s: own[t].append(n)
dup = {t: v for t, v in own.items() if len(v) > 1}
print(f'{len(fl.routes)} routes, covered {len(own)}, targets on >1 route: {len(dup)}')
print(f'sum J_i now {sum(RI.cost(r["tank"]) for r in fl.routes.values()):.3f}  (cover optimum 12.355)')
# assign each duplicated target to the route with the smallest tank; the others drop it
drop = collections.defaultdict(list)
for t, v in dup.items():
    keep = min(v, key=lambda n: fl.routes[n]['tank'])
    for n in v:
        if n != keep: drop[n].append(t)
jobs = [(n, fl.routes[n]['st'], drop[n]) for n in sorted(fl.routes) if drop[n]]
print('dropping per route: ' + ' '.join(f'{n}:-{len(drop[n])}' for n, _, _ in jobs))
with mp.get_context('fork').Pool(8) as pl: res = pl.map(RI.w_remove, jobs)
tot = 0.0; nfb = 0
for n in sorted(fl.routes):
    r = [x for x in res if x['name'] == n]
    if r and r[0].get('ok'):
        t0 = fl.routes[n]['tank']; t1 = r[0]['tank']; k = len(r[0]['st']['asts'])
        print(f'  {n}: {len(A[n])} -> {k} fb, {t0:.0f} -> {t1:.0f} kg, J_i {RI.cost(t0):.3f} -> {RI.cost(t1):.3f}')
        tot += RI.cost(t1); nfb += k
    else:
        tot += RI.cost(fl.routes[n]['tank']); nfb += len(A[n])
        if drop[n]: print(f'  {n}: prune FAILED, kept as is')
print(f'\npruned fleet: {nfb} flybys, sum J_i {tot:.3f} -> J {tot+2:.3f}   (need < 14.153 today)')
