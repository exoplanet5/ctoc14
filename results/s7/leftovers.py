import sys, json, glob, pathlib, collections
sys.path.insert(0,'.'); sys.path.insert(0,'tools')
import numpy as np, run_ialns as RI
from ctoc14.kepler import load_mea
from ctoc14.search import UNREACHABLE
ALL=set(range(1,301))-set(UNREACHABLE)
runs={}
for d in ('results/s7/F1/fleet','results/s7/F1p/fleet','results/s7/R1/fleet','results/s7/T1/fleet'):
    if pathlib.Path(d,'fleet.json').exists(): runs[d.split('/')[2]]=ALL-set(RI.IFleet(d).coverage())
# stage 6 best fill
for d in ('results/n8s6/PART/best','results/n8/fill2'):
    if pathlib.Path(d,'fleet.json').exists(): runs['s6:'+d.split('/')[-1]]=ALL-set(RI.IFleet(d).coverage())
cnt=collections.Counter()
for v in runs.values(): cnt.update(v)
print('runs:', {k:len(v) for k,v in runs.items()})
n=len(runs); per=sorted(cnt.items(), key=lambda kv:-kv[1])
always=[t for t,c in per if c==n]
print(f'left in ALL {n} runs: {len(always)}')
print(f'left in >= {n-1}: {len([t for t,c in per if c>=n-1])}')
ids,el=load_mea()
E=dict(zip(ids.tolist(), el))
def stats(S,name):
    if not S: return
    A=np.array([E[t] for t in S])
    q=A[:,0]*(1-A[:,1]); Q=A[:,0]*(1+A[:,1])
    print(f'{name} (n={len(S)}): a {np.median(A[:,0]):.2f} e {np.median(A[:,1]):.2f} i {np.median(A[:,2]):.1f} deg | '
          f'q med {np.median(q):.2f} Q med {np.median(Q):.2f} | i>10deg {int((A[:,2]>10).sum())} ({100*(A[:,2]>10).mean():.0f}%), '
          f'a>1.6 {int((A[:,0]>1.6).sum())} ({100*(A[:,0]>1.6).mean():.0f}%), q>1.05 {int((q>1.05).sum())}, Q<0.95 {int((Q<0.95).sum())}')
stats(always,'always left'); stats(sorted(ALL-set(always)),'rest'); stats(sorted(ALL),'all')
json.dump(dict(always=always, counts={str(k):v for k,v in cnt.items()}, nruns=n), open('results/s7/leftovers.json','w'), indent=1)
print('always:',always)
