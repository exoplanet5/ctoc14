"""Balanced partition by MERGING an existing disjoint fleet, then re-beaming each merged group.

A greedy cover cannot be balanced (46, 43, 38, 33, ... and the rest unchainable), but J_i is convex in the route
length, so the fleet we want is N equal groups of ~298/N targets. `ifleet_t10d` already IS a disjoint partition of
all 298 into 10 routes; merge it down to N groups (greedily, by how close each route's targets pass the other's
trajectory) and run the beam inside each group with the new levers (m0 1600, vinf 4.0, prize 1 J per target).
Usage: merge_beam.py src_fleet out_dir [--n 7] [--m0 1600] [--nproc 8]
"""
import os
for _v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'): os.environ.setdefault(_v,'1')
import sys, json, time, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/'tools'))
import numpy as np
import run_ialns as RI
from ctoc14.search import Params, beam_search, tour_json, UNREACHABLE
from ctoc14.colgen import CostModel
from ctoc14.impulsive import settle_tour
from ctoc14.constants import DAY, AU

ALL=[t for t in range(1,301) if t not in UNREACHABLE]; CM=CostModel(s=0.70,reserve=2.0)

class DeepCollect(list):
    def __init__(self,minn): super().__init__(); self.minn=minn
    def extend(self,states):
        for s in states:
            if len(s.seq)>=self.minn: list.append(self,s)

class BP(Params):
    def __getstate__(self):
        d=dict(self.__dict__); d.pop('collect',None); return d

ap=argparse.ArgumentParser(); ap.add_argument('src'); ap.add_argument('out'); ap.add_argument('--n',type=int,default=7)
ap.add_argument('--m0',type=float,default=1600.0); ap.add_argument('--vinf',type=float,default=4.0)
ap.add_argument('--beam',type=int,default=300); ap.add_argument('--nproc',type=int,default=8)
ap.add_argument('--pad',type=int,default=12,help='extra nearest targets added to each group')
ap.add_argument('--cap',type=float,default=1.35,help='size cap as a multiple of 298/n')
ap.add_argument('--subset',type=int,default=60,help='pad every group to at least this many targets')
ap.add_argument('--min-keep',type=int,default=20); ap.add_argument('--columns',default='')
a=ap.parse_args()
out=pathlib.Path(a.out); out.mkdir(parents=True,exist_ok=True); say=RI.logger(out/'log.txt')
fleet=RI.IFleet(a.src); RI.eph()
say(f'source: {fleet.summary()}')
names=sorted(fleet.routes)
with mp.get_context('fork').Pool(a.nproc) as pool:
    res=pool.map(RI.w_cands,[(n,fleet.routes[n]['st'],ALL,1.5) for n in names])
dist={}
for n,lst in zip(names,res):
    d={}
    for c in lst: d[c['ast']]=min(d.get(c['ast'],9e9),c['dist'])
    dist[n]=d
grp={n:set(int(x) for x in fleet.routes[n]['st']['asts']) for n in names}
def aff(x,y):                       # mean distance of y's targets to x's trajectory (and vice versa)
    dx=[dist[x].get(t,2.0) for t in grp[y]]; dy=[dist[y].get(t,2.0) for t in grp[x]]
    return 0.5*(np.mean(dx)+np.mean(dy))
TGT=len(ALL)/a.n; CAP=a.cap*TGT
while len(grp)>a.n:
    # balanced agglomeration: the SMALLEST group joins its most compatible neighbour that stays under the size cap
    y=min(sorted(grp), key=lambda n: len(grp[n]))
    opts=[(aff(x,y),x) for x in sorted(grp) if x!=y and len(grp[x])+len(grp[y])<=CAP]
    if not opts:
        opts=[(aff(x,y)*(1+len(grp[x])/TGT),x) for x in sorted(grp) if x!=y]
    opts.sort(); _,x=opts[0]
    grp[x]=grp[x]|grp[y]; dist[x]={t:min(dist[x].get(t,9e9),dist[y].get(t,9e9)) for t in set(dist[x])|set(dist[y])}
    del grp[y]; del dist[y]
    say(f'merged {y} into {x}: now {len(grp[x])} targets ({len(grp)} groups)')
say('groups: '+' '.join(f'{n}:{len(g)}' for n,g in sorted(grp.items())))
json.dump({n:sorted(g) for n,g in grp.items()},open(out/'groups.json','w'))
wf=CM.w_fuel(0.5*(a.m0-640),a.m0)
P=BP(beam=a.beam,w_fuel=wf,m_margin=40.0,dv_max=1.2,lin_tofs=np.arange(20,600+1e-9,10)*DAY,
     lin_drmax=0.15*AU,vinf_cap=a.vinf,tof_refine=True)
cf=open(a.columns,'a') if a.columns else None
newf=RI.IFleet()
for k,(n,g) in enumerate(sorted(grp.items()),1):
    sub=set(g); want=max(len(g)+a.pad, a.subset)
    for t in sorted(dist[n], key=lambda t: dist[n][t]):     # pad with the nearest outsiders
        if len(sub)>=want: break
        sub.add(t)
    prize=np.zeros(300)
    for t in sub: prize[t-1]=1.0
    P.prize=prize; P.collect=DeepCollect(a.min_keep); tic=time.time()
    best,beam=beam_search(RI.eph(),excluded=sorted(set(ALL)-sub),m0=a.m0,
                          t_launch_grid=np.arange(0,800+1e-9,20)*DAY,P=P,n_proc=a.nproc,verbose=False)
    states=list(P.collect)+list(beam)
    uniq={}
    for s in states:
        key=frozenset(x[0] for x in s.seq)
        if key not in uniq or s.fuel<uniq[key].fuel: uniq[key]=s
    cand=sorted(uniq.values(),key=lambda s: float(CM.cost(s.fuel,s.m0))-len(s.seq))
    if cf:
        for s in cand[:300]:
            r=tour_json(s); r['tag']=f'merge{n}'; r['targets']=sorted(x[0] for x in s.seq)
            r['cost_planner']=float(CM.cost(s.fuel,s.m0)); cf.write(json.dumps(r)+'\n')
        cf.flush()
    took=None
    for s in cand[:4]:
        tg=sorted(x[0] for x in s.seq)
        try:
            ip,miss,lag=settle_tour(RI.eph(),tour_json(s),lambda ip: RI.settle(ip,100))
        except Exception as e:
            say(f'   {len(tg)}: settle raised {type(e).__name__}'); continue
        if ip is None or miss>150.0: say(f'   {len(tg)}: settle failed ({miss:.0f} km)'); continue
        tank=float(ip.tank())
        say(f'   {len(tg)} flybys: twin tank {tank:.0f} (J_i {RI.cost(tank):.3f}, {ip.dv()/len(tg):.3f}/flyby)')
        if took is None or RI.cost(tank)-len(tg) < took[3]: took=(ip,tank,tg,RI.cost(tank)-len(tg))
    if took is None: say(f'group {n}: no settle'); continue
    ip,tank,tg,_=took
    newf.routes[f'{k:02d}']=dict(st=RI.ist(ip),tank=tank); newf.save(out,note=f'group {n}')
    say(f'group {n}: {len(g)} targets (+{a.pad} pad) -> route {len(tg)} flybys, tank {tank:.0f}; '
        f'fleet {newf.summary()}  ({time.time()-tic:.0f} s)')
left=sorted(set(ALL)-set(newf.coverage()))
say(f'final: {newf.summary()}'); say(f'uncovered ({len(left)}): {left}')
json.dump(dict(uncovered=left,covered=sorted(newf.coverage()),J=newf.J()),open(out/'result.json','w'))
