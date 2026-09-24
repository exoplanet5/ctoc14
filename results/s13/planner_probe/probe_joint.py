"""Stage-13 joint-planner probe: K craft planned JOINTLY by ctoc14.jointsearch.joint_beam_search over one shared pool
(one visited mask, the assignment of targets to craft decided inside the beam), same leg model and settings as
probe_beam.py.  prize 1 on every pool target and P.w_rare = 1 (so the prize enters the joint score, as in
tools/run_lns.py pair moves); every non-pool target is EXCLUDED.  Reports the coverage/cost front over all collected
joint states, beam-diversity diagnostics, and settles every craft of the best state.
Usage: probe_joint.py out_stem --pool routes:08,09 --K 2 --beam 100 --nproc 1"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, signal, pathlib, argparse, warnings
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from greedy_cover import BP, ALL, CM, _timeout
from ctoc14.search import tour_json
import ctoc14.jointsearch as JS
import ctoc14.impulsive as IM
from ctoc14.impulsive import settle_tour
from ctoc14.constants import DAY, AU
warnings.filterwarnings('ignore')

ap = argparse.ArgumentParser(); ap.add_argument('out')
ap.add_argument('--pool', default='routes:08,09'); ap.add_argument('--K', type=int, default=2)
ap.add_argument('--beam', type=int, default=100); ap.add_argument('--nproc', type=int, default=1)
ap.add_argument('--m0', type=float, default=1600.0); ap.add_argument('--dvmax', type=float, default=1.2)
ap.add_argument('--launch-max', type=float, default=800.0); ap.add_argument('--wait', type=float, default=60.0)
ap.add_argument('--wrare', type=float, default=1.0); ap.add_argument('--parent-cap', type=int, default=0)
ap.add_argument('--settle', action='store_true'); ap.add_argument('--canon', action='store_true')
a = ap.parse_args()
IM.LIN_FALLBACK = 2.5
F = RI.IFleet(ROOT / 'results/s12/cover_clean/fleet')
names = a.pool[7:].split(',') if a.pool.startswith('routes:') else []
if a.pool == 'all':
    pool = list(ALL)
elif names:
    pool = sorted(set().union(*[set(int(x) for x in F.routes[n]['st']['asts']) for n in names]))
else:
    pool = sorted(set(json.load(open(a.pool))))
ref = {n: (len(F.routes[n]['st']['asts']), round(F.routes[n]['tank'], 1), round(RI.cost(F.routes[n]['tank']), 4)) for n in names}
pset = set(pool); prize = np.zeros(300)
for t in pool: prize[t - 1] = 1.0
wf = CM.w_fuel(0.5 * (a.m0 - 640.0), a.m0)
P = BP(beam=a.beam, w_fuel=wf, m_margin=40.0, dv_max=a.dvmax, lin_tofs=np.arange(20, 600 + 1e-9, 10) * DAY,
       tofs=np.arange(15, 400 + 1e-9, 5) * DAY, lin_drmax=0.15 * AU, vinf_cap=4.0, tof_refine=True, prize=prize)
P.w_rare = a.wrare; P.wait = a.wait * DAY; P.parent_cap = a.parent_cap; P.collect_joint = []
t0 = time.time()
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent)); import jointsym
STATS = []
best = jointsym.joint_beam_search(RI.eph(), [a.m0] * a.K, np.arange(0, a.launch_max + 1e-9, 20) * DAY, P=P,
                                  excluded=sorted(set(ALL) - pset), n_proc=a.nproc, verbose=True, canon=a.canon, stats=STATS)
t_run = time.time() - t0
cands = P.collect_joint + [best]
def ev(j):
    rs = [s for s in j.craft if s is not None and s.n() > 0]
    cov = set().union(*[{x[0] for x in s.seq} for s in rs]) if rs else set()
    tanks = [float(CM.tank(s.fuel, s.m0)) for s in rs]
    return len(cov), sum(RI.cost(t) for t in tanks), rs, tanks
front = {}
for j in cands:
    c, J, rs, tanks = ev(j)
    if c not in front or J < front[c][0]: front[c] = (J, [s.n() for s in rs], [round(t) for t in tanks], j)
cmax = max(front)
# beam diversity at the deepest level collected
depth_of = lambda j: j.n()
dmax = max(depth_of(j) for j in P.collect_joint) if P.collect_joint else 0
last = [j for j in P.collect_joint if depth_of(j) == dmax]
splits = set(tuple(s.n() if s is not None else 0 for s in j.craft) for j in last)
per_craft = [len(set(frozenset(x[0] for x in j.craft[i].seq) if j.craft[i] is not None else frozenset() for j in last)) for i in range(a.K)]
res = dict(pool=a.pool, pool_size=len(pool), K=a.K, beam=a.beam, nproc=a.nproc, m0=a.m0, t_run_s=round(t_run, 1),
           n_joint_states=len(cands), t10d_reference=ref, t10d_sumJ=round(sum(v[2] for v in ref.values()), 4),
           max_covered=cmax, front={int(c): dict(sumJ_planner=round(v[0], 4), depths=v[1], tanks_planner=v[2]) for c, v in sorted(front.items()) if c >= cmax - 6},
           last_level=dict(depth=dmax, n_states=len(last), distinct_splits=len(splits), distinct_sets_per_craft=per_craft),
           settles=[], canon=a.canon, wait_d=a.wait,
           perm_free_frac=dict(mean=round(float(np.mean([x['perm_free'] / x['kept'] for x in STATS])), 3),
                               min=round(float(np.min([x['perm_free'] / x['kept'] for x in STATS])), 3)),
           levels=STATS)
print(json.dumps({k: v for k, v in res.items() if k not in ('front', 'levels')}), flush=True)
for c in sorted(front, reverse=True)[:6]:
    print(f'  covered {c}: planner sum J_i {front[c][0]:.3f}, depths {front[c][1]}, tanks {front[c][2]}', flush=True)
J, depths, tanks, j = front[cmax]
tours = [tour_json(s) for s in j.craft if s is not None and s.n() > 0]
json.dump(tours, open(f'{a.out}_tours.json', 'w'), indent=1)
res['waits_per_craft'] = [sum(1 for k in range(1, len(t['legs'])) if t['legs'][k]['t_flyby'] - t['legs'][k-1]['t_flyby'] > t['legs'][k]['tof'] + DAY) for t in tours]
print('  waits per craft', res['waits_per_craft'], flush=True)
if a.settle:
    for i, s in enumerate([s for s in j.craft if s is not None and s.n() > 0]):
        tour = tour_json(s); t1 = time.time(); rec = dict(craft=i, n=s.n(), tank_planner=round(float(CM.tank(s.fuel, s.m0)), 1))
        try:
            signal.signal(signal.SIGALRM, _timeout); signal.setitimer(signal.ITIMER_REAL, 240.0)
            ip, miss, lag = settle_tour(RI.eph(), tour, lambda ip: RI.settle(ip, 100))
        except Exception as e:
            ip, miss = None, float('nan'); rec['error'] = type(e).__name__
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
        rec['t_settle_s'] = round(time.time() - t1, 1)
        if ip is not None and miss <= 150:
            tank = float(ip.tank()); rec.update(ok=True, tank_twin=round(tank, 1), J_i=round(RI.cost(tank), 4),
                                                kg_per_fb=round((tank - 600) / s.n(), 2), dv_per_fb=round(ip.dv() / s.n(), 4))
            np.savez(f'{a.out}_c{i}.npz', **RI.ist(ip))
        else:
            rec.update(ok=False, miss_km=float(miss) if np.isfinite(miss) else None)
        res['settles'].append(rec); print('  settle', rec, flush=True)
json.dump(res, open(f'{a.out}.json', 'w'), indent=1, default=str)
