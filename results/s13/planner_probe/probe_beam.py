"""Stage-13 planner probe: ONE route planned freely by ctoc14.search.beam_search inside a prescribed target pool,
with the production settings of tools/skel_fill.py / tools/greedy_cover.py (m0 1600, vinf 4, dvmax 1.2, Lambert
15-400 d + linear 20-600 d legs, lin_drmax 0.15 AU, tof_refine, m_margin 40, launch grid 0-800 d / 20 d,
w_fuel = CostModel(0.70, 2).w_fuel(480, 1600)).  The Pareto front (flybys -> min planner fuel) is settled into the
impulsive twin (impulsive.settle_tour + run_ialns.settle), deepest first.

  --pool all | routes:00,01 (target sets of results/s12/cover_clean/fleet = t10d) | file.json (list of ids)
  --mode strict : every non-pool target is EXCLUDED (never a leg end)
         soft   : non-pool targets stay legal leg ends at --soft-prize (0 = stepping stones worth nothing)
Writes <out>.json.  Usage: probe_beam.py out_stem --pool routes:08,09 --mode strict --beam 300 --nproc 2 --settle 3
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, signal, pathlib, argparse, warnings, resource
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from greedy_cover import BP, ALL, CM, DeepCollect, _timeout
from ctoc14.search import beam_search, tour_json
import ctoc14.search as S
import ctoc14.impulsive as IM
from ctoc14.impulsive import settle_tour
from ctoc14.constants import DAY, AU
warnings.filterwarnings('ignore')

ap = argparse.ArgumentParser(); ap.add_argument('out')
ap.add_argument('--pool', default='all'); ap.add_argument('--mode', default='strict', choices=['strict', 'soft'])
ap.add_argument('--soft-prize', type=float, default=0.0)
ap.add_argument('--beam', type=int, default=300); ap.add_argument('--nproc', type=int, default=2)
ap.add_argument('--m0', type=float, default=1600.0); ap.add_argument('--dvmax', type=float, default=1.2)
ap.add_argument('--drmax', type=float, default=0.15); ap.add_argument('--max-depth', type=int, default=60)
ap.add_argument('--grid-max', type=float, default=800.0); ap.add_argument('--settle', type=int, default=3)
ap.add_argument('--settle-timeout', type=float, default=240.0); ap.add_argument('--min-keep', type=int, default=10)
a = ap.parse_args()
IM.LIN_FALLBACK = 2.5

fleet_dir = ROOT / 'results/s12/cover_clean/fleet'
if a.pool == 'all':
    pool = list(ALL)
elif a.pool.startswith('routes:'):
    F = RI.IFleet(fleet_dir); pool = sorted(set().union(*[set(int(x) for x in F.routes[n]['st']['asts']) for n in a.pool[7:].split(',')]))
else:
    pool = sorted(set(json.load(open(a.pool))))
pset = set(pool)
prize = np.zeros(300)
for t in ALL:
    prize[t - 1] = 1.0 if t in pset else a.soft_prize
excluded = sorted(set(ALL) - pset) if a.mode == 'strict' else []

# count expansions
_n_exp = [0]; _orig = S.expand
def _cnt(eph, s, P):
    _n_exp[0] += 1; return _orig(eph, s, P)
if a.nproc == 1:
    S.expand = _cnt

wf = CM.w_fuel(0.5 * (a.m0 - 640.0), a.m0); grid = np.arange(0, a.grid_max + 1e-9, 20) * DAY
P = BP(beam=a.beam, w_fuel=wf, m_margin=40.0, dv_max=a.dvmax, lin_tofs=np.arange(20, 600 + 1e-9, 10) * DAY,
       tofs=np.arange(15, 400 + 1e-9, 5) * DAY, lin_drmax=a.drmax * AU, vinf_cap=4.0, tof_refine=True,
       max_depth=a.max_depth)
P.prize = prize; P.collect = DeepCollect(a.min_keep)
t0 = time.time(); _r0 = resource.getrusage(resource.RUSAGE_SELF); c0 = _r0.ru_utime + _r0.ru_stime
best, beam = beam_search(RI.eph(), excluded=excluded, m0=a.m0, t_launch_grid=grid, P=P, n_proc=a.nproc, verbose=False)
t_beam = time.time() - t0
_ru = resource.getrusage(resource.RUSAGE_SELF); _rc = resource.getrusage(resource.RUSAGE_CHILDREN)
cpu_beam = (_ru.ru_utime + _ru.ru_stime + _rc.ru_utime + _rc.ru_stime) - c0
states = list(P.collect) + list(beam)
npool = lambda s: sum(1 for x in s.seq if x[0] in pset)
byk = {}; bytot = {}
for s in states:
    k = npool(s)
    if k not in byk or s.fuel < byk[k].fuel: byk[k] = s
    n = len(s.seq)
    if n not in bytot or s.fuel < bytot[n].fuel: bytot[n] = s
deep = max(len(s.seq) for s in states) if states else 0
res = dict(pool=a.pool, pool_size=len(pool), mode=a.mode, soft_prize=a.soft_prize, beam=a.beam, nproc=a.nproc, m0=a.m0,
           t_beam_s=round(t_beam, 1), cpu_beam_s=round(cpu_beam, 1), n_expand=_n_exp[0] or None, n_states=len(states), deepest_total=deep,
           deepest_pool=max(byk) if byk else 0,
           pareto_pool={int(k): dict(n_tot=len(byk[k].seq), fuel=round(byk[k].fuel, 1), tank_planner=round(float(CM.tank(byk[k].fuel, a.m0)), 1))
                        for k in sorted(byk)},
           settles=[])
print(f'[{a.out}] beam {a.beam} nproc {a.nproc} pool {len(pool)} {a.mode}: {t_beam:.0f} s wall / {cpu_beam:.0f} s cpu, deepest {deep} total / '
      f'{res["deepest_pool"]} pool; Pareto(pool) ' + ' '.join(f'{k}:{v["tank_planner"]:.0f}' for k, v in res['pareto_pool'].items() if k >= deep - 8), flush=True)
for k in sorted(byk, reverse=True)[:a.settle]:
    s = byk[k]; tour = tour_json(s); t1 = time.time()
    rec = dict(k_pool=int(k), n_tot=len(s.seq), tank_planner=round(float(CM.tank(s.fuel, a.m0)), 1))
    try:
        signal.signal(signal.SIGALRM, _timeout); signal.setitimer(signal.ITIMER_REAL, a.settle_timeout)
        ip, miss, lag = settle_tour(RI.eph(), tour, lambda ip: RI.settle(ip, 100))
    except Exception as e:
        ip, miss = None, float('nan'); rec['error'] = type(e).__name__
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    rec['t_settle_s'] = round(time.time() - t1, 1)
    if ip is None or not (miss <= 150.0):
        rec.update(ok=False, miss_km=float(miss) if np.isfinite(miss) else None)
    else:
        tank = float(ip.tank()); n = len(ip.asts)
        rec.update(ok=True, miss_km=round(float(miss), 1), lag_d=round(lag / DAY, 2), tank_twin=round(tank, 1), J_i=round(RI.cost(tank), 4),
                   kg_per_fb=round((tank - 600.0) / n, 2), dv_per_fb=round(ip.dv() / n, 4),
                   t_end_yr=round(float(np.max(ip.tf)) / DAY / 365.25, 2), t_launch_d=round(float(ip.tL) / DAY, 0))
        np.savez(f'{a.out}_k{k}.npz', **RI.ist(ip))
    res['settles'].append(rec)
    print(f'   settle {rec}', flush=True)
json.dump(res, open(f'{a.out}.json', 'w'), indent=1)
