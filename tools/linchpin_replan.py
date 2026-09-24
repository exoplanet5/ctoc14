"""Stage-13 linchpin, step 3: RE-PLAN one route FROM SCRATCH, free choice inside a prescribed target pool, and settle
the result in the impulsive twin.

Pool (--cond):  own  = the route's own target set in the source fleet (control: can the builder re-fly the route?)
                pool = pools.json[route] (own targets + targets of dissolved routes, tools/linchpin_pools.py)
Every target outside the pool is EXCLUDED (never a leg end), so every flyby counts for the pool.

Builders (--builder):
  beam  ctoc14.search.beam_search from Earth, prize 1 per pool target (no skeleton, no windows): the purest re-plan.
  seg   tools/skel_fill.fill_segmented (the stage-5/6 builder): the pool's HARD targets (prizes_H >= 1.5) are waypoints
        with windows +-half around an epoch -- the source route's own flyby epoch for its own hard targets, the
        closest-approach epoch to the source trajectory (linchpin_dist.json) for a hard target it did not fly --
        and the easy pool targets are packed between them by waypoint-to-waypoint beams (K states kept per waypoint).
Production settings (skel_fill / greedy_cover): m0 1600, vinf 4, dvmax 1.2 (waypoint legs 2.5), Lambert 15-400 d +
linear 20-600 d legs, tof_refine, m_margin 40, launch grid 0-800 d / 20 d, w_fuel = CM.w_fuel(480, 1600).

The Pareto front (pool targets flown -> min planner fuel, 2 states per depth) is settled deepest first with
impulsive.settle_tour + run_ialns.settle (miss <= 150 km).  Writes <out>/result.json, cand_k<k>.npz per settled
candidate, and route_<name>.npz = the settled candidate with the best value J_i - lam * k.
Usage: linchpin_replan.py src_fleet pools.json route out_dir --cond pool --builder beam [--beam 1000] [--nproc 1]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, signal, pathlib, argparse, warnings
from argparse import Namespace
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from greedy_cover import BP, ALL, CM, DeepCollect, _timeout
from ctoc14.search import beam_search, tour_json
import ctoc14.impulsive as IM
from ctoc14.impulsive import settle_tour
from ctoc14.constants import DAY, AU
import skel_fill as SF
warnings.filterwarnings('ignore')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('src'); ap.add_argument('pools'); ap.add_argument('name')
    ap.add_argument('out')
    ap.add_argument('--cond', default='pool', choices=['own', 'pool']); ap.add_argument('--builder', default='beam',
                                                                                     choices=['beam', 'seg'])
    ap.add_argument('--beam', type=int, default=1000); ap.add_argument('--drmax', type=float, default=0.15)
    ap.add_argument('--half', type=float, default=45.0); ap.add_argument('--seg-k', type=int, default=50)
    ap.add_argument('--seg-depth', type=int, default=20); ap.add_argument('--max-depth', type=int, default=60)
    ap.add_argument('--m0', type=float, default=1600.0); ap.add_argument('--vinf', type=float, default=4.0)
    ap.add_argument('--dvmax', type=float, default=1.2); ap.add_argument('--wp-dvmax', type=float, default=2.5)
    ap.add_argument('--wp-prize', type=float, default=2.0); ap.add_argument('--nproc', type=int, default=1)
    ap.add_argument('--grid', default='0,800,20'); ap.add_argument('--min-keep', type=int, default=10)
    ap.add_argument('--try', dest='ntry', type=int, default=6, help='Pareto depths to settle (deepest first)')
    ap.add_argument('--per-depth', type=int, default=2, help='states tried per depth if the cheapest fails')
    ap.add_argument('--settle-timeout', type=float, default=300.0); ap.add_argument('--iters', type=int, default=100)
    ap.add_argument('--lam', type=float, default=0.15, help='J value of one pool target flown (route choice)')
    ap.add_argument('--hard', default=str(ROOT / 'results/n8/prizes_H.json'))
    ap.add_argument('--dist', default=str(ROOT / 'results/s13/linchpin/dist.json'))
    ap.add_argument('--lin-fallback', type=float, default=2.5)
    ap.add_argument('--soft-prize', type=float, default=-1.0, help='>= 0: SOFT pool -- targets outside the pool stay legal '
                    'leg ends (stepping stones) at this prize instead of being excluded; default strict')
    a = ap.parse_args()
    IM.LIN_FALLBACK = a.lin_fallback if a.lin_fallback > 0 else None
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt')
    E = RI.eph()
    src = RI.IFleet(a.src); st0 = src.routes[a.name]['st']
    own = [int(x) for x in st0['asts']]; own_t = {int(x): float(t) for x, t in zip(st0['asts'], st0['tf'])}
    tank0 = src.routes[a.name]['tank']
    pool = sorted(own) if a.cond == 'own' else sorted(int(x) for x in json.load(open(a.pools))['pools'][a.name])
    pset = set(pool); extra = sorted(pset - set(own))
    H = np.asarray(json.load(open(a.hard)), float); hard = set(t for t in ALL if H[t - 1] >= 1.5)
    say(f'=== {a.name} cond {a.cond} builder {a.builder}: pool {len(pool)} = own {len(pset & set(own))} + extra '
        f'{len(extra)} {extra}; hard in pool {len(pset & hard)}; source route {len(own)} fb @ {tank0:.1f} kg '
        f'({(tank0 - 600) / len(own):.2f} kg/fb)')
    wf = CM.w_fuel(0.5 * (a.m0 - 640.0), a.m0)
    g = [float(x) for x in a.grid.split(',')]; grid = np.arange(g[0], g[1] + 1e-9, g[2]) * DAY
    t0 = time.time()
    if a.builder == 'beam':
        prize = np.zeros(300)
        for t in pool: prize[t - 1] = 1.0
        soft = a.soft_prize >= 0.0
        if soft:
            for t in ALL:
                if t not in pset: prize[t - 1] = a.soft_prize
        P = BP(beam=a.beam, w_fuel=wf, m_margin=40.0, dv_max=a.dvmax, lin_tofs=np.arange(20, 600 + 1e-9, 10) * DAY,
               tofs=np.arange(15, 400 + 1e-9, 5) * DAY, lin_drmax=a.drmax * AU, vinf_cap=a.vinf, tof_refine=True,
               max_depth=a.max_depth)
        P.prize = prize; P.collect = DeepCollect(a.min_keep)
        best, beam = beam_search(E, excluded=[] if soft else sorted(set(ALL) - pset), m0=a.m0, t_launch_grid=grid, P=P,
                                 n_proc=a.nproc, verbose=False)
        states = list(P.collect) + list(beam)
        wps = []
    else:
        D = json.load(open(a.dist))['minima'][a.name]
        wps = []
        for x in sorted(pset & hard):
            if x in own_t:
                wps.append((own_t[x], x))
            else:
                lst = D[str(x)]
                if not lst:
                    say(f'  hard extra {x}: no approach to the source trajectory -- left as an easy pool target'); continue
                t, d = min(lst, key=lambda z: z[1]); wps.append((float(t), x))
                say(f'  hard extra {x}: waypoint at {t / DAY:.0f} d (closest approach {d:.3f} AU to the source route)')
        wps.sort()
        wpset = set(x for _, x in wps); easy = sorted(pset - wpset)
        prize = np.zeros(300)
        for t in easy: prize[t - 1] = 1.0
        for x in wpset: prize[x - 1] = a.wp_prize
        if a.soft_prize >= 0.0:
            stones = sorted(set(ALL) - pset)
            for t in stones: prize[t - 1] = a.soft_prize
            easy = sorted(set(easy) | set(stones))
        ns = Namespace(dvmax=a.dvmax, wp_dvmax=a.wp_dvmax, vinf=a.vinf, nproc=a.nproc, wp_prize=a.wp_prize, m0=a.m0)
        say(f'  skeleton: {len(wps)} waypoints ' + ' '.join(f'{x}@{t / DAY:.0f}' for t, x in wps) + f'; easy {len(easy)}')
        states = SF.fill_segmented(a.name, wps, easy, ns, wf, grid, say, a.half, a.beam, a.drmax, K=a.seg_k,
                                   seg_depth=a.seg_depth, prize=prize)
    t_build = time.time() - t0
    npool = lambda s: sum(1 for x in s.seq if x[0] in pset)
    byk = {}
    for s in states:
        k = npool(s); byk.setdefault(k, []).append(s)
    for k in byk:
        uniq = {}
        for s in sorted(byk[k], key=lambda s: s.fuel):
            key = frozenset(x[0] for x in s.seq)
            if key not in uniq: uniq[key] = s
        byk[k] = list(uniq.values())[:a.per_depth]
    deep = max(byk) if byk else 0
    par = {int(k): round(float(CM.tank(byk[k][0].fuel, a.m0)), 1) for k in sorted(byk)}
    say(f'  build {t_build:.0f} s: {len(states)} states, deepest {deep} of pool {len(pool)}; Pareto ' +
        ' '.join(f'{k}:{v:.0f}' for k, v in par.items() if k >= deep - 10))
    res = dict(name=a.name, cond=a.cond, builder=a.builder, soft_prize=a.soft_prize, pool=pool, pool_size=len(pool), own=sorted(own),
               extra=extra, source_fb=len(own), source_tank=tank0, beam=a.beam, drmax=a.drmax, half=a.half,
               seg_k=a.seg_k, t_build_s=round(t_build, 1), n_states=len(states), deepest=deep, pareto_planner=par,
               waypoints=[[x, t] for t, x in wps], settled=[])
    json.dump(res, open(out / 'result.json', 'w'), indent=1)
    best = None
    for k in sorted(byk, reverse=True)[:a.ntry]:
        for j, s in enumerate(byk[k]):
            tour = tour_json(s); t1 = time.time()
            rec = dict(k=int(k), n=len(s.seq), alt=j, tank_planner=round(float(CM.tank(s.fuel, a.m0)), 1))
            try:
                signal.signal(signal.SIGALRM, _timeout); signal.setitimer(signal.ITIMER_REAL, a.settle_timeout)
                ip, miss, lag = settle_tour(E, tour, lambda ip: RI.settle(ip, a.iters))
            except Exception as e:
                ip, miss = None, float('nan'); rec['error'] = type(e).__name__
            finally:
                signal.setitimer(signal.ITIMER_REAL, 0)
            rec['t_settle_s'] = round(time.time() - t1, 1)
            if ip is None or not (miss <= 150.0):
                rec.update(ok=False, miss_km=float(miss) if np.isfinite(miss) else None)
                say(f'    {k} pool fb (alt {j}): settle FAILED (miss {miss}, {rec["t_settle_s"]:.0f} s)')
                res['settled'].append(rec); continue
            tank = float(ip.tank()); n = len(ip.asts); tg = sorted(int(x) for x in ip.asts)
            rec.update(ok=True, miss_km=round(float(miss), 1), tank=round(tank, 2), J_i=round(RI.cost(tank), 5),
                       kg_per_fb=round((tank - 600.0) / n, 3), kg_per_pool=round((tank - 600.0) / max(k, 1), 3),
                       stones=int(n - k), dv_per_fb=round(ip.dv() / n, 4),
                       J_per_fb=round(RI.cost(tank) / n, 5), t_launch_d=round(float(ip.tL) / DAY, 1),
                       t_end_d=round(float(np.max(ip.tf)) / DAY, 1), targets=tg,
                       own_flown=len(set(tg) & set(own)), extra_flown=len(set(tg) & set(extra)),
                       value=round(RI.cost(tank) - a.lam * k, 5))
            np.savez(out / f'cand_k{k}_a{j}.npz', **RI.ist(ip))
            say(f'    {k} pool fb (alt {j}): planner {rec["tank_planner"]:.0f} -> twin {tank:.1f} kg, '
                f'{rec["kg_per_fb"]:.2f} kg/fb, {rec["dv_per_fb"]:.3f} km/s/fb, J_i {rec["J_i"]:.4f}, own '
                f'{rec["own_flown"]}/{len(set(own) & pset)} extra {rec["extra_flown"]}/{len(extra)} '
                f'(miss {miss:.0f} km, {rec["t_settle_s"]:.0f} s)')
            res['settled'].append(rec)
            if best is None or rec['value'] < best[0]['value']:
                best = (rec, ip)
            break                                    # a settled state at this depth: go to the next depth
        json.dump(res, open(out / 'result.json', 'w'), indent=1)
    if best is not None:
        rec, ip = best
        np.savez(out / f'route_{a.name}.npz', **RI.ist(ip))
        res['best'] = {k: v for k, v in rec.items() if k != 'targets'}
        missing = sorted(pset - set(rec['targets']))
        res['best']['orphans'] = missing
        say(f'  BEST {a.name}: {rec["k"]}/{len(pool)} pool targets @ {rec["tank"]:.1f} kg ({rec["kg_per_fb"]:.2f} kg/fb; '
            f'source {len(own)} @ {tank0:.1f} = {(tank0 - 600) / len(own):.2f} kg/fb); orphans {missing}')
    else:
        say(f'  {a.name}: NOTHING settled')
    res['t_total_s'] = round(time.time() - t0, 1)
    json.dump(res, open(out / 'result.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
