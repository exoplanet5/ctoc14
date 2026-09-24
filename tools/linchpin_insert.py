"""Stage-13 linchpin, step 5: the INSERTION contrast on identical targets.  Insert the extra pool targets of one survivor
(pools.json extra[route]) into its ORIGINAL settled route, greedy cheapest-first with full re-evaluation every round:
each round, every remaining target is tried at its --cands closest approaches to the CURRENT route (local minima,
run_ialns.w_cands, up to --dmax AU) with run_ialns.w_insert (aim-point homotopy + full settle, miss <= 150 km); the
cheapest successful insertion is applied.  A target with no successful candidate in --fails consecutive rounds is an
orphan.  Writes <out>/result.json (per-step kg, J, depth) and route_<name>.npz (final route).
Usage: linchpin_insert.py src_fleet pools.json route out_dir [--nproc 3] [--cands 3] [--dmax 0.6]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from ctoc14.constants import DAY


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('src'); ap.add_argument('pools'); ap.add_argument('name')
    ap.add_argument('out'); ap.add_argument('--nproc', type=int, default=3); ap.add_argument('--cands', type=int, default=3)
    ap.add_argument('--dmax', type=float, default=0.6); ap.add_argument('--fails', type=int, default=2)
    ap.add_argument('--targets', default='', help='comma list overriding pools.json extra[route]')
    a = ap.parse_args()
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt')
    RI.eph()
    F = RI.IFleet(a.src); st = F.routes[a.name]['st']; tank0 = tank = F.routes[a.name]['tank']
    n0 = len(st['asts'])
    todo = ([int(x) for x in a.targets.split(',')] if a.targets else
            [int(e[0]) for e in json.load(open(a.pools))['extra'][a.name]])
    say(f'=== insertion contrast {a.name}: {n0} fb @ {tank0:.1f} kg; inserting {len(todo)} targets {todo}')
    steps = []; fails = {x: 0 for x in todo}; orphans = []; t_all = time.time(); ntrial = 0
    with mp.get_context('fork').Pool(a.nproc, maxtasksperchild=20) as pool:
        while todo:
            cj = pool.map(RI.w_cands, [(a.name, st, [x], a.dmax) for x in todo])
            jobs = []
            for x, lst in zip(todo, cj):
                for c in sorted(lst, key=lambda c: c['dist'])[:a.cands]:
                    jobs.append((a.name, st, x, c['t']))
            dist = {(x, round(c['t'])): c['dist'] for x, lst in zip(todo, cj) for c in lst}
            t0 = time.time(); res = pool.map(RI.w_insert, jobs) if jobs else []; ntrial += len(jobs)
            ok = [r for r in res if r['ok']]
            by = {}
            for r in ok:
                if r['ast'] not in by or r['tank'] < by[r['ast']]['tank']: by[r['ast']] = r
            say(f'  round {len(steps) + len(orphans)}: {len(todo)} targets, {len(jobs)} trials, {len(ok)} settled '
                f'({time.time() - t0:.0f} s); per-target best dkg: ' +
                ', '.join(f'{x}:{by[x]["tank"] - tank:+.1f}' if x in by else f'{x}:FAIL' for x in todo))
            for x in todo:
                fails[x] = 0 if x in by else fails[x] + 1
            for x in [x for x in todo if fails[x] >= a.fails]:
                orphans.append(x); todo.remove(x); say(f'  {x}: no successful insertion in {a.fails} rounds -> ORPHAN')
            if not by:
                continue
            r = min(by.values(), key=lambda r: r['tank'])
            d = dist.get((r['ast'], round(r['t'])), float('nan'))
            steps.append(dict(ast=r['ast'], t_d=round(r['t'] / DAY, 1), dist_au=round(d, 4), dkg=round(r['tank'] - tank, 2),
                              dJ=round(RI.cost(r['tank']) - RI.cost(tank), 5), depth_before=len(st['asts']),
                              tank_after=round(r['tank'], 2), miss_km=round(r['miss'], 1)))
            say(f'  INSERT {r["ast"]} at {r["t"] / DAY:.0f} d (closest approach {d:.3f} AU): {tank:.1f} -> {r["tank"]:.1f} kg '
                f'(+{r["tank"] - tank:.1f} kg, dJ {RI.cost(r["tank"]) - RI.cost(tank):+.4f}), depth {len(st["asts"]) + 1}')
            st = r['st']; tank = r['tank']; todo.remove(r['ast'])
            np.savez(out / f'route_{a.name}.npz', **st)
    n = len(st['asts'])
    res = dict(name=a.name, source_fb=n0, source_tank=tank0, source_kg_per_fb=round((tank0 - 600) / n0, 3),
               inserted=len(steps), orphans=orphans, final_fb=n, final_tank=round(tank, 2),
               final_kg_per_fb=round((tank - 600) / n, 3), dkg_total=round(tank - tank0, 2),
               dkg_per_inserted=round((tank - tank0) / max(len(steps), 1), 2),
               dJ_total=round(RI.cost(tank) - RI.cost(tank0), 5), trials=ntrial, t_s=round(time.time() - t_all), steps=steps)
    json.dump(res, open(out / 'result.json', 'w'), indent=1)
    say(f'  DONE {a.name}: {n0} fb @ {tank0:.1f} -> {n} fb @ {tank:.1f} kg; +{len(steps)} inserted for {tank - tank0:+.1f} kg '
        f'({res["dkg_per_inserted"]:.1f} kg each), orphans {orphans}; kg/fb {res["source_kg_per_fb"]:.2f} -> '
        f'{res["final_kg_per_fb"]:.2f}; {ntrial} trials, {res["t_s"]} s')


if __name__ == '__main__':
    main()
