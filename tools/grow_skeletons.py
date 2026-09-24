"""Grow a fleet of skeleton routes (impulsive twin states, route_*.npz) into a full cover by cheapest insertion.
1. De-duplicate: a target on several routes stays where removing it saves the least (w_remove trials elsewhere).
2. Phase 1: insert the uncovered HARD targets (given list) by cheapest tank increase over all (route, event) candidates;
   phase 2: the remaining targets. One accepted insertion per route per round (the other trials of that route are
   stale), candidates = closest approaches of the target to each route (run_ialns.candidates, --dmax AU).
Fleet is saved after every round (run_ialns.IFleet format) so `run_ialns.py search/export` can take over.
Usage: grow_skeletons.py skel_dir out_dir --hard 8,26,... [--dmax 0.3] [--per-target 3] [--nproc 8] [--max-trials 240]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse, multiprocessing as mp
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1])); sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import numpy as np
import run_ialns as RI
from ctoc14.constants import VE
from ctoc14.search import UNREACHABLE

ap = argparse.ArgumentParser(); ap.add_argument('src'); ap.add_argument('dst'); ap.add_argument('--hard', default='')
ap.add_argument('--dmax', type=float, default=0.3); ap.add_argument('--dmax-wide', type=float, default=0.6)
ap.add_argument('--per-target', type=int, default=3); ap.add_argument('--nproc', type=int, default=8)
ap.add_argument('--max-trials', type=int, default=240); ap.add_argument('--max-rounds', type=int, default=40)
ap.add_argument('--max-tank', type=float, default=1150.0, help='do not grow a route past this tank')
ap.add_argument('--max-dtank', type=float, default=1e9, help='reject an insertion that costs more than this (kg)')
ap.add_argument('--hard-only', action='store_true', help='stop after the hard phase (skeleton merge: no easy growth)')
a = ap.parse_args()
out = pathlib.Path(a.dst); out.mkdir(parents=True, exist_ok=True); say = RI.logger(out.with_suffix('.log'))
fleet = RI.IFleet(a.src)
hard = set(int(x) for x in a.hard.split(',') if x)
say(f'skeletons: {fleet.summary()}')
RI.eph()


def dedupe(pool):
    cov = fleet.coverage()
    dup = {t: sorted(v) for t, v in cov.items() if len(v) > 1}
    if not dup: return
    jobs = [(n, fleet.routes[n]['st'], [t]) for t, names in dup.items() for n in names]
    res = [r for r in pool.map(RI.w_remove, jobs) if r['ok']]
    # for each duplicate target: remove it from every route except the one whose removal saves the least
    by = {}
    for r in res: by.setdefault(r['asts'][0], []).append(r)
    plan = {}
    for t, lst in by.items():
        lst.sort(key=lambda r: fleet.routes[r['name']]['tank'] - r['tank'])       # saving ascending
        keep = lst[0]['name']
        for r in lst[1:]:
            plan.setdefault(r['name'], []).append(t)
    jobs = [(n, fleet.routes[n]['st'], asts) for n, asts in plan.items()]
    for r in pool.map(RI.w_remove, jobs):
        if r['ok']:
            say(f'  dedupe {r["name"]}: removed {r["asts"]}: tank {fleet.routes[r["name"]]["tank"]:.1f} -> {r["tank"]:.1f}')
            fleet.routes[r['name']] = dict(st=r['st'], tank=r['tank'])
        else:
            say(f'  dedupe {r["name"]}: removal of {r["asts"]} FAILED (kept)')


def grow(pool, targets, label):
    dmax = a.dmax
    for rnd in range(a.max_rounds):
        cov = fleet.coverage(); todo = [t for t in targets if t not in cov]
        if not todo: return
        tic = time.time()
        by = RI.candidates(fleet, pool, todo, dmax, a.per_target)
        cands = [c for lst in by.values() for c in lst]
        if not cands:
            if dmax < a.dmax_wide:
                dmax = a.dmax_wide; say(f'  {label} round {rnd}: no candidates, widening to {dmax} AU'); continue
            say(f'  {label} round {rnd}: no candidates for {todo}'); return
        cands.sort(key=lambda c: c['dist']); cands = cands[:a.max_trials]
        jobs = [(c['host'], fleet.routes[c['host']]['st'], c['ast'], c['t']) for c in cands
                if fleet.routes[c['host']]['tank'] < a.max_tank]
        res = [r for r in pool.map(RI.w_insert, jobs) if r['ok']]
        for r in res: r['dtank'] = r['tank'] - fleet.routes[r['name']]['tank']
        res = [r for r in res if r['dtank'] <= a.max_dtank]
        res.sort(key=lambda r: r['dtank'])
        used_host = set(); used_ast = set(); acc = []
        for r in res:
            if r['name'] in used_host or r['ast'] in used_ast: continue
            used_host.add(r['name']); used_ast.add(r['ast']); acc.append(r)
            fleet.routes[r['name']] = dict(st=r['st'], tank=r['tank'])
        fleet.save(out, note=f'{label} round {rnd}')
        say(f'  {label} round {rnd}: {len(todo)} to place, {len(jobs)} trials, {len(res)} feasible, accepted {len(acc)}: ' +
            ' '.join(f'{r["ast"]}->{r["name"]}(+{r["dtank"]:.0f}kg)' for r in acc) + f'; J {fleet.J():.4f} ({time.time()-tic:.0f} s)')
        if not acc:
            if dmax < a.dmax_wide:
                dmax = a.dmax_wide; say(f'  {label}: widening to {dmax} AU'); continue
            say(f'  {label}: stuck with {len(todo)} unplaced: {todo}'); return


with mp.get_context('fork').Pool(a.nproc) as pool:
    dedupe(pool)
    fleet.save(out, note='deduped'); say(f'after dedupe: {fleet.summary()}')
    allt = [t for t in range(1, 301) if t not in UNREACHABLE]
    grow(pool, sorted(hard), 'hard')
    if not a.hard_only:
        grow(pool, allt, 'all')
fleet.save(out, note='grown'); say(f'final: {fleet.summary()}')
