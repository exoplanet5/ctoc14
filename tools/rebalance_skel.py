"""Rebalance the hard-target load of a skeleton fleet: move flybys from the routes with more than --max-hard targets
to the routes with the fewest, cheapest twin cost first (candidates = closest approaches, run_ialns.candidates;
insertion by insert_homotopy, removal by remove_flybys, both re-settled).  Balance matters for the windowed fill:
a skeleton with 12 waypoints fills to ~23 flybys, one with 6-7 to 33-34 (docs/stage5_skeleton_fill.md).
Usage: rebalance_skel.py src_dir dst_dir [--max-hard 9] [--dmax 0.5] [--per-target 4] [--nproc 8]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, time, pathlib, argparse, multiprocessing as mp
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1])); sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import numpy as np
import run_ialns as RI

ap = argparse.ArgumentParser(); ap.add_argument('src'); ap.add_argument('dst')
ap.add_argument('--max-hard', type=int, default=9); ap.add_argument('--dmax', type=float, default=0.5)
ap.add_argument('--per-target', type=int, default=4); ap.add_argument('--nproc', type=int, default=8)
ap.add_argument('--max-moves', type=int, default=12)
a = ap.parse_args()
out = pathlib.Path(a.dst); out.mkdir(parents=True, exist_ok=True); say = RI.logger(out.with_suffix('.log'))
fleet = RI.IFleet(a.src); RI.eph()
nfb = lambda n: len(fleet.routes[n]['st']['asts'])
say(f'start: {fleet.summary()}')
with mp.get_context('fork').Pool(a.nproc) as pool:
    for mv in range(a.max_moves):
        src = max(fleet.routes, key=nfb)
        if nfb(src) <= a.max_hard:
            say('balanced'); break
        light = [n for n in fleet.routes if nfb(n) < a.max_hard]
        targets = [int(t) for t in fleet.routes[src]['st']['asts']]
        by = RI.candidates(fleet, pool, targets, a.dmax, a.per_target, exclude=[n for n in fleet.routes if n not in light])
        cands = [c for lst in by.values() for c in lst]
        if not cands:
            say(f'{src}: no candidates on {light}'); break
        tic = time.time()
        res = [r for r in pool.map(RI.w_insert, [(c['host'], fleet.routes[c['host']]['st'], c['ast'], c['t']) for c in cands]) if r['ok']]
        if not res:
            say(f'{src}: {len(cands)} insertion trials, none feasible'); break
        for r in res: r['dtank'] = r['tank'] - fleet.routes[r['name']]['tank']
        res.sort(key=lambda r: r['dtank'])
        done = False
        for r in res[:4]:
            rm = pool.map(RI.w_remove, [(src, fleet.routes[src]['st'], [r['ast']])])[0]
            if not rm['ok']:
                continue
            saved = fleet.routes[src]['tank'] - rm['tank']
            fleet.routes[r['name']] = dict(st=r['st'], tank=r['tank']); fleet.routes[src] = dict(st=rm['st'], tank=rm['tank'])
            say(f'move {mv}: {r["ast"]} {src}->{r["name"]} (+{r["dtank"]:.0f} kg there, -{saved:.0f} kg here); '
                f'{fleet.summary()}  ({time.time()-tic:.0f} s)')
            fleet.save(out, note=f'move {mv}'); done = True; break
        if not done:
            say(f'{src}: removal failed for the best {min(4, len(res))} candidates'); break
fleet.save(out, note='rebalanced'); say(f'final: {fleet.summary()}')
