"""Parallel regret-style growth of a whole carrier fleet (stage 7).

Sequential building (tools/carrier_fleet.py) repeats the old failure: craft 0 takes 40 targets, craft 7 gets 20, and
the last 50 targets cost 60-70 kg each.  Here all routes grow TOGETHER from their sparse carrier bases.  Each round:
closest approaches of every uncovered target to every route; targets are ranked by how FEW routes can see them and by
the gap between their best and second-best approach (regret proxy), so a target with one home is placed before that
home fills up with targets that had alternatives.  Half of each route's trials go to high-regret targets, half to its
nearest ones.  Accept per route the cheapest successes (<= --max-step kg, --sep days apart), one joint re-settle.

Usage: carrier_regret.py base_fleet out_dir [--per-route 6] [--max-step 45] [--max-tank 1150]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, time, pathlib, argparse, collections, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from ctoc14.search import UNREACHABLE
from ctoc14.constants import DAY
from carrier_grow import w_multi


def w_multi_named(job):
    n, st, items = job
    return n, w_multi((st, items))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('src'); ap.add_argument('out')
    ap.add_argument('--per-route', type=int, default=6); ap.add_argument('--naccept', type=int, default=3)
    ap.add_argument('--max-step', type=float, default=45.0); ap.add_argument('--step-hi', type=float, default=90.0)
    ap.add_argument('--max-tank', type=float, default=1150.0); ap.add_argument('--dmax', type=float, default=0.2)
    ap.add_argument('--sep', type=float, default=200.0); ap.add_argument('--gap', type=float, default=12.0)
    ap.add_argument('--pools', default='', help='json {route: [targets]}: each route sees only its pool until growth stalls')
    ap.add_argument('--nproc', type=int, default=8); ap.add_argument('--rounds', type=int, default=80)
    a = ap.parse_args()
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt')
    fleet = RI.IFleet(a.src); names = sorted(fleet.routes)
    say('start: ' + ' '.join(f'{n}:{len(fleet.routes[n]["st"]["asts"])}@{fleet.routes[n]["tank"]:.0f}' for n in names))
    import json
    pools = {n: set(v) for n, v in json.load(open(a.pools)).items()} if a.pools else None
    pool = mp.get_context('fork').Pool(a.nproc); t0 = time.time(); step = a.max_step
    for rnd in range(1, a.rounds + 1):
        cov = set(fleet.coverage()); left = [t for t in range(1, 301) if t not in cov and t not in UNREACHABLE]
        if not left: break
        res = pool.map(RI.w_cands, [(n, fleet.routes[n]['st'], left, a.dmax) for n in names])
        per = collections.defaultdict(dict)                      # target -> route -> best candidate
        for lst in res:
            for c in lst:
                tf = np.asarray(fleet.routes[c['host']]['st']['tf'])
                if np.abs(tf - c['t']).min() <= a.gap * DAY: continue
                if pools is not None and c['ast'] not in pools.get(c['host'], ()): continue
                if c['host'] not in per[c['ast']] or c['dist'] < per[c['ast']][c['host']]['dist']:
                    per[c['ast']][c['host']] = c
        def regret(t):
            d = sorted(c['dist'] for c in per[t].values())
            return (len(d), -(d[1] - d[0]) if len(d) > 1 else -9.0)
        jobs = []; used = collections.Counter()
        half = a.per_route // 2
        for t in sorted(per, key=regret):                        # high regret first, to its nearest route
            c = min(per[t].values(), key=lambda c: c['dist'])
            if used[c['host']] < half and fleet.routes[c['host']]['tank'] < a.max_tank - 10:
                jobs.append(c); used[c['host']] += 1
        tried = set((c['host'], c['ast']) for c in jobs)
        for n in names:                                          # then each route's nearest
            mine = sorted((c for t in per for h, c in per[t].items() if h == n and (n, t) not in tried), key=lambda c: c['dist'])
            for c in mine[:a.per_route - used[n]]: jobs.append(c)
        if not jobs and pools is not None:
            pools = None; say(f'round {rnd}: pools exhausted -> open growth'); continue
        if not jobs: say('no candidates'); break
        out_r = pool.map(RI.w_insert, [(c['host'], fleet.routes[c['host']]['st'], c['ast'], c['t']) for c in jobs], chunksize=1)
        ok = [r for r in out_r if r.get('ok') and r['tank'] <= a.max_tank and r['tank'] - fleet.routes[r['name']]['tank'] <= step]
        ok.sort(key=lambda r: r['tank'] - fleet.routes[r['name']]['tank'])
        picks = collections.defaultdict(list); taken = set()
        for r in ok:
            if r['ast'] in taken or len(picks[r['name']]) >= a.naccept: continue
            if all(abs(r['t'] - q['t']) > a.sep * DAY for q in picks[r['name']]):
                picks[r['name']].append(r); taken.add(r['ast'])
        if not picks and pools is not None:
            pools = None; say(f'round {rnd}: pools exhausted -> open growth'); continue
        if not picks:
            if step < a.step_hi: step = a.step_hi; say(f'round {rnd}: nothing under +{a.max_step:.0f} kg -> step limit {step:.0f}'); continue
            say(f'round {rnd}: nothing acceptable'); break
        multi = dict(pool.map(w_multi_named, [(n, fleet.routes[n]['st'], [(r['ast'], r['t']) for r in p])
                                              for n, p in picks.items() if len(p) > 1]))
        added = 0
        for n, p in picks.items():
            j = multi.get(n) if len(p) > 1 else None
            if j is not None and j[1] <= a.max_tank: st, tk = j; k = len(p)
            else: st, tk = p[0]['st'], p[0]['tank']; k = 1
            fleet.routes[n] = dict(st=st, tank=tk); added += k
        fleet.save(out / 'fleet', note=f'regret growth round {rnd}')
        say(f'round {rnd}: {len(ok)}/{len(jobs)} usable, +{added} -> covered {len(fleet.coverage())}, '
            f'sum J_i {sum(RI.cost(r["tank"]) for r in fleet.routes.values()):.3f} | '
            + ' '.join(f'{len(fleet.routes[n]["st"]["asts"])}@{fleet.routes[n]["tank"]:.0f}' for n in names) + f'  ({time.time()-t0:.0f} s)')
    left = sorted(set(range(1, 301)) - set(fleet.coverage()) - set(UNREACHABLE))
    say(f'done: covered {len(fleet.coverage())}, J {fleet.J():.3f}; left ({len(left)}): {left}')


if __name__ == '__main__':
    main()
