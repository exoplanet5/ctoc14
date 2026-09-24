"""Dissolve a craft by redistributing its unique targets (stage 7, section 10).

J = sum_i J_i + N_miss, so removing a craft saves its whole J_i (1.2-1.5) and the fleet must absorb its unique targets.
The 2026-09-16 campaign measured redistribution at 0.10-0.25 J per target against a budget of ~0.06, and gave up.  Two
things changed: routes are now DEEP (42-48 flybys) and grown by the same insertion operator, and insertion into a route
built by insertion is cheap (4-50 kg), so the budget per target is met if the hosts have room.

For each candidate victim: its unique targets are offered to every other route at their closest approaches (cheapest
host first, one target at a time, hosts re-settled as they grow); the trial is accepted if the fleet J drops.

Usage: fleet_dissolve.py fleet_dir out_dir [--victims all] [--max-tank 1400]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, time, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from ctoc14.constants import DAY
_A = None


def w_try(job):
    """insert target X into route n at time t; returns the J increase."""
    n, st, tank, X, t = job
    r = RI.w_insert((n, st, X, t))
    if not r.get('ok') or r['tank'] > _A.max_tank: return n, X, None
    return n, X, (r['st'], float(r['tank']), RI.cost(r['tank']) - RI.cost(tank))


def dissolve(fl, victim, pl, say):
    routes = {n: dict(st=r['st'], tank=r['tank']) for n, r in fl.routes.items() if n != victim}
    own = set(int(x) for x in fl.routes[victim]['st']['asts'])
    others = set(int(x) for n, r in routes.items() for x in r['st']['asts'])
    todo = sorted(own - others)
    saved = RI.cost(fl.routes[victim]['tank']); spent = 0.0; placed = []
    say(f'  dissolve {victim}: saves J_i {saved:.3f}, must place {len(todo)} targets '
        f'(budget {saved/max(1,len(todo)):.3f} J each)')
    while todo:
        jobs = []
        for X in todo:
            for n, r in routes.items():
                c = [c for c in RI.w_cands((n, r['st'], [X], _A.dmax))
                     if np.abs(np.asarray(r['st']['tf']) - c['t']).min() > _A.gap * DAY]
                for cc in sorted(c, key=lambda c: c['dist'])[:_A.per]:
                    jobs.append((n, r['st'], r['tank'], X, cc['t']))
        if not jobs:
            say(f'    no candidate for {todo}'); return None
        res = [r for r in pl.map(w_try, jobs, chunksize=1) if r[2] is not None]
        if not res:
            say(f'    nothing insertable, {len(todo)} left: {todo}'); return None
        best = {}
        for n, X, v in res:
            if X not in best or v[2] < best[X][1][2]: best[X] = (n, v)
        X = min(best, key=lambda X: best[X][1][2]); n, (st, tk, dj) = best[X]
        routes[n] = dict(st=st, tank=tk); spent += dj; placed.append((X, n, dj))
        todo.remove(X)
        if spent >= saved:
            say(f'    over budget after {len(placed)} of {len(placed)+len(todo)} (spent {spent:.3f} >= {saved:.3f})')
            return None
    say(f'    placed all {len(placed)} for {spent:.3f} J (saved {saved:.3f}) -> net {saved - spent:+.3f}')
    return routes, saved - spent


def main():
    global _A
    ap = argparse.ArgumentParser(); ap.add_argument('src'); ap.add_argument('out')
    ap.add_argument('--victims', default='all'); ap.add_argument('--max-tank', type=float, default=1400.0)
    ap.add_argument('--dmax', type=float, default=0.35); ap.add_argument('--per', type=int, default=2)
    ap.add_argument('--gap', type=float, default=10.0); ap.add_argument('--nproc', type=int, default=8)
    a = ap.parse_args(); _A = a
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt')
    fl = RI.IFleet(a.src)
    J0 = sum(RI.cost(r['tank']) for r in fl.routes.values()) + 300 - len(fl.coverage())
    say(f'{a.src}: {len(fl.routes)} craft, covered {len(fl.coverage())}, J {J0:.3f}')
    vic = sorted(fl.routes) if a.victims == 'all' else a.victims.split(',')
    vic.sort(key=lambda n: RI.cost(fl.routes[n]['tank']) / max(1, len(set(int(x) for x in fl.routes[n]['st']['asts']))))
    best = None
    with mp.get_context('fork').Pool(a.nproc) as pl:
        for v in vic[::-1]:
            t0 = time.time()
            r = dissolve(fl, v, pl, say)
            say(f'  ({time.time()-t0:.0f} s)')
            if r and (best is None or r[1] > best[1]):
                best = r; nf = RI.IFleet(); nf.routes = {k: dict(st=x['st'], tank=x['tank']) for k, x in r[0].items()}
                nf.save(out / 'fleet', note=f'dissolved {v}')
                for k, x in r[0].items(): np.savez(out / f'route_x{k}.npz', **x['st'])
                say(f'  SAVED: {len(nf.routes)} craft, covered {len(nf.coverage())}, J {nf.J():.3f}')
    if best is None: say('no dissolve succeeded')


if __name__ == '__main__':
    main()
