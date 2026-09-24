"""Deepen every route of a fleet so the set cover gains SLACK (stage 7, section 10).

The cover over 471 routes picks 10 craft whose union is exactly 298 targets: no redundancy, so no craft can ever be
dropped.  Deepening changes that.  A flyby added to a route is worth it whenever its J cost is below the fleet's
current J per flyby (~0.0415), i.e. ~34 kg of tank; every target thus added is one more route that covers it, and the
MILP can then drop the craft that was only kept for it.

Each route is grown in its own process (nearest-first insertion trials, accept the first within --max-step) on every
target it does not already fly.  Output: the grown routes as a pool, to be handed straight back to route_cover.py.

Usage: fleet_deepen.py fleet_dir out_dir [--max-step 34] [--max-tank 1200] [--extra 12]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, time, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from ctoc14.search import UNREACHABLE
from ctoc14.constants import DAY
ALL = [t for t in range(1, 301) if t not in UNREACHABLE]
_A = None


def w_deep(job):
    n, st, tank = job; a = _A
    t0 = time.time(); added = []; nfail = 0
    while len(added) < a.extra and time.time() - t0 < a.tmax:
        have = set(int(x) for x in st['asts'])
        left = [t for t in ALL if t not in have]
        if not left: break
        cands = RI.w_cands((n, st, left, a.dmax))
        tf = np.asarray(st['tf']); best = {}
        for x in sorted(cands, key=lambda x: x['dist']):
            if np.abs(tf - x['t']).min() > a.gap * DAY: best.setdefault(x['ast'], x)
        trial = sorted(best.values(), key=lambda x: x['dist'])[:a.ntrial]
        if not trial: break
        got = False
        for x in trial:
            r = RI.w_insert((n, st, x['ast'], x['t']))
            if r.get('ok') and r['tank'] <= a.max_tank and r['tank'] - tank <= a.max_step:
                st, tank = r['st'], r['tank']; added.append(x['ast']); got = True; break
        if not got:
            nfail += 1
            if nfail >= a.nfail: break
        else:
            nfail = 0
    return n, st, float(tank), added, time.time() - t0


def main():
    global _A
    ap = argparse.ArgumentParser(); ap.add_argument('src'); ap.add_argument('out')
    ap.add_argument('--max-step', type=float, default=34.0); ap.add_argument('--max-tank', type=float, default=1200.0)
    ap.add_argument('--extra', type=int, default=12); ap.add_argument('--dmax', type=float, default=0.2)
    ap.add_argument('--ntrial', type=int, default=10); ap.add_argument('--gap', type=float, default=12.0)
    ap.add_argument('--nfail', type=int, default=3); ap.add_argument('--tmax', type=float, default=1500.0)
    ap.add_argument('--nproc', type=int, default=8)
    a = ap.parse_args(); _A = a
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt')
    fl = RI.IFleet(a.src)
    say(f'deepening {len(fl.routes)} routes of {a.src} (J {fl.J():.3f}, covered {len(fl.coverage())})')
    jobs = [(n, r['st'], r['tank']) for n, r in sorted(fl.routes.items())]
    new = RI.IFleet()
    with mp.get_context('fork').Pool(a.nproc) as pl:
        for n, st, tank, added, dt in pl.imap_unordered(w_deep, jobs):
            n0 = len(jobs[[j[0] for j in jobs].index(n)][1]['asts']); t0 = jobs[[j[0] for j in jobs].index(n)][2]
            new.routes[n] = dict(st=st, tank=tank)
            np.savez(out / f'route_d{n}.npz', **st)
            say(f'  {n}: {n0} -> {len(st["asts"])} fb, {t0:.0f} -> {tank:.0f} kg, '
                f'{RI.cost(t0)/n0:.4f} -> {RI.cost(tank)/len(st["asts"]):.4f} J/fb  (+{added})  ({dt:.0f} s)')
    new.save(out / 'fleet', note='deepened')
    say(f'deepened fleet: covered {len(new.coverage())}, sum J_i {sum(RI.cost(r["tank"]) for r in new.routes.values()):.3f}, '
        f'total flybys {sum(len(r["st"]["asts"]) for r in new.routes.values())} (surplus '
        f'{sum(len(r["st"]["asts"]) for r in new.routes.values()) - len(new.coverage())})')


if __name__ == '__main__':
    main()
