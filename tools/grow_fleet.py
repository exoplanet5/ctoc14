"""Grow a whole fleet from scratch with the impulsive whole-trajectory model (simultaneous cheapest insertion).

K routes launched at --tl (days, one value or K values) with tangential launch v_inf spread uniformly over
[-vt, +vt] km/s (so their phases drift apart), no flybys. Every route keeps a cache of evaluated insertion trials
(nearest unassigned targets within --dmax AU of its trajectory). Each step applies, for every route, its cheapest cached
insertion whose target is still unassigned (one per route, conflicts resolved by cheapest dJ first), provided
dJ <= --lam; the modified routes' caches are recomputed in parallel. Stops when no route has an insertion <= lam or all
targets are covered. Output: <out>/route_<k>.npz (run_ialns fleet format) + grow.log; continue with
`run_ialns.py search <out> ... --relocate` (tail targets are placed by relocate/dissolve logic or stay missed).
Usage: grow_fleet.py outdir --k 11 [--tl 0] [--vt 1.5] [--lam 0.05] [--dmax 0.06] [--trials 12] [--nproc 8]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse, multiprocessing as mp
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1])); sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import numpy as np
from run_ialns import IFleet, ist, ipr, w_insert, w_cands, cost, logger, eph
from ctoc14.impulsive import ImpulsiveProblem
from ctoc14.constants import DAY, T_MISSION
from ctoc14.search import UNREACHABLE


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out')
    ap.add_argument('--k', type=int, default=11); ap.add_argument('--tl', default='0'); ap.add_argument('--vt', type=float, default=1.5)
    ap.add_argument('--lam', type=float, default=0.05); ap.add_argument('--dmax', type=float, default=0.06)
    ap.add_argument('--trials', type=int, default=12); ap.add_argument('--nproc', type=int, default=8)
    a = ap.parse_args()
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True); say = logger(out / 'grow.log')
    E = eph()
    tls = [float(x) for x in a.tl.split(',')]
    tls = tls if len(tls) == a.k else [tls[0]] * a.k
    vts = np.linspace(-a.vt, a.vt, a.k) if a.k > 1 else np.array([0.0])
    fleet = IFleet()
    for k in range(a.k):
        tL = tls[k] * DAY
        rE, vE = E.earth_state(tL)
        vinf = vE / np.linalg.norm(vE) * vts[k]
        ts = np.arange(tL + 10 * DAY, T_MISSION - 2 * DAY, 20 * DAY)
        ip = ImpulsiveProblem(E, tL, vinf, ts, np.zeros((len(ts), 3)), np.zeros(0), [])
        fleet.routes[f'{k + 1:02d}'] = dict(st=ist(ip), tank=float(ip.tank()))
    U = set(range(1, 301)) - set(UNREACHABLE)
    say(f'grow fleet: {a.k} routes, launch {tls}, tangential v_inf {np.round(vts, 2).tolist()}, lambda {a.lam}, dmax {a.dmax}')
    cache = {n: [] for n in fleet.routes}; stale = set(fleet.routes)
    step = 0
    with mp.get_context('fork').Pool(a.nproc, maxtasksperchild=40) as pool:
        while U:
            step += 1; tic = time.time()
            jobs = []
            for n in sorted(stale):
                st = fleet.routes[n]['st']
                cands = sorted(w_cands((n, st, sorted(U), a.dmax)), key=lambda c: c['dist'])[:a.trials]
                jobs += [(n, st, c['ast'], c['t']) for c in cands]
                cache[n] = []
            for r in pool.map(w_insert, jobs):
                if r['ok']:
                    r['dJ'] = cost(r['tank']) - fleet.routes[r['name']]['tank'] * 0 - cost(fleet.routes[r['name']]['tank'])
                    cache[r['name']].append(r)
            offers = []
            for n, lst in cache.items():
                lst = [r for r in lst if r['ast'] in U]
                cache[n] = lst
                if lst:
                    b = min(lst, key=lambda r: r['dJ'])
                    if b['dJ'] <= a.lam:
                        offers.append(b)
            stale = set()
            taken = set()
            for b in sorted(offers, key=lambda r: r['dJ']):
                if b['ast'] in taken:
                    continue
                taken.add(b['ast']); U.discard(b['ast'])
                fleet.routes[b['name']] = dict(st=b['st'], tank=b['tank']); stale.add(b['name'])
            if not stale:
                # routes whose cache is empty (no candidates within dmax) get no offers; stop when nothing is offered
                say(f'step {step}: no insertion <= lambda; stop with {len(U)} unassigned'); break
            fleet.save(out, note=f'grow step {step}')
            say(f'step {step}: +{len(taken)} ({", ".join(f"{b["ast"]}->{b["name"]} {b["dJ"]:+.3f}" for b in sorted(offers, key=lambda r: r["dJ"]) if b["ast"] in taken)}); '
                f'covered {298 - len(U)}, sum J_i {sum(cost(r["tank"]) for r in fleet.routes.values()):.3f} ({time.time() - tic:.0f} s, {len(jobs)} trials)')
    fleet.save(out, note='grow final')
    say(f'final: {fleet.summary()}; unassigned {sorted(U)}')


if __name__ == '__main__':
    main()
