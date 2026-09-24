"""Strip the REDUNDANT flybys out of a fleet: the one column-generation mechanism the set cover cannot do itself.

The measurement that motivates this.  Two fleets, both covering all 298 reachable targets:

    results/s7/deep1/fleet   10 routes, 365 flybys, 298 covered   J/fb 0.0368   sum J_i 13.430   J 15.430
    set-cover optimum        10 routes, 298 flybys, 298 covered   J/fb 0.0415   sum J_i 12.355   J 14.355

The first is made of the best routes we can build -- 0.0368 J per flyby, at the 0.0369 the contest target needs -- and
it loses anyway, because 67 of its 365 flybys are DUPLICATES: two craft visiting the same asteroid, which the objective
pays for twice and counts once.  The second wins by being perfectly disjoint at a mediocre per-flyby rate.  J < 13 needs
both properties at once, and no route we know how to build has them.

The set cover cannot fix this: it may only pick whole routes, never edit one.  Column generation cannot either -- the
duals are degenerate here (a column priced at reduced cost +2.3 moved the LP bound by zero).  But a REMOVAL is a new,
cheaper column for the same coverage, generated without any pricing problem: drop a duplicated target from every host
but the one that carries it most cheaply, re-settle, and the tank falls by roughly what the flyby cost to insert.

Greedy, exactly as the objective is written.  Every candidate removal is a real `run_ialns.w_remove` (drop the flyby,
re-settle the whole route) and is kept only if the settled tank drops and the target is still covered by someone.

Usage: fleet_dedupe.py fleet_dir out_dir [--rounds 4] [--nproc 8]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse, collections, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI


def report(fl):
    cov = fl.coverage(); s = sum(RI.cost(r['tank']) for r in fl.routes.values())
    tot = sum(len(set(int(x) for x in r['st']['asts'])) for r in fl.routes.values())
    return s, len(cov), tot, s + 298 - len(cov) + 2


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('src'); ap.add_argument('out')
    ap.add_argument('--rounds', type=int, default=4); ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--per-round', type=int, default=40, help='max removals trialled per round')
    ap.add_argument('--bulk', action='store_true', help='one pass: strip every profitable duplicate at once')
    a = ap.parse_args()
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt')
    fl = RI.IFleet(a.src)
    s, c, tot, J = report(fl)
    say(f'{a.src}: {len(fl.routes)} routes, {tot} flybys, {c} covered ({tot - c} duplicate flybys), '
        f'sum J_i {s:.3f}, J {J:.3f}, J/fb {s/tot:.4f}')
    with mp.get_context('fork').Pool(a.nproc) as pl:
        if a.bulk:
            # one pass: price every (duplicated target, host) removal, keep each target on the host that carries it
            # most cheaply, and strip it from all the others in ONE w_remove per route (states never go stale).
            cov = fl.coverage(); dup = {t: sorted(h) for t, h in cov.items() if len(h) > 1}
            jobs = [(h, fl.routes[h]['st'], [t]) for t, H in dup.items() for h in H]
            tic = time.time(); res = pl.map(RI.w_remove, jobs)
            g = {}
            for r in res:
                if r.get('ok'): g[(r['name'], int(r['asts'][0]))] = fl.routes[r['name']]['tank'] - r['tank']
            say(f'bulk: {len(dup)} duplicated targets, {len(jobs)} trials, {len(g)} settled ({time.time()-tic:.0f} s)')
            drop = collections.defaultdict(list)
            for t, H in dup.items():
                pr = [(g.get((h, t), -1e9), h) for h in H]
                pr.sort()                                   # keep the host with the SMALLEST saving
                for _, h in pr[1:]:
                    if g.get((h, t), 0.0) > 0.5: drop[h].append(t)
            say('  drops per route: ' + ' '.join(f'{h}:{len(v)}' for h, v in sorted(drop.items())))
            tic = time.time()
            out2 = pl.map(RI.w_remove, [(h, fl.routes[h]['st'], v) for h, v in sorted(drop.items())])
            for r in out2:
                h = r['name']
                if not r.get('ok'):
                    say(f'  {h}: bulk removal of {len(r["asts"])} targets did NOT settle -- kept as is'); continue
                say(f'  {h}: -{len(r["asts"])} flybys, tank {fl.routes[h]["tank"]:.0f} -> {r["tank"]:.0f} kg')
                fl.routes[h] = dict(st=r['st'], tank=r['tank'])
            s, c, tot, J = report(fl)
            say(f'bulk pass ({time.time()-tic:.0f} s): {tot} flybys, {c} covered, sum J_i {s:.3f}, J {J:.3f}, J/fb {s/tot:.4f}')
        for rd in range(0 if a.bulk else a.rounds):
            cov = fl.coverage()
            dup = {t: sorted(h) for t, h in cov.items() if len(h) > 1}
            if not dup: say('no duplicates left'); break
            # trial: remove target t from host h, for every (t, h) with t covered more than once
            jobs = []
            for t, H in dup.items():
                for h in H: jobs.append((h, fl.routes[h]['st'], [t]))
            jobs = jobs[:a.per_round * len(fl.routes)]
            tic = time.time()
            res = pl.map(RI.w_remove, jobs)
            gains = []
            for r in res:
                if not r.get('ok'): continue
                h = r['name']; t = int(r['asts'][0])
                g = fl.routes[h]['tank'] - r['tank']
                if g <= 0.5: continue
                gains.append((g, h, t, r['st'], float(r['tank'])))
            gains.sort(key=lambda x: -x[0])
            say(f'round {rd}: {len(dup)} duplicated targets, {len(jobs)} removal trials, '
                f'{len(gains)} profitable ({time.time()-tic:.0f} s); best saves {gains[0][0]:.0f} kg' if gains
                else f'round {rd}: {len(dup)} duplicated targets, no profitable removal')
            if not gains: break
            # accept greedily: one removal per route per round (the states go stale after a change),
            # and never drop a target's last host
            left = {t: len(H) for t, H in cov.items()}
            used = set(); nacc = 0
            for g, h, t, st, tank in gains:
                if h in used or left[t] <= 1: continue
                fl.routes[h] = dict(st=st, tank=tank); used.add(h); left[t] -= 1; nacc += 1
            s, c, tot, J = report(fl)
            say(f'  accepted {nacc} removals -> {tot} flybys, {c} covered, sum J_i {s:.3f}, J {J:.3f}, J/fb {s/tot:.4f}')
    fl.save(out / 'fleet', note='redundancy-stripped')
    for n, r in fl.routes.items(): np.savez(out / f'route_{n}.npz', **r['st'])
    s, c, tot, J = report(fl)
    say(f'FINAL: {len(fl.routes)} routes, {tot} flybys, {c} covered, sum J_i {s:.3f}, J {J:.3f}, J/fb {s/tot:.4f}')


if __name__ == '__main__':
    main()
