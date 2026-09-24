"""Measurement: SEQUENTIAL twin closing of a settled pass fleet (the plan's Gate-2 quantity), then relocate/swap polish.
Round = candidates (RI.candidates, per-target radius widening 0.08 -> 0.15 -> 0.30 AU on failure) -> RI.w_insert trials
(parallel) -> apply the cheapest insertion per host (one per host per round, re-priced every round, so escalation is
measured) if dJ < CAP.  A target that fails 3 rounds is retried once with a longer homotopy (stages 20, iters 10);
then it is unplaceable.  After closing: relocate passes (removal saving vs insertion elsewhere, swaps) to saturation.
usage: close_probe.py SRC_FLEET OUT_DIR [CAP] [NPROC] [RELOCATE_PASSES]"""
import os, sys, json, time, pathlib, argparse, multiprocessing as mp
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
os.chdir(ROOT)
import numpy as np
import run_ialns as RI
from ctoc14.globalopt import insert_homotopy
from greedy_cover import ALL

def w_insert_long(job):
    name, st, ast, t = job
    ip = RI.ipr(st)
    try:
        d = insert_homotopy(ip, ast, t, stages=20, iters_stage=10, tol_km=100, log=RI.QUIET)
        if d.max() > 1e4:
            return dict(name=name, ast=ast, t=t, ok=False, miss=float(d.max()))
        miss = RI.settle(ip, 100)
    except Exception as e:
        return dict(name=name, ast=ast, t=t, ok=False, miss=np.inf, error=repr(e))
    return dict(name=name, ast=ast, t=t, ok=miss <= 150, st=RI.ist(ip), tank=float(ip.tank()), miss=miss)

def main():
    src, out = sys.argv[1], pathlib.Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
    CAP = float(sys.argv[3]) if len(sys.argv) > 3 else 0.25
    NP = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    NREL = int(sys.argv[5]) if len(sys.argv) > 5 else 3
    say = RI.logger(out / 'log.txt'); fl = RI.IFleet(src)
    base = {n: len(r['st']['asts']) for n, r in fl.routes.items()}
    sumJ = lambda f: sum(RI.cost(r['tank']) for r in f.routes.values())
    U = sorted(set(ALL) - set(fl.coverage())); J0 = sumJ(fl)
    say(f'start {fl.summary()}; sum J_i {J0:.4f}; {len(U)} leftovers; CAP {CAP}')
    fails = {u: 0 for u in U}; long_tried = set(); placed = []; dead = []
    radii = (0.08, 0.15, 0.30)
    with mp.get_context('fork').Pool(NP, maxtasksperchild=20) as pool:
        rnd = 0
        while True:
            rem = [u for u in U if u not in {p['ast'] for p in placed} and u not in dead]
            if not rem: break
            rnd += 1; tic = time.time()
            by = {}
            for r in sorted(set(radii[min(fails[u], 2)] for u in rem)):
                grp = [u for u in rem if radii[min(fails[u], 2)] == r]
                by.update(RI.candidates(fl, pool, grp, r, 3))
            jobs = []; longj = []
            for u in rem:
                for c in by.get(u, []):
                    j = (c['host'], fl.routes[c['host']]['st'], c['ast'], c['t'])
                    (longj if fails[u] >= 3 else jobs).append(j)
            res = (pool.map(RI.w_insert, jobs) if jobs else []) + (pool.map(w_insert_long, longj) if longj else [])
            best = {}
            for r in res:
                if not r.get('ok'): continue
                r['dJ'] = RI.cost(r['tank']) - RI.cost(fl.routes[r['name']]['tank'])
                if r['ast'] not in best or r['dJ'] < best[r['ast']]['dJ']: best[r['ast']] = r
            used = set(); nplace = 0
            for r in sorted(best.values(), key=lambda r: r['dJ']):
                if r['name'] in used or r['dJ'] >= CAP: continue
                used.add(r['name']); kg = r['tank'] - fl.routes[r['name']]['tank']
                fl.routes[r['name']] = dict(st=r['st'], tank=r['tank']); nplace += 1
                cands = by.get(r['ast'], []); d = min((c['dist'] for c in cands if c['host'] == r['name']), default=None)
                placed.append(dict(ast=r['ast'], host=r['name'], rnd=rnd, dJ=round(r['dJ'], 4), kg=round(kg, 1),
                                   d=None if d is None else round(d, 3), depth=len(r['st']['asts']), tank=round(r['tank'], 1),
                                   long=r['ast'] in {j[2] for j in longj}))
                say(f'  r{rnd}: placed {r["ast"]} on {r["name"]} (d {d if d is None else round(d,3)} AU, host now {len(r["st"]["asts"])} fb) '
                    f'dJ {r["dJ"]:+.4f} ({kg:+.1f} kg)')
            for u in rem:
                if u in {p['ast'] for p in placed}: continue
                if u in best and best[u]['dJ'] >= CAP and fails[u] < 3:
                    fails[u] += 1; continue          # priced but too dear this round: try again (host may change)
                if u not in best:
                    fails[u] += 1
                if fails[u] >= 4 or (fails[u] >= 3 and u in long_tried):
                    dead.append(u); say(f'  r{rnd}: {u} unplaceable (fails {fails[u]}, cheapest {best[u]["dJ"]:.3f})' if u in best else f'  r{rnd}: {u} unplaceable')
                if fails[u] >= 3: long_tried.add(u)
            say(f'round {rnd}: {len(jobs)}+{len(longj)} trials, {sum(1 for r in res if r.get("ok"))} settled, placed {nplace}; '
                f'total placed {len(placed)}, dead {len(dead)}, left {len(rem) - nplace}; sum J_i {sumJ(fl):.4f} ({time.time()-tic:.0f} s)')
            json.dump(dict(placed=placed, dead=dead), open(out / 'closing.json', 'w'), indent=1)
            fl.save(out / 'closed', note=f'closing round {rnd}')
            if nplace == 0 and all(fails[u] >= 3 for u in rem if u not in dead) and not any(fails[u] < 4 and u not in long_tried for u in rem if u not in dead):
                break
            if rnd >= 30: break
        J1 = sumJ(fl)
        dj = [p['dJ'] for p in placed]
        say(f'CLOSING: placed {len(placed)}/{len(U)} for sum dJ {sum(dj):.4f} (mean {np.mean(dj) if dj else 0:.4f}, '
            f'median {np.median(dj) if dj else 0:.4f}); unplaceable {len(dead)} {dead}; sum J_i {J0:.4f} -> {J1:.4f}; '
            f'covered {len(fl.coverage())}')
        a = argparse.Namespace(min_saving=0.005, relocate_top=40, dmax=0.3, per_target=6, retime=False, swap=True, swap_max=40)
        for k in range(NREL):
            tic = time.time()
            if not RI.relocate_pass(fl, pool, a, say): say(f'relocate pass {k+1}: no move'); break
            say(f'relocate pass {k+1}: sum J_i {sumJ(fl):.4f} ({time.time()-tic:.0f} s)')
            fl.save(out / 'relocated', note=f'relocate pass {k+1}')
        J2 = sumJ(fl)
        say(f'FINAL: covered {len(fl.coverage())}, sum J_i {J0:.4f} -> closing {J1:.4f} -> relocate {J2:.4f}; '
            f'J = {J2 + 2 + (298 - len(fl.coverage())):.4f}')
        json.dump(dict(src=src, J0=J0, J_closed=J1, J_relocated=J2, placed=placed, dead=dead, covered=len(fl.coverage()),
                       depth_before=base, depth_after={n: len(r['st']['asts']) for n, r in fl.routes.items()}),
                  open(out / 'closing.json', 'w'), indent=1)

if __name__ == '__main__':
    main()
