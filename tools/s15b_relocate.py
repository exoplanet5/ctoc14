"""s15b: RELOCATE expensive targets of a fleet's routes into cheaper hosts (twin level).

1. removal savings (run_ialns.w_remove: remove + settle) of every target of the --from routes (default: every route with
   kg/flyby > --kgfb), kept when dJ >= --min-save and the target is covered only by that route;
2. for the top --top targets: approaches to every OTHER route within --dmax AU, lin_price screened (s15b_close._price_host),
   twin insertion attempts (s15b_close._insert) where lin dJ < --screen x saving, at most --per-target hosts each;
3. a move is (saving - insertion dJ) > --min-gain; non-conflicting moves (each route in at most one move per round) are
   applied best first; rounds repeat on the changed routes until nothing improves or --rounds.
Every settled removal / insertion is saved as a column under OUT/cols (the selector can recombine them).

usage: s15b_relocate.py OUT FLEET [--from t1,s1] [--kgfb 11] [--top 20] [--dmax 0.2] [--nproc 4] [--rounds 2]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)
import numpy as np
import run_ialns as RI
import s15b_fleet as FA
import s15b_close as CL
from greedy_cover import ALL
from ctoc14.constants import DAY


def savecol(out, tag, st):
    n = len(st['asts']); p = out / 'cols' / f'route_{tag}_{n}.npz'; i = 1
    while p.exists():
        p = out / 'cols' / f'route_{tag}_{n}_{i}.npz'; i += 1
    np.savez(p, **st)
    return str(p.relative_to(ROOT))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('fleet')
    ap.add_argument('--from', dest='frm', default=''); ap.add_argument('--kgfb', type=float, default=11.0)
    ap.add_argument('--top', type=int, default=20); ap.add_argument('--dmax', type=float, default=0.2)
    ap.add_argument('--nproc', type=int, default=4); ap.add_argument('--rounds', type=int, default=2)
    ap.add_argument('--min-save', type=float, default=0.008); ap.add_argument('--min-gain', type=float, default=0.002)
    ap.add_argument('--screen', type=float, default=1.3); ap.add_argument('--per-target', type=int, default=3)
    ap.add_argument('--timeout', type=float, default=300.0)
    a = ap.parse_args()
    if a.nproc > 8:
        raise SystemExit('CPU budget: --nproc <= 8')
    out = ROOT / a.out; (out / 'cols').mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt'); RI.eph()
    files = FA.fleet_files(a.fleet)
    fl = {k: dict(st=FA.load_st(f), f=f) for k, f in files.items()}
    for r in fl.values():
        r['tank'] = float(RI.ipr(r['st']).tank())
    sj0 = sum(RI.cost(r['tank']) for r in fl.values())
    say(f'relocate {a.fleet}: sum J_i {sj0:.4f}')
    frm = [x for x in a.frm.split(',') if x] or [k for k, r in fl.items() if (r['tank'] - 600) / len(r['st']['asts']) > a.kgfb]
    moves_done = []; dirty = set(frm); rem_cache = {}
    with mp.get_context('fork').Pool(a.nproc, maxtasksperchild=4) as pool:
        for rnd in range(1, a.rounds + 1):
            cov = {}
            for k, r in fl.items():
                for x in r['st']['asts']:
                    cov.setdefault(int(x), set()).add(k)
            jobs = [(k, fl[k]['st'], [int(x)]) for k in sorted(dirty) if k in frm for x in fl[k]['st']['asts']
                    if len(cov[int(x)]) == 1]
            if jobs:
                R = pool.map(RI.w_remove, jobs, chunksize=1)
                for r in R:
                    if r['ok']:
                        rem_cache[(r['name'], r['asts'][0])] = r
            dirty = set()
            sav = []
            for (k, x), r in rem_cache.items():
                if k not in fl or x not in set(int(y) for y in fl[k]['st']['asts']) or len(cov.get(x, ())) != 1:
                    continue
                s = RI.cost(fl[k]['tank']) - RI.cost(r['tank'])
                if s >= a.min_save:
                    sav.append((s, k, x))
            sav.sort(reverse=True); top = sav[:a.top]
            say(f' round {rnd}: {len(sav)} removals >= {a.min_save}; top ' + ', '.join(f'{x}@{k} {s:.3f}' for s, k, x in top[:12]))
            if not top:
                break
            # price insertions of the top targets into every other route
            X_by_host = {}
            for s, k, x in top:
                for h in fl:
                    if h != k:
                        X_by_host.setdefault(h, set()).add(x)
            P = pool.map(CL._price_host, [(h, fl[h]['st'], sorted(xs), a.dmax) for h, xs in X_by_host.items()], chunksize=1)
            priced = [c for lst in P for c in lst]
            src = {x: (s, k) for s, k, x in top}
            att = []
            for x, (s, k) in src.items():
                cs = sorted([c for c in priced if c['ast'] == x and c['host'] != k and c['lin_dJ'] < a.screen * s],
                            key=lambda c: c['lin_dJ'])
                seen_h = set()
                for c in cs:
                    if c['host'] in seen_h:
                        continue
                    seen_h.add(c['host']); att.append(c)
                    if len(seen_h) >= a.per_target:
                        break
            say(f'  {len(priced)} priced approaches; {len(att)} insertion attempts: ' +
                ', '.join(f'{c["ast"]}->{c["host"]} lin {c["lin_dJ"]:.3f} (save {src[c["ast"]][0]:.3f})' for c in att))
            if not att:
                break
            res = pool.map(CL._insert, [(c['host'], fl[c['host']]['st'], c['ast'], c['t'], c['dist'], a.timeout,
                                         fl[c['host']]['tank']) for c in att], chunksize=1)
            cand = []
            for c, r in zip(att, res):
                s, k = src[c['ast']]
                if r['ok']:
                    r['col'] = savecol(out, f'{r["host"]}x{r["ast"]}', r['st'])
                    g = s - r['dJ']; cand.append((g, c['ast'], k, r))
                say(f'    {c["ast"]}: {k} -> {c["host"]} ({c["dist"]:.3f} AU, lin {c["lin_dJ"]:.3f}): ' +
                    (f'OK +{r["dkg"]:.1f} kg dJ {r["dJ"]:.4f}; gain {s - r["dJ"]:+.4f}' if r['ok'] else f'fail {r["err"]}') +
                    f' ({r["sec"]:.0f} s)')
            used = set(); nm = 0
            for g, x, k, r in sorted(cand, key=lambda z: -z[0]):
                if g < a.min_gain or k in used or r['host'] in used:
                    continue
                rm = rem_cache[(k, x)]
                fl[k].update(st=rm['st'], tank=rm['tank'], f=savecol(out, f'{k}m{x}', rm['st']))
                fl[r['host']].update(st=r['st'], tank=r['tank'], f=r['col'])
                used |= {k, r['host']}; dirty |= {k, r['host']}; nm += 1
                moves_done.append(dict(ast=x, src=k, dst=r['host'], gain=round(g, 4), round=rnd))
                say(f'  MOVE {x}: {k} -> {r["host"]} gain {g:.4f}')
            for k in used:                                     # removal prices of changed routes are stale
                for key in [kk for kk in rem_cache if kk[0] == k]:
                    del rem_cache[key]
            if not nm:
                break
    F = RI.IFleet()
    for k, r in fl.items():
        F.routes[k] = dict(st=r['st'], tank=r['tank'])
    F.save(out / 'fleet', note=f's15b_relocate of {a.fleet}')
    sj = sum(RI.cost(r['tank']) for r in fl.values())
    rep = dict(fleet=a.fleet, sumJi0=round(sj0, 4), sumJi=round(sj, 4), covered=len(F.coverage()), moves=moves_done,
               files={k: r['f'] for k, r in fl.items()})
    json.dump(rep, open(out / 'relocate.json', 'w'), indent=1)
    say(f'relocate done: sum J_i {sj0:.4f} -> {sj:.4f} ({len(moves_done)} moves), covered {rep["covered"]}')


if __name__ == '__main__':
    main()
