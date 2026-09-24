"""s15b: rank near-complete 9-craft fleets by their ESTIMATED CLOSED cost.

Selection (cached column pool, s15b_select.load_cached): K=298 (if feasible: done), then --alt alternatives at K=297 and
K=296 (no-good cuts on each found miss set).  For every alternative fleet: the misses' approaches to every route within
--dmax AU, lin_price screened (s15b_close._price_host); estimated closed sum J_i = sum J_i + for each miss the cheapest
lin dJ (hosts distinct, greedy), 'inf' if a miss has no priced approach.  Writes OUT/eval.json (fleets sorted by the
estimate) and saves every alternative fleet as OUT/fleets/<K>_<i>/ (IFleet).

usage: s15b_eval.py OUT [--alt 4] [--Ks 297,296] [--dmax 0.2] [--nproc 1] [--ftlim 240]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)
import numpy as np
import run_ialns as RI
import s15_select as SS
import s15b_select as BS
import s15b_fleet as FA
import s15b_close as CL


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--alt', type=int, default=4)
    ap.add_argument('--Ks', default='297,296'); ap.add_argument('--dmax', type=float, default=0.2)
    ap.add_argument('--nproc', type=int, default=1); ap.add_argument('--ftlim', type=float, default=240.0)
    ap.add_argument('--force', default='')
    a = ap.parse_args()
    out = ROOT / a.out; out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt'); RI.eph()
    cols = BS.load_cached(BS.PATS, a.nproc, say)
    force = [int(x) for x in a.force.split(',') if x]
    pick, msg = BS.solve(cols, 9, 298, force, [], a.ftlim)
    rows = []
    if pick is not None:
        J, cv, sj = SS.evaluate(cols, pick)
        say(f'K=298 FEASIBLE: sum J_i {sj:.4f}')
        rows.append(dict(K=298, i=0, covered=cv, sumJi=round(sj, 4), misses=[], est=round(sj, 4),
                         picks=[cols[j]['f'] for j in pick]))
    else:
        say(f'K=298: {msg[:50]}')
    cuts = []
    for K in [int(x) for x in a.Ks.split(',')]:
        for i in range(a.alt):
            pick, msg = BS.solve(cols, 9, K, force, cuts, a.ftlim)
            if pick is None:
                say(f'  K>={K} alt {i}: {msg[:50]}'); break
            J, cv, sj = SS.evaluate(cols, pick)
            miss = sorted(set(SS.ALL) - set().union(*[cols[j]['set'] for j in pick]))
            cuts.append(miss)
            rows.append(dict(K=K, i=i, covered=cv, sumJi=round(sj, 4), misses=miss, picks=[cols[j]['f'] for j in pick]))
            say(f'  K>={K} alt {i}: {cv} covered @ {sj:.4f}, misses {miss}')
    # price the misses of every alternative
    for r in rows:
        files = {f'{j + 1:02d}': f for j, f in enumerate(r['picks'])}
        sts = {k: FA.load_st(f) for k, f in files.items()}
        if not r['misses']:
            F = RI.IFleet()
            for k, st in sts.items():
                F.routes[k] = dict(st=st, tank=float(RI.ipr(st).tank()))
            F.save(out / 'fleets' / f'{r["K"]}_{r["i"]}', note=f'eval K {r["K"]} #{r["i"]} (closed)')
            r['dir'] = str((out / 'fleets' / f'{r["K"]}_{r["i"]}').relative_to(ROOT)); r['plan'] = []
            continue
        with mp.get_context('fork').Pool(a.nproc) as pool:
            P = pool.map(CL._price_host, [(k, st, r['misses'], a.dmax) for k, st in sts.items()], chunksize=1)
        pr = sorted([c for lst in P for c in lst], key=lambda c: c['lin_dJ'])
        used = set(); est = r['sumJi']; plan = []
        for X in r['misses']:
            cs = [c for c in pr if c['ast'] == X and c['host'] not in used]
            if not cs:
                est = float('inf'); plan.append((X, None, None)); continue
            c = cs[0]; used.add(c['host']); est += c['lin_dJ']; plan.append((X, c['host'], round(c['lin_dJ'], 4)))
        r['est'] = round(est, 4) if np.isfinite(est) else None; r['plan'] = plan
        r['priced'] = [dict(ast=c['ast'], host=c['host'], dist=round(c['dist'], 3), lin_dJ=round(c['lin_dJ'], 4)) for c in pr[:12]]
        say(f'  K={r["K"]} alt {r["i"]}: {r["covered"]} @ {r["sumJi"]:.4f} misses {r["misses"]} -> est closed '
            f'{"inf" if r["est"] is None else round(r["est"], 4)} plan {plan}')
        F = RI.IFleet()
        for k, st in sts.items():
            F.routes[k] = dict(st=st, tank=float(RI.ipr(st).tank()))
        F.save(out / 'fleets' / f'{r["K"]}_{r["i"]}', note=f'eval alt K {r["K"]} #{r["i"]}')
        r['dir'] = str((out / 'fleets' / f'{r["K"]}_{r["i"]}').relative_to(ROOT))
    rows.sort(key=lambda r: (r.get('est') is None, r.get('est') or 99))
    json.dump(rows, open(out / 'eval.json', 'w'), indent=1)
    say('ranking: ' + ' | '.join(f'{r["K"]}#{r["i"]} {r["sumJi"]:.4f} -> {r.get("est")}' for r in rows))


if __name__ == '__main__':
    main()
