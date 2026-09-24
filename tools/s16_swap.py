"""Stage 16: EXCHANGE (2-opt) columns for a 9-craft fleet: A' = A - x + y and B' = B - y + x.

Single-target moves (s15b_relocate, s16_iter) cannot find exchanges whose halves are each unprofitable.  Here:
  1. removal columns (A - x) come from the s16_iter cache (results/s16/iter/*/cache.jsonl, keyed by route content) or are
     settled here (run_ialns.w_remove);
  2. lin_price of every target y of route B into every other route A (s15b_close._price_host, approaches <= --dmax);
  3. pair estimate  g = save(A,x) + save(B,y) - lin(y -> A) - lin(x -> B)  over pairs where x passes B and y passes A;
     the top --pairs pairs are built: twin insertion (s15b_close._insert) of y into the settled (A - x) state and of x into
     (B - y); every settled result is saved as a column OUT/cols/route_<hashA8>m<x>x<y>_<n>.npz.
The exact selector (tools/s16_select.py / s16_iter) then combines them with everything else.

usage: s16_swap.py OUT FLEET [--nproc 6] [--pairs 60] [--dmax 0.2] [--min-save 0.004] [--cache results/s16/iter/A/cache.jsonl]"""
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
from s16_colgen import w_rem
from s16_iter import chash


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('fleet')
    ap.add_argument('--nproc', type=int, default=6); ap.add_argument('--pairs', type=int, default=60)
    ap.add_argument('--dmax', type=float, default=0.2); ap.add_argument('--min-save', type=float, default=0.004)
    ap.add_argument('--cache', default='results/s16/iter/A/cache.jsonl'); ap.add_argument('--timeout', type=float, default=300.0)
    a = ap.parse_args()
    out = ROOT / a.out; (out / 'cols').mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt'); RI.eph()
    files = FA.fleet_files(a.fleet)
    fl = {k: dict(st=FA.load_st(f), f=f) for k, f in files.items()}
    for r in fl.values():
        r['tank'] = float(RI.ipr(r['st']).tank()); r['h'] = chash(r['st'])
    cache = {}
    for c in a.cache.split(','):
        if c and (ROOT / c).exists():
            for l in open(ROOT / c):
                try:
                    m = json.loads(l)
                    if m['key'][0] == 'rem':
                        cache[(m['key'][1], int(m['key'][2]))] = m
                except Exception:
                    pass
    sj0 = sum(RI.cost(r['tank']) for r in fl.values())
    say(f'swap {a.fleet}: sum J_i {sj0:.4f}; {len(cache)} cached removals')
    tic = time.time()
    with mp.get_context('fork').Pool(a.nproc, maxtasksperchild=8) as pool:
        # ---- 1. removals
        jobs = [(k, r['st'], int(x), a.timeout) for k, r in fl.items() for x in r['st']['asts'] if (r['h'], int(x)) not in cache]
        for r in pool.imap_unordered(w_rem, jobs, chunksize=1):
            k = r['name']; x = int(r['asts'][0]); m = dict(ok=bool(r.get('ok')))
            if m['ok']:
                p = out / 'cols' / f'route_{fl[k]["h"][:8]}m{x}_{len(r["st"]["asts"])}.npz'; np.savez(p, **r['st'])
                m.update(tank=float(r['tank']), col=str(p.relative_to(ROOT)))
            cache[(fl[k]['h'], x)] = m
        rem = {}
        for k, r in fl.items():
            for x in r['st']['asts']:
                m = cache.get((r['h'], int(x)))
                if m and m['ok']:
                    rem[(k, int(x))] = dict(save=RI.cost(r['tank']) - RI.cost(m['tank']), col=m['col'], tank=m['tank'])
        say(f'  removals: {len(jobs)} new, {len(rem)} settled ({time.time() - tic:.0f} s)')
        # ---- 2. lin prices of every target into every other route
        own = {k: [int(x) for x in r['st']['asts']] for k, r in fl.items()}
        jobs = [(h, fl[h]['st'], sorted(x for k in fl if k != h for x in own[k]), a.dmax) for h in fl]
        P = pool.map(CL._price_host, jobs, chunksize=1)
        lin = {}
        for lst in P:
            for c in lst:
                key = (c['host'], c['ast'])
                if key not in lin or c['lin_dJ'] < lin[key]['lin_dJ']:
                    lin[key] = c
        say(f'  {len(lin)} (host, target) lin prices ({time.time() - tic:.0f} s)')
        # ---- 3. pairs
        where = {x: k for k, xs in own.items() for x in xs}
        pairs = []
        for (A, y), cy in lin.items():                   # y (of route B) into A
            B = where.get(y)
            if B is None or B == A or (B, y) not in rem:
                continue
            for x in own[A]:
                cx = lin.get((B, x))
                if cx is None or (A, x) not in rem:
                    continue
                g = rem[(A, x)]['save'] + rem[(B, y)]['save'] - cy['lin_dJ'] - cx['lin_dJ']
                pairs.append((g, A, x, B, y, cy, cx))
        pairs.sort(key=lambda z: -z[0])
        seen = set(); top = []
        for p in pairs:
            key = frozenset([(p[1], p[2]), (p[3], p[4])])
            if key in seen:
                continue
            seen.add(key); top.append(p)
            if len(top) >= a.pairs:
                break
        say(f'  {len(pairs)} pairs; top estimates ' + ', '.join(f'{x}:{A}<->{y}:{B} {g:+.4f}' for g, A, x, B, y, _, _ in top[:12]))
        # ---- 4. build the halves: y into (A - x), x into (B - y)
        halves = {}
        for g, A, x, B, y, cy, cx in top:
            halves.setdefault((A, x, y), (cy['t'], cy['dist'])); halves.setdefault((B, y, x), (cx['t'], cx['dist']))
        jobs = []
        for (A, x, y), (t, d) in halves.items():
            st = FA.load_st(rem[(A, x)]['col'])
            jobs.append((f'{A}|{x}|{y}', st, y, t, d, a.timeout, rem[(A, x)]['tank']))
        say(f'  {len(jobs)} half-exchange insertions')
        res = {}
        for r in pool.imap_unordered(CL._insert, jobs, chunksize=1):
            A, x, y = r['host'].split('|'); x = int(x); y = int(y)
            if r['ok']:
                p = out / 'cols' / f'route_{fl[A]["h"][:8]}m{x}x{y}_{len(r["st"]["asts"])}.npz'; np.savez(p, **r['st'])
                res[(A, x, y)] = dict(tank=float(r['tank']), col=str(p.relative_to(ROOT)))
        rows = []
        for g, A, x, B, y, _, _ in top:
            ra = res.get((A, x, y)); rb = res.get((B, y, x))
            if ra and rb:
                gain = (RI.cost(fl[A]['tank']) - RI.cost(ra['tank'])) + (RI.cost(fl[B]['tank']) - RI.cost(rb['tank']))
                rows.append(dict(A=A, x=x, B=B, y=y, est=round(g, 4), gain=round(gain, 4)))
        rows.sort(key=lambda r: -r['gain'])
        say(f'  exchanges built: {len(rows)} of {len(top)}; best ' + ', '.join(f'{r["x"]}:{r["A"]}<->{r["y"]}:{r["B"]} {r["gain"]:+.4f} (est {r["est"]:+.4f})' for r in rows[:12])
            + f' ({time.time() - tic:.0f} s)')
    json.dump(dict(fleet=a.fleet, sumJi=round(sj0, 5), exchanges=rows), open(out / 'swap.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
