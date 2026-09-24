"""Stage 16: ITERATED MILP-coordinated relocation (move columns + exact selection), with a content-hash cache.

Each iteration, for the current 9-craft fleet F:
  1. removals   every (route, target): run_ialns.w_remove (remove + settle) -> column   [cached by route content]
  2. insertions targets whose removal saves >= --min-save, into every other route within --dmax AU, lin_price screen
                (lin dJ < --screen x saving, <= --per-target hosts), twin insertion s15b_close._insert -> column
                [cached by (host content, target)]
  3. re-settle  every route of F that is new this iteration (settle iters 150; lighter state kept as a column)
  4. select     exact K=298 N<=9 MILP (tools/s16_select.solve) over the honest pool + every stage-16 column
                (results/s16/**/cols) -> F'; stop when F' is not lighter by --min-gain.
Columns: OUT/cols/route_<hash8>m<x>_<n>.npz, route_<hash8>x<x>_<n>.npz, route_<hash8>rs_<n>.npz; cache OUT/cache.jsonl;
fleets OUT/it<i>/fleet (RI.IFleet) + OUT/iter.json.  Resumable (cache + iteration dirs).

usage: s16_iter.py OUT FLEET [--nproc 8] [--iters 8] [--min-save 0.004] [--dmax 0.25] [--screen 3.0] [--per-target 6]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, signal, hashlib, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)
import numpy as np
import run_ialns as RI
import s15b_fleet as FA
import s15b_close as CL
import s16_select as SEL
from s16_colgen import w_rem
from greedy_cover import _timeout, ALL


def chash(st):
    h = hashlib.sha1()
    for k in ('tL', 'vinf', 'ts', 'Ts', 'tf', 'asts'):
        h.update(np.asarray(st[k], float).round(9).tobytes())
    return h.hexdigest()[:12]


def w_resettle(job):
    k, st, timeout = job
    tic = time.time()
    try:
        signal.signal(signal.SIGALRM, _timeout); signal.setitimer(signal.ITIMER_REAL, timeout)
        ip = RI.ipr(st); miss = float(RI.settle(ip, 150))
        return dict(name=k, ok=miss <= 150.0, tank=float(ip.tank()), miss=miss, st=RI.ist(ip), sec=round(time.time() - tic))
    except Exception as e:
        return dict(name=k, ok=False, err=type(e).__name__, sec=round(time.time() - tic))
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('fleet')
    ap.add_argument('--nproc', type=int, default=8); ap.add_argument('--iters', type=int, default=8)
    ap.add_argument('--min-save', type=float, default=0.004); ap.add_argument('--dmax', type=float, default=0.25)
    ap.add_argument('--screen', type=float, default=3.0); ap.add_argument('--per-target', type=int, default=6)
    ap.add_argument('--timeout', type=float, default=300.0); ap.add_argument('--min-gain', type=float, default=5e-4)
    ap.add_argument('--ftlim', type=float, default=600.0)
    ap.add_argument('--seed-colgen', default='', help='s16_colgen OUT dir whose moves.jsonl seeds the cache (its fleet = FLEET dir)')
    ap.add_argument('--seed-fleet', default='', help='the fleet dir that colgen run was made on')
    a = ap.parse_args()
    out = ROOT / a.out; (out / 'cols').mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt'); RI.eph()
    cpath = out / 'cache.jsonl'; cache = {}
    if cpath.exists():
        for l in open(cpath):
            try:
                m = json.loads(l); cache[tuple(m['key'])] = m
            except Exception:
                pass
    def put(m, fo):
        cache[tuple(m['key'])] = m; fo.write(json.dumps(m) + '\n'); fo.flush()
    if a.seed_colgen and not cpath.exists():
        sf = FA.fleet_files(a.seed_fleet); hs = {k: chash(FA.load_st(f)) for k, f in sf.items()}; n = 0
        with open(cpath, 'a') as fo:
            for l in open(ROOT / a.seed_colgen / 'moves.jsonl'):
                m = json.loads(l)
                if m['route'] not in hs:
                    continue
                e = dict(key=[m['kind'], hs[m['route']], int(m['ast'])], ok=bool(m['ok']))
                if m['ok']:
                    e.update(tank=m['tank'], col=m['col'])
                put(e, fo); n += 1
        say(f'seeded {n} cache entries from {a.seed_colgen}')
    hist = json.load(open(out / 'iter.json')) if (out / 'iter.json').exists() else []
    spec = hist[-1]['fleet'] if hist else a.fleet
    seen_rs = set(m['key'][1] for m in cache.values() if m['key'][0] == 'rs')
    for it in range(len(hist), a.iters):
        tic = time.time()
        files = FA.fleet_files(spec)
        fl = {k: dict(st=FA.load_st(f), f=f) for k, f in files.items()}
        for r in fl.values():
            r['tank'] = float(RI.ipr(r['st']).tank()); r['h'] = chash(r['st'])
        sj0 = sum(RI.cost(r['tank']) for r in fl.values())
        say(f'iter {it}: {spec}: sum J_i {sj0:.4f}, covered {len(set(int(x) for r in fl.values() for x in r["st"]["asts"]))}; '
            + ' '.join(f'{k}:{len(r["st"]["asts"])}@{r["tank"]:.1f}' for k, r in sorted(fl.items())))
        with mp.get_context('fork').Pool(a.nproc, maxtasksperchild=8) as pool, open(cpath, 'a') as fo:
            # ---- 3 (first): re-settle routes never re-settled
            jobs = [(k, r['st'], a.timeout) for k, r in fl.items() if r['h'] not in seen_rs]
            for r in pool.imap_unordered(w_resettle, jobs, chunksize=1):
                k = r['name']; m = dict(key=['rs', fl[k]['h']], ok=bool(r.get('ok')), sec=r['sec'])
                if m['ok'] and r['tank'] < fl[k]['tank'] - 0.01:
                    p = out / 'cols' / f'route_{fl[k]["h"][:8]}rs_{len(r["st"]["asts"])}.npz'; np.savez(p, **r['st'])
                    m.update(tank=r['tank'], col=str(p.relative_to(ROOT)))
                    say(f'  re-settle {k}: {fl[k]["tank"]:.2f} -> {r["tank"]:.2f} kg')
                put(m, fo); seen_rs.add(fl[k]['h'])
                if m['ok'] and r['tank'] < fl[k]['tank'] - 0.01:       # work on the lighter state from here on
                    st2 = FA.load_st(m['col']); fl[k].update(st=st2, tank=float(RI.ipr(st2).tank()), f=m['col'], h=chash(st2))
                    seen_rs.add(fl[k]['h']); put(dict(key=['rs', fl[k]['h']], ok=False, sec=0), fo)
            # ---- 1. removals
            jobs = [(k, r['st'], int(x), a.timeout) for k, r in fl.items() for x in r['st']['asts']
                    if ('rem', r['h'], int(x)) not in cache]
            by_h = {r['h']: k for k, r in fl.items()}
            t1 = time.time()
            for r in pool.imap_unordered(w_rem, jobs, chunksize=1):
                k = r['name']; x = int(r['asts'][0]); h = fl[k]['h']
                m = dict(key=['rem', h, x], ok=bool(r.get('ok')))
                if m['ok']:
                    p = out / 'cols' / f'route_{h[:8]}m{x}_{len(r["st"]["asts"])}.npz'; np.savez(p, **r['st'])
                    m.update(tank=float(r['tank']), col=str(p.relative_to(ROOT)))
                put(m, fo)
            sav = []
            for k, r in fl.items():
                for x in r['st']['asts']:
                    m = cache.get(('rem', r['h'], int(x)))
                    if m and m['ok']:
                        s = RI.cost(r['tank']) - RI.cost(m['tank'])
                        if s >= a.min_save:
                            sav.append((s, k, int(x)))
            sav.sort(reverse=True)
            say(f'  removals: {len(jobs)} new ({time.time() - t1:.0f} s); {len(sav)} save >= {a.min_save}; top ' +
                ', '.join(f'{x}@{k} {s:.3f}' for s, k, x in sav[:10]))
            # ---- 2. insertions
            src = {x: (s, k) for s, k, x in sav}
            by_host = {}
            for s, k, x in sav:
                for h in fl:
                    if h != k:
                        by_host.setdefault(h, set()).add(x)
            P = pool.map(CL._price_host, [(h, fl[h]['st'], sorted(xs), a.dmax) for h, xs in by_host.items()], chunksize=1)
            priced = [c for lst in P for c in lst]
            att = []
            for x, (s, k) in src.items():
                cs = sorted([c for c in priced if c['ast'] == x and c['host'] != k and c['lin_dJ'] < a.screen * s],
                            key=lambda c: c['lin_dJ'])
                hs = set()
                for c in cs:
                    if c['host'] in hs:
                        continue
                    hs.add(c['host'])
                    if ('ins', fl[c['host']]['h'], x) not in cache:
                        att.append(c)
                    if len(hs) >= a.per_target:
                        break
            t2 = time.time()
            jobs = [(c['host'], fl[c['host']]['st'], c['ast'], c['t'], c['dist'], a.timeout, fl[c['host']]['tank']) for c in att]
            gains = []
            for r in pool.imap_unordered(CL._insert, jobs, chunksize=1):
                h = fl[r['host']]['h']; x = int(r['ast'])
                m = dict(key=['ins', h, x], ok=bool(r['ok']), dist=r['dist'])
                if r['ok']:
                    p = out / 'cols' / f'route_{h[:8]}x{x}_{len(r["st"]["asts"])}.npz'; np.savez(p, **r['st'])
                    m.update(tank=float(r['tank']), col=str(p.relative_to(ROOT)))
                put(m, fo)
            for x, (s, k) in src.items():
                for hname, r in fl.items():
                    m = cache.get(('ins', r['h'], x))
                    if m and m['ok']:
                        gains.append((s - (RI.cost(m['tank']) - RI.cost(r['tank'])), x, k, hname))
            gains.sort(reverse=True)
            say(f'  insertions: {len(att)} new attempts ({time.time() - t2:.0f} s); single-move gains ' +
                ', '.join(f'{x}:{k}->{h} {g:+.4f}' for g, x, k, h in gains[:10]))
        # ---- 4. select
        cols = SEL.load_cached(SEL.PATS, a.nproc, say)
        pick, msg = SEL.solve(cols, 9, 298, [], [], a.ftlim)
        if pick is None:
            say(f'  select failed: {msg[:60]}'); break
        sj = sum(cols[j]['cost'] for j in pick)
        F = RI.IFleet()
        for i, j in enumerate(sorted(pick, key=lambda j: (-len(cols[j]['asts']), cols[j]['tank']))):
            st = FA.load_st(cols[j]['f']); F.routes[f'r{i + 1}'] = dict(st=st, tank=float(RI.ipr(st).tank()))
        d = out / f'it{it + 1}' / 'fleet'; F.save(d, note=f's16_iter it{it + 1}')
        cov = len(set(F.coverage()) - {131, 144})
        hist.append(dict(it=it + 1, sumJi0=round(sj0, 5), sumJi=round(sj, 5), covered=cov, fleet=str(d.relative_to(ROOT)),
                         picks=[cols[j]['f'] for j in pick], wall=round(time.time() - tic)))
        json.dump(hist, open(out / 'iter.json', 'w'), indent=1)
        say(f'  SELECT it{it + 1}: sum J_i {sj0:.4f} -> {sj:.4f}, covered {cov} ({time.time() - tic:.0f} s); picks '
            + ', '.join(pathlib.Path(cols[j]['f']).stem for j in pick))
        spec = hist[-1]['fleet']
        if sj > sj0 - a.min_gain:
            say('  converged'); break


if __name__ == '__main__':
    main()
