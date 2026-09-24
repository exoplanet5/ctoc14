"""Stage 15: close CANDIDATE 9-craft fleets built from chains, not only the selector's J-optimal pick.

Round 1 measured that no 9 columns of the pool cover more than 290 targets (s15_select frontier: coverage >= 291 is
infeasible), while the chains on residual pools cover ALL isolated targets at 285-288 (5-6 kept heads of the best
fleet + the chain's taken routes), leaving 10-13 ORDINARY targets.  This tool assembles such fleets, then offers every
miss (nearest first) to every fleet route passing within --close-d AU (s15_isogen.run_close: lin_price screen, twin
insertion, sequential greedy absorption per host, every intermediate route saved as a column), all hosts of all fleets
in one fork Pool(--nproc <= 8).  The selector then sees the absorbed columns.

usage: s15_fleetclose.py OUT --round 1 [--chains C00h6,C02h5,C03h5] [--close-d 0.2] [--nproc 8] [--max-abs 4]
Fleets: OUT/fleets/<chain>/ (RI.IFleet dirs); closer columns: results/s15/cols/r<RR>/fclose/<chain>_<host>/;
report: OUT/fleetclose.json."""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, glob, time, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)
import numpy as np
import run_ialns as RI
import s15_isogen as IG
from greedy_cover import ALL


def load_st(f):
    z = np.load(ROOT / f); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL'])
    return st


def chain_fleet(cdir):
    """(fleet routes {name: (file, st)}, meta) of a finished chain job: its kept heads + the taken route of every step."""
    meta = json.load(open(cdir / 'meta.json'))
    job = meta['job']; C = meta['chain']
    files = {}
    for i, h in enumerate(job.get('heads', [])):
        files[f'h{i + 1}'] = str(pathlib.Path(h).relative_to(ROOT)) if h.startswith('/') else h
    for t in C.get('taken', []):
        tg = sorted(t['targets'])
        m = [c for c in C['cols'] if sorted(c['targets']) == tg]
        if not m:
            continue
        c = min(m, key=lambda c: abs(c['tank'] - t['tank']))
        files[f's{t["k"]}'] = c['f']
    return files, meta


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--round', type=int, default=1)
    ap.add_argument('--chains', default=''); ap.add_argument('--close-d', type=float, default=0.20)
    ap.add_argument('--nproc', type=int, default=8); ap.add_argument('--max-abs', type=int, default=4)
    ap.add_argument('--min-cov', type=int, default=280)
    a = ap.parse_args()
    if a.nproc > 8:
        raise SystemExit('CPU budget: --nproc <= 8')
    out = ROOT / a.out; out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt')
    RI.eph()
    base = ROOT / f'results/s15/cols/r{a.round:02d}'
    names = [x for x in a.chains.split(',') if x] or sorted(p.name for p in list(base.glob('C*')) + list(base.glob('K*'))
                                                         if (p / 'meta.json').exists())
    fleets = {}
    for nm in names:
        files, meta = chain_fleet(base / nm)
        F = RI.IFleet()
        for k, f in files.items():
            st = load_st(f); F.routes[k] = dict(st=st, tank=float(RI.ipr(st).tank()))
        cov = set(F.coverage()); M = sorted(set(ALL) - cov)
        sj = sum(RI.cost(r['tank']) for r in F.routes.values())
        dup = sum(len(h) - 1 for h in F.coverage().values())
        say(f'fleet {nm}: {len(F.routes)} craft, covered {len(cov)} (dup {dup}), sum J_i {sj:.4f}; misses {len(M)} {M}; '
            f'iso missing {sorted(set(M) & set(IG.ISO))}')
        if len(F.routes) > 9 or len(cov) < a.min_cov:
            say(f'  skipped (craft {len(F.routes)}, covered {len(cov)})'); continue
        F.save(out / 'fleets' / nm, note=f'chain {nm} fleet')
        fleets[nm] = dict(files=files, covered=len(cov), sumJi=sj, misses=M)
    # approaches of each fleet's misses to its routes
    jobs = []; seen = {}
    for nm, fl in fleets.items():
        scan = [(k, load_st(f), fl['misses'], a.close_d) for k, f in fl['files'].items()]
        with mp.get_context('fork').Pool(min(a.nproc, len(scan))) as pool:
            R = pool.map(RI.w_cands, scan)
        fl['approach'] = {}
        for (k, st, _, _), c in zip(scan, R):
            tf = np.asarray(st['tf']); d = {}
            for x in c:
                if np.abs(tf - x['t']).min() > 10 * 86400.0 and (x['ast'] not in d or x['dist'] < d[x['ast']]):
                    d[x['ast']] = x['dist']
            fl['approach'][k] = {int(X): round(v, 4) for X, v in d.items()}
            Xs = tuple(X for X, _ in sorted(d.items(), key=lambda kv: kv[1])[:a.max_abs + 1])
            if not Xs:
                continue
            key = (fl['files'][k], Xs)
            if key in seen:
                seen[key]['fleets'].append(nm); continue
            tag = f'{nm}_{k}'
            j = dict(kind='close', tag=tag, out=f'results/s15/cols/r{a.round:02d}/fclose/{tag}',
                     items=[(fl['files'][k], list(Xs))], absorb_dmax=a.close_d, absorb_kg=150.0, item_trials=3,
                     max_abs=a.max_abs, wall_budget=2400.0, fleets=[nm])
            seen[key] = j; jobs.append(j)
        near = sorted({X for k in fl['approach'] for X in fl['approach'][k]})
        say(f'fleet {nm}: misses with an approach <= {a.close_d} AU: {len(near)} of {len(fl["misses"])} '
            f'(none: {sorted(set(fl["misses"]) - set(near))})')
    say(f'{len(jobs)} closer host jobs over {len(fleets)} fleets')
    tic = time.time()
    with mp.get_context('fork').Pool(a.nproc, maxtasksperchild=1) as pool:
        res = pool.map(IG.run_job, [{k: v for k, v in j.items() if k != 'fleets'} for j in jobs], chunksize=1)
    rep = dict(fleets={}, jobs=len(jobs), wall_s=None)
    for j, r in zip(jobs, res):
        for h in r.get('hosts', []):
            for nm in j['fleets']:
                rep['fleets'].setdefault(nm, []).append(dict(host=j['tag'].split('_')[-1], f=h['f'], tried=h['tried'],
                                                             got=h['got'], dkg=h['dkg']))
        say(f'  {j["tag"]}: ' + '; '.join(f'tried {h["tried"]} got {h["got"]} (+{h["dkg"]} kg)' for h in r.get('hosts', []))
            + (f' ERROR {r.get("error", "")[-120:]}' if r.get('error') else ''))
    for nm, fl in fleets.items():
        got = {}
        for h in rep['fleets'].get(nm, []):
            for X in h['got']:
                got.setdefault(X, []).append(h['host'])
        left = sorted(set(fl['misses']) - set(got))
        say(f'fleet {nm}: {fl["covered"]} covered @ {fl["sumJi"]:.4f}; closer placed {len(got)} of {len(fl["misses"])} '
            f'somewhere ({sorted(got)}); never placed {left}')
        fl['placed'] = sorted(got); fl['unplaced'] = left
    rep['summary'] = {nm: {k: v for k, v in fl.items() if k in ('covered', 'sumJi', 'misses', 'placed', 'unplaced', 'approach')}
                      for nm, fl in fleets.items()}
    rep['wall_s'] = round(time.time() - tic)
    IG.jdump(rep, out / 'fleetclose.json')
    say(f'fleetclose done ({rep["wall_s"]} s)')


if __name__ == '__main__':
    main()
