"""P1 (stage 17 redteam): lin-priced insertion ARCS for every column of the honest pool (+ s16 best fleet routes).

Columns: results/s16/honest (status != rejected; all regular or regridded = honest, law 8) + results/s16/best/fleet
route_r1..r9 (honest: 0 irregular nodes, window cap <= 0.99, best.json).  Deduped by target set (lightest tank kept),
n >= --nmin.  For each column: approach minima (run_ialns.w_cands, 2-d sampled track) to every reachable target it
does not fly, dist <= --dmax AU, not within 10 d of its own flybys, <= 2 nearest per target; s14_twinbeam.lin_price
at t (and t +- 4 d when the t price is <= 0.15 J), as s15b_close._price_host.  Arcs with lin dJ <= --keep are stored.
Output: p1/cols.jsonl (one row per column: key, file, n, tank, asts, arcs[{ast, t_day, dist, lin_kg, lin_dJ, res_km}])
Resumable (skips keys already in cols.jsonl).  usage: p1_price.py [--nproc 4] [--nmin 20] [--dmax 0.10] [--keep 0.10]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, signal, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)
import numpy as np
import run_ialns as RI
from ctoc14.constants import DAY, VE

OUT = ROOT / 'results/s17/redteam/p1'
REACH = [x for x in range(1, 301) if x not in (131, 144)]


def columns(nmin):
    cols = {}
    for l in open(ROOT / 'results/s16/honest/index.jsonl'):
        r = json.loads(l)
        if r['status'] == 'rejected':
            continue
        f = f'results/s16/honest/route_{r["key"]}.npz'
        z = np.load(ROOT / f); a = tuple(sorted(int(x) for x in z['asts']))
        cols.setdefault(a, []).append((float(r['tank_out']), r['key'], f))
    for f in sorted((ROOT / 'results/s16/best/fleet').glob('route_*.npz')):
        z = np.load(f); a = tuple(sorted(int(x) for x in z['asts']))
        st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL'])
        cols.setdefault(a, []).append((float(RI.ipr(st).tank()) - 1e-6, 'fleet_' + f.stem[6:], str(f.relative_to(ROOT))))
    out = []
    for a, lst in cols.items():
        lst.sort()
        if len(a) >= nmin:
            out.append(dict(key=lst[0][1], file=lst[0][2], n=len(a)))
    return out


def _alarm(*a):
    raise TimeoutError


def job(args):
    import s14_twinbeam as TB
    col, dmax, keep, timeout = args
    tic = time.time()
    z = np.load(ROOT / col['file']); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL'])
    own = set(int(a) for a in st['asts'])
    ip = RI.ipr(st); tank = float(ip.tank())
    row = dict(col, tank=tank, asts=sorted(own), arcs=[], ok=True)
    try:
        signal.signal(signal.SIGALRM, _alarm); signal.setitimer(signal.ITIMER_REAL, timeout)
        apps = RI.w_cands((col['key'], st, [x for x in REACH if x not in own], dmax))
        tf = np.asarray(st['tf'], float)
        apps = [a for a in apps if np.abs(tf - a['t']).min() > 10 * DAY]
        by = {}
        for a in sorted(apps, key=lambda a: a['dist']):
            by.setdefault(a['ast'], [])
            if len(by[a['ast']]) < 2:
                by[a['ast']].append(a)
        par = dict(st=st, nreg=len(TB.centres(float(st['tL']), float(np.max(st['tf'])))), dv=float(ip.dv()))
        def price(X, t):
            try:
                dvm, lin, cm = TB.lin_price(par, X, t)
            except Exception:
                return None
            if not np.isfinite(dvm):
                return None
            kg = float(tank * (np.exp(max(dvm, 0.0) / VE) - 1.0))
            return dict(t=t, kg=kg, dJ=float(RI.cost(tank + kg) - RI.cost(tank)), res=float(lin['res_km']))
        for X, lst in by.items():
            best = None
            for a in lst:
                p = price(X, a['t'])
                if p is not None and p['dJ'] <= 0.15:
                    for dt in (-4.0, 4.0):
                        q = price(X, a['t'] + dt * DAY)
                        if q is not None and q['dJ'] < p['dJ']:
                            p = q
                if p is not None and (best is None or p['dJ'] < best[0]['dJ']):
                    best = (p, a['dist'])
            if best is not None and best[0]['dJ'] <= keep:
                p, d = best
                row['arcs'].append(dict(ast=int(X), t_day=round(p['t'] / DAY, 2), dist=round(float(d), 4), lin_kg=round(p['kg'], 2),
                                        lin_dJ=round(p['dJ'], 5), res_km=round(p['res'])))
        row['napp'] = len(apps)
    except TimeoutError:
        row['ok'] = False
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    row['sec'] = round(time.time() - tic, 1)
    return row


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--nproc', type=int, default=4); ap.add_argument('--nmin', type=int, default=20)
    ap.add_argument('--dmax', type=float, default=0.10); ap.add_argument('--keep', type=float, default=0.10)
    ap.add_argument('--timeout', type=float, default=240.0); ap.add_argument('--limit', type=int, default=0)
    a = ap.parse_args(); RI.eph()
    cols = columns(a.nmin)
    # fleet routes and deep columns first (most relevant for N = 8)
    cols.sort(key=lambda c: (not c['key'].startswith('fleet_'), -c['n']))
    done = set()
    fo_path = OUT / 'cols.jsonl'
    if fo_path.exists():
        for l in open(fo_path):
            done.add(json.loads(l)['key'])
    todo = [c for c in cols if c['key'] not in done]
    if a.limit:
        todo = todo[:a.limit]
    print(f'{len(cols)} columns (n >= {a.nmin}), {len(done)} done, {len(todo)} to price', flush=True)
    tic = time.time()
    with open(fo_path, 'a') as fo, mp.get_context('fork').Pool(a.nproc, maxtasksperchild=20) as pool:
        for i, r in enumerate(pool.imap_unordered(job, [(c, a.dmax, a.keep, a.timeout) for c in todo], chunksize=1)):
            fo.write(json.dumps(r) + '\n'); fo.flush()
            if i % 20 == 0:
                print(f'[{time.time() - tic:6.0f} s] {i + 1}/{len(todo)} {r["key"]} n {r["n"]} arcs {len(r["arcs"])} '
                      f'({r["sec"]} s, ok {r["ok"]})', flush=True)
    print('done', round(time.time() - tic), 's', flush=True)
