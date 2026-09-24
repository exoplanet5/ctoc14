"""Stage 10 -- measure the REAL cost of adding one target to a finished (settled) route.

Phase 1 (`cands`): for every host route in a list of impulsive fleets, sample the route trajectory,
find local minima of the distance to every asteroid not already on the route (run_ialns.w_cands),
and stratify the resulting (host, target, epoch, closest-approach) candidates over distance buckets.

Phase 2 (`run`): for each candidate call run_ialns.w_insert -- aim-point homotopy + SCP + full settle --
and record the REAL post-settle tank, i.e. the true kg cost of that insertion, or the failure.

Rows written to a JSONL: host dir/name, depth (flybys before insertion), tank0, dist (AU), t (s),
margin_d (days to the nearest existing flyby of the host), ok, tank1, dkg, miss, sec.
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, glob, pathlib, argparse, random, multiprocessing as mp
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import numpy as np
import run_ialns as RI
from ctoc14.constants import DAY

UNREACH = {131, 144}
BUCKETS = [(0.0, 0.02), (0.02, 0.04), (0.04, 0.06), (0.06, 0.09), (0.09, 0.13),
           (0.13, 0.18), (0.18, 0.25), (0.25, 0.35)]


def host_list(spec):
    """spec = 'dir' or 'dir:route,route' or 'dir:*N' (N random routes)."""
    out = []
    for s in spec:
        if ':' in s:
            d, sel = s.rsplit(':', 1)
        else:
            d, sel = s, ''
        files = sorted(glob.glob(os.path.join(d, 'route_*.npz')))
        names = [os.path.basename(f)[6:-4] for f in files]
        if sel.startswith('*'):
            k = int(sel[1:]); random.shuffle(names); names = sorted(names[:k])
        elif sel:
            names = [n for n in names if n in sel.split(',')]
        for n in names:
            out.append((d, n))
    return out


def w_c(job):
    d, n, dmax = job
    z = np.load(os.path.join(d, f'route_{n}.npz'))
    st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL'])
    have = {int(a) for a in st['asts']}
    T = [x for x in range(1, 301) if x not in have and x not in UNREACH]
    cs = RI.w_cands((n, st, T, dmax))
    tf = np.asarray(st['tf'], float)
    for c in cs:
        c['dir'] = d; c['depth'] = len(have)
        c['margin_d'] = float(np.abs(tf - c['t']).min() / DAY)
    return cs


def cmd_cands(a):
    hosts = host_list(a.hosts)
    print(f'{len(hosts)} hosts', flush=True)
    with mp.get_context('fork').Pool(a.nproc) as pool:
        res = pool.map(w_c, [(d, n, a.dmax) for d, n in hosts])
    jobs = []
    rng = random.Random(a.seed)
    for cs in res:
        if not cs:
            continue
        # one candidate per target per host: the closest approach
        best = {}
        for c in sorted(cs, key=lambda c: c['dist']):
            best.setdefault(c['ast'], c)
        by = {b: [] for b in range(len(BUCKETS))}
        for c in best.values():
            for bi, (lo, hi) in enumerate(BUCKETS):
                if lo <= c['dist'] < hi:
                    by[bi].append(c); break
        for bi, lst in by.items():
            rng.shuffle(lst)
            jobs += lst[:a.per_bucket]
    rng.shuffle(jobs)
    json.dump(jobs, open(a.out, 'w'))
    d0 = np.array([j['dist'] for j in jobs]); dp = np.array([j['depth'] for j in jobs])
    print(f'{len(jobs)} candidate insertions; dist {d0.min():.3f}-{d0.max():.3f}, depth {dp.min()}-{dp.max()}')
    for lo, hi in BUCKETS:
        print(f'  {lo:.2f}-{hi:.2f}: {int(((d0 >= lo) & (d0 < hi)).sum())}')


def w_run(job):
    c = job
    d, n = c['dir'], c['host']
    z = np.load(os.path.join(d, f'route_{n}.npz'))
    st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL'])
    tank0 = float(RI.ipr(st).tank())
    tic = time.time()
    try:
        r = RI.w_insert((n, st, c['ast'], c['t']))
    except Exception as e:
        r = dict(ok=False, miss=float('inf'), error=repr(e)[:80])
    row = dict(dir=d, host=n, ast=int(c['ast']), t=float(c['t']), dist=float(c['dist']),
               depth=int(c['depth']), margin_d=float(c['margin_d']), tank0=tank0,
               ok=bool(r.get('ok')), miss=float(min(r.get('miss', np.inf), 1e12)),
               sec=time.time() - tic)
    if r.get('ok'):
        row['tank1'] = float(r['tank']); row['dkg'] = float(r['tank']) - tank0
    else:
        row['err'] = r.get('error', '')
    return row


def cmd_run(a):
    jobs = json.load(open(a.cands))
    done = set()
    if os.path.exists(a.out) and a.resume:
        for line in open(a.out):
            try:
                r = json.loads(line); done.add((r['dir'], r['host'], r['ast'], round(r['t'], 1)))
            except Exception:
                pass
    jobs = [c for c in jobs if (c['dir'], c['host'], c['ast'], round(c['t'], 1)) not in done]
    if a.limit:
        jobs = jobs[:a.limit]
    print(f'{len(jobs)} trials ({len(done)} already done)', flush=True)
    f = open(a.out, 'a')
    nok = 0; tic = time.time()
    with mp.get_context('fork').Pool(a.nproc, maxtasksperchild=8) as pool:
        for i, row in enumerate(pool.imap_unordered(w_run, jobs), 1):
            f.write(json.dumps(row) + '\n'); f.flush()
            nok += row['ok']
            if i % 10 == 0:
                print(f'  {i}/{len(jobs)} ok {nok} ({nok / i:.0%}) {time.time() - tic:.0f}s', flush=True)
    print(f'done: {len(jobs)} trials, {nok} ok, {time.time() - tic:.0f}s')


def main():
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('cands'); p.add_argument('hosts', nargs='+'); p.add_argument('--out', required=True)
    p.add_argument('--dmax', type=float, default=0.35); p.add_argument('--per-bucket', type=int, default=3)
    p.add_argument('--nproc', type=int, default=10); p.add_argument('--seed', type=int, default=0)
    p.set_defaults(fn=cmd_cands)
    p = sub.add_parser('run'); p.add_argument('cands'); p.add_argument('--out', required=True)
    p.add_argument('--nproc', type=int, default=10); p.add_argument('--limit', type=int, default=0)
    p.add_argument('--resume', action='store_true', default=True)
    p.set_defaults(fn=cmd_run)
    a = ap.parse_args(); a.fn(a)


if __name__ == '__main__':
    main()
