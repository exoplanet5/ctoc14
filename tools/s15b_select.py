"""s15b: frontier selection with FORCED / BANNED targets and a no-good list (s15_select.frontier variant).

min sum J_i  s.t.  sum y <= N,  coverage >= K,  every target of --force covered, (--ban-miss) for every listed miss set
M: not all of M uncovered ... in practice used as: find the best K-fleets whose misses differ from the known ones, so the
closer can try misses that pass near some route.  Several solutions per K with --alt (each new solution adds a no-good
cut on its miss set: sum_{t in miss set} c_t >= 1).

usage: s15b_select.py OUT PAT [PAT ...] --K 297 [--force 64] [--alt 3] [--N 9] [--ftlim 300] [--nproc 2]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
_CWD0 = os.getcwd()
os.chdir(ROOT)
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import csr_matrix, hstack, identity
import run_ialns as RI
import s15_select as SS

ALL = SS.ALL; TI = SS.TI; NT = SS.NT
PATS = ['results/s15/cols/**/route_*.npz', 'results/s13/**/route_*.npz', 'results/s14/**/route_*.npz',
        'results/newgen/**/route_*.npz', 'results/s12/cover_clean/fleet/route_*.npz', 'results/s15b/cols/**/route_*.npz',
        'results/s15b/close/**/cols/route_*.npz']


CACHE = ROOT / 'results/s15b/colcache.pkl'


def load_cached(pats, nproc, say, max_tank=1400.0):
    """s15_select.load_columns with a (path, mtime, size) cache of s13_master.load results (the miss check integrates
    every route: ~0.1 s per file)."""
    import glob, fnmatch, pickle, collections, multiprocessing as mp
    import s13_master as SM
    files = sorted(set(f for p in pats for f in glob.glob(p, recursive=True) if not fnmatch.fnmatch(f, '*/relocated/*')))
    try:
        cache = pickle.load(open(CACHE, 'rb'))
    except Exception:
        cache = {}
    def key(f):
        s = os.stat(f); return (f, s.st_mtime, s.st_size)
    todo = [f for f in files if key(f) not in cache]
    if todo:
        with mp.get_context('fork').Pool(nproc) as pool:
            res = pool.map(SM.load, [(f, max_tank) for f in todo], chunksize=8)
        for f, r in zip(todo, res):
            cache[key(f)] = r
        tmp = CACHE.with_suffix('.tmp'); pickle.dump(cache, open(tmp, 'wb')); tmp.replace(CACHE)
    res = [cache[key(f)] for f in files]
    why = collections.Counter(r for _, r in res)
    best = {}
    for c, _ in res:
        if c is None:
            continue
        k = frozenset(c['asts'])
        if k not in best or c['tank'] < best[k]['tank']:
            best[k] = dict(c)
    cols = sorted(best.values(), key=lambda x: (-len(x['asts']), x['tank']))
    for c in cols:
        c['cost'] = RI.cost(c['tank']); c['set'] = frozenset(c['asts'])
    by_t = collections.defaultdict(list)
    for i, c in enumerate(cols):
        for t in c['asts']:
            by_t[t].append(i)
    keep = []; dom = 0
    for j, c2 in enumerate(cols):
        t0 = min(c2['asts'], key=lambda t: len(by_t[t]))
        if any(i != j and cols[i]['cost'] <= c2['cost'] + 1e-12 and len(cols[i]['set']) >= len(c2['set'])
               and c2['set'] <= cols[i]['set'] and (cols[i]['cost'] < c2['cost'] - 1e-12 or len(cols[i]['set']) > len(c2['set']) or i < j)
               for i in by_t[t0]):
            dom += 1; continue
        keep.append(c2)
    say(f'load: {len(files)} files ({len(todo)} new; {dict(why)}), {len(cols)} distinct, {dom} dominated -> {len(keep)} columns')
    return keep


def solve(cols, N, K, force, cuts, tlim):
    A, c = SS.matrices(cols); nr = len(cols)
    Mcov = hstack([A, -identity(NT, format='csr')]).tocsr()
    cons = [LinearConstraint(Mcov, np.zeros(NT), np.full(NT, np.inf)),
            LinearConstraint(csr_matrix(np.concatenate([np.ones(nr), np.zeros(NT)])[None, :]), -np.inf, N),
            LinearConstraint(csr_matrix(np.concatenate([np.zeros(nr), np.ones(NT)])[None, :]), K, np.inf)]
    lo = np.zeros(nr + NT); hi = np.ones(nr + NT)
    for t in force:
        lo[nr + TI[t]] = 1.0
    for ms in cuts:                                   # at least one of this miss set covered
        row = np.zeros(nr + NT)
        for t in ms:
            row[nr + TI[t]] = 1.0
        cons.append(LinearConstraint(csr_matrix(row[None, :]), 1, np.inf))
    r = milp(np.concatenate([c, np.zeros(NT)]), constraints=cons, integrality=np.ones(nr + NT), bounds=Bounds(lo, hi),
             options=dict(time_limit=tlim, mip_rel_gap=1e-4))
    if r.x is None:
        return None, r.message
    return [j for j in range(nr) if r.x[j] > 0.5], r.message


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('pats', nargs='*')
    ap.add_argument('--K', default='297'); ap.add_argument('--force', default='')
    ap.add_argument('--alt', type=int, default=1); ap.add_argument('--N', type=int, default=9)
    ap.add_argument('--ftlim', type=float, default=300.0); ap.add_argument('--nproc', type=int, default=2)
    a = ap.parse_args()
    out = pathlib.Path(a.out) if os.path.isabs(a.out) else pathlib.Path(_CWD0) / a.out
    out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt'); RI.eph()
    cols = load_cached(a.pats or PATS, a.nproc, say)
    force = [int(x) for x in a.force.split(',') if x]
    res = {}
    for K in [int(x) for x in a.K.split(',')]:
        cuts = []; res[K] = []
        for i in range(a.alt):
            pick, msg = solve(cols, a.N, K, force, cuts, a.ftlim)
            if pick is None:
                say(f'  K>={K} force {force} alt {i}: {msg[:60]}'); break
            J, cv, sj = SS.evaluate(cols, pick)
            miss = sorted(set(ALL) - set().union(*[cols[j]['set'] for j in pick]))
            say(f'  K>={K} force {force} alt {i}: {cv} covered sum J_i {sj:.4f}; misses {miss}; depths '
                f'{sorted([len(cols[j]["asts"]) for j in pick], reverse=True)}')
            res[K].append(dict(covered=cv, sumJi=round(sj, 4), misses=miss, picks=[cols[j]['f'] for j in pick]))
            if not miss:
                break
            cuts.append(miss)
    json.dump({str(k): v for k, v in res.items()}, open(out / 'alts.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
