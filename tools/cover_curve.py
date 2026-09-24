"""J-vs-pool-size curve for the exact set cover (stage 11).

Loads every route npz listed in a file (or matched by globs), dedupes by target set, then solves the exact set cover
on random subsamples of growing size.  The SHAPE of J(pool size) is the deliverable: still falling -> buy compute;
plateaued -> the pool is not the binding constraint.

Usage: cover_curve.py --list pool.txt --out results/s11/curve.json [--sizes 100,200,...] [--reps 3]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, glob, time, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from ctoc14.search import UNREACHABLE
import route_cover as RC

ALL = RC.ALL


def load_cols(files, nproc=8, max_tank=1400.0, cache=None):
    if cache and pathlib.Path(cache).exists():
        old = {x['f']: x for x in json.load(open(cache))}
    else:
        old = {}
    need = [f for f in files if f not in old]
    out = [old[f] for f in files if f in old]
    if need:
        with mp.get_context('fork').Pool(nproc) as pl:
            out += [x for x in pl.map(RC.w_load, need, chunksize=8) if x]
    out = [x for x in out if x['tank'] <= max_tank and x['asts']]
    if cache:
        json.dump(out, open(cache, 'w'))
    return out


def dedupe(cols):
    best = {}
    for x in cols:
        k = frozenset(x['asts'])
        if k not in best or x['tank'] < best[k]['tank']:
            best[k] = x
    return sorted(best.values(), key=lambda x: -len(x['asts']))


def solve_quiet(cols, tlim=120.0, gap=0.0005):
    from scipy.optimize import milp, LinearConstraint, Bounds, linprog
    from scipy.sparse import coo_matrix
    nr = len(cols); ti = {t: i for i, t in enumerate(ALL)}; nt = len(ALL)
    c = np.concatenate([[RI.cost(x['tank']) for x in cols], np.ones(nt)])
    R = []; C = []
    for j, x in enumerate(cols):
        for t in x['asts']:
            R.append(ti[t]); C.append(j)
    R += list(range(nt)); C += list(range(nr, nr + nt))
    A = coo_matrix((np.ones(len(R)), (R, C)), shape=(nt, nr + nt)).tocsr()
    lp = linprog(c, A_ub=-A, b_ub=-np.ones(nt), bounds=(0, 1), method='highs')
    res = milp(c, constraints=LinearConstraint(A, np.ones(nt), np.full(nt, np.inf)),
               integrality=np.ones(nr + nt), bounds=Bounds(0, 1),
               options=dict(time_limit=tlim, mip_rel_gap=gap))
    if res.x is None:
        return None
    pick = [j for j in range(nr) if res.x[j] > 0.5]
    miss = int(round(res.x[nr:].sum()))
    return dict(J=float(res.fun) + 2, lp=float(lp.fun) + 2, n=len(pick), miss=miss,
                sumJ=float(sum(RI.cost(cols[j]['tank']) for j in pick)),
                covered=nt - miss, pick=[cols[j]['f'] for j in pick])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('pats', nargs='*')
    ap.add_argument('--list', default='')
    ap.add_argument('--out', default='results/s11/curve.json')
    ap.add_argument('--cache', default='results/s11/cols.json')
    ap.add_argument('--sizes', default='')
    ap.add_argument('--reps', type=int, default=3)
    ap.add_argument('--tlim', type=float, default=120.0)
    ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--max-tank', type=float, default=1400.0)
    a = ap.parse_args()
    files = []
    if a.list:
        files += [l.strip() for l in open(a.list) if l.strip()]
    files += sorted(set(f for p in a.pats for f in glob.glob(p, recursive=True)))
    files = sorted(set(files))
    out = pathlib.Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    cols = load_cols(files, a.nproc, a.max_tank, a.cache)
    cols = dedupe(cols)
    print(f'{len(files)} files -> {len(cols)} distinct routes in {time.time()-t0:.0f} s; '
          f'union {len(set(t for x in cols for t in x["asts"]))}/{len(ALL)}', flush=True)
    sizes = [int(s) for s in a.sizes.split(',')] if a.sizes else \
        sorted(set([100, 200, 300, 400, 500, 600, 800, 1000, 1200, len(cols)]))
    sizes = [s for s in sizes if s <= len(cols)]
    rows = []
    rng = np.random.default_rng(0)
    for s in sizes:
        for rep in range(1 if s >= len(cols) else a.reps):
            sub = cols if s >= len(cols) else [cols[i] for i in rng.choice(len(cols), s, replace=False)]
            if len(set(t for x in sub for t in x['asts'])) < len(ALL):
                pass
            t1 = time.time()
            r = solve_quiet(sub, a.tlim)
            if r is None:
                continue
            r.update(size=s, rep=rep, secs=time.time() - t1)
            rows.append(r)
            print(f'  pool {s:5d} rep{rep}: J {r["J"]:.3f}  LP {r["lp"]:.3f}  craft {r["n"]:2d}  '
                  f'covered {r["covered"]}  sumJ {r["sumJ"]:.3f}  ({r["secs"]:.0f} s)', flush=True)
            json.dump(rows, open(out, 'w'), indent=1)
    print(f'-> {out}')


if __name__ == '__main__':
    main()
