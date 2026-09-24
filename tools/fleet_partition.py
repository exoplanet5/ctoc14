"""Build a DISJOINT fleet by partitioning the catalogue first and deepening second.

The measured wall of this problem is not route quality but the PARTITION: a target a route chooses costs
0.032-0.042 J, a target forced onto a finished route costs 0.115-0.148.  Every pool-and-cover attempt has
foundered on redundancy -- 28 guided-deepened routes with 33 flybys each covered only 220 distinct targets
(4.2x redundancy), so the set cover still returns the same 10 craft at J 14.355.

This tool removes the redundancy by construction:

  1. every settled base route is sampled against all 298 catalogue orbits (`run_ialns.w_cands`) to get
     d[route][target] = closest approach of that target to that route's trajectory;
  2. K host routes are chosen by k-medoids on -min(d) so their trajectories are spread;
  3. each target is assigned to the host that passes CLOSEST to it -- a Voronoi partition in
     trajectory-distance space, which is exactly the coordinate the insertion surrogate prices in;
  4. each host is then deepened (tools/deep_guided) with `avoid` = every target that is not its own,
     so the resulting fleet is disjoint by construction and its sum J_i is the fleet objective directly.

The number this produces is honest in both directions: targets no host can actually reach are simply missed and
cost +1 each, so the printed J already contains the partition's failures.

Usage: fleet_partition.py base_fleet_dir out_dir [--k 9] [--dmax 0.35] [--tmax 300] [--nproc 8]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from ctoc14.search import UNREACHABLE
ALL = [t for t in range(1, 301) if t not in UNREACHABLE]
TI = {t: i for i, t in enumerate(ALL)}
BIG = 9.0


def w_dist(job):
    """row of closest approaches (AU) of every catalogue target to one settled route, BIG where unseen."""
    name, st, dmax = job
    row = np.full(len(ALL), BIG)
    have = set(int(x) for x in st['asts'])
    for t in have:
        if t in TI: row[TI[t]] = 0.0
    for c in RI.w_cands((name, st, [t for t in ALL if t not in have], dmax)):
        i = TI[c['ast']]
        if c['dist'] < row[i]: row[i] = c['dist']
    return name, row


def kmedoids(D, k, rng, iters=60):
    """choose k rows so that sum_t min_k D[k, t] is small (targets close to SOME host)."""
    n = len(D)
    best = None
    for _ in range(iters):
        S = [int(rng.integers(n))]
        while len(S) < k:                                # k-means++ style spread
            m = D[S].min(0)
            gain = np.array([np.minimum(D[j], m).sum() for j in range(n)])
            gain[S] = np.inf
            S.append(int(np.argmin(gain)))
        for _ in range(12):                              # local swap improvement
            cur = np.minimum.reduce(D[S]).sum(); moved = False
            for pos in range(k):
                for j in range(n):
                    if j in S: continue
                    T = list(S); T[pos] = j
                    v = np.minimum.reduce(D[T]).sum()
                    if v < cur - 1e-9: cur = v; S = T; moved = True
            if not moved: break
        v = np.minimum.reduce(D[S]).sum()
        if best is None or v < best[0]: best = (v, list(S))
    return best[1], best[0]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('src'); ap.add_argument('out')
    ap.add_argument('--k', type=int, default=9); ap.add_argument('--dmax', type=float, default=0.35)
    ap.add_argument('--tmax', type=float, default=300.0); ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--seed', type=int, default=0); ap.add_argument('--assign-max', type=float, default=BIG)
    ap.add_argument('--hosts', default='', help='comma separated route names, overrides k-medoids')
    ap.add_argument('--deep-args', default='--extra 40 --nfail 10 --recand 2 --pmin 0.25')
    a = ap.parse_args()
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt')
    fl = RI.IFleet(a.src); names = sorted(fl.routes)
    say(f'{len(names)} base routes from {a.src}; sampling closest approaches to {len(ALL)} targets')
    tic = time.time()
    with mp.get_context('fork').Pool(a.nproc) as pl:
        rows = dict(pl.map(w_dist, [(n, fl.routes[n]['st'], a.dmax) for n in names]))
    D = np.array([rows[n] for n in names])
    say(f'distance matrix in {time.time()-tic:.0f} s; median min-distance over hosts '
        f'{np.median(D.min(0)):.4f} AU, targets with some host inside 0.035 AU: {(D.min(0)<0.035).sum()}')
    if a.hosts:
        S = [names.index(x) for x in a.hosts.split(',')]
    else:
        S, v = kmedoids(D, a.k, np.random.default_rng(a.seed))
        say(f'k-medoids k={a.k}: sum of min distance {v:.2f} AU')
    hosts = [names[i] for i in S]
    sub = D[S]
    own = sub.argmin(0); dmin = sub.min(0)
    assign = {h: [] for h in hosts}
    for i, t in enumerate(ALL):
        if dmin[i] <= a.assign_max: assign[hosts[own[i]]].append(t)
    say('hosts and their Voronoi cells:')
    for j, h in enumerate(hosts):
        c = assign[h]; dd = sub[j][[TI[t] for t in c]]
        say(f'  {h:12s} base {len(fl.routes[h]["st"]["asts"]):2d} fb @{fl.routes[h]["tank"]:5.0f} kg  cell {len(c):3d} targets, '
            f'{(dd<0.035).sum():3d} inside 0.035 AU, {(dd<0.06).sum():3d} inside 0.06')
    AV = {h: sorted(set(ALL) - set(assign[h])) for h in hosts}
    json.dump(AV, open(out / 'avoid.json', 'w'))
    json.dump({h: assign[h] for h in hosts}, open(out / 'cells.json', 'w'))
    say(f'-> {out}/avoid.json; now run:\n  deep_guided.py {a.src} {out}/fleet_deep --avoid {out}/avoid.json '
        f'--only {",".join(hosts)} --tmax {a.tmax:.0f} --nproc {a.nproc} {a.deep_args}')


if __name__ == '__main__':
    main()
