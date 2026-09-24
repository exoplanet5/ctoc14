"""Precomputed LAMBERT TRANSFER GRAPH over encounters (stage 9).

This is the layer every GTOC winner builds and we never did.  JPL, GTOC9: "databases of transfers between all bodies
on a fine time grid ... an easy-to-compute yet accurate estimate of the transfer Delta-V".  Tsinghua, GTOC11: a
database of rendezvous opportunities for 83 453 asteroids, then beam search over 3-8 flyby FRAGMENTS, then a genetic
algorithm to pick the ten motherships.  NUDT, GTOC9 (2nd, and the problem is ours exactly -- partition N objects
among missions, cost = sum of a fixed-plus-variable cost per mission): ant colony optimisation builds the WHOLE
partition in one constructive pass, with pheromone learning across passes.

Our search re-solves Lambert inside every beam node, so a fleet costs hours and the fleet layer only ever sees ~500
finished routes.  With the graph precomputed, a chain is pure table lookup and the compute goes where it belongs --
into the combinatorics of WHICH targets go to WHICH craft.

Graph.  Nodes = encounters (tools/encounter_db.py: local |z| minima of each near-1-AU pass, 11 per target).
Edges p->q = the Lambert arcs from r_p at t_p to r_q at t_q, all branches up to nrev_max; each edge stores the
DEPARTURE and ARRIVAL velocity.  A flyby imposes no velocity constraint, so the cost of a chain is the sum of
junction impulses  |v_dep(p->q) - v_arr(o->p)| -- exact, and O(1) per junction once the table exists.  The search
state is therefore the last EDGE, not the last node (a line graph), which is what makes it a shortest-path problem.

Usage: transfer_db.py enc.npz out.npz [--dtmin 15 --dtmax 500 --dvmax 1.2 --nrev 1]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, time, pathlib, argparse
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
from ctoc14.constants import DAY, AU
from ctoc14.lambert import lambert_all


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('enc'); ap.add_argument('out')
    ap.add_argument('--dtmin', type=float, default=15.0); ap.add_argument('--dtmax', type=float, default=500.0)
    ap.add_argument('--nrev', type=int, default=1); ap.add_argument('--vmax', type=float, default=45.0)
    ap.add_argument('--chunk', type=int, default=400000); ap.add_argument('--per-target', type=int, default=3)
    a = ap.parse_args()
    z = np.load(a.enc); A = z['A']
    T = A[:, 0]; R = A[:, 4:7]; n = len(A)
    TG = A[:, 1].astype(np.int64)
    src = []; dst = []
    for i in range(n):                                     # encounters are time sorted
        lo = np.searchsorted(T, T[i] + a.dtmin * DAY); hi = np.searchsorted(T, T[i] + a.dtmax * DAY)
        j = np.arange(lo, hi)
        j = j[TG[j] != TG[i]]
        if a.per_target and len(j):
            # encounters of one target within a window are near-duplicates: keep the --per-target earliest of each,
            # spread over the window, so the out-degree stays bounded as the time grid gets finer
            o = np.lexsort((T[j], TG[j])); j = j[o]
            g = TG[j]; first = np.concatenate([[True], g[1:] != g[:-1]])
            grp = np.cumsum(first) - 1
            rank = np.arange(len(j)) - np.repeat(np.nonzero(first)[0], np.bincount(grp))
            cntg = np.bincount(grp)[grp]
            stride = np.maximum(1, cntg // a.per_target)
            j = j[(rank % stride == 0) & (rank // np.maximum(stride, 1) < a.per_target)]
        src.append(np.full(len(j), i)); dst.append(j)
    src = np.concatenate(src); dst = np.concatenate(dst)
    print(f'{n} encounters -> {len(src)} candidate edges ({len(src)/n:.0f} per node)', flush=True)
    tic = time.time()
    ES = []; ED = []; E1 = []; E2 = []
    for c0 in range(0, len(src), a.chunk):
        sl = slice(c0, min(c0 + a.chunk, len(src)))
        r1 = R[src[sl]]; r2 = R[dst[sl]]; tof = T[dst[sl]] - T[src[sl]]
        v1, v2, nr, br = lambert_all(r1, r2, tof, nrev_max=a.nrev)
        sp = np.linalg.norm(v1, axis=-1)
        bad = ~np.isfinite(sp) | (sp > a.vmax) | (np.linalg.norm(v2, axis=-1) > a.vmax)
        # keep EVERY revolution branch as its own edge: which one is cheapest depends on the velocity the craft
        # arrives with, which is not known when the table is built (picking one branch here cost ~8 flybys per route)
        for k in range(v1.shape[0]):
            m = ~bad[k]
            if not m.any(): continue
            ES.append(src[sl][m]); ED.append(dst[sl][m]); E1.append(v1[k][m]); E2.append(v2[k][m])
        print(f'  {sl.stop}/{len(src)}  ({time.time()-tic:.0f} s)', flush=True)
    src = np.concatenate(ES); dst = np.concatenate(ED)
    V1 = np.concatenate(E1); V2 = np.concatenate(E2)
    o = np.lexsort((dst, src)); src, dst, V1, V2 = src[o], dst[o], V1[o], V2[o]
    ok = np.ones(len(src), bool)
    ptr = np.searchsorted(src, np.arange(n + 1))
    np.savez_compressed(a.out, src=src.astype(np.int32), dst=dst.astype(np.int32),
                        v1=V1.astype(np.float32), v2=V2.astype(np.float32), ptr=ptr.astype(np.int64),
                        A=A, ids=z['ids'])
    print(f'{len(src)} feasible edges ({100.0*len(src)/len(ok):.0f} %), out-degree median '
          f'{np.median(np.diff(ptr)):.0f} -> {a.out}  ({time.time()-tic:.0f} s)')


if __name__ == '__main__':
    main()
