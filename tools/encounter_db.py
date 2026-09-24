"""Encounter database + element-space transfer cost (stage 9, the layer every GTOC winner has and we lack).

JPL's GTOC9 write-up: "databases of transfers between all bodies on a fine time grid are made, containing an
easy-to-compute yet accurate estimate of the transfer Delta-V"; Tsinghua's GTOC11: a database of rendezvous times
"for all the 83453 asteroids at different phase angles".  Every winner turns trajectory design into TABLE LOOKUP,
then spends the compute on the combinatorics.  We re-solve Lambert inside every beam node, which is why our fleet
layer only ever sees ~500 finished routes instead of millions of partitions.

This builds the two tables:

  ENCOUNTERS  every epoch at which a target can be met cheaply by a craft near 1 AU -- the local minima of |z| over
              each pass of the reachable ring (r in [0.75, 1.35] AU).  ~10-25 per target over 15 years.
  COST        the element-space metric of ctoc14/phasemodel (calibrated to 2 % on the 10-craft fleet):
              dv = K_TOTAL [ K_DRIFT |d drift| + K_PLANE |d i_vec| ], i.e. what it costs a craft that has just met
              encounter p to meet encounter q next.  O(1) per pair, no Lambert, no propagation.

Usage: encounter_db.py out.npz [--rmin 0.75 --rmax 1.35 --step 2]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, time, pathlib, argparse
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
from ctoc14.constants import AU, DAY, MU, T_MISSION
from ctoc14.kepler import Ephemeris, load_mea
from ctoc14.search import UNREACHABLE
from ctoc14.phasemodel import K_DRIFT, K_PLANE, K_TOTAL, drift as drift_of


def build(a):
    eph = Ephemeris(); ids, elem = load_mea()
    keep = np.array([int(x) not in UNREACHABLE for x in ids]); ids = ids[keep]
    idx = np.nonzero(keep)[0]
    tt = np.arange(0.0, T_MISSION, a.step * DAY)
    E = []
    for k, j in enumerate(idx):
        R, V = eph.ast_states_at(np.full(len(tt), j), tt)
        r = np.linalg.norm(R, axis=1) / AU
        ok = (r > a.rmin) & (r < a.rmax)
        z = np.abs(R[:, 2]) / AU
        good = ok & (z < a.zmax)
        # SAMPLE each pass, do not take only its |z| minimum: a flyby epoch 10 d away is a different problem, and
        # against a real 10-craft fleet the minima-only database held only 20 % of the actual flybys within 2 d.
        idxs = np.nonzero(good)[0]
        if len(idxs) == 0: continue
        brk = np.nonzero(np.diff(idxs) > 1)[0]
        for seg in np.split(idxs, brk + 1):
            if len(seg) == 0: continue
            take = seg[::max(1, int(round(a.every / a.step)))]
            if seg[np.argmin(z[seg])] not in take: take = np.append(take, seg[np.argmin(z[seg])])
            for i in take:
                E.append((tt[i], int(ids[k]), r[i], R[i, 2] / AU, *R[i], *V[i]))
    E.sort()
    A = np.array(E, float)
    print(f'{len(A)} encounters for {len(ids)} targets '
          f'({len(A)/len(ids):.1f} each), |z| < {a.zmax} AU, r in [{a.rmin}, {a.rmax}] AU')
    return A, ids


def frames(A):
    """per-encounter craft-orbit signature: the craft must match position and velocity direction loosely, so we take
    the ASTEROID's osculating drift rate and orbit normal at the encounter as the state a craft needs there."""
    R = A[:, 4:7]; V = A[:, 7:10]
    h = np.cross(R, V); hn = np.linalg.norm(h, axis=1)
    n = h / hn[:, None]
    rn = np.linalg.norm(R, axis=1)
    sma = 1.0 / (2.0 / rn - np.einsum('ij,ij->i', V, V) / MU)
    sma = np.where(sma > 0, sma, np.nan)          # phasemodel.drift wants a in km, not AU
    return drift_of(sma), n


def cost_matrix_rows(A, dr, nv, i, dtmin, dtmax):
    """element-space dv from encounter i to every later encounter (vectorised)."""
    dt = A[:, 0] - A[i, 0]
    ok = (dt > dtmin * DAY) & (dt < dtmax * DAY) & (A[:, 1] != A[i, 1])
    dv = K_TOTAL * (K_DRIFT * np.abs(dr - dr[i]) + K_PLANE * np.linalg.norm(nv - nv[i], axis=1))
    dv[~ok] = np.inf
    return dv


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('out')
    ap.add_argument('--rmin', type=float, default=0.75); ap.add_argument('--rmax', type=float, default=1.35)
    ap.add_argument('--zmax', type=float, default=0.12); ap.add_argument('--step', type=float, default=2.0); ap.add_argument('--every', type=float, default=6.0)
    a = ap.parse_args()
    tic = time.time()
    A, ids = build(a)
    dr, nv = frames(A)
    np.savez_compressed(a.out, A=A, drift=dr, nvec=nv, ids=ids)
    print(f'-> {a.out}  ({time.time()-tic:.0f} s)')
    # sanity: distribution of the cheapest onward hop
    n = len(A); rng = np.random.default_rng(0); s = rng.choice(n, size=min(400, n), replace=False)
    best = []
    for i in s:
        d = cost_matrix_rows(A, dr, nv, int(i), 15.0, 400.0)
        if np.isfinite(d).any(): best.append(np.sort(d[np.isfinite(d)])[:5].mean())
    print(f'mean of the 5 cheapest onward hops: median {np.median(best):.3f} km/s '
          f'(t10d fleet average 0.49 km/s per flyby)')
