"""Carrier x target closest-approach table, UNCAPPED.

`results/s7/libopt.pkl` stores each carrier's `neigh` only out to 0.15 AU, but real farm routes take a third of
their extra targets from beyond that.  This builds the full 240 x 298 table of

    d[c,t]   minimum |r_carrier(t) - r_ast(t)| over the mission, in AU   (ballistic carrier, no steering)
    e[c,t]   the epoch of that minimum [s]

on a common 3 d grid, so the fleet layer can price ANY (carrier, target) pair.

Usage: carrier_reach.py [--lib results/s7/libopt.pkl] [--out results/s10/reach.npz] [--step 3]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, time, pickle, pathlib, argparse
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
from ctoc14.constants import AU, DAY, T_MISSION
from ctoc14.kepler import Ephemeris, propagate_batch
from ctoc14.search import UNREACHABLE

ALL = np.array([t for t in range(1, 301) if t not in UNREACHABLE])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--lib', default=str(ROOT / 'results/s7/libopt.pkl'))
    ap.add_argument('--out', default=str(ROOT / 'results/s10/reach.npz'))
    ap.add_argument('--step', type=float, default=3.0)
    a = ap.parse_args()
    eph = Ephemeris(); lib = pickle.load(open(a.lib, 'rb'))
    tt = np.arange(0.0, T_MISSION - DAY, a.step * DAY)
    tic = time.time()
    RA = np.zeros((len(ALL), len(tt), 3))
    for i, x in enumerate(ALL):
        RA[i], _ = eph.ast_states_at(np.full(len(tt), int(x) - 1), tt)
    print(f'asteroid grid {RA.shape} in {time.time()-tic:.0f} s', flush=True)
    D = np.full((len(lib), len(ALL)), 9.9); E = np.zeros((len(lib), len(ALL)))
    for c, car in enumerate(lib):
        tL = float(car['tL'])
        rE, vE = eph.earth_state(tL)
        r0 = np.array(rE, float); v0 = np.array(vE, float) + np.array(car['vinf'], float)
        m = tt > tL + 20 * DAY
        dt = tt[m] - tL
        rc, _ = propagate_batch(np.repeat(r0[None], len(dt), 0), np.repeat(v0[None], len(dt), 0), dt)
        d = np.linalg.norm(RA[:, m, :] - rc[None], axis=2) / AU     # (nast, nt)
        k = d.argmin(1)
        D[c] = d[np.arange(len(ALL)), k]; E[c] = tt[m][k]
        if c % 40 == 0: print(f'  carrier {c} ({time.time()-tic:.0f} s)', flush=True)
    pathlib.Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.out, D=D, E=E, ids=ALL, tL=np.array([c['tL'] for c in lib]))
    print(f'{a.out}: D {D.shape}, median min-d {np.median(D):.3f} AU, '
          f'per carrier <0.05 {np.median((D<0.05).sum(1)):.0f}, <0.10 {np.median((D<0.10).sum(1)):.0f}, '
          f'<0.20 {np.median((D<0.20).sum(1)):.0f}  ({time.time()-tic:.0f} s)')


if __name__ == '__main__':
    main()
