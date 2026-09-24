"""Stage-13 linchpin, step 1: target-to-trajectory proximity for every route of a settled fleet.

For every route (impulsive twin, run_ialns.trajectory, Kepler arcs sampled every 2 d) and every reachable target,
all LOCAL MINIMA of the craft-asteroid distance below --dmax AU (parabolic refinement, exactly run_ialns.w_cands),
including the route's own targets (whose minimum is ~0 at their flyby).  Output JSON:
    {"routes": {name: {"n": flybys, "tank": kg, "asts": [...]}}, "minima": {name: {target: [[t_s, d_au], ...]}}}
Usage: linchpin_dist.py fleet_dir out.json [--dmax 1.0] [--nproc 8]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from ctoc14.constants import AU
from ctoc14.search import UNREACHABLE

ALL = [t for t in range(1, 301) if t not in UNREACHABLE]


def minima(job):
    name, st, dmax = job
    E = RI.eph(); ip = RI.ipr(st); tt, r = RI.trajectory(ip)
    out = {}
    for X in ALL:
        ra, _ = E.ast_states_at(np.full(len(tt), X - 1), tt)
        d = np.linalg.norm(r - ra, axis=1) / AU
        i = np.where((d[1:-1] <= d[:-2]) & (d[1:-1] < d[2:]) & (d[1:-1] < dmax))[0] + 1
        lst = []
        for k in i:
            a, b, c = d[k - 1], d[k], d[k + 1]; den = a - 2 * b + c
            off = 0.5 * (a - c) / den if den > 0 else 0.0
            dm = b - 0.125 * (a - c) ** 2 / den if den > 0 else b
            lst.append([float(tt[k] + off * (tt[1] - tt[0])), float(max(dm, 0.0))])
        out[X] = lst
    return name, out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('fleet'); ap.add_argument('out')
    ap.add_argument('--dmax', type=float, default=1.0); ap.add_argument('--nproc', type=int, default=8)
    a = ap.parse_args()
    F = RI.IFleet(a.fleet); RI.eph()
    jobs = [(n, r['st'], a.dmax) for n, r in sorted(F.routes.items())]
    with mp.get_context('fork').Pool(min(a.nproc, len(jobs))) as pool:
        res = dict(pool.map(minima, jobs))
    json.dump(dict(fleet=a.fleet, dmax=a.dmax,
                   routes={n: dict(n=len(r['st']['asts']), tank=r['tank'], asts=[int(x) for x in r['st']['asts']],
                                   tf=[float(x) for x in r['st']['tf']]) for n, r in F.routes.items()},
                   minima={n: {str(k): v for k, v in m.items()} for n, m in res.items()}), open(a.out, 'w'))
    print(f'wrote {a.out}: {len(res)} routes x {len(ALL)} targets')


if __name__ == '__main__':
    main()
