"""How many targets does a settled route pass CLOSE to?  The quantity that actually caps a route's depth.

A deepened route stops growing when it runs out of untaken targets whose closest approach to its trajectory is
inside the affordable radius.  With the fitted surrogate the affordable radius at the contest breakeven
(0.0369 J per flyby = 34 kg at 950 kg) is only d ~ 0.03-0.04 AU, so a route's REACHABLE DEPTH is essentially
'#targets within 0.035 AU of its trajectory, at an epoch at least `gap` days from an existing flyby'.

This measures exactly that for every settled route of one or more fleets, so the drifting-carrier question
('does sweeping through more of the catalogue's orbital space buy candidate density?') gets a direct answer
instead of an inference from flyby counts.

Usage: cand_density.py fleet_dir [fleet_dir ...] [--dmax 0.12] [--gap 10] [--nproc 3]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, pathlib, argparse, collections, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from insert_model import predict_kg, predict_ok, predict_dJ
from ctoc14.search import UNREACHABLE
from ctoc14.constants import DAY
ALL = [t for t in range(1, 301) if t not in UNREACHABLE]
_A = None
BINS = (0.02, 0.035, 0.06, 0.10)


def w(job):
    tag, name, st, tank = job
    have = set(int(x) for x in st['asts'])
    left = [t for t in ALL if t not in have]
    C = RI.w_cands((name, st, left, _A.dmax))
    tf = np.asarray(st['tf'], float); depth = len(have)
    best = {}
    for c in C:
        m = float(np.abs(tf - c['t']).min() / DAY)
        if m < _A.gap: continue
        cur = best.get(c['ast'])
        if cur is None or c['dist'] < cur[0]: best[c['ast']] = (c['dist'], m)
    n = [sum(1 for d, m in best.values() if d < b) for b in BINS]
    # affordable set: predicted kg under the contest breakeven AND p_ok above pmin
    x = (tank - 600.0) / 1400.0; cap = _A.breakeven * 1400.0 / (1.0 + 2.0 * x)
    aff = [t for t, (d, m) in best.items()
           if predict_kg(d, depth, tank, m) <= cap and predict_ok(d, depth, tank, m) >= _A.pmin]
    eaff = sum(predict_ok(best[t][0], depth, tank, best[t][1]) for t in aff)
    return dict(tag=tag, name=name, depth=depth, tank=tank, Jfb=RI.cost(tank) / depth,
                n002=n[0], n0035=n[1], n006=n[2], n010=n[3], n_aff=len(aff), e_aff=float(eaff),
                span=float((tf.max() - tf.min()) / DAY), cap_kg=float(cap))


def main():
    global _A
    ap = argparse.ArgumentParser(); ap.add_argument('dirs', nargs='+')
    ap.add_argument('--dmax', type=float, default=0.12); ap.add_argument('--gap', type=float, default=10.0)
    ap.add_argument('--breakeven', type=float, default=0.0369); ap.add_argument('--pmin', type=float, default=0.30)
    ap.add_argument('--nproc', type=int, default=3); ap.add_argument('--out', default='')
    _A = a = ap.parse_args()
    jobs = []
    for d in a.dirs:
        fl = RI.IFleet(d); tag = pathlib.Path(d).parts[-2] if pathlib.Path(d).name == 'fleet' else pathlib.Path(d).name
        for n, r in sorted(fl.routes.items()): jobs.append((tag, n, r['st'], r['tank']))
    with mp.get_context('fork').Pool(a.nproc) as pl:
        rows = pl.map(w, jobs)
    by = collections.defaultdict(list)
    for r in rows: by[r['tag']].append(r)
    print(f'{"fleet":10s} {"n":>3s} {"depth":>6s} {"tank":>6s} {"J/fb":>7s} {"<.02":>6s} {"<.035":>6s} '
          f'{"<.06":>6s} {"<.10":>6s} {"afford":>7s} {"E[aff]":>7s} {"span d":>7s}')
    for tag, R in sorted(by.items()):
        f = lambda k: np.median([x[k] for x in R])
        print(f'{tag:10s} {len(R):3d} {f("depth"):6.1f} {f("tank"):6.0f} {f("Jfb"):7.4f} {f("n002"):6.1f} '
              f'{f("n0035"):6.1f} {f("n006"):6.1f} {f("n010"):6.1f} {f("n_aff"):7.1f} {f("e_aff"):7.1f} {f("span"):7.0f}')
    print()
    for r in sorted(rows, key=lambda x: (x['tag'], -x['n_aff'])):
        print(f'  {r["tag"]:9s} {r["name"]:6s} depth {r["depth"]:3d} tank {r["tank"]:5.0f} J/fb {r["Jfb"]:.4f}  '
              f'close {r["n002"]:3d}/{r["n0035"]:3d}/{r["n006"]:3d}/{r["n010"]:3d}  affordable {r["n_aff"]:3d} '
              f'(E {r["e_aff"]:.1f}) cap {r["cap_kg"]:.0f} kg')
    if a.out: json.dump(rows, open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()
