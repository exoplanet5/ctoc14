"""Re-grid an impulsive twin onto REGULAR dtau bins before export (stage 15b export fix).

ImpulsiveProblem.tcap is what 0.43 N delivers over ONE dtau bin, but it is enforced PER IMPULSE.  Routes born from the
planner (impulsive.from_tour / settle_tour) carry extra irregular nodes -- junction and linear-leg impulses a fraction of
a day after other nodes -- so two impulses inside one bin can each sit at the cap: twice the continuous-thrust
capability.  Their twin tanks are optimistic and impulsive.to_exact (which maps each node's impulse onto its own bin)
starts 1e7-3e8 km off; the long routes then diverge in export (stage 15b: h01-h04, h06 of the 9-craft fleet).
Routes imported from a thrust problem (t10d) only ever had regular nodes, which is why they converted.

regrid: sum every impulse into the regular bin that contains it, then run_ialns.settle (restore + optimise + enforce_cap
+ restore), now with the cap meaning one impulse per window.  Output: honest twin tank, clean input for to_exact.

Usage: s15b_regrid.py fleet_dir out_dir [--names h01,h02] [--dtau 20] [--nproc 5]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, time, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from ctoc14.impulsive import ImpulsiveProblem, to_exact
from ctoc14.constants import DAY


def regrid(ip, dtau=20 * DAY):
    t_end = max(float(np.max(ip.tf)), float(ip.ts[-1])) + dtau
    edges = np.arange(ip.tL, t_end + dtau, dtau); ts = 0.5 * (edges[:-1] + edges[1:])
    b = np.clip(np.searchsorted(edges, ip.ts, side='right') - 1, 0, len(ts) - 1)
    Ts = np.zeros((len(ts), 3)); np.add.at(Ts, b, ip.Ts)
    return ImpulsiveProblem(ip.eph, ip.tL, ip.vinf, ts, Ts, ip.tf, ip.asts, margin=ip.margin, dtau=dtau)


def job(args):
    name, st, dtau, iters = args
    t0 = time.time(); ip = RI.ipr(st); tank0 = float(ip.tank())
    g = regrid(ip, dtau * DAY)
    over0 = float((np.linalg.norm(g.Ts, axis=1) / g.tcap).max())
    miss = RI.settle(g, iters)
    over1 = float((np.linalg.norm(g.Ts, axis=1) / g.tcap).max())
    gm = to_exact(g, RI.eph()); Yg, _ = gm.integrate(dense_samples=True); dg, _ = gm.misses(Yg)
    return name, RI.ist(g), tank0, float(g.tank()), float(miss), over0, over1, float(np.linalg.norm(dg, axis=1).max()), time.time() - t0


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('src'); ap.add_argument('out')
    ap.add_argument('--names', default=''); ap.add_argument('--dtau', type=float, default=20.0)
    ap.add_argument('--iters', type=int, default=150); ap.add_argument('--nproc', type=int, default=5)
    a = ap.parse_args()
    fl = RI.IFleet(a.src); names = a.names.split(',') if a.names else sorted(fl.routes)
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    with mp.get_context('fork').Pool(a.nproc) as pl:
        res = pl.map(job, [(n, fl.routes[n]['st'], a.dtau, a.iters) for n in names])
    J = lambda m: 1 + (m - 600) / 1400 + ((m - 600) / 1400) ** 2
    for n, st, t0, t1, miss, o0, o1, dex, dt in res:
        np.savez(out / f'route_{n}.npz', **st)
        print(f'{n}: twin {t0:.1f} -> regridded {t1:.1f} kg ({t1-t0:+.1f}, dJ {J(t1)-J(t0):+.4f}); miss {miss:.0f} km; '
              f'max imp/cap {o0:.2f} -> {o1:.2f}; to_exact start miss {dex:.2e} km ({dt:.0f} s)', flush=True)


if __name__ == '__main__':
    main()
