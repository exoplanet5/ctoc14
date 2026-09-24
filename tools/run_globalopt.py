"""Whole-trajectory minimum-propellant re-optimisation of converted craft (fixed flyby sequence), then tank tightening.

For every selected craft of a submission file: SCP (ctoc14.globalopt.optimise) -> feasibility restoration -> tank
m0 = 600 + fuel + margin -> restoration at the new mass -> second SCP stage -> restoration -> tighten again -> rows with the
validator integrator + validation. Every stage is checkpointed (<out>/sc<k>_<stage>.npz); the final fragment is
<out>/go_sc<k>.txt (+ .json summary). Craft whose result does not validate keep no fragment.
Usage: run_globalopt.py submission.txt outdir [--sc 1,2,...] [--iters1 120] [--iters2 60] [--margin 1.5] [--nproc 8]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse, multiprocessing as mp
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from ctoc14.kepler import Ephemeris
from ctoc14.globalopt import load_craft, GlobalProblem, optimise, restore, write_craft, tighten as tighten_gp
from ctoc14.constants import DAY, M_DRY


def run_one(job):
    src, sc, out, a = job
    out = pathlib.Path(out); logf = open(out / f'sc{sc}.log', 'a')
    def log(s):
        logf.write(f'[{time.strftime("%H:%M:%S")}] {s}\n'); logf.flush()
    eph = Ephemeris(); c = load_craft(src, sc)
    gp = GlobalProblem(eph, c, hstep=a.hstep * DAY)
    fuel0 = c['m0'] - c['m_end']; nfb = len(c['flybys'])
    log(f'start {src} sc{sc}: m0 {c["m0"]:.2f}, fuel {fuel0:.2f}, flybys {nfb}')
    def ckpt(stage):
        np.savez(out / f'sc{sc}_{stage}.npz', Ts=gp.Ts, tf=gp.tf, vinf=gp.vinf, m0=gp.m0, ts=gp.ts, asts=np.array(gp.asts), tL=gp.tL)
    def load(stage):
        f = out / f'sc{sc}_{stage}.npz'
        if not f.exists(): return False
        z = np.load(f); gp.Ts, gp.tf, gp.vinf, gp.m0 = z['Ts'], z['tf'], z['vinf'], float(z['m0']); return True
    def tighten():
        return tighten_gp(gp, a.margin, log=log, verbose=True)
    t0 = time.time()
    if not load('s1'):
        optimise(gp, iters=a.iters1, rho0=a.rho0, log=log); restore(gp, log=log); ckpt('s1')
    log(f'stage 1: fuel {gp.fuel():.2f} kg ({time.time() - t0:.0f} s)')
    if not load('s2'):
        tighten(); ckpt('s2')
    log(f'stage 2 (m0 {gp.m0:.2f}): fuel {gp.fuel():.2f} kg')
    if a.iters2 > 0 and not load('s3'):
        optimise(gp, iters=a.iters2, rho0=a.rho0, log=log); restore(gp, log=log); d = tighten(); ckpt('s3')
    log(f'stage 3 (m0 {gp.m0:.2f}): fuel {gp.fuel():.2f} kg ({time.time() - t0:.0f} s)')
    frag = out / f'go_sc{sc}.txt'
    rep, fuel, pk = write_craft(gp, str(frag), eph=eph)
    ok = bool(rep.ok) and len(rep.flybys) >= nfb
    res = dict(src=src, sc=sc, ok=ok, errors=rep.errors[:5], m0=gp.m0, fuel=fuel, fuel0=fuel0, m0_0=c['m0'], flybys=len(rep.flybys),
               flybys0=nfb, peak_thrust=pk, J_i=float(1 + (gp.m0 - 600) / 1400 + ((gp.m0 - 600) / 1400) ** 2),
               J_i0=float(1 + (c['m0'] - 600) / 1400 + ((c['m0'] - 600) / 1400) ** 2), runtime=time.time() - t0)
    json.dump(res, open(out / f'go_sc{sc}.json', 'w'), indent=1)
    if not ok:
        frag.rename(out / f'go_sc{sc}_INVALID.txt')
    log(f'RESULT {res}')
    return res


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('submission'); ap.add_argument('outdir')
    ap.add_argument('--sc', default=None); ap.add_argument('--iters1', type=int, default=120); ap.add_argument('--iters2', type=int, default=60)
    ap.add_argument('--margin', type=float, default=1.5); ap.add_argument('--tol', type=float, default=50.0)
    ap.add_argument('--rho0', type=float, default=100.0); ap.add_argument('--hstep', type=float, default=0.25)
    ap.add_argument('--nproc', type=int, default=8)
    a = ap.parse_args()
    out = pathlib.Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    from ctoc14.validator import parse
    scs = sorted({r.sc for r in parse(a.submission)}) if a.sc is None else [int(x) for x in a.sc.split(',')]
    jobs = [(a.submission, sc, str(out), a) for sc in scs]
    with mp.get_context('fork').Pool(a.nproc, maxtasksperchild=1) as pool:
        for r in pool.imap_unordered(run_one, jobs):
            print(f"sc{r['sc']}: ok={r['ok']} flybys {r['flybys']}/{r['flybys0']} m0 {r['m0_0']:.1f} -> {r['m0']:.1f}, fuel {r['fuel0']:.1f} -> {r['fuel']:.1f}, "
                  f"J_i {r['J_i0']:.4f} -> {r['J_i']:.4f} ({r['runtime']:.0f} s)", flush=True)


if __name__ == '__main__':
    main()
