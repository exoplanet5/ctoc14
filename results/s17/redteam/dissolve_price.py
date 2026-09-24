"""s17 redteam: price dissolving route V of the s16 best twin fleet by insertion into the other 8 routes.

For every target X of V and every other host H: approach minima of H's trajectory to X (run_ialns.w_cands, dmax),
skip approaches within 10 d of H's own flybys, keep <= 3 nearest per (H, X), and price each with the linearised
whole-route twin s14_twinbeam.lin_price at t, t-4 d, t+4 d (as s15b_close._price_host).  Output: one JSON row per
(V, X, H, approach) with dist, lin_dv, lin_kg, lin_dJ, res_km (nonlinear residual of the linear step).

usage: dissolve_price.py V [V ...] [--dmax 0.25] [--nproc 2]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)
import numpy as np
import run_ialns as RI
from ctoc14.constants import DAY, VE

FLEET = ROOT / 'results/s16/best/fleet'
OUT = ROOT / 'results/s17/redteam'


def load():
    F = {}
    for f in sorted(FLEET.glob('route_*.npz')):
        z = np.load(f); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL']); F[f.stem[6:]] = st
    return F


def job(args):
    import s14_twinbeam as TB
    V, host, st, targets, dmax = args
    tic = time.time()
    apps = RI.w_cands((host, st, targets, dmax))
    tf = np.asarray(st['tf'], float)
    apps = [a for a in apps if np.abs(tf - a['t']).min() > 10 * DAY]
    by = {}
    for a in sorted(apps, key=lambda a: a['dist']):
        by.setdefault(a['ast'], [])
        if len(by[a['ast']]) < 3:
            by[a['ast']].append(a)
    ip = RI.ipr(st); tank = float(ip.tank())
    par = dict(st=st, nreg=len(TB.centres(float(st['tL']), float(np.max(st['tf'])))), dv=float(ip.dv()))
    rows = []
    for X, lst in by.items():
        for a in lst:
            best = None
            for dt in (0.0, -4.0, 4.0):
                try:
                    dvm, lin, cm = TB.lin_price(par, X, a['t'] + dt * DAY)
                except Exception:
                    continue
                if np.isfinite(dvm) and (best is None or dvm < best[0]):
                    best = (float(dvm), a['t'] + dt * DAY, float(lin['res_km']))
            if best is None:
                continue
            kg = float(tank * (np.exp(max(best[0], 0.0) / VE) - 1.0))
            rows.append(dict(V=V, host=host, ast=int(X), t_day=best[1] / DAY, dist=float(a['dist']), lin_dv=best[0],
                             lin_kg=kg, lin_dJ=float(RI.cost(tank + kg) - RI.cost(tank)), res_km=best[2],
                             host_n=len(st['asts']), host_tank=tank))
    return dict(V=V, host=host, rows=rows, sec=time.time() - tic, napp=len(apps))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('V', nargs='+'); ap.add_argument('--dmax', type=float, default=0.25)
    ap.add_argument('--nproc', type=int, default=2); ap.add_argument('--hosts', default='')
    a = ap.parse_args()
    RI.eph(); F = load()
    jobs = []
    for V in a.V:
        targets = [int(x) for x in F[V]['asts']]
        for h in F:
            if h == V or (a.hosts and h not in a.hosts.split(',')):
                continue
            jobs.append((V, h, F[h], targets, a.dmax))
    out = OUT / f'dissolve_{"_".join(a.V)}.jsonl'
    with open(out, 'w') as fo, mp.get_context('fork').Pool(a.nproc) as pool:
        for r in pool.imap_unordered(job, jobs):
            for row in r['rows']:
                fo.write(json.dumps(row) + '\n')
            fo.flush()
            print(f'{r["V"]} host {r["host"]}: {r["napp"]} approaches, {len(r["rows"])} priced, {r["sec"]:.0f} s', flush=True)
