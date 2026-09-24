"""Settle tools/drift_dp.py routes in the impulsive twin and report J per flyby (the only verdict that counts).

A phase-continuous drifting route is ONE smooth trajectory, so it settles like any tour: Lambert chain seed with the
linearised-leg fallback (ctoc14.impulsive.settle_tour), then the SCP (run_ialns.settle).  A route only counts if the
twin converges to miss <= 150 km.

Usage: drift_settle2.py routes.json out_dir [--iters 120] [--nproc 8]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, pathlib, argparse, time, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np, run_ialns as RI
from ctoc14.impulsive import settle_tour
from ctoc14.constants import VE, DAY

_A = None


def w(job):
    k, r = job
    ev = sorted([(float(t), int(x)) for x, t in r['legs']])
    if len(ev) < 4: return k, ('short', 0.0)
    tour = dict(t_launch=float(r['tL']), vinf=[float(q) for q in r['vinf']],
                legs=[dict(ast=x, t_flyby=t) for t, x in ev])
    try:
        ip, miss, lag = settle_tour(RI.eph(), tour, lambda ip: RI.settle(ip, _A.iters))
    except Exception as e:
        return k, ('error', repr(e)[:70])
    if ip is None: return k, ('miss', float(miss))
    return k, (RI.ist(ip), float(ip.tank()), len(ip.asts))


def main():
    global _A
    ap = argparse.ArgumentParser(); ap.add_argument('routes'); ap.add_argument('out')
    ap.add_argument('--iters', type=int, default=120); ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--tag', default='d')
    _A = a = ap.parse_args()
    R = json.load(open(a.routes)); out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    say = RI.logger(out / 'log.txt')
    tic = time.time()
    with mp.get_context('fork').Pool(a.nproc) as pl:
        res = dict(pl.map(w, list(enumerate(R))))
    fl = RI.IFleet(); rows = []
    for k in sorted(res):
        v = res[k]; r = R[k]
        if v[0] in ('miss', 'error', 'short'):
            say(f'{k:02d}: seed {r.get("seed","?")} {r["n"]} fb / {r.get("ncar",1)} carriers -> NOT settled ({v[0]} {v[1]})')
            continue
        st, tk, nf = v; J = RI.cost(tk)
        fl.routes[f'{a.tag}{k:02d}'] = dict(st=st, tank=tk)
        np.savez(out / f'route_{a.tag}{k:02d}.npz', **st)
        say(f'{k:02d}: seed {r.get("seed","?")} {nf} fb / {r.get("ncar",1)} carriers, tank {tk:.0f} kg, '
            f'dv {VE*np.log(tk/601.5):.2f} km/s, J_i {J:.3f}, {J/nf:.4f} J/fb  '
            f'[model dv {r.get("dv_phase",0)+r.get("dv_gap",0)+r.get("dv_switch",0):.2f}]')
        rows.append(dict(k=k, seed=r.get('seed'), n=nf, ncar=r.get('ncar', 1), tank=tk, Jfb=J / nf,
                         nswitch=r.get('nswitch', 0), dv_model=r.get('dv_phase', 0) + r.get('dv_gap', 0) + r.get('dv_switch', 0),
                         dv_true=float(VE * np.log(tk / 601.5))))
    if rows:
        arr = np.array([[x['n'], x['tank'], x['Jfb']] for x in rows])
        say(f'\nsettled {len(rows)}/{len(R)} in {time.time()-tic:.0f} s; flybys median {np.median(arr[:,0]):.0f} max {arr[:,0].max():.0f}; '
            f'J/fb best {arr[:,2].min():.4f} median {np.median(arr[:,2]):.4f}')
    if fl.routes:
        fl.save(out / 'fleet', note='phase-continuous drifting carrier routes')
        say(f'fleet: {len(fl.routes)} routes, {len(fl.coverage())} covered, sum J_i '
            f'{sum(RI.cost(r["tank"]) for r in fl.routes.values()):.3f}')
    json.dump(rows, open(out / 'settled.json', 'w'), indent=1)
    say('bar: best fixed-carrier deep route 0.0339 J/fb, fleet target 0.0369, t10d fleet 0.0415')


if __name__ == '__main__':
    main()
