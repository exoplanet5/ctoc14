"""Robust settle of drifting-carrier routes: settle_tour with a longer lag ladder, and a PRUNE fallback.

Why.  A phase-continuous drifting route (tools/drift_dp.py) has more flybys and more plane change than a fixed-carrier
route, so its Lambert-chain seed is more likely to contain one degenerate leg -- and one bad leg fails the whole tour
(13 of 16 failed at 120 iterations, two of them at a miss of only 171 and 2123 km).  This script keeps the verdict
honest -- a route still only counts if the twin converges to miss <= 150 km -- but does not throw the route away for one
bad leg: on failure it prices every leg's Lambert junction, drops the target at the end of the worst one, and retries.
The reported flyby count is always the count that actually settled.

Usage: drift_settle3.py routes.json out_dir [--iters 200] [--prune 5] [--nproc 8] [--tag d]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, pathlib, argparse, time, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np, run_ialns as RI
from ctoc14.impulsive import settle_tour
from ctoc14.lambert import lambert_best
from ctoc14.constants import VE, DAY

_A = None
LADDER = (2.0 * DAY, 1.0 * DAY, 0.5 * DAY, 0.25 * DAY, 0.05 * DAY)


def leg_dv(eph, tL, vinf, ev):
    """junction |dv| of every leg of the Lambert chain (ev = sorted [(t, ast)])."""
    asts = [x for _, x in ev]; tf = np.array([t for t, _ in ev])
    rE, vE = eph.earth_state(tL)
    R, _ = eph.ast_states_at(np.array(asts) - 1, tf)
    r_prev = np.array(rE, float); v_prev = np.array(vE, float) + np.asarray(vinf, float); t_prev = float(tL)
    out = []
    for k in range(len(asts)):
        v1, v2, _, _ = lambert_best(r_prev[None], R[k][None], np.array([tf[k] - t_prev]), v_prev[None], nrev_max=2)
        v1 = v1[0]; v2 = v2[0]
        if not np.all(np.isfinite(v1)):
            out.append(1e9); return np.array(out + [1e9] * (len(asts) - len(out)))
        out.append(0.0 if k == 0 else float(np.linalg.norm(v1 - v_prev)))
        r_prev = R[k]; v_prev = v2; t_prev = tf[k]
    return np.array(out)


def w(job):
    k, r = job; a = _A
    ev = sorted([(float(t), int(x)) for x, t in r['legs']])
    n0 = len(ev); dropped = []
    for attempt in range(a.prune + 1):
        if len(ev) < 4: return k, ('short', float(len(ev)))
        tour = dict(t_launch=float(r['tL']), vinf=[float(q) for q in r['vinf']],
                    legs=[dict(ast=x, t_flyby=t) for t, x in ev])
        try:
            ip, miss, lag = settle_tour(RI.eph(), tour, lambda ip: RI.settle(ip, a.iters), lags=LADDER)
        except Exception as e:
            ip, miss = None, np.inf
        if ip is not None:
            return k, (RI.ist(ip), float(ip.tank()), len(ip.asts), n0, dropped)
        if attempt == a.prune: return k, ('miss', float(miss), n0, dropped)
        try:
            d = leg_dv(RI.eph(), r['tL'], r['vinf'], ev)
        except Exception:
            return k, ('miss', float(miss), n0, dropped)
        j = int(np.argmax(d))
        dropped.append(ev[j][1]); ev = ev[:j] + ev[j + 1:]
    return k, ('miss', np.inf, n0, dropped)


def main():
    global _A
    ap = argparse.ArgumentParser(); ap.add_argument('routes'); ap.add_argument('out')
    ap.add_argument('--iters', type=int, default=200); ap.add_argument('--prune', type=int, default=5)
    ap.add_argument('--nproc', type=int, default=8); ap.add_argument('--tag', default='d')
    _A = a = ap.parse_args()
    R = json.load(open(a.routes)); out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    say = RI.logger(out / 'log.txt'); tic = time.time()
    with mp.get_context('fork').Pool(a.nproc) as pl:
        res = dict(pl.map(w, list(enumerate(R)), chunksize=1))
    fl = RI.IFleet(); rows = []
    for k in sorted(res):
        v = res[k]; r = R[k]
        if v[0] in ('miss', 'error', 'short'):
            say(f'{k:02d}: seed {r.get("seed","?")} {r["n"]} fb / {r.get("ncar",1)} car -> NOT settled ({v[0]} {v[1]:.0f}, dropped {v[3] if len(v)>3 else []})')
            continue
        st, tk, nf, n0, dr = v; J = RI.cost(tk)
        fl.routes[f'{a.tag}{k:02d}'] = dict(st=st, tank=tk)
        np.savez(out / f'route_{a.tag}{k:02d}.npz', **st)
        say(f'{k:02d}: seed {r.get("seed","?")} {nf} fb (of {n0}, dropped {dr}) / {r.get("ncar",1)} car, tank {tk:.0f} kg, '
            f'dv {VE*np.log(tk/601.5):.2f}, J_i {J:.3f}, {J/nf:.4f} J/fb  [model {r.get("dv_phase",0)+r.get("dv_gap",0)+r.get("dv_switch",0):.2f}]')
        rows.append(dict(k=k, seed=r.get('seed'), n=nf, n0=n0, ncar=r.get('ncar', 1), tank=tk, Jfb=J / nf,
                         dropped=dr, dv_true=float(VE * np.log(tk / 601.5)),
                         dv_model=r.get('dv_phase', 0) + r.get('dv_gap', 0) + r.get('dv_switch', 0)))
    if rows:
        arr = np.array([[x['n'], x['tank'], x['Jfb']] for x in rows])
        say(f'\nsettled {len(rows)}/{len(R)} in {time.time()-tic:.0f} s; flybys median {np.median(arr[:,0]):.0f} max {arr[:,0].max():.0f}; '
            f'J/fb best {arr[:,2].min():.4f} median {np.median(arr[:,2]):.4f}')
    if fl.routes:
        fl.save(out / 'fleet', note='phase-continuous drifting carrier routes (pruned settle)')
        say(f'fleet: {len(fl.routes)} routes, {len(fl.coverage())} covered, sum J_i {sum(RI.cost(x["tank"]) for x in fl.routes.values()):.3f}')
    json.dump(rows, open(out / 'settled.json', 'w'), indent=1)
    say('bar: best fixed-carrier BASE route 0.0514 J/fb (this harness), deep fixed route 0.0339, fleet target 0.0369')


if __name__ == '__main__':
    main()


