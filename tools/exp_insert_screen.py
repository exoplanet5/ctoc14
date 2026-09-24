"""Screen insertions of a victim craft's unique targets into host trajectories with the whole-trajectory linear model.

For every host craft H (all craft except the victims) and every target X unique to a victim: close approaches of H's
current trajectory to X (local minima of the distance on the sample grid, < --dmax AU); for each candidate time the
linearised minimum-propellant problem of H with the extra flyby constraint is solved by IRLS; the estimate is
dF = fuel_lin(with X) - fuel_lin(without). Output: JSON list sorted by dF.
Usage: exp_insert_screen.py submission.txt --victims 11 [--dmax 0.3] [--out screen.json] [--nproc 8]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse, multiprocessing as mp
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from ctoc14.kepler import Ephemeris, propagate_batch
from ctoc14.globalopt import load_craft, GlobalProblem, irls_min_fuel
from ctoc14.constants import DAY, AU, T_MISSION
from ctoc14.validator import parse

EPS_FAST = (0.02, 0.005, 0.001)


def close_approaches(eph, gp, Ys, ast, dmax_km, t_min):
    """Local minima of |r_sc - r_ast| on the sample grid below dmax_km, refined by a parabola. Returns [(t, d)]."""
    k = ast - 1
    ra, _ = eph.ast_states_at(np.full(len(gp.ts), k), gp.ts)
    d = np.linalg.norm(Ys[:, 0:3] - ra, axis=1)
    out = []
    for i in range(1, len(d) - 1):
        if d[i] <= d[i - 1] and d[i] < d[i + 1] and d[i] < dmax_km and gp.ts[i] >= t_min:
            a, b, c = d[i - 1], d[i], d[i + 1]; den = a - 2 * b + c
            off = 0.5 * (a - c) / den if den > 0 else 0.0
            out.append((float(gp.ts[i] + off * (gp.ts[1] - gp.ts[0])), float(b), i))
    return out


def screen_host(job):
    src, sc, targets, a = job
    eph = Ephemeris(); c = load_craft(src, sc); gp = GlobalProblem(eph, c, hstep=0.25 * DAY)
    Yf, Ys = gp.integrate(dense_samples=True); dm, vr = gp.misses(Yf)
    ch = gp.chain(Yf, Ys); kind = ch['kind']
    node_of_f = {int(-2 - kind[n]): n for n in range(len(kind)) if kind[n] <= -2}
    B, Lv = gp.rows_at(ch, gp.tf, nodes=[node_of_f[j] for j in range(len(gp.tf))])
    T0, _, _ = irls_min_fuel(B, Lv, vr, dm, gp.Ts, gp.vinf, gp.tcap, rho=0.0, eps_sched=EPS_FAST, inner=3)
    f0 = gp.fuel(T0)
    res = []
    for X in targets:
        if X in gp.asts:
            continue
        for (tc, dist, i) in close_approaches(eph, gp, Ys, X, a.dmax * AU, gp.tL + 20 * DAY):
            if tc > gp.ts[-1] - 2 * DAY or tc > T_MISSION - 2 * DAY:
                continue
            rs, vs = propagate_batch(Ys[i, 0:3][None], Ys[i, 3:6][None], np.array([tc - gp.ts[i]]))
            ra, va = eph.ast_states_at(np.array([X - 1]), np.array([tc]))
            Bc, Lvc = gp.rows_at(ch, [tc])
            B2 = np.concatenate([B, Bc]); Lv2 = np.concatenate([Lv, Lvc])
            dm2 = np.concatenate([dm, rs - ra]); vr2 = np.concatenate([vr, vs - va])
            T1, dt, dv = irls_min_fuel(B2, Lv2, vr2, dm2, gp.Ts, gp.vinf, gp.tcap, rho=0.0, eps_sched=EPS_FAST, inner=3)
            f1 = gp.fuel(T1)
            res.append(dict(host=sc, ast=int(X), t=tc / DAY, dist_au=dist / AU, dF=f1 - f0, fuel_lin=f1, fuel_base_lin=f0,
                            fuel_cur=gp.fuel(), vrel=float(np.linalg.norm(vs - va)), n_host=len(gp.asts)))
    return res


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('submission'); ap.add_argument('--victims', required=True)
    ap.add_argument('--hosts', default=None); ap.add_argument('--dmax', type=float, default=0.3)
    ap.add_argument('--out', default=None); ap.add_argument('--nproc', type=int, default=8)
    a = ap.parse_args()
    rows = parse(a.submission)
    cov = {}
    for r in rows:
        if r.event == 3: cov.setdefault(r.ast, set()).add(r.sc)
    victims = [int(x) for x in a.victims.split(',')]
    scs = sorted({r.sc for r in rows})
    hosts = [s for s in scs if s not in victims] if a.hosts is None else [int(x) for x in a.hosts.split(',')]
    targets = sorted(x for x, S in cov.items() if S <= set(victims))
    print(f'victims {victims}: {len(targets)} targets only they cover: {targets}', flush=True)
    jobs = [(a.submission, h, targets, a) for h in hosts]
    allres = []
    with mp.get_context('fork').Pool(a.nproc) as pool:
        for r in pool.imap_unordered(screen_host, jobs):
            allres += r
            if r: print(f"host {r[0]['host']}: {len(r)} candidates, best {min(x['dF'] for x in r):.1f} kg", flush=True)
    allres.sort(key=lambda x: x['dF'])
    out = a.out or f'results/newgen/screen_v{a.victims.replace(",", "_")}.json'
    json.dump(allres, open(out, 'w'), indent=1)
    best = {}
    for x in allres:
        best.setdefault(x['ast'], x)
    for X in targets:
        b = best.get(X)
        print(f'target {X}: ' + (f"best host {b['host']} t={b['t']:.0f} d dist {b['dist_au']:.3f} AU vrel {b['vrel']:.1f} dF_lin {b['dF']:.1f} kg"
                                 if b else 'no candidate'))


if __name__ == '__main__':
    main()
