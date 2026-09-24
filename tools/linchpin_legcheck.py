"""Stage-13 linchpin diagnostic: can the beam planner's leg model represent the legs of a settled fleet's routes?

For every leg (previous flyby -> next flyby; leg 0 = launch -> first flyby) of every route, using the TWIN's actual
state at the departure flyby (best case for the planner, which only knows its own Lambert/linear arrival velocity):
  tof, twin dv spent on the leg (impulse magnitudes between the two flybys), coast miss (ballistic propagation of the
  departure state to the arrival epoch vs the target position, AU), single-rev Lambert junction dv, best 0-2 rev
  Lambert junction dv, linear continuous-thrust cost (ctoc14/linleg.LinLeg, ub 0.8, a = 0.5 N / --mass).
Planner admissibility (ctoc14/search.expand at production settings): Lambert leg if 15 <= tof <= 400 d and junction dv
<= min(dv_max, 0.6 a tof); linear leg if 20 <= tof <= 600 d, coast miss <= lin_drmax and cost <= dv_max.
Usage: linchpin_legcheck.py fleet_dir out.json [--dvmax 1.2] [--drmax 0.15] [--mass 1300]
"""
import sys, json, pathlib, argparse
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from ctoc14.kepler import propagate_batch
from ctoc14.lambert import lambert, lambert_best
from ctoc14.linleg import LinLeg
from ctoc14.constants import DAY, AU, TMAX


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('fleet'); ap.add_argument('out')
    ap.add_argument('--dvmax', type=float, default=1.2); ap.add_argument('--drmax', type=float, default=0.15)
    ap.add_argument('--mass', type=float, default=1300.0)
    a = ap.parse_args()
    E = RI.eph(); F = RI.IFleet(a.fleet); acc = TMAX / a.mass * 1e-3
    out = {}
    tot = dict(n=0, lam=0, lin=0, any=0, any_lax=0)
    for name, r in sorted(F.routes.items()):
        ip = RI.ipr(r['st']); Yf, Ys = ip.integrate()
        o = np.argsort(ip.tf); tf = ip.tf[o]; asts = [ip.asts[k] for k in o]; Y = Yf[o]
        imp = np.linalg.norm(ip.Ts, axis=1)
        legs = []
        for k in range(len(tf)):
            if k == 0:
                t1 = ip.tL; r1 = ip.rE; v1 = ip.vE + ip.vinf
            else:
                t1 = tf[k - 1]; r1 = Y[k - 1, :3]; v1 = Y[k - 1, 3:6]
            t2 = tf[k]; tof = t2 - t1
            R2, _ = E.ast_states_at(np.array([asts[k] - 1]), np.array([t2])); R2 = R2[0]
            dv_twin = float(imp[(ip.ts > t1) & (ip.ts <= t2)].sum())
            rc, _ = propagate_batch(r1[None], v1[None], np.array([tof])); miss = float(np.linalg.norm(rc[0] - R2) / AU)
            vl, _ = lambert(r1[None], R2[None], np.array([tof])); dv1 = float(np.linalg.norm(vl[0] - v1))
            vb, _, nrev, _ = lambert_best(r1[None], R2[None], np.array([tof]), v1[None], nrev_max=2)
            dvb = float(np.linalg.norm(vb[0] - v1))
            L = LinLeg(r1, v1, np.array([tof]), dtau=10 * DAY)
            c, _, _ = L.solve(0, R2[None], acc, ub=0.8); dvl = float(c[0])
            td = tof / DAY
            ok_lam = (15 <= td <= 400) and np.isfinite(dv1) and dv1 <= min(a.dvmax, 0.6 * acc * tof)
            ok_lin = (20 <= td <= 600) and miss <= a.drmax and np.isfinite(dvl) and dvl <= a.dvmax
            ok_lax = (np.isfinite(dvb) and dvb <= 2.5) or (td <= 1000 and miss <= 0.30 and np.isfinite(dvl) and dvl <= 2.5)
            legs.append(dict(k=k, ast=int(asts[k]), tof_d=round(td, 1), dv_twin=round(dv_twin, 3), coast_miss_au=round(miss, 4),
                             lam1_dv=round(dv1, 3), lam_best_dv=round(dvb, 3), lam_best_nrev=int(nrev[0]),
                             lin_dv=round(dvl, 3) if np.isfinite(dvl) else None, ok_lambert=bool(ok_lam), ok_linear=bool(ok_lin),
                             ok_planner=bool(ok_lam or ok_lin), ok_relaxed=bool(ok_lax)))
        n = len(legs); nl = sum(l['ok_lambert'] for l in legs); nn = sum(l['ok_linear'] for l in legs)
        na = sum(l['ok_planner'] for l in legs); nx = sum(l['ok_relaxed'] for l in legs)
        runs = []; cur = 0
        for l in legs:
            cur = cur + 1 if l['ok_planner'] else 0; runs.append(cur)
        out[name] = dict(n=n, ok_lambert=nl, ok_linear=nn, ok_planner=na, ok_relaxed=nx, longest_admissible_run=max(runs),
                         tof_gt300=sum(1 for l in legs if l['tof_d'] > 300), miss_gt015=sum(1 for l in legs if l['coast_miss_au'] > 0.15),
                         legs=legs)
        for kk, v in (('n', n), ('lam', nl), ('lin', nn), ('any', na), ('any_lax', nx)): tot[kk] += v
        print(f'{name}: {n} legs; planner-admissible {na} (Lambert {nl}, linear {nn}); relaxed {nx}; longest admissible '
              f'chain {max(runs)}; legs >300 d {out[name]["tof_gt300"]}, coast miss >0.15 AU {out[name]["miss_gt015"]}; '
              f'median tof {np.median([l["tof_d"] for l in legs]):.0f} d')
    print(f'ALL: {tot["n"]} legs, planner-admissible {tot["any"]} ({tot["any"] / tot["n"]:.0%}), Lambert {tot["lam"]}, '
          f'linear {tot["lin"]}, relaxed {tot["any_lax"]}')
    json.dump(dict(fleet=a.fleet, dvmax=a.dvmax, drmax=a.drmax, mass=a.mass, total=tot, routes=out), open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()
