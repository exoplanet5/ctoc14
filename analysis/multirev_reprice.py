"""Re-price planned tours with multi-revolution Lambert legs (analysis for analysis/multirev_lambert.md).

For each tour (results/phasing/camp{1,2}/milp_n12/tour_sc*.json) the Lambert chain of insertion.build_tour is rebuilt,
but every leg may use any of the 5 solutions (0,R), (1,L), (1,R), (2,L), (2,R) of lambert_all. Junction i couples the
arrival velocity of leg i-1 with the departure velocity of leg i, so the choice is a chain problem: it is solved
EXACTLY by dynamic programming over the 5 options per leg (Viterbi), minimising the planner's propellant proxy
sum_i dm_i(dv_i; m_i, tof_i) at the baseline mass profile (plus a large penalty per junction violating dv <= eta a tof).
Two other selections are reported for comparison: 'greedy' (sequential min junction dv given the previous leg's choice,
i.e. what lambert_best does inside a forward search) and 'pairwise' (per junction min over both adjacent legs' options
independently - a non-realisable lower bound). The exact propellant of the chosen chain is re-accumulated with
insertion.mass_chain.

Usage: python analysis/multirev_reprice.py [--out analysis/multirev_reprice.json]
"""
import sys, glob, json, time, argparse, pathlib
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ctoc14.constants import MU, DAY, VE, TMAX, VINF_MAX, AU
from ctoc14.kepler import Ephemeris
from ctoc14.lambert import lambert, lambert_all, lambert_tof_min, BRANCH_LEFT
from ctoc14.insertion import tour_from_json, mass_chain, SLACK_TOL

ETA = 0.6
PEN = 1e4      # kg-equivalent penalty per infeasible junction in the DP objective


def dm_of(dv, m, tof):
    """Planner propellant of a junction dv delivered with mass m over a leg of duration tof (search.py rule)."""
    with np.errstate(all='ignore'):
        a = TMAX / m * 1e-3
        kappa = 1.0 + dv / (2.0 * a * tof)
        return m * (1.0 - np.exp(-kappa * dv / VE))


def junction_costs(V1, V2, vE, m_start, tof):
    """Cost tensors of the chain. c0[o]: launch junction with option o of leg 0; c[i][p, o]: junction i between option p
    of leg i-1 and option o of leg i. Returns dv arrays (same shapes) and the DP cost (propellant proxy + penalty)."""
    K, n, _ = V1.shape
    with np.errstate(invalid='ignore'):
        dv0 = np.maximum(0.0, np.linalg.norm(V1[:, 0] - vE, axis=-1) - VINF_MAX)              # (K,)
        dvi = np.linalg.norm(V1[:, 1:][None, :, :, :] - V2[:, :-1][:, None, :, :], axis=-1)    # (K_p, K_o, n-1)
    lim = ETA * TMAX / m_start * 1e-3 * tof                                                    # (n,)
    c0 = dm_of(dv0, m_start[0], tof[0]) + PEN * (dv0 > lim[0] + SLACK_TOL)
    ci = dm_of(dvi, m_start[1:], tof[1:]) + PEN * (dvi > lim[1:] + SLACK_TOL)
    c0 = np.where(np.isfinite(c0), c0, np.inf); ci = np.where(np.isfinite(ci), ci, np.inf)
    return dv0, dvi, c0, ci, lim


def dp_chain(c0, ci):
    """Viterbi over the option chain. Returns the option index per leg."""
    K = c0.shape[0]; n = ci.shape[2] + 1
    best = c0.copy(); arg = np.zeros((n, K), int)
    for i in range(1, n):
        tot = best[:, None] + ci[:, :, i - 1]          # (p, o)
        arg[i] = np.argmin(tot, axis=0); best = tot[arg[i], np.arange(K)]
    o = np.empty(n, int); o[-1] = int(np.argmin(best))
    for i in range(n - 1, 0, -1):
        o[i - 1] = arg[i, o[i]]
    return o, float(best.min())


def greedy_chain(dv0, dvi):
    n = dvi.shape[2] + 1
    o = np.empty(n, int); o[0] = int(np.nanargmin(np.where(np.isfinite(dv0), dv0, np.inf)))
    for i in range(1, n):
        row = np.where(np.isfinite(dvi[o[i - 1], :, i - 1]), dvi[o[i - 1], :, i - 1], np.inf)
        o[i] = int(np.argmin(row))
    return o


def chain_dv(dv0, dvi, o):
    n = len(o)
    dv = np.empty(n); dv[0] = dv0[o[0]]
    for i in range(1, n):
        dv[i] = dvi[o[i - 1], o[i], i - 1]
    return dv


def reprice(eph, path, nrev_max=2):
    tour = json.load(open(path))
    TM = tour_from_json(eph, tour, name=pathlib.Path(path).stem, eta=ETA)
    n = TM.n
    r1 = np.vstack([TM.rE[None, :], TM.r[:-1]]); r2 = TM.r; tof = TM.tof
    V1, V2, nrevs, branches = lambert_all(r1, r2, tof, nrev_max=nrev_max)          # (K, n, 3)
    tmin1 = lambert_tof_min(r1, r2, nrev=1)
    dv0, dvi, c0, ci, lim = junction_costs(V1, V2, TM.vE, TM.m_start, tof)
    K = V1.shape[0]
    # sanity: option 0 chain == baseline
    base_dv = chain_dv(dv0, dvi, np.zeros(n, int))
    assert np.nanmax(np.abs(base_dv - TM.dv)) < 1e-7, (path, np.nanmax(np.abs(base_dv - TM.dv)))
    o_dp, _ = dp_chain(c0, ci)
    o_gr = greedy_chain(dv0, dvi)
    dv_dp = chain_dv(dv0, dvi, o_dp); dv_gr = chain_dv(dv0, dvi, o_gr)
    # pairwise lower bound: per junction min over (p, o)
    dv_pw = np.empty(n); dv_pw[0] = np.nanmin(dv0)
    for i in range(1, n):
        dv_pw[i] = np.nanmin(dvi[:, :, i - 1])
    out = {}
    for lbl, dv, o in (('base', TM.dv, np.zeros(n, int)), ('dp', dv_dp, o_dp), ('greedy', dv_gr, o_gr), ('pairwise', dv_pw, None)):
        dm, ms, slack, m_end = mass_chain(TM.m0, dv, tof, ETA)
        out[lbl] = dict(sum_dv=float(np.sum(dv)), fuel=float(TM.m0 - m_end), n_infeasible=int(np.sum(slack < -SLACK_TOL)),
                        dv=dv.tolist(), n_multirev=int(np.sum(nrevs[o] >= 1)) if o is not None else None,
                        options=[[int(nrevs[k]), 'L' if branches[k] == BRANCH_LEFT else 'R'] for k in o] if o is not None else None)
    # per-leg single-switch saving: leg i switched alone (neighbours at baseline); cost = junction i + junction i+1
    single = []
    for i in range(n):
        opts = []
        for o in range(K):
            if not np.isfinite(V1[o, i, 0]):
                continue
            d_in = dv0[o] if i == 0 else dvi[0, o, i - 1]
            d_out = dvi[o, 0, i] if i < n - 1 else 0.0
            dm_in = dm_of(d_in, TM.m_start[i], tof[i]); dm_out = dm_of(d_out, TM.m_start[i + 1], tof[i + 1]) if i < n - 1 else 0.0
            feas = (d_in <= lim[i] + SLACK_TOL) and (i == n - 1 or d_out <= lim[i + 1] + SLACK_TOL)
            opts.append(dict(o=o, nrev=int(nrevs[o]), br='L' if branches[o] == BRANCH_LEFT else 'R', dv_in=float(d_in), dv_out=float(d_out),
                             dv2=float(d_in + d_out), dm2=float(dm_in + dm_out), feas=bool(feas)))
        base = opts[0]
        cands = [c for c in opts[1:] if c['feas']] or []
        bestc = min(cands, key=lambda c: c['dm2']) if cands else None
        single.append(dict(leg=i, ast=int(TM.asts[i]), tof_d=float(tof[i] / DAY), tmin1_d=float(tmin1[i] / DAY),
                           has_multirev=bool(any(c['nrev'] >= 1 for c in opts)),
                           dv_base=float(TM.dv[i]), dv2_base=base['dv2'], dm2_base=base['dm2'],
                           save_dv=float(base['dv2'] - bestc['dv2']) if bestc else 0.0,
                           save_dm=float(base['dm2'] - bestc['dm2']) if bestc else 0.0,
                           best=(f"{bestc['nrev']}{bestc['br']}" if bestc else None),
                           n_opts=len(opts), n_feas_multirev=len(cands)))
    return dict(name=TM.name, path=path, n=n, m0=float(TM.m0), t_launch_d=float(TM.t_launch / DAY), t_end_d=float(TM.t_end / DAY),
                results=out, legs=single)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='analysis/multirev_reprice.json')
    ap.add_argument('--nrev-max', type=int, default=2)
    args = ap.parse_args()
    eph = Ephemeris()
    paths = sorted(glob.glob('results/phasing/camp1/milp_n12/tour_sc*.json')) + sorted(glob.glob('results/phasing/camp2/milp_n12/tour_sc*.json'))
    res = []
    tic = time.time()
    for p in paths:
        r = reprice(eph, p, nrev_max=args.nrev_max)
        r['camp'] = p.split('/')[2]
        res.append(r)
        b = r['results']['base']; d = r['results']['dp']; g = r['results']['greedy']; w = r['results']['pairwise']
        print(f"{r['camp']}/{r['name']:9s} n {r['n']:2d} m0 {r['m0']:7.1f}  sum dv base {b['sum_dv']:6.3f} -> dp {d['sum_dv']:6.3f} (greedy {g['sum_dv']:6.3f}, pairwise LB {w['sum_dv']:6.3f}) km/s;"
              f"  fuel {b['fuel']:6.1f} -> {d['fuel']:6.1f} kg (greedy {g['fuel']:6.1f});  infeasible {b['n_infeasible']} -> {d['n_infeasible']};  legs nrev>=1: dp {d['n_multirev']}, greedy {g['n_multirev']}"
              f";  legs with a multirev solution {sum(l['has_multirev'] for l in r['legs'])}")
    print(f'{len(paths)} tours re-priced in {time.time() - tic:.1f} s')
    # ---- aggregate
    legs = [l for r in res for l in r['legs']]
    tofs = np.array([l['tof_d'] for l in legs]); has = np.array([l['has_multirev'] for l in legs])
    sdv = np.array([l['save_dv'] for l in legs]); sdm = np.array([l['save_dm'] for l in legs])
    print(f'\nlegs total {len(legs)}; > 150 d: {np.sum(tofs > 150)}; > 250 d: {np.sum(tofs > 250)}; > 400 d: {np.sum(tofs > 400)}')
    print(f'legs with a nrev>=1 solution: {has.sum()} (tof range {tofs[has].min() if has.any() else 0:.0f}-{tofs[has].max() if has.any() else 0:.0f} d); '
          f'tof_min(nrev=1) median over all legs {np.median([l["tmin1_d"] for l in legs]):.0f} d, min {np.min([l["tmin1_d"] for l in legs]):.0f} d')
    for lbl, sel in (('legs > 150 d', tofs > 150), ('legs > 250 d', tofs > 250), ('legs with multirev option', has)):
        if not sel.any():
            print(f'{lbl}: none'); continue
        s = sdv[sel]; t = sdm[sel]
        print(f'{lbl}: {sel.sum()} legs; single-switch saving > 0 in {np.sum(s > 1e-9)} legs; '
              f'dv saving quantiles [50,90,max] {np.percentile(s, 50):.3f} {np.percentile(s, 90):.3f} {s.max():.3f} km/s, sum {s.sum():.3f}; '
              f'propellant saving sum {t.sum():.2f} kg, max {t.max():.2f} kg')
    B = sum(r['results']['base']['fuel'] for r in res); D = sum(r['results']['dp']['fuel'] for r in res); Gs = sum(r['results']['greedy']['fuel'] for r in res)
    BD = sum(r['results']['base']['sum_dv'] for r in res); DD = sum(r['results']['dp']['sum_dv'] for r in res)
    print(f'\nTOTAL over {len(res)} tours: sum dv {BD:.2f} -> {DD:.2f} km/s ({100 * (1 - DD / BD):.2f}%); planner propellant {B:.1f} -> {D:.1f} kg (DP), {Gs:.1f} kg (greedy); '
          f'legs choosing nrev>=1 (DP) {sum(r["results"]["dp"]["n_multirev"] for r in res)} / {len(legs)}')
    json.dump(res, open(args.out, 'w'), indent=1)
    print('written', args.out)


if __name__ == '__main__':
    main()
