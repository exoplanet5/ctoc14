"""Calibrate the linearised leg model against (a) the planner's Lambert/kappa estimate and (b) the converter's actual fuel per leg.
Usage: exp_legcalib.py tour.json info.json"""
import sys, json, pathlib, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ctoc14.kepler import Ephemeris
from ctoc14.lambert import lambert
from ctoc14.linleg import LinLeg
from ctoc14.constants import DAY, TMAX, VE
eph = Ephemeris(); tour = json.load(open(sys.argv[1])); info = json.load(open(sys.argv[2]))
act = {l['ast']: l['fuel_leg'] for l in info['legs']}
legs = tour['legs']; m = tour['m0']
rE, vE = eph.earth_state(tour['t_launch']); r_prev, v_prev, t_prev = rE, None, tour['t_launch']
rows = []
for k, l in enumerate(legs):
    R = eph.ast_state(l['ast'] - 1, l['t_flyby'])[0]; tof = l['tof']
    v1, v2 = lambert(r_prev[None], R[None], np.array([tof])); v1, v2 = v1[0], v2[0]
    a = TMAX / m * 1e-3
    if k > 0:
        dv = np.linalg.norm(v1 - v_prev); kappa = 1 + dv / (2 * a * tof); dm_plan = m * (1 - np.exp(-kappa * dv / VE))
        L = LinLeg(r_prev, v_prev, np.array([tof]))
        cost, varr, un = L.solve(0, R[None], a, ub=0.8)
        # also unbounded-ish (ub=1.0) and check non-linear miss for the ub=0.8 profile
        K = L.K[0]; Mk = L.Prv[0, :K]
        miss = np.nan
        if np.isfinite(cost[0]):
            # recover U by re-solving (solve returns norms only) -> quick re-run keeping U
            w = np.ones((1, K)); MMk = np.einsum('kij,klj->kil', Mk, Mk); dR = (R - L.rc[0])[None]
            for it in range(5):
                G = np.einsum('tk,kij->tij', 1 / w, MMk); y = np.linalg.solve(G, dR[..., None])[..., 0]
                U = np.einsum('tk,kji,tj->tki', 1 / w, Mk, y); w = np.maximum(np.linalg.norm(U, axis=-1), 1e-9)
            miss = L.check(0, U[0], R) / np.linalg.norm(dR[0])
        dm_lin = m * (1 - np.exp(-cost[0] / VE)) if np.isfinite(cost[0]) else np.nan
        rows.append((l['ast'], tof / DAY, dv, dm_plan, cost[0], dm_lin, act.get(l['ast'], np.nan), miss, np.linalg.norm(R - L.rc[0]) / 1.496e8))
        m -= dm_plan
    else:
        m -= 0.0
    r_prev, v_prev, t_prev = R, v2, l['t_flyby']
print('ast   tof[d]  dv_lam  fuel_plan  dv_lin  fuel_lin  fuel_actual  nonlin_miss/|dR|  |dR|[AU]')
for r in rows: print('%4d %7.0f %7.2f %9.1f %7.2f %9.1f %11.1f %10.3f %9.3f' % r)
A = np.array(rows, float); ok = np.isfinite(A[:, 4]) & np.isfinite(A[:, 6])
print(f'sum fuel: plan {np.nansum(A[:,3]):.0f}  linear {np.nansum(A[ok,5]):.0f} (feasible {ok.sum()}/{len(A)})  actual {np.nansum(A[:,6]):.0f}  (actual over same legs {np.nansum(A[ok,6]):.0f})')
print(f'median ratio actual/plan {np.nanmedian(A[:,6]/A[:,3]):.2f}, actual/linear {np.nanmedian(A[ok,6]/A[ok,5]):.2f}; median nonlinear miss/|dR| {np.nanmedian(A[:,7]):.3f}')
