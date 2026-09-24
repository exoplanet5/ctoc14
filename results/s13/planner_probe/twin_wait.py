"""Wait-aware impulsive-twin seed for JOINT-planner tours (ctoc14.jointsearch lets a craft with no feasible leg coast
P.wait days and try again; the coast is not in the tour record except as t_flyby[k] - t_flyby[k-1] > tof[k]).
ctoc14.impulsive.from_tour spans each leg from flyby to flyby, i.e. it replaces 'coast W days + Lambert arc of tof'
by ONE Lambert arc of W + tof -- a different (often multi-rev / degenerate) transfer, so the seed is garbage and the
settle fails.  This copy (the shared module is not modified) propagates the coast first and starts the Lambert arc
(or the LIN_FALLBACK linear leg) at the end of the coast, exactly as the planner did."""
import numpy as np
import ctoc14.impulsive as IM
from ctoc14.impulsive import ImpulsiveProblem
from ctoc14.kepler import propagate_twobody
from ctoc14.constants import DAY, VINF_MAX, TMAX


def from_tour_wait(eph, tour, dtau=20 * DAY, margin=1.5, lag=1.0 * DAY):
    from ctoc14.lambert import lambert_best
    from ctoc14.linleg import LinLeg
    tL = float(tour['t_launch'])
    asts = [int(l['ast']) for l in tour['legs']]; tf = np.array([float(l['t_flyby']) for l in tour['legs']])
    tofs = np.array([float(l['tof']) for l in tour['legs']])
    rE, vE = eph.earth_state(tL)
    R, _ = eph.ast_states_at(np.array(asts) - 1, tf)
    r_prev = rE; v_prev = vE + np.asarray(tour['vinf'], float); t_prev = tL
    imp_t = []; imp_dv = []; vinf = None
    for k in range(len(asts)):
        t_dep = tf[k] - tofs[k]
        if k > 0 and t_dep > t_prev + 1.0 * DAY:                     # a planner WAIT: coast first
            r2, v2c = propagate_twobody(r_prev, v_prev, np.array([t_dep - t_prev]))
            r_prev = r2[0]; v_prev = v2c[0]; t_prev = t_dep
        v1, v2, _, _ = lambert_best(r_prev[None], R[k][None], np.array([tf[k] - t_prev]), v_prev[None], nrev_max=2)
        v1 = v1[0]; v2 = v2[0]
        if not np.all(np.isfinite(v1)):
            return None
        if k == 0:
            vinf = v1 - vE; nv = np.linalg.norm(vinf)
            if nv > VINF_MAX:
                imp_t.append(tL + lag); imp_dv.append(vinf * (1 - VINF_MAX / nv)); vinf = vinf * (VINF_MAX / nv)
        else:
            dvj = float(np.linalg.norm(v1 - v_prev))
            if IM.LIN_FALLBACK is not None and (not np.isfinite(dvj) or dvj > IM.LIN_FALLBACK):
                tof = tf[k] - t_prev
                L = LinLeg(r_prev, v_prev, np.array([tof]), dtau=10 * DAY)
                cost, varr, _, U = L.solve(0, R[k][None], TMAX / 1000.0 * 1e-3, ub=0.9, return_U=True)
                if np.isfinite(cost[0]) and cost[0] < dvj and L.check(0, U[0], R[k]) < IM.LIN_MISS_MAX:
                    for kk in range(L.K[0]):
                        if np.linalg.norm(U[0][kk]) * L.dtau > 1e-6:
                            imp_t.append(t_prev + L.taus[kk]); imp_dv.append(U[0][kk] * L.dtau)
                    r_prev = R[k]; v_prev = varr[0]; t_prev = tf[k]
                    continue
            imp_t.append(t_prev + lag); imp_dv.append(v1 - v_prev)
        r_prev = R[k]; v_prev = v2; t_prev = tf[k]
    edges = np.arange(tL, tf.max() + dtau, dtau); ts = 0.5 * (edges[:-1] + edges[1:]); Ts = np.zeros((len(ts), 3))
    if imp_t:
        ts = np.concatenate([ts, np.array(imp_t)]); Ts = np.concatenate([Ts, np.array(imp_dv).reshape(-1, 3)])
        o = np.argsort(ts); ts = ts[o]; Ts = Ts[o]
    return ImpulsiveProblem(eph, tL, vinf, ts, Ts, tf, asts, margin=margin, dtau=dtau)


def settle_tour_wait(eph, tour, settle, lags=(1.0 * DAY, 0.25 * DAY, 0.05 * DAY), tol_km=150.0, screen_tank=2000.0):
    best = np.inf
    for lag in lags:
        ip = from_tour_wait(eph, tour, lag=lag)
        if ip is None:
            return None, np.inf, None
        if screen_tank and float(ip.tank()) > screen_tank:
            return None, float('inf'), None
        miss = float(settle(ip))
        if miss <= tol_km:
            return ip, miss, lag
        best = min(best, miss)
    return None, best, None
