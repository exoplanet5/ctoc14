"""Convert a Lambert-planned tour (from search.py) into a continuous-thrust trajectory that satisfies the validator model.

Method (sequential 2-leg sliding window):
  * thrust is parametrised by samples on a uniform grid (h >= 8640 s) with the validator's sliding-window cubic Lagrange
    interpolation (one long arc per thrust segment; zero samples produce exactly zero thrust where the 4-sample window is zero);
  * for legs (i, i+1) the free variables are the samples inside those legs; constraints are the 2-D B-plane miss vectors at
    the closest approach to each target (flyby time free); objective ~ sum |T_k| via iteratively reweighted least-norm updates;
  * sensitivities by forward finite differences on a vectorised fixed-step RK4 integrator (thrust model identical to the validator);
  * converged legs are frozen; the final rows are produced with the validator's own integrator (SCTrajectory) and checked.

Failure handling of a window (all behind keyword options of convert_tour, see there):
  * adaptive_iter  - keep iterating past max_iter while the miss is still shrinking geometrically (slow-but-converging windows);
  * widen_on_fail  - un-freeze the previous leg and re-solve a 3-leg window (i-1, i, i+1): more thrust capacity before the
                     tight junction;
  * smart_replan   - before the full beam search, try to keep the remaining sequence order and only re-time it (with up to two
                     skips) from the actual state; never accept a candidate identical to the window that just failed;
  * fuel_guard     - propellant-aware guard: if the planned propellant of the remaining legs (calibrated by the actual/planned
                     ratio so far) exceeds the remaining tank, drop the most expensive planned legs first, not the last ones.
"""
import numpy as np
from .constants import MU, VE, TMAX, DAY, M_DRY, T_MISSION, D_FLYBY, VINF_MAX
from .kepler import propagate_twobody
from .thrust import thrust_interp, arc_max_thrust
from .submission import SCTrajectory

# ------------------------------------------------------------------ thrust weights on the uniform sample grid
def _window_weights(ts, t, t_sel=None):
    """For scalar t return (l, w[4]) such that T(t) = sum_q w[q] * Ts[l+q] (sliding-window cubic). The interval (and hence
    the window) is selected with t_sel (default t): integrators pass the step midpoint so that a step never straddles the
    arc on/off discontinuity or a window change."""
    n = len(ts)
    if t_sel is None:
        t_sel = t
    if t_sel < ts[0] or t_sel > ts[-1]:
        return None, None
    j = min(max(int(np.searchsorted(ts, t_sel, side='right') - 1), 0), n - 2)
    if n < 4:
        l = 0; k = n
    else:
        l = min(max(j - 1, 0), n - 4); k = 4
    tw = ts[l:l + k]
    w = np.ones(k)
    for a in range(k):
        for b in range(k):
            if a != b:
                w[a] *= (t - tw[b]) / (tw[a] - tw[b])
    return l, w

def _thrust_batch(ts, Ts, t, t_sel=None):
    """Ts: (P, n, 3) -> T(t): (P, 3)."""
    l, w = _window_weights(ts, t, t_sel)
    if l is None:
        return np.zeros((Ts.shape[0], 3))
    return np.einsum('q,pqc->pc', w, Ts[:, l:l + len(w), :])

def _rhs(t, Y, ts, Ts, t_sel=None):
    """Y: (P,7) -> dY/dt (P,7)."""
    r = Y[:, 0:3]; v = Y[:, 3:6]; m = Y[:, 6:7]
    T = _thrust_batch(ts, Ts, t, t_sel)
    Tn = np.linalg.norm(T, axis=1, keepdims=True)
    rn = np.linalg.norm(r, axis=1, keepdims=True)
    acc = -MU * r / rn ** 3 + T / m * 1e-3
    return np.concatenate([v, acc, -Tn / (VE * 1e3)], axis=1)

def integrate_rk4(Y0, t0, t_out, ts, Ts, hstep=0.1 * DAY, dense=False):
    """Vectorised RK4 from t0 to each time in t_out (sorted, > t0). Y0: (P,7). Returns list of Y at t_out
    (and, if dense, the full (times, states[steps,P,7]) history)."""
    Y = Y0.copy(); t = t0; out = []; hist_t = [t0]; hist_Y = [Y0.copy()]
    for te in t_out:
        while t < te - 1e-9:
            h = min(hstep, te - t)
            tm = t + h / 2
            k1 = _rhs(t, Y, ts, Ts, tm); k2 = _rhs(tm, Y + h / 2 * k1, ts, Ts, tm)
            k3 = _rhs(tm, Y + h / 2 * k2, ts, Ts, tm); k4 = _rhs(t + h, Y + h * k3, ts, Ts, tm)
            Y = Y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4); t += h
            if dense:
                hist_t.append(t); hist_Y.append(Y.copy())
        out.append(Y.copy())
    if dense:
        return out, (np.array(hist_t), np.array(hist_Y))
    return out

# ------------------------------------------------------------------ closest approach on the reference trajectory
def closest_approach(eph, ast_idx, hist_t, hist_Y, t_nom, search=20 * DAY, ts=None, Ts=None):
    """Find the closest-approach time to asteroid ast_idx near t_nom using the dense reference (P=1) history:
    coarse minimum on the history grid inside +-search, then a bounded 1-D minimisation of the true distance (states obtained
    by integrating from the nearest earlier history point). Returns (t_star, Y_star[7], d_vec[3], v_rel[3], dist)."""
    from scipy.optimize import minimize_scalar
    body = eph.ast[ast_idx]
    for _ in range(4):   # widen the coarse window if the minimum sits on its edge
        sel = np.where((hist_t >= t_nom - search) & (hist_t <= t_nom + search))[0]
        if len(sel) == 0:
            sel = np.array([int(np.argmin(np.abs(hist_t - t_nom)))])
        tt = hist_t[sel]; YY = hist_Y[sel, 0, :]
        ra, va = body.state(tt)
        d = np.linalg.norm(YY[:, 0:3] - ra, axis=1)
        k = int(np.argmin(d)); t_c = tt[k]
        if (0 < k < len(sel) - 1) or search >= 200 * DAY:
            break
        search *= 3
    def state_at(t):
        k0 = int(np.searchsorted(hist_t, t, side='right') - 1)
        k0 = min(max(k0, 0), len(hist_t) - 1)
        if t <= hist_t[k0] + 1e-9:
            return hist_Y[k0][0]
        return integrate_rk4(hist_Y[k0], hist_t[k0], [t], ts, Ts)[0][0]
    def dist(t):
        Y = state_at(t); r1, _ = body.state(t)
        return float(np.linalg.norm(Y[0:3] - r1))
    dtg = (tt[1] - tt[0]) if len(tt) > 1 else 0.1 * DAY
    lo, hi = t_c - 2 * dtg, t_c + 2 * dtg
    for _ in range(8):   # expand the bracket while the minimum lies on a bound (slow encounters)
        res = minimize_scalar(dist, bounds=(lo, hi), method='bounded', options=dict(xatol=1e-3, maxiter=200))
        t_star = float(res.x); w = hi - lo
        if t_star - lo < 0.02 * w:
            lo -= w
        elif hi - t_star < 0.02 * w:
            hi += w
        else:
            break
    # one Newton polish on the along-track component
    Y = state_at(t_star); r1, v1 = body.state(t_star)
    dvec = Y[0:3] - r1; vrel = Y[3:6] - v1
    dt = -np.dot(dvec, vrel) / np.dot(vrel, vrel)
    if abs(dt) < dtg:
        t_star += dt; Y = state_at(t_star); r1, v1 = body.state(t_star)
        dvec = Y[0:3] - r1; vrel = Y[3:6] - v1
    return t_star, Y, dvec, vrel, float(np.linalg.norm(dvec))

def _bplane_basis(vrel):
    vh = vrel / np.linalg.norm(vrel)
    z = np.array([0.0, 0.0, 1.0])
    e1 = np.cross(vh, z)
    if np.linalg.norm(e1) < 1e-6:
        e1 = np.cross(vh, np.array([1.0, 0.0, 0.0]))
    e1 /= np.linalg.norm(e1); e2 = np.cross(vh, e1)
    return e1, e2

def _exact_state(Y0, t0, t1, ts, Ts):
    """Integrate with the validator integrator (DOP853 per sample interval, exact coast outside samples)."""
    from .validator import integrate_segment
    r, v, m = Y0[0:3].copy(), Y0[3:6].copy(), float(Y0[6]); t = t0
    # sample boundaries between t0 and t1
    ks = np.where((ts > t0) & (ts < t1))[0]
    for k in ks:
        r, v, m = integrate_segment(r, v, m, t, ts[k], ts, Ts); t = ts[k]
    r, v, m = integrate_segment(r, v, m, t, t1, ts, Ts)
    return np.concatenate([r, v, [m]])

# ------------------------------------------------------------------ window solver (Newton / IRLS on samples + flyby times)
def _solve_window(eph, Y_start, t_start, ts, Ts, free, win, tcap, miss_tol, max_iter, fd_step, hstep, verbose, i,
                  adaptive_iter=False, max_iter_ext=36):
    """Iterate on the free samples Ts[free] and the flyby times win[q]['t_nom'] (both mutated in place) until every miss is
    below miss_tol. Returns (converged, dists, n_iter, m_min) with m_min the lowest model mass at the window's flybys.
    With adaptive_iter the loop continues past max_iter (up to max_iter_ext) as long as the largest miss still halves every
    two iterations."""
    n_s = len(ts)
    converged = False
    dists = [1e12] * len(win); m_min = float(Y_start[6])
    if len(free) == 0:
        return converged, dists, 0, m_min
    hist = []; it = 0
    while it < max_iter or (adaptive_iter and it < max_iter_ext and len(hist) >= 3 and hist[-1] < 0.5 * hist[-3]):
        # states at the current flyby-time variables t_f (no closest-approach search: t_f is a Newton variable)
        tfs = [leg['t_nom'] for leg in win]
        order = np.argsort(tfs); t_sorted = [tfs[o] for o in order]
        nf = len(free); P = 3 * nf; nw = len(win)
        Tb = np.repeat(Ts[None, :, :], P + 1, axis=0)
        for q in range(nf):
            for c in range(3):
                Tb[1 + 3 * q + c, free[q], c] += fd_step
        Y0b = np.repeat(Y_start[None, :], P + 1, axis=0)
        outs = integrate_rk4(Y0b, t_start, t_sorted, ts, Tb, hstep)
        m_min = float(min(o[0, 6] for o in outs))
        miss = np.zeros(3 * nw); G = np.zeros((3 * nw, P + nw)); dists = [0.0] * nw; rows = [None] * nw
        for oi, o in enumerate(order):
            k = win[o]['ast'] - 1; t_f = tfs[o]
            Yb = outs[oi]
            ra, va = eph.ast[k].state(t_f)
            dvecs = Yb[:, 0:3] - ra
            vrel = Yb[0, 3:6] - va
            miss[3 * o:3 * o + 3] = dvecs[0]; dists[o] = float(np.linalg.norm(dvecs[0])); rows[o] = (k, t_f, vrel)
            G[3 * o:3 * o + 3, :P] = ((dvecs[1:] - dvecs[0]) / fd_step).T
            G[3 * o:3 * o + 3, P + o] = vrel * DAY          # flyby-time shift variable in days
        if verbose:
            print(f'  leg {i:3d} it {it:2d}: miss [km] ' + ' '.join(f'{d:9.1f}' for d in dists)
                  + f'  |vrel| ' + ' '.join(f'{np.linalg.norm(rows[q][2]):5.2f}' for q in range(nw))
                  + f'  t_f ' + ' '.join(f'{rows[q][1]/DAY:8.2f}' for q in range(nw)) + f'  max|T|={np.linalg.norm(Ts[free],axis=1).max():.3f}', flush=True)
        hist.append(max(dists))
        if max(dists) <= miss_tol:
            converged = True; break
        # IRLS weights: L1 objective on sample norms; heavier for samples at the cap; time shifts are cheap
        Tn = np.linalg.norm(Ts[free], axis=1)
        w = 1.0 / np.sqrt(Tn ** 2 + 0.02 ** 2)
        w = np.where(Tn >= tcap - 1e-9, w * 50.0, w)
        Winv = np.concatenate([np.repeat(1.0 / w, 3), np.full(nw, 1.0 / 1e-4)])
        GW = G * Winv[None, :]
        A = GW @ G.T + 1e-9 * np.eye(G.shape[0])
        if not np.all(np.isfinite(miss)) or not np.all(np.isfinite(A)):
            break   # degenerate window -> treated as not converged
        try:
            lam = np.linalg.solve(A, -miss)
        except np.linalg.LinAlgError:
            lam = np.linalg.lstsq(A, -miss, rcond=None)[0]
        dx = Winv * (G.T @ lam)
        dT = dx[:P].reshape(nf, 3); dtau = dx[P:]
        # damp very large steps (linear model validity): limit thrust change per sample to 0.3 N and time shift to 60 d
        sc = min(1.0, 0.3 / max(1e-12, np.linalg.norm(dT, axis=1).max()), 60.0 / max(1e-12, np.abs(dtau).max()))
        dT *= sc; dtau *= sc
        Tnew = Ts[free] + dT
        nn = np.linalg.norm(Tnew, axis=1)
        over = nn > tcap
        Tnew[over] *= (tcap / nn[over])[:, None]
        Ts[free] = Tnew
        # interpolation-overshoot guard: the cubic interpolant may exceed the sample cap between samples;
        # scale down samples around any interval whose interpolant peaks above 0.49 N (the loop then re-converges)
        for _g in range(3):
            lo_i = max(0, free[0] - 3); hi_i = min(n_s, free[-1] + 4)
            pk, tpk = arc_max_thrust(ts[lo_i:hi_i], Ts[lo_i:hi_i], n_sub=100)
            if pk <= 0.49:
                break
            j = int(np.clip(np.searchsorted(ts, tpk) - 1, 0, n_s - 1))
            idx = [q for q in range(max(0, j - 2), min(n_s, j + 3)) if q in set(free.tolist())]
            for q in idx:
                Ts[q] *= 0.9
        for q, leg in enumerate(win):
            leg['t_nom'] = min(leg['t_nom'] + float(dtau[q]) * DAY, T_MISSION - 600.0)
        it += 1
    return converged, dists, it, m_min

# ------------------------------------------------------------------ planning-model helpers (Lambert chain, re-timing)
def _plan_propellant(m, legs):
    """Planning-model propellant of each leg (uses leg['dv_est'], leg['tof']) starting from mass m; sets leg['dm_plan']."""
    for l in legs:
        dv = float(l.get('dv_est', 0.0)); tof = float(l.get('tof', 0.0))
        a = TMAX / m * 1e-3
        kappa = 1 + dv / (2 * a * tof) if (dv > 0 and tof > 0) else 1.0
        dm = m * (1 - np.exp(-kappa * dv / VE))
        l['dm_plan'] = float(dm); m -= dm
    return m

def _retime_sequence(eph, t, r, v, m, remaining, m0, tL, vinf, P, max_skip=0, beam=60, per_pos=6):
    """Order-preserving re-planning with the Lambert-junction model of search.py: re-time the targets of `remaining` in their
    original order from the actual state, allowing up to max_skip consecutive targets to be skipped at each step. Beam search
    over the planner's TOF grid, stratified by the position reached in `remaining` (top `per_pos` states per position by
    (targets visited, score), `beam` positions at most) so that paths skipping an expensive leg early survive to the end.
    Returns the State with the most targets (ties: best score) or None."""
    from .search import State, expand, mask_of
    if not remaining:
        return None
    # keep every feasible TOF per target: from an off-plan state the arrival velocity (which depends on the TOF) decides
    # whether the NEXT leg is feasible, and the planner's default of 2 TOFs per target discards those options
    try:
        import dataclasses
        P = dataclasses.replace(P, n_per_target=len(P.tofs), max_children=max(P.max_children, 4 * len(P.tofs)))
    except Exception:
        pass
    idx = {a: k for k, a in enumerate(remaining)}
    all_mask = (1 << 300) - 1
    s0 = State(t=float(t), r=np.asarray(r, float), v=np.asarray(v, float), m=float(m), visited=0, seq=(),
               fuel=float(m0 - m), t_launch=float(tL), vinf=np.asarray(vinf, float), m0=float(m0))
    states = [s0]; best = None
    for _depth in range(len(remaining)):
        children = []
        for s in states:
            p = idx[s.seq[-1][0]] + 1 if s.seq else 0
            allowed = remaining[p:p + 1 + max_skip]
            if not allowed:
                continue
            s.visited = all_mask & ~mask_of(allowed)          # expand() only targets the allowed next asteroids
            children += expand(eph, s, P)
        if not children:
            break
        groups = {}
        for c in children:
            c.visited = mask_of([x[0] for x in c.seq])
            groups.setdefault(idx[c.seq[-1][0]], []).append(c)
        nb = []
        for p in sorted(groups)[:beam]:
            g = sorted(groups[p], key=lambda s: (-s.n(), s.score(P)))
            seen = set(); kept = []
            for c in g:
                key = (c.visited, int(c.t // (5 * DAY)))
                if key in seen:
                    continue
                seen.add(key); kept.append(c)
                if len(kept) >= per_pos:
                    break
            nb += kept
        states = nb
        top = max(states, key=lambda s: (s.n(), -s.score(P)))
        if best is None or top.n() > best.n() or (top.n() == best.n() and top.score(P) < best.score(P)):
            best = top
    return best

def _identical_to_failed(seq, failed, tol=10 * DAY):
    """seq: list of (ast, t). failed: list of (asts, t_noms_initial, t_noms_final) of window attempts that failed."""
    for asts, tn0, tn1 in failed:
        k = min(len(asts), len(seq))
        if k == 0 or [s[0] for s in seq[:k]] != list(asts[:k]):
            continue
        if all(min(abs(seq[q][1] - tn0[q]), abs(seq[q][1] - tn1[q])) < tol for q in range(k)):
            return True
    return False

def _smart_replan(eph, t, Y, remaining, visited, m0, tL, vinf, RP, failed, verbose):
    """Candidates: (1) original order re-timed without skips (accepted at once if complete and not identical to a failed
    window), else (2) original order re-timed with <= 2 skips per step and (3) the full beam search from the actual state.
    Returns (new_legs or None, kind)."""
    from .search import beam_search_from_state
    def legs_of(s):
        return [dict(ast=int(x[0]), t_nom=float(x[1]), dv_est=float(x[2]), tof=float(x[3])) for x in s.seq]
    def ok(s):
        return s is not None and s.n() > 0 and not _identical_to_failed([(x[0], x[1]) for x in s.seq], failed)
    r, v, m = Y[0:3], Y[3:6], float(Y[6])
    # the order-preserving candidates use the planner's own feasibility limit (eta >= 0.6): the remaining legs were planned
    # with it, and a candidate identical to the window that just failed is rejected anyway
    try:
        import dataclasses
        RPr = dataclasses.replace(RP, eta=max(RP.eta, 0.6))
    except Exception:
        RPr = RP
    try:
        s1 = _retime_sequence(eph, t, r, v, m, remaining, m0, tL, vinf, RPr, max_skip=0)
    except Exception:
        s1 = None
    if s1 is not None and s1.n() == len(remaining) and ok(s1):
        return legs_of(s1), 'retime'
    cands = []
    try:
        s2 = _retime_sequence(eph, t, r, v, m, remaining, m0, tL, vinf, RPr, max_skip=2)
        if ok(s2):
            cands.append(('retime_skip', s2))
    except Exception:
        pass
    best, _ = beam_search_from_state(eph, t, r, v, m, visited, m0, tL, vinf, P=RP, n_proc=1, verbose=False)
    if ok(best):
        cands.append(('beam', best))
    if not cands:
        return None, 'none'
    cands.sort(key=lambda ks: (-ks[1].n(), ks[1].fuel))     # most targets; ties keep the order-preserving candidate first
    kind, s = cands[0]
    return legs_of(s), kind

# ------------------------------------------------------------------ main conversion
class ConversionError(Exception):
    pass

def convert_tour(eph, tour, h=1.0 * DAY, tcap=0.45, miss_tol=300.0, max_iter=12, fd_step=0.02, verbose=True,
                 hstep=0.1 * DAY, window_legs=2, drop_on_fail=True, replan=True, global_excluded=(), allowed_extra=(),
                 replan_params=None, max_replans=30,
                 adaptive_iter=True, max_iter_ext=36, widen_on_fail=True, max_widen=8, smart_replan=False,
                 fuel_guard=False, fuel_guard_factor='auto', fuel_guard_max_avail=None, expensive_drop=True):
    """tour: dict from search.tour_json. Returns (SCTrajectory, info dict).

    Failure-handling options (see module docstring): adaptive_iter / max_iter_ext, widen_on_fail / max_widen (total 3-leg
    retries per tour), smart_replan, fuel_guard / fuel_guard_factor ('auto' = actual/planned propellant ratio so far,
    clipped to [1, 1.3] once >= 300 kg were planned; or a number; fuel_guard_max_avail restricts the guard to the terminal
    phase, remaining tank <= that many kg). Setting all four flags False reproduces the original
    2-leg window + beam re-plan + drop behaviour exactly."""
    tL = float(tour['t_launch']); m0 = float(tour['m0'])
    vinf = np.array(tour['vinf'], float)
    if np.linalg.norm(vinf) > VINF_MAX - 1e-6:
        vinf = vinf * ((VINF_MAX - 1e-6) / np.linalg.norm(vinf))
    rE, vE = eph.earth_state(tL)
    legs = [dict(ast=int(l['ast']), t_nom=float(l['t_flyby']), dv_est=float(l.get('dv_est', 0.0)), tof=float(l.get('tof', 0.0)))
            for l in tour['legs']]
    _plan_propellant(m0, legs)
    t_last = legs[-1]['t_nom'] + 5 * DAY
    n_s = int(np.floor((min(t_last, T_MISSION - DAY) - tL) / h))
    ts = tL + h * (np.arange(n_s) + 1)              # first sample one step after launch
    Ts = np.zeros((n_s, 3))
    # state at the start of the current window (after the last frozen flyby)
    Y_start = np.concatenate([rE, vE + vinf, [m0]]); t_start = tL
    frozen_upto = 0                                  # samples with index < frozen_upto are frozen
    flybys = []                                      # (t_star, ast_id, dist)
    info = dict(legs=[], dropped=[], events=[], n_replans=0, n_widen=0, n_widen_ok=0, n_guard_drops=0)
    i = 0
    fuel_reserve = 3.0   # kg kept above dry mass
    last_snap = None     # state at the start of the last frozen leg's window (for widen_on_fail)
    fail_state = None    # state at the failure that triggered a widened retry (restored if the retry fails too)
    pending_nwin = None  # window length of the next iteration (window_legs + 1 after a backtrack)
    widen_at = {}
    failed_windows = {}  # leg index -> [(asts, t_noms at start, t_noms at end)] of failed window attempts
    plan_used = 0.0      # planning-model propellant of the frozen legs
    guard_at = set()     # leg indices where the fuel guard already acted
    while i < len(legs):
        if fuel_guard and i < len(legs) - 1 and i not in guard_at and pending_nwin is None:
            # (c) propellant-aware guard on the planning model: the planned propellant of the remaining legs (the plan's own
            # per-leg numbers - a Lambert chain re-evaluated from the actual state is hypersensitive to small flyby-time
            # shifts and useless here), scaled by the actual/planned ratio of the legs flown so far. If it does not fit the
            # remaining tank, drop the fewest legs: the most expensive planned legs first, without overshooting the deficit
            # (never the imminent leg i). The gaps are left to the window solver and the failure fallbacks.
            if fuel_guard_factor == 'auto':
                pu = sum(l['dm_plan'] for l in info['legs'][1:]); au = sum(l['fuel_leg'] for l in info['legs'][1:])
                factor = 1.0 if pu < 300.0 else float(np.clip(au / pu, 1.0, 1.3))
            else:
                factor = float(fuel_guard_factor)
            avail = Y_start[6] - M_DRY - fuel_reserve
            need = factor * sum(l['dm_plan'] for l in legs[i:])
            if need > avail + 10.0 and (fuel_guard_max_avail is None or avail <= fuel_guard_max_avail):
                guard_at.add(i)
                deficit = need - avail
                cand = list(legs[i + 1:]); drop = []
                while deficit > 0 and cand:
                    over = [l for l in cand if factor * l['dm_plan'] >= deficit]
                    pick = min(over, key=lambda l: l['dm_plan']) if over else max(cand, key=lambda l: l['dm_plan'])
                    drop.append(pick); cand.remove(pick); deficit -= factor * pick['dm_plan']
                lost = [l['ast'] for l in drop]
                legs[i + 1:] = [l for l in legs[i + 1:] if l['ast'] not in set(lost)]
                info['dropped'] += lost; info['n_guard_drops'] += len(lost)
                info['events'].append(dict(kind='fuel_guard', i=i, lost=lost, factor=factor, need=need, avail=avail,
                                           dm_lost=[l['dm_plan'] for l in drop]))
                if verbose: print(f'  !! fuel guard at leg {i}: planned {need:.0f} kg (x{factor:.2f}) > avail {avail:.0f}: '
                                  f'dropped {lost} ({[round(l["dm_plan"]) for l in drop]} kg)', flush=True)
        if Y_start[6] < M_DRY + fuel_reserve + 1.0:
            info['dropped'] += [l['ast'] for l in legs[i:]]
            if verbose: print(f'  !! out of propellant at leg {i}: ending tour', flush=True)
            break
        nwin = pending_nwin or window_legs; pending_nwin = None
        win = legs[i:i + nwin]
        t_nom0 = [l['t_nom'] for l in win]
        t_end = win[-1]['t_nom'] + 40 * DAY
        # free samples: index >= frozen_upto and time <= t_end (+ a little)
        free = np.where((np.arange(n_s) >= max(frozen_upto, 1)) & (ts <= min(t_end, ts[-1])) & (np.arange(n_s) < n_s - 1))[0]
        import os
        if os.environ.get('DUMP_LEG') == str(i):
            np.savez(os.environ.get('DUMP_PATH', 'results/dump_leg.npz'), Y_start=Y_start, t_start=t_start, ts=ts, Ts=Ts, free=free,
                     asts=np.array([l['ast'] for l in win]), t_noms=np.array([l['t_nom'] for l in win]), t_end=t_end)
            print(f'  [dumped window {i}]', flush=True)
        converged, dists, n_it, m_min = _solve_window(eph, Y_start, t_start, ts, Ts, free, win, tcap, miss_tol, max_iter, fd_step,
                                                      hstep, verbose, i, adaptive_iter=adaptive_iter, max_iter_ext=max_iter_ext)
        if converged and nwin > window_legs and m_min < M_DRY + fuel_reserve:
            # a widened window that closes its flybys only by draining the tank below dry mass is a failure, not a solution
            converged = False
            if verbose: print(f'  !! 3-leg window converged but exhausts propellant (m={m_min:.1f}): treated as failed', flush=True)
        if not converged:
            Ts[free] = 0.0
            if nwin > window_legs and fail_state is not None:
                # the widened (backtracked) window failed too: return to the state at the original failure
                fs = fail_state; fail_state = None
                Y_start = fs['Y_start']; t_start = fs['t_start']; frozen_upto = fs['frozen_upto']; Ts[:] = fs['Ts']
                flybys[:] = fs['flybys']; info['legs'][:] = fs['legs_info']; plan_used = fs['plan_used']
                for l, tn in zip(legs, fs['t_noms']): l['t_nom'] = tn
                i = fs['i']; win = legs[i:i + window_legs]; dists = fs['dists']; t_nom0 = fs['t_nom0']
                info['events'].append(dict(kind='widen_fail', i=i, asts=[l['ast'] for l in win], n_iter=n_it))
                if verbose: print(f'  !! 3-leg retry failed; back to leg {i}', flush=True)
            elif (widen_on_fail and nwin == window_legs and last_snap is not None and i >= 1 and n_it > 0
                  and info['n_widen'] < max_widen and widen_at.get((i, tuple(l['ast'] for l in win)), 0) < 1):
                # (a) un-freeze the previous leg and re-solve legs (i-1, i, i+1) together
                widen_at[(i, tuple(l['ast'] for l in win))] = 1; info['n_widen'] += 1
                fail_state = dict(Y_start=Y_start.copy(), t_start=t_start, frozen_upto=frozen_upto, Ts=Ts.copy(),
                                  flybys=list(flybys), legs_info=list(info['legs']), t_noms=[l['t_nom'] for l in legs],
                                  i=i, dists=list(dists), plan_used=plan_used, t_nom0=list(t_nom0))
                for l, tn in zip(win, t_nom0): l['t_nom'] = tn
                Y_start = last_snap['Y_start'].copy(); t_start = last_snap['t_start']; frozen_upto = last_snap['frozen_upto']
                Ts[:] = last_snap['Ts']; plan_used = last_snap['plan_used']
                flybys.pop(); info['legs'].pop()
                i -= 1; pending_nwin = window_legs + 1
                info['events'].append(dict(kind='widen', i=i, asts=[l['ast'] for l in legs[i:i + pending_nwin]]))
                if verbose: print(f'  !! window failed (miss {max(dists):.0f} km): retry legs {i}..{i + pending_nwin - 1} as a '
                                  f'{pending_nwin}-leg window', flush=True)
                continue
            failed_windows.setdefault(i, []).append(([l['ast'] for l in win], list(t_nom0), [l['t_nom'] for l in win]))
            if replan and info['n_replans'] < max_replans:
                # re-plan the rest of the tour from the actual state with the Lambert-junction search
                from .search import beam_search_from_state, Params
                info['n_replans'] += 1
                flown = {f[1] for f in flybys}
                remaining = [l['ast'] for l in legs[i:]]
                candidates = set(remaining) | set(allowed_extra)
                visited = (set(range(1, 301)) - candidates) | flown | set(global_excluded)
                RP = replan_params or Params(beam=150, eta=0.5, w_fuel=0.7, tofs=np.arange(15, 401, 5) * DAY)
                new_legs = None; kind = 'beam'
                ra_ = info.setdefault('replans_at', {})
                key = remaining[0] if remaining else i     # count re-plans per TARGET, not per leg index (avoids cascades after drops)
                if smart_replan:
                    ra_[key] = ra_.get(key, 0) + 1
                    new_legs, kind = _smart_replan(eph, t_start, Y_start, remaining, visited, m0, tL, vinf, RP,
                                                   failed_windows.get(i, []), verbose)
                else:
                    best, _ = beam_search_from_state(eph, t_start, Y_start[0:3], Y_start[3:6], Y_start[6], visited, m0, tL, vinf, P=RP, n_proc=1, verbose=False)
                    same = best is not None and best.n() > 0 and [x[0] for x in best.seq][:2] == remaining[:2] \
                        and ra_.get(key, 0) >= 1
                    ra_[key] = ra_.get(key, 0) + 1
                    if best is not None and best.n() > 0 and not same:
                        new_legs = [dict(ast=int(x[0]), t_nom=float(x[1]), dv_est=float(x[2]), tof=float(x[3])) for x in best.seq]
                if new_legs is not None and ra_[key] <= 3:
                    _plan_propellant(float(Y_start[6]), new_legs)
                    lost = [a for a in remaining if a not in {l['ast'] for l in new_legs}]
                    info['dropped'] += lost
                    info['events'].append(dict(kind='replan', how=kind, i=i, n_new=len(new_legs), n_old=len(remaining), lost=lost))
                    if verbose: print(f'  !! replanned from leg {i} ({kind}): {len(new_legs)} legs (was {len(remaining)}), lost {lost}', flush=True)
                    legs[i:] = new_legs
                    continue
            if drop_on_fail and len(win) >= 1:
                worst = int(np.argmax(dists))
                dropped = win[worst]['ast']
                info['dropped'].append(dropped)
                info['events'].append(dict(kind='drop', i=i, ast=dropped, miss=float(dists[worst])))
                if verbose: print(f'  !! dropping asteroid {dropped} (miss {dists[worst]:.0f} km)')
                legs.pop(i + worst)
                continue
            raise ConversionError(f'leg {i} did not converge')
        # freeze leg i: polish its flyby time with a closest-approach search on the dense reference, then freeze
        if nwin > window_legs:
            info['n_widen_ok'] += 1
            info['events'].append(dict(kind='widen_ok', i=i, asts=[l['ast'] for l in win], n_iter=n_it))
            if verbose: print(f'  !! 3-leg window converged at leg {i}', flush=True)
        fail_state = None
        last_snap = dict(Y_start=Y_start.copy(), t_start=t_start, frozen_upto=frozen_upto, Ts=Ts.copy(), plan_used=plan_used)
        t_star_i = win[0]['t_nom']
        try:
            (Yend,), (ht, hY) = integrate_rk4(Y_start[None, :], t_start, [t_star_i + 5 * DAY], ts, Ts[None, :, :], hstep, dense=True)
            t_p, Y_p, dvec_p, vrel_p, dist_p = closest_approach(eph, win[0]['ast'] - 1, ht, hY, t_star_i, search=3 * DAY, ts=ts, Ts=Ts[None])
            if dist_p < dists[0] and abs(t_p - t_star_i) < 3 * DAY:
                t_star_i = t_p; win[0]['t_nom'] = t_p
        except Exception:
            pass
        # zero out negligible samples (< 1e-4 N) in the frozen part for clean arcs
        fz = free[ts[free] <= t_star_i]
        small = np.linalg.norm(Ts[fz], axis=1) < 1e-4
        Ts[fz[small]] = 0.0
        # recompute the state at t_star_i with frozen samples using the EXACT (validator) integrator
        Y_i = _exact_state(Y_start, t_start, t_star_i, ts, Ts)[None, :]
        # verify the flyby distance with the cleaned samples
        ra, va = eph.ast[win[0]['ast'] - 1].state(t_star_i)
        dist_i = float(np.linalg.norm(Y_i[0, 0:3] - ra))
        if dist_i > D_FLYBY * 0.8:
            # zeroing small samples moved it: undo zeroing and recompute
            Ts[fz[small]] = Ts[fz[small]]  # (kept zero) -> fallback: re-run window without cleaning
            info.setdefault('warn', []).append(f'leg {i}: cleaned miss {dist_i:.0f} km')
        fuel_leg_i = float(Y_start[6] - Y_i[0, 6]); dm_plan_i = float(win[0].get('dm_plan', 0.0))
        if expensive_drop and fuel_leg_i > 80.0 and fuel_leg_i > max(3.0 * dm_plan_i, dm_plan_i + 60.0):
            # the converted leg costs far more than planned: a miss costs 1 unit, this propellant would buy several flybys
            info['dropped'].append(win[0]['ast'])
            info['events'].append(dict(kind='drop_expensive', i=i, ast=win[0]['ast'], fuel=fuel_leg_i, planned=dm_plan_i))
            if verbose: print(f'  !! dropping asteroid {win[0]["ast"]}: leg would cost {fuel_leg_i:.0f} kg vs {dm_plan_i:.0f} planned', flush=True)
            Ts[free] = 0.0
            legs.pop(i)
            continue
        if Y_i[0, 6] < M_DRY + fuel_reserve:
            info['dropped'] += [l['ast'] for l in legs[i:]]
            if verbose: print(f'  !! leg {i} would exhaust propellant (m={Y_i[0,6]:.1f}); ending tour before it', flush=True)
            Ts[free] = 0.0
            break
        flybys.append((t_star_i, win[0]['ast'], dist_i))
        info['legs'].append(dict(ast=win[0]['ast'], t_flyby=t_star_i, dist=dist_i, m=float(Y_i[0, 6]),
                                 fuel_leg=float(Y_start[6] - Y_i[0, 6]), dm_plan=float(win[0].get('dm_plan', 0.0))))
        plan_used += float(win[0].get('dm_plan', 0.0))
        if verbose:
            print(f'leg {i:3d} frozen: ast {win[0]["ast"]:3d} t={t_star_i/DAY:8.1f} d dist={dist_i:6.1f} km m={Y_i[0,6]:.1f} kg', flush=True)
        Y_start = Y_i[0]; t_start = t_star_i
        # freeze two extra samples past the flyby: the 4-sample Lagrange window couples them to the thrust before t_star_i
        frozen_upto = min(n_s, int(np.searchsorted(ts, t_star_i, side='right')) + 2)
        i += 1
    # ---- build rows with the exact integrator
    traj = build_trajectory(tL, rE, vE + vinf, m0, ts, Ts, flybys, verbose=verbose)
    info['n_flybys'] = len(flybys); info['fuel'] = m0 - traj.m
    return traj, info

def build_trajectory(tL, r0, v0, m0, ts, Ts, flybys, verbose=False):
    """Create SCTrajectory rows: ONE thrust arc containing every sample (identical to the optimisation model), flyby rows
    inserted inside the arc (nudged by 1 s if they coincide with a sample time), end after the last flyby."""
    fb = sorted(flybys)
    t_end = fb[-1][0] if fb else ts[-1]
    # samples after the last flyby are unnecessary: keep up to the first sample after t_end (thrust there is ~0)
    n_keep = int(np.searchsorted(ts, t_end, side='right')) + 1
    n_keep = min(max(n_keep, 4), len(ts))
    arc_ts = ts[:n_keep]; arc_Ts = Ts[:n_keep].copy()
    inside = []
    for (t, a, d) in fb:
        if np.min(np.abs(arc_ts - t)) < 1.0:
            t = t + 1.0
        inside.append((t, a))
    traj = SCTrajectory(tL, r0, v0, m0)
    traj.thrust_arc(arc_ts, arc_Ts, flybys=[f for f in inside if f[0] <= arc_ts[-1]])
    for f in inside:
        if f[0] > arc_ts[-1]:
            traj.coast_to(f[0]); traj.flyby(f[1])
    dt_end = min(60.0, max(1.0, T_MISSION - 1.0 - traj.t))   # end row strictly later than the last row, inside the window
    traj.coast(dt_end)
    traj.end()
    return traj
