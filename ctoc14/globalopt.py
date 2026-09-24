"""Whole-trajectory minimum-propellant optimisation of a converted craft with a FIXED flyby sequence (sequential convex
programming on the validator's thrust model).

The converter (lowthrust.convert_tour) closes flybys with a feasibility Newton on 2-leg windows and stops at 300 km, so
propellant is never minimised over the trajectory and the lead time of a correction is at most two legs. Here every
thrust sample of the arc, every flyby time and the launch v_inf direction are free at once:

  min  sum_s |T_s|                        (propellant ~ h / ve * sum_s |T_s|)
  s.t. r(t_j) = r_ast_j(t_j)  for every flyby j (linearised: STM chain along the current trajectory, flyby-time shift
       along the relative velocity), |T_s| <= tcap, |v_inf| <= 4 km/s

solved by IRLS (weighted least norm, closed form, rows = 3 x flybys) with a proximal term; each iterate is re-integrated
with the exact thrust model (RK4 on the sliding-window cubic) and re-linearised."""
import numpy as np
from .constants import MU, VE, TMAX, DAY, M_DRY, T_MISSION, VINF_MAX
from .kepler import propagate_batch
from .lowthrust import integrate_rk4


def load_craft(path, sc=None):
    """Rows of one craft -> dict(tL, r0, v0, m0, ts, Ts, flybys[(t, ast)]). One continuous thrust arc expected."""
    from .validator import parse
    rows = parse(path)
    if sc is None:
        sc = rows[0].sc
    rr = [r for r in rows if r.sc == sc]
    L = rr[0]
    ts = np.array([r.t for r in rr if r.event == 1]); Ts = np.array([r.T for r in rr if r.event == 1])
    fb = [(r.t, r.ast) for r in rr if r.event == 3]
    return dict(tL=L.t, r0=L.r.copy(), v0=L.v.copy(), m0=L.m, ts=ts, Ts=Ts, flybys=fb, m_end=rr[-1].m)


def stm_chain(times, Y):
    """Kepler STMs between consecutive states (central differences, one batched propagation) and their cumulative
    products Phi_k = Phi(t_k, t_0). times (K,), Y (K,7). Returns (K,6,6)."""
    K = len(times); dt = np.diff(times)
    er, ev = 1.0, 1e-6
    base_r = np.repeat(Y[:-1, 0:3], 12, axis=0); base_v = np.repeat(Y[:-1, 3:6], 12, axis=0)
    pert = np.zeros((12, 6))
    for i in range(6):
        pert[2 * i, i] = er if i < 3 else ev; pert[2 * i + 1, i] = -(er if i < 3 else ev)
    P = np.tile(pert, (K - 1, 1))
    rp, vp = propagate_batch(base_r + P[:, 0:3], base_v + P[:, 3:6], np.repeat(dt, 12))
    X = np.concatenate([rp, vp], axis=1).reshape(K - 1, 12, 6)
    step = np.array([er, er, er, ev, ev, ev])
    S = ((X[:, 0::2, :] - X[:, 1::2, :]) / (2 * step)[None, :, None]).transpose(0, 2, 1)   # (K-1, 6, 6): d x_{k+1} / d x_k
    Phi = np.zeros((K, 6, 6)); Phi[0] = np.eye(6)
    for k in range(K - 1):
        Phi[k + 1] = S[k] @ Phi[k]
    return Phi

def _lagrange_w(tw, t):
    """tw (N,4) window nodes, t (N,) -> (N,4) cubic Lagrange weights."""
    L = np.ones(tw.shape)
    for a in range(4):
        for b in range(4):
            if a != b:
                L[:, a] *= (t - tw[:, b]) / (tw[:, a] - tw[:, b])
    return L


def rk4_rows(Y, t0, dt, tw, W, hstep, massless=False):
    """Batched fixed-step RK4 of the validator dynamics, one row per segment with its own start time, duration, window
    nodes tw (N,4) and window samples W (N,4,3) (zeros = coast). Equal sub-steps ceil(dt / hstep)."""
    Y = Y.copy(); nsub = np.maximum(1, np.ceil(dt / hstep - 1e-9)).astype(int)
    for ns_val in np.unique(nsub):
        sel = nsub == ns_val
        y = Y[sel]; h = (dt[sel] / ns_val)[:, None]; tws = tw[sel]; Ws = W[sel]; tt = t0[sel]
        def f(t, y):
            T = np.einsum('nq,nqc->nc', _lagrange_w(tws, t), Ws)
            r = y[:, 0:3]; rn = np.linalg.norm(r, axis=1, keepdims=True)
            acc = -MU * r / rn ** 3 + T / y[:, 6:7] * 1e-3
            mdot = np.zeros((len(y), 1)) if massless else -np.linalg.norm(T, axis=1, keepdims=True) / (VE * 1e3)
            return np.concatenate([y[:, 3:6], acc, mdot], axis=1)
        for i in range(ns_val):
            t = tt + i * h[:, 0]
            k1 = f(t, y); k2 = f(t + h[:, 0] / 2, y + h / 2 * k1); k3 = f(t + h[:, 0] / 2, y + h / 2 * k2); k4 = f(t + h[:, 0], y + h * k3)
            y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        Y[sel] = y
    return Y


def integrate_rk4_ml(Y0, t0, t_out, ts, Ts, hstep=0.1 * DAY, massless=False):
    """lowthrust.integrate_rk4 with an optional mass-free mode (mass derivative zero: the samples act as accelerations
    Ts / m0 with m0 the initial mass)."""
    if not massless:
        return integrate_rk4(Y0, t0, t_out, ts, Ts, hstep)
    from .lowthrust import _rhs
    def rhs(t, Y, tm):
        dY = _rhs(t, Y, ts, Ts, tm); dY[:, 6] = 0.0; return dY
    Y = Y0.copy(); t = t0; out = []
    for te in t_out:
        while t < te - 1e-9:
            h = min(hstep, te - t); tm = t + h / 2
            k1 = rhs(t, Y, tm); k2 = rhs(tm, Y + h / 2 * k1, tm); k3 = rhs(tm, Y + h / 2 * k2, tm); k4 = rhs(t + h, Y + h * k3, tm)
            Y = Y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4); t += h
        out.append(Y.copy())
    return out


class GlobalProblem:
    def __init__(self, eph, craft, tcap=0.45, hstep=0.1 * DAY):
        self.eph = eph; self.c = craft; self.tcap = tcap; self.hstep = hstep
        self.tL = craft['tL']; self.m0 = craft['m0']
        self.ts = craft['ts'].copy(); self.Ts = craft['Ts'].copy()
        self.asts = [a for _, a in craft['flybys']]; self.tf = np.array([t for t, _ in craft['flybys']])
        rE, vE = eph.earth_state(self.tL); self.rE = rE; self.vE = vE
        self.vinf = craft['v0'] - vE
        self.h = np.median(np.diff(self.ts))
        self.massless = False

    def integrate(self, Ts=None, tf=None, vinf=None, dense_samples=False):
        """Exact-model RK4 states at the flyby times (and at every sample time if dense_samples)."""
        Ts = self.Ts if Ts is None else Ts; tf = self.tf if tf is None else tf; vinf = self.vinf if vinf is None else vinf
        Y0 = np.concatenate([self.rE, self.vE + vinf, [self.m0]])[None, :]
        if dense_samples:
            tout = np.concatenate([self.ts, tf]); order = np.argsort(tout, kind='stable'); tout_s = tout[order]
        else:
            order = np.argsort(tf, kind='stable'); tout_s = tf[order]
        outs = integrate_rk4_ml(Y0, self.tL, list(tout_s), self.ts, Ts[None], self.hstep, massless=self.massless)
        Yo = np.array([o[0] for o in outs]); Y = np.empty_like(Yo); Y[order] = Yo
        if dense_samples:
            return Y[len(self.ts):], Y[:len(self.ts)]
        return Y, None

    def misses(self, Yf, tf=None):
        """Relative position / velocity at the flyby times; an aim-point offset (homotopy of an inserted flyby,
        self.offsets[ast]) is subtracted from the relative position."""
        tf = self.tf if tf is None else tf
        idx = np.array(self.asts) - 1
        ra, va = self.eph.ast_states_at(idx, tf)
        dm = Yf[:, 0:3] - ra
        offs = getattr(self, 'offsets', None)
        if offs:
            for j, a in enumerate(self.asts):
                if a in offs:
                    dm[j] = dm[j] - offs[a]
        return dm, Yf[:, 3:6] - va

    def chain(self, Yf, Ys, Ts=None, vinf=None):
        """Exact discrete sensitivities of the RK4 model (mass coupling and interpolation spread included).
        Segments = consecutive nodes of {tL} u ts u tf; per segment a batched central-difference RK4 on the 7 state
        components and the 12 window-sample components; cumulative 7x7 products Pi(node) = d x(node) / d x(launch)."""
        Ts = self.Ts if Ts is None else Ts; vinf = self.vinf if vinf is None else vinf
        ns = len(self.ts); nf = len(self.tf); ts = self.ts
        Y0 = np.concatenate([self.rE, self.vE + vinf, [self.m0]])
        tn = np.concatenate([[self.tL], ts, self.tf]); Yn = np.concatenate([Y0[None], Ys, Yf], axis=0)
        kind = np.concatenate([[-1], np.arange(ns), -2 - np.arange(nf)])        # -1 launch, s >= 0 sample, -2-j flyby j
        order = np.argsort(tn, kind='stable'); tn = tn[order]; Yn = Yn[order]; kind = kind[order]
        G = len(tn) - 1
        t0 = tn[:-1]; dt = tn[1:] - tn[:-1]; tm = t0 + dt / 2
        inside = (tm > ts[0]) & (tm < ts[-1])
        j_int = np.clip(np.searchsorted(ts, tm, side='right') - 1, 0, ns - 2)
        l = np.clip(j_int - 1, 0, ns - 4)
        idx = l[:, None] + np.arange(4)[None, :]
        tw = ts[idx]; W = Ts[idx] * inside[:, None, None]
        er, ev, em, eT = 1.0, 1e-6, 1e-3, 1e-4
        steps = [er, er, er, ev, ev, ev, em]
        nP = 2 * 7 + 2 * 12
        Yb = np.repeat(Yn[:-1], nP, axis=0); Wb = np.repeat(W, nP, axis=0)
        base = np.arange(G) * nP
        for i in range(7):
            Yb[base + 2 * i, i] += steps[i]; Yb[base + 2 * i + 1, i] -= steps[i]
        for q in range(4):
            for c in range(3):
                o = 14 + 2 * (3 * q + c)
                Wb[base + o, q, c] += eT; Wb[base + o + 1, q, c] -= eT
        Wb[np.repeat(~inside, nP)] = 0.0
        Out = rk4_rows(Yb, np.repeat(t0, nP), np.repeat(dt, nP), np.repeat(tw, nP, axis=0), Wb, self.hstep,
                       massless=self.massless).reshape(G, nP, 7)
        S = np.zeros((G, 7, 7)); C = np.zeros((G, 7, 4, 3))
        for i in range(7):
            S[:, :, i] = (Out[:, 2 * i] - Out[:, 2 * i + 1]) / (2 * steps[i])
        for q in range(4):
            for c in range(3):
                o = 14 + 2 * (3 * q + c)
                C[:, :, q, c] = (Out[:, o] - Out[:, o + 1]) / (2 * eT) * inside[:, None]
        Pi = np.zeros((G + 1, 7, 7)); Pi[0] = np.eye(7)
        for g in range(G):
            Pi[g + 1] = S[g] @ Pi[g]
        return dict(Pi=Pi, Pinv=np.linalg.inv(Pi), C=C, idx=idx, tn=tn, Yn=Yn, kind=kind, ns=ns)

    def rows_at(self, ch, times, nodes=None):
        """Sensitivity rows of r(t) for arbitrary times (node = last chain node <= t, plus a straight-line shift inside the
        partial segment): B (n, ns, 3, 3), Lv (n, 3, 3)."""
        tn = ch['tn']; Pi = ch['Pi']; Pinv = ch['Pinv']; C = ch['C']; idx = ch['idx']; ns = ch['ns']
        n = len(times); B = np.zeros((n, ns, 3, 3)); Lv = np.zeros((n, 3, 3))
        for j, t in enumerate(times):
            nj = int(np.searchsorted(tn, t + 1e-6, side='right') - 1) if nodes is None else nodes[j]
            P = np.zeros((3, 7)); P[:, 0:3] = np.eye(3); P[:, 3:6] = np.eye(3) * (t - tn[nj])
            Pr = P @ Pi[nj]
            Lv[j] = Pr[:, 3:6]
            if nj == 0:
                continue
            M = np.einsum('ab,gbc->gac', Pr, Pinv[1:nj + 1])
            contrib = np.einsum('gab,gbqc->gqac', M, C[:nj])
            np.add.at(B[j], idx[:nj].reshape(-1), contrib.reshape(-1, 3, 3))
        return B, Lv

    def linearise(self, Yf, Ys, Ts=None, vinf=None):
        """B (nf, ns, 3, 3) d r(t_j) / d T_s [km/N], Lv (nf, 3, 3) d r(t_j) / d v_inf."""
        ch = self.chain(Yf, Ys, Ts, vinf)
        kind = ch['kind']; nf = len(self.tf)
        node_of_f = {int(-2 - kind[n]): n for n in range(len(kind)) if kind[n] <= -2}
        return self.rows_at(ch, self.tf, nodes=[node_of_f[j] for j in range(nf)])

    def fuel(self, Ts=None):
        """Propellant [kg]: thrust mode h/ve sum|T_s|; mass-free mode m0 (1 - exp(-dv/ve)) with dv = h sum|T_s| / m0."""
        Ts = self.Ts if Ts is None else Ts
        if self.massless:
            return self.m0 * (1.0 - np.exp(-self.dv(Ts) / VE))
        return self.h / (VE * 1e3) * np.linalg.norm(Ts, axis=1).sum()

    def dv(self, Ts=None):
        """Delta-v [km/s] of the mass-free acceleration samples Ts / m0."""
        Ts = self.Ts if Ts is None else Ts
        return self.h * np.linalg.norm(Ts, axis=1).sum() / self.m0 * 1e-3


def _bplane(vrel):
    """(n, 2, 3) orthonormal bases perpendicular to each relative velocity."""
    vh = vrel / np.linalg.norm(vrel, axis=1, keepdims=True)
    z = np.tile([0.0, 0.0, 1.0], (len(vh), 1)); x = np.tile([1.0, 0.0, 0.0], (len(vh), 1))
    e1 = np.cross(vh, z); bad = np.linalg.norm(e1, axis=1) < 1e-6; e1[bad] = np.cross(vh[bad], x[bad])
    e1 /= np.linalg.norm(e1, axis=1, keepdims=True); e2 = np.cross(vh, e1)
    return np.stack([e1, e2], axis=1)


def irls_min_fuel(B, Lv, vrel, miss, T_ref, vinf, tcap, rho=0.0, eps_sched=(0.05, 0.02, 0.01, 0.005, 0.002, 0.001),
                  inner=4, w_t=None, w_v=1.0, cap_v=VINF_MAX - 1e-6, free_mask=None):
    """Linear min sum|T_s| over the full sample set. The flyby time of j is free, so only the B-plane components
    (perpendicular to vrel_j) of the linearised miss must vanish:
         E_j [ sum_s B_js (T_s - T_ref_s) + Lv_j dv + miss_j ] = 0
    and the time shift follows from the along-vrel component. Weighted least norm per IRLS pass (weights 1/|T_s| plus
    the proximal rho; v_inf change weighted w_v), Jacobi-scaled normal equations (2 rows per flyby).
    Returns (T_new, dt_days, dv). w_t is accepted for compatibility and ignored."""
    nf, ns = B.shape[0], B.shape[1]
    free = np.ones(ns, bool) if free_mask is None else free_mask
    E = _bplane(vrel)                                                   # (nf, 2, 3)
    BT = B.transpose(0, 2, 1, 3).reshape(nf, 3, ns * 3)                 # (nf, 3, 3ns)
    A_T = np.einsum('jpa,jac->jpc', E, BT).reshape(2 * nf, ns * 3)
    A_v = np.einsum('jpa,jab->jpb', E, Lv).reshape(2 * nf, 3)
    c = -np.einsum('jpa,ja->jp', E, miss).reshape(-1) + A_T @ T_ref.reshape(-1)
    T = T_ref.copy(); dv = np.zeros(3)
    capw = np.ones(ns)
    for eps in eps_sched:
        for _ in range(inner):
            Tn = np.linalg.norm(T, axis=1)
            d = 1.0 / np.maximum(Tn, eps) * capw + rho
            d = np.where(free, d, 1e12)
            Dinv_T = np.repeat(1.0 / d, 3)
            q_T = rho * T_ref.reshape(-1)
            G = (A_T * Dinv_T[None, :]) @ A_T.T + (A_v / w_v) @ A_v.T
            rhs = c - A_T @ (Dinv_T * q_T)
            sc = 1.0 / np.sqrt(np.maximum(np.diag(G), 1e-300))
            Gs = G * sc[:, None] * sc[None, :]
            y = np.linalg.solve(Gs + 1e-12 * np.eye(len(Gs)), rhs * sc)
            lam = y * sc
            T = (Dinv_T * (q_T + A_T.T @ lam)).reshape(ns, 3)
            dv = (A_v.T @ lam) / w_v
            Tn = np.linalg.norm(T, axis=1)
            over = Tn > tcap
            capw = np.where(over, capw * 4.0, np.maximum(1.0, capw * 0.8))
    Tn = np.linalg.norm(T, axis=1); over = Tn > tcap
    T[over] *= (tcap / Tn[over])[:, None]
    # time shifts from the along-vrel component of the predicted relative position
    dr = miss + np.einsum('jsac,sc->ja', B, T - T_ref) + Lv @ dv
    dt = -np.einsum('ja,ja->j', dr, vrel) / np.einsum('ja,ja->j', vrel, vrel) / DAY
    return T, dt, dv


def optimise(gp, iters=60, tol_km=100.0, rho0=1.0, mu_km=2e4, verbose=True, min_rel_gain=1e-4, log=print,
             eps_sched=(0.01, 0.003, 0.001), inner=2):
    """SCP on gp (mutates gp.Ts, gp.tf, gp.vinf). Merit = fuel [kg] + sum_j max(0, |miss_j| - tol) / mu_km.
    Each iteration: linearise at the current iterate, IRLS min-fuel step with proximal weight rho, backtracking on the
    step fraction; rho grows on rejection and shrinks on success. Stops when feasible and the fuel gain stalls."""
    import time
    def evaluate(Ts, tf, vinf, dense=False):
        Yf, Ys = gp.integrate(Ts=Ts, tf=tf, vinf=vinf, dense_samples=dense)
        dm, vr = gp.misses(Yf, tf)
        d = np.linalg.norm(dm, axis=1)
        fuel = gp.fuel(Ts)
        merit = fuel + np.maximum(0.0, d - tol_km).sum() / mu_km
        return dict(Yf=Yf, Ys=Ys, dm=dm, vr=vr, d=d, fuel=fuel, merit=merit)
    cur = evaluate(gp.Ts, gp.tf, gp.vinf, dense=True)
    rho = rho0; hist = []; alpha_last = 0.25
    for it in range(iters):
        tic = time.time()
        B, Lv = gp.linearise(cur["Yf"], cur["Ys"])
        T_new, dt, dv = irls_min_fuel(B, Lv, cur['vr'], cur['dm'], gp.Ts, gp.vinf, gp.tcap, rho=rho, eps_sched=eps_sched, inner=inner)
        vinf_new = gp.vinf + dv
        nv = np.linalg.norm(vinf_new)
        if nv > VINF_MAX - 1e-6:
            vinf_new *= (VINF_MAX - 1e-6) / nv
        accepted = False
        a_try = min(1.0, 2.0 * alpha_last)
        for alpha in (a_try, a_try / 3, a_try / 10):
            Ts_a = gp.Ts + alpha * (T_new - gp.Ts); tf_a = gp.tf + alpha * dt * DAY; vinf_a = gp.vinf + alpha * (vinf_new - gp.vinf)
            tf_a = np.minimum(tf_a, T_MISSION - 2 * DAY)
            cand = evaluate(Ts_a, tf_a, vinf_a, dense=True)
            if cand['merit'] < cur['merit'] - 1e-6:
                accepted = True; break
        if accepted:
            gp.Ts, gp.tf, gp.vinf = Ts_a, tf_a, vinf_a
            cur = cand
            rho = max(rho * 0.7, 1e-3) if alpha >= a_try - 1e-12 else rho * 1.5
            alpha_last = alpha
        else:
            alpha_last = max(alpha_last / 10, 0.01)
            rho *= 4.0
        hist.append((it, cur['fuel'], float(cur['d'].max()), rho, accepted))
        if verbose:
            log(f'  it {it:3d}: fuel {cur["fuel"]:7.2f} kg, max miss {cur["d"].max():10.1f} km, median {np.median(cur["d"]):8.1f}, '
                f'rho {rho:8.3g}, {"acc a=%.2f" % alpha if accepted else "REJ"} ({time.time() - tic:.1f} s)')
        if len(hist) > 6 and cur['d'].max() <= tol_km:
            f_old = [h[1] for h in hist[-6:]]
            if max(f_old) - min(f_old) < min_rel_gain * max(f_old) * 6 and all(h[2] <= tol_km for h in hist[-6:]):
                break
        if rho > 1e6:
            break
    return cur, hist


def restore(gp, tol_km=150.0, iters=12, rho=1e3, log=print, verbose=True):
    """Feasibility restoration: IRLS steps with a strong proximal weight (also on flyby-time shifts and v_inf, growing with
    rho) and a line search on the step fraction, accepted when the summed excess miss decreases."""
    import time
    def excess(d):
        return float(np.maximum(0.0, d - tol_km).sum())
    Yf, Ys = gp.integrate(dense_samples=True); dm, vr = gp.misses(Yf); d = np.linalg.norm(dm, axis=1)
    for it in range(iters):
        if d.max() <= tol_km:
            break
        tic = time.time()
        B, Lv = gp.linearise(Yf, Ys)
        T_new, dt, dv = irls_min_fuel(B, Lv, vr, dm, gp.Ts, gp.vinf, gp.tcap, rho=rho, w_v=1.0 + 1e-4 * rho)
        vinf_new = gp.vinf + dv; nv = np.linalg.norm(vinf_new)
        if nv > VINF_MAX - 1e-6:
            vinf_new *= (VINF_MAX - 1e-6) / nv
        acc = None
        for alpha in (1.0, 0.5, 0.25, 0.1):
            Ts_a = gp.Ts + alpha * (T_new - gp.Ts); vinf_a = gp.vinf + alpha * (vinf_new - gp.vinf)
            tf_a = np.minimum(gp.tf + alpha * dt * DAY, T_MISSION - 2 * DAY)
            Yf2, Ys2 = gp.integrate(Ts=Ts_a, tf=tf_a, vinf=vinf_a, dense_samples=True); dm2, vr2 = gp.misses(Yf2, tf_a)
            d2 = np.linalg.norm(dm2, axis=1)
            if excess(d2) < excess(d):
                acc = alpha; break
        if acc is not None:
            gp.Ts, gp.tf, gp.vinf = Ts_a, tf_a, vinf_a
            Yf, Ys, dm, vr, d = Yf2, Ys2, dm2, vr2, d2
            rho = max(rho / 2, 10.0) if acc == 1.0 else rho * 2
        else:
            rho *= 4
        if verbose:
            log(f'  restore {it:2d}: fuel {gp.fuel():7.2f} kg, max miss {d.max():9.1f} km, rho {rho:8.3g}, '
                f'{"a=%.2f" % acc if acc is not None else "REJ"} ({time.time() - tic:.1f} s)')
    return d


def write_craft(gp, out, eph=None, header=None):
    """Rows with the validator integrator (lowthrust.build_trajectory) + validation. Returns (report, fuel, flyby dists)."""
    from .lowthrust import build_trajectory
    from .submission import write_submission
    from .validator import validate
    from .thrust import arc_max_thrust
    pk, _ = arc_max_thrust(gp.ts, gp.Ts, n_sub=40)
    fb = [(float(t), int(a), 0.0) for t, a in zip(gp.tf, gp.asts)]
    traj = build_trajectory(gp.tL, gp.rE, gp.vE + gp.vinf, gp.m0, gp.ts, gp.Ts, fb)
    write_submission(out, [traj], header=header or f'globalopt m0={gp.m0}')
    rep = validate(out, eph=eph or gp.eph, verbose=False)
    return rep, gp.m0 - traj.m, pk


def add_flyby(gp, ast, t):
    """Insert a flyby constraint (ast at time t) into the problem, keeping the flyby list sorted by time."""
    tf = np.append(gp.tf, t); asts = list(gp.asts) + [int(ast)]
    o = np.argsort(tf, kind='stable')
    gp.tf = tf[o]; gp.asts = [asts[i] for i in o]


def remove_flybys(gp, asts):
    keep = [i for i, a in enumerate(gp.asts) if a not in set(asts)]
    gp.tf = gp.tf[keep]; gp.asts = [gp.asts[i] for i in keep]


def state_of(gp):
    return dict(tL=float(gp.tL), m0=float(gp.m0), ts=gp.ts.copy(), Ts=gp.Ts.copy(), tf=np.array(gp.tf, float), asts=np.array(gp.asts, int),
                vinf=np.array(gp.vinf, float))


def problem_from_state(eph, st, hstep=0.25 * DAY, tcap=0.45):
    rE, vE = eph.earth_state(float(st['tL']))
    craft = dict(tL=float(st['tL']), r0=rE, v0=vE + np.asarray(st['vinf'], float), m0=float(st['m0']), ts=np.asarray(st['ts'], float),
                 Ts=np.asarray(st['Ts'], float), flybys=[(float(t), int(a)) for t, a in zip(st['tf'], st['asts'])], m_end=np.nan)
    return GlobalProblem(eph, craft, tcap=tcap, hstep=hstep)


def save_state(st, path):
    np.savez(path, **st)


def load_state(path):
    z = np.load(path)
    return {k: z[k] for k in z.files}


def tighten(gp, margin=1.5, log=print, verbose=False):
    """Exact tank resizing: with the acceleration history fixed, m(t) = m0 exp(-int|a|/ve) is proportional to m0, so
    scaling every thrust sample by k = m0_new / m0_old leaves the trajectory unchanged. k is chosen so that the model end
    mass is 600 + margin: m0_new = (600 + margin) / (1 - F / m0_old). A short restoration absorbs the sample-norm effects."""
    F = gp.fuel(); k = (M_DRY + margin) / (1.0 - F / gp.m0) / gp.m0
    gp.Ts = gp.Ts * k; gp.m0 = gp.m0 * k
    return restore(gp, log=log, verbose=verbose)


def extend_arc(gp, t_end):
    """Extend the sample grid (zero thrust) up to t_end so flybys after the current last sample can be targeted."""
    h = gp.h; t_last = gp.ts[-1]
    if t_end <= t_last:
        return
    n = int(np.ceil((min(t_end, T_MISSION - DAY) - t_last) / h))
    if n <= 0:
        return
    gp.ts = np.concatenate([gp.ts, t_last + h * (np.arange(n) + 1)]); gp.Ts = np.concatenate([gp.Ts, np.zeros((n, 3))])


def to_massless(gp):
    """Thrust problem -> mass-free problem with the same trajectory: samples T_s m0 / m(t_s) act as accelerations / m0."""
    import copy
    Yf, Ys = gp.integrate(dense_samples=True)
    g2 = copy.copy(gp); g2.Ts = gp.Ts * (gp.m0 / Ys[:, 6])[:, None]; g2.massless = True
    return g2


def to_thrust(gm, margin=1.5):
    """Mass-free problem -> thrust problem with the tank sized for the delta-v: m0 = (600 + margin) exp(dv / ve),
    T_s = a_s m(t_s), m(t_s) = m0 exp(-dv(t_s) / ve) (cumulative sample delta-v, half-sample centred)."""
    import copy
    A = gm.Ts / gm.m0 * 1e-3                                           # km/s^2
    an = np.linalg.norm(A, axis=1) * gm.h
    dv_s = np.cumsum(an) - 0.5 * an
    M0 = (M_DRY + margin) * np.exp(an.sum() / VE)
    g2 = copy.copy(gm); g2.massless = False; g2.m0 = M0
    g2.Ts = gm.Ts / gm.m0 * (M0 * np.exp(-dv_s / VE))[:, None]
    return g2


def polish_ml(gp, iters=60, margin=1.5, log=print, verbose=False, restore_first=True):
    """Thrust problem -> mass-free SCP (restoration, L1 optimisation, restoration) -> thrust problem with the tank sized
    for the delta-v -> thrust-space restoration -> exact tightening. Returns (thrust problem, final misses)."""
    gm = to_massless(gp)
    if restore_first:
        restore(gm, tol_km=300, iters=40, rho=1e3, log=log, verbose=verbose)
    if iters > 0:
        optimise(gm, iters=iters, rho0=100.0, log=log, verbose=verbose)
    restore(gm, tol_km=150, iters=20, log=log, verbose=verbose)
    gt = to_thrust(gm, margin)
    restore(gt, tol_km=150, iters=20, log=log, verbose=verbose)
    d = tighten(gt, margin, log=log, verbose=verbose)
    return gt, d



def _bplane_part(gp, items):
    """Aim-point offsets d0[ast] = B-plane part of the current miss at each inserted flyby (items = [(ast, t), ...])."""
    Yf, _ = gp.integrate(dense_samples=True)
    dm, vr = gp.misses(Yf)
    d0 = {}
    for ast, t in items:
        j = [i for i, a in enumerate(gp.asts) if a == ast and abs(gp.tf[i] - t) < 1e-6][0]
        vh = vr[j] / np.linalg.norm(vr[j])
        d0[int(ast)] = dm[j] - np.dot(dm[j], vh) * vh
    return d0


def insert_block(gp, items, tol_km=100.0, rho=100.0, ds0=0.15, ds_min=0.01, iters_stage=6, log=print, verbose=False,
                 tol_stage=None):
    """Insert one or several flybys at once with an ADAPTIVE joint aim-point homotopy.

    items = [(ast, t), ...]: every aim point starts at the craft's current B-plane position and all of them move to
    their asteroids together, parameterised by s in [0, 1]. After each trial step the flybys are restored; a step that
    leaves a miss above the stage tolerance is REVERTED (thrust, flyby times, v_inf) and retried with half the step, so
    the homotopy backtracks instead of diverging (the fixed-stage `insert_homotopy` diverges irrecoverably once one
    stage fails). Returns the final misses; on failure the problem is left at the last feasible s < 1 and the returned
    misses are large.
    """
    items = sorted([(int(a), float(t)) for a, t in items], key=lambda it: it[1])
    extend_arc(gp, max(max(t for _, t in items) + 5 * DAY, gp.ts[-1]))
    for ast, t in items:
        add_flyby(gp, ast, t)
    d0 = _bplane_part(gp, items)
    dmax0 = max(np.linalg.norm(v) for v in d0.values())
    tol_stage = tol_stage if tol_stage is not None else max(10.0 * tol_km, 0.02 * dmax0)
    gp.offsets = dict(getattr(gp, 'offsets', {}) or {})
    save = lambda: (gp.Ts.copy(), gp.tf.copy(), np.array(gp.vinf, float))
    keep = save(); s = 0.0; ds = ds0; d = np.zeros(1)
    while s < 1.0 - 1e-9:
        st = min(1.0, s + ds); last = st >= 1.0 - 1e-9
        for ast, v in d0.items():
            gp.offsets[ast] = (1.0 - st) * v
        d = restore(gp, tol_km=tol_km if last else tol_stage, iters=4 * iters_stage if last else iters_stage,
                    rho=rho, log=log, verbose=verbose)
        if d.max() <= (tol_km if last else tol_stage) or (last and d.max() < 1e4):
            s = st; keep = save(); ds = min(ds * 1.5, 0.3)
            if verbose:
                log(f'  block homotopy s {s:.3f} (+{ds:.3f}): max miss {d.max():.0f} km, fuel {gp.fuel():.1f}')
        else:
            gp.Ts, gp.tf, gp.vinf = keep
            ds *= 0.5
            if verbose:
                log(f'  block homotopy s {s:.3f}: step to {st:.3f} rejected (miss {d.max():.0f} km), ds -> {ds:.3f}')
            if ds < ds_min:
                for ast in d0:
                    gp.offsets.pop(ast, None)
                for ast, v in d0.items():
                    gp.offsets[ast] = (1.0 - s) * v
                Yf, _ = gp.integrate(dense_samples=True)
                dm, _ = gp.misses(Yf)
                for ast in d0:
                    gp.offsets.pop(ast, None)
                Yf, _ = gp.integrate(dense_samples=True)
                dm, _ = gp.misses(Yf)
                return np.linalg.norm(dm, axis=1)
    for ast in d0:
        gp.offsets.pop(ast, None)
    d = restore(gp, tol_km=tol_km, iters=30, rho=rho, log=log, verbose=verbose)
    return d


def insert_homotopy(gp, ast, t, stages=10, iters_stage=4, rho=100.0, tol_km=300.0, log=print, verbose=False):
    """Add a flyby of ast near time t with an aim-point homotopy: the target point moves from the current relative position
    to the asteroid in `stages` equal steps (B-plane part only), a short restoration after each. Returns final misses."""
    extend_arc(gp, max(t + 5 * DAY, gp.ts[-1]))
    add_flyby(gp, ast, t)
    Yf, _ = gp.integrate(dense_samples=True)
    j = [i for i, a in enumerate(gp.asts) if a == ast and abs(gp.tf[i] - t) < 1e-6][0]
    dm, vr = gp.misses(Yf)
    vh = vr[j] / np.linalg.norm(vr[j]); d0 = dm[j] - np.dot(dm[j], vh) * vh
    gp.offsets = dict(getattr(gp, 'offsets', {}) or {})
    stages = max(stages, int(np.ceil(np.linalg.norm(d0) / 1.5e6)))
    for k in range(1, stages + 1):
        gp.offsets[ast] = (1.0 - k / stages) * d0
        tol = tol_km if k == stages else max(tol_km, 0.02 * np.linalg.norm(d0) / stages)
        d = restore(gp, tol_km=tol, iters=iters_stage if k < stages else 4 * iters_stage, rho=rho, log=log, verbose=verbose)
        if verbose:
            log(f'  homotopy stage {k}/{stages}: max miss {d.max():.0f} km, fuel {gp.fuel():.1f}')
    del gp.offsets[ast]
    if tol_km < d.max() < 1e5:
        d = restore(gp, tol_km=tol_km, iters=30, rho=rho, log=log, verbose=verbose)
    return d
