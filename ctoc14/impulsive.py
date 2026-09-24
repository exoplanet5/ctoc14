"""Fast impulsive twin of the whole-trajectory model (ctoc14/globalopt.py).

A craft = launch (Earth state + v_inf) + impulses Dv_m at fixed nodes tau_m (every `dtau` days) + Kepler arcs between
nodes; flybys at free times t_j. Same interface as globalopt.GlobalProblem (Ts = impulses [km/s], tf, vinf, tcap,
integrate, misses, linearise, fuel), so globalopt.optimise / restore / insert_homotopy / add_flyby / remove_flybys work
unchanged. Cost: fuel(Ts) = tank-sized propellant (600 + margin) (exp(sum|Dv| / ve) - 1) [kg]; J_i follows from
m0 = (600 + margin) exp(sum|Dv| / ve). Linearisation: central-difference Kepler STMs per node segment (one batched
propagation), cumulative products, rows d r(t_j) / d Dv_m = Phi(t_j, tau_m)[r, v]."""
import numpy as np
from .constants import VE, DAY, M_DRY, T_MISSION, VINF_MAX, TMAX
from .kepler import propagate_batch


class ImpulsiveProblem:
    massless = False

    def __init__(self, eph, tL, vinf, ts, Ts, tf, asts, margin=1.5, tcap=1.0, dtau=20 * DAY):
        self.eph = eph; self.tL = float(tL); self.vinf = np.asarray(vinf, float)
        self.ts = np.asarray(ts, float); self.Ts = np.asarray(Ts, float)
        self.tf = np.asarray(tf, float); self.asts = [int(a) for a in asts]
        self.margin = margin; self.tcap_fixed = tcap; self.h = dtau
        rE, vE = eph.earth_state(self.tL); self.rE = rE; self.vE = vE
        self.m0 = M_DRY + margin     # nominal (fuel() is tank-sized)

    @property
    def tcap(self):
        """Impulse cap [km/s]: what 0.43 N delivers over one bin at the current tank mass (continuous-thrust capability),
        never above the fixed cap."""
        return min(self.tcap_fixed, 0.43 * self.h / (self.tank() * 1e3))

    # ---------------------------------------------------------------- dynamics
    def _nodes(self, tf):
        """Sorted node list: launch, impulse nodes, flybys. Returns (times, kinds) with kind -1 launch, s >= 0 impulse,
        -2-j flyby j."""
        tn = np.concatenate([[self.tL], self.ts, tf]); kind = np.concatenate([[-1], np.arange(len(self.ts)), -2 - np.arange(len(tf))])
        o = np.argsort(tn, kind='stable')
        return tn[o], kind[o]

    def integrate(self, Ts=None, tf=None, vinf=None, dense_samples=False):
        """States (6 + dummy mass column) at the flyby times and just BEFORE every impulse node."""
        Ts = self.Ts if Ts is None else Ts; tf = self.tf if tf is None else tf; vinf = self.vinf if vinf is None else vinf
        tn, kind = self._nodes(tf)
        r = self.rE.copy(); v = self.vE + vinf; t = self.tL
        Yf = np.zeros((len(tf), 7)); Ys = np.zeros((len(self.ts), 7))
        # propagate node to node (vectorising is not needed: ~300 short arcs)
        for n in range(1, len(tn)):
            dt = tn[n] - t
            if dt > 0:
                rr, vv = propagate_batch(r[None], v[None], np.array([dt])); r, v = rr[0], vv[0]; t = tn[n]
            k = kind[n]
            if k >= 0:
                Ys[k, :3] = r; Ys[k, 3:6] = v; Ys[k, 6] = 1.0
                v = v + Ts[k]
            else:
                j = -2 - k; Yf[j, :3] = r; Yf[j, 3:6] = v; Yf[j, 6] = 1.0
        return Yf, Ys

    def misses(self, Yf, tf=None):
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

    def dv(self, Ts=None):
        Ts = self.Ts if Ts is None else Ts
        return float(np.linalg.norm(Ts, axis=1).sum())

    def fuel(self, Ts=None):
        return (M_DRY + self.margin) * (np.exp(self.dv(Ts) / VE) - 1.0)

    def tank(self, Ts=None):
        return (M_DRY + self.margin) * np.exp(self.dv(Ts) / VE)

    # ---------------------------------------------------------------- linearisation
    def chain(self, Yf, Ys, Ts=None, vinf=None):
        Ts = self.Ts if Ts is None else Ts; vinf = self.vinf if vinf is None else vinf
        tn, kind = self._nodes(self.tf)
        # state just AFTER each node (impulse applied), launch state after v_inf
        X = np.zeros((len(tn), 6)); X[0, :3] = self.rE; X[0, 3:] = self.vE + vinf
        for n in range(1, len(tn)):
            k = kind[n]
            if k >= 0:
                X[n, :3] = Ys[k, :3]; X[n, 3:] = Ys[k, 3:6] + Ts[k]
            else:
                X[n] = Yf[-2 - k, :6]
        G = len(tn) - 1; dt = np.diff(tn)
        er, ev = 100.0, 1e-4
        pert = np.zeros((12, 6))
        for i in range(6):
            e = er if i < 3 else ev; pert[2 * i, i] = e; pert[2 * i + 1, i] = -e
        B0 = np.repeat(X[:-1], 12, axis=0) + np.tile(pert, (G, 1))
        rp, vp = propagate_batch(B0[:, :3], B0[:, 3:], np.repeat(dt, 12))
        O = np.concatenate([rp, vp], axis=1).reshape(G, 12, 6)
        step = np.array([er, er, er, ev, ev, ev])
        S = ((O[:, 0::2, :] - O[:, 1::2, :]) / (2 * step)[None, :, None]).transpose(0, 2, 1)
        Pi = np.zeros((G + 1, 6, 6)); Pi[0] = np.eye(6)
        for g in range(G):
            Pi[g + 1] = S[g] @ Pi[g]
        return dict(Pi=Pi, Pinv=np.linalg.inv(Pi), tn=tn, kind=kind)

    def rows_at(self, ch, times, nodes=None):
        """B (n, M, 3, 3): d r(t) / d Dv_m; Lv (n, 3, 3): d r(t) / d v_inf. Impulse m acts after its node, so
        d x(node_n) / d Dv_m = Pi[n] Pi[node_m]^-1 [0; I] for node_m < n."""
        tn = ch['tn']; kind = ch['kind']; Pi = ch['Pi']; Pinv = ch['Pinv']; M = len(self.ts)
        node_of_imp = np.zeros(M, int)
        for n, k in enumerate(kind):
            if k >= 0: node_of_imp[k] = n
        n_t = len(times); B = np.zeros((n_t, M, 3, 3)); Lv = np.zeros((n_t, 3, 3))
        for j, t in enumerate(times):
            nj = int(np.searchsorted(tn, t + 1e-6, side='right') - 1) if nodes is None else nodes[j]
            P = np.zeros((3, 6)); P[:, :3] = np.eye(3); P[:, 3:] = np.eye(3) * (t - tn[nj])
            Pr = P @ Pi[nj]
            Lv[j] = Pr[:, 3:6]
            act = node_of_imp < nj
            if np.any(act):
                B[j, act] = np.einsum('ab,mbc->mac', Pr, Pinv[node_of_imp[act]][:, :, 3:6])
        return B, Lv

    def linearise(self, Yf, Ys, Ts=None, vinf=None):
        ch = self.chain(Yf, Ys, Ts, vinf)
        kind = ch['kind']
        node_of_f = {int(-2 - kind[n]): n for n in range(len(kind)) if kind[n] <= -2}
        return self.rows_at(ch, self.tf, nodes=[node_of_f[j] for j in range(len(self.tf))])


def from_thrust_problem(gp, dtau=20 * DAY, margin=1.5):
    """Impulsive twin of a (thrust or mass-free) globalopt problem: impulses = acceleration integrated over dtau bins."""
    Yf, Ys = gp.integrate(dense_samples=True)
    m = Ys[:, 6] if not gp.massless else np.full(len(gp.ts), gp.m0)
    acc = gp.Ts / m[:, None] * 1e-3                         # km/s^2 at the samples
    t_end = max(gp.tf.max(), gp.ts[-1])
    edges = np.arange(gp.tL, t_end + dtau, dtau)
    ts = 0.5 * (edges[:-1] + edges[1:])
    b = np.clip(np.searchsorted(edges, gp.ts, side='right') - 1, 0, len(ts) - 1)
    Ts = np.zeros((len(ts), 3)); np.add.at(Ts, b, acc * gp.h)
    return ImpulsiveProblem(gp.eph, gp.tL, gp.vinf, ts, Ts, gp.tf, gp.asts, margin=margin, dtau=dtau)


LIN_FALLBACK = 2.5        # km/s: a Lambert junction above this seeds the leg with the linear thrust profile (None = off)
LIN_MISS_MAX = 3.0e6      # km: accept the linear seed if its non-linear end-point miss is below this (restore fixes it)


def from_tour(eph, tour, dtau=20 * DAY, margin=1.5, lag=1.0 * DAY):
    """Impulsive twin of a planner tour (search.tour_json record): the legs are rebuilt as a Lambert chain (best of
    0-2 revolutions per leg, judged by the junction dv), the junction impulses sit at exact extra nodes 1 s after each
    departure (`lag` after it: a node closer than the flyby-time shifts of the SCP makes the linearisation invalid as soon
    as a flyby moves past its own departure impulse), the launch v_inf is the first leg's (capped at VINF_MAX with the excess as an impulse); regular dtau nodes
    carry zero. Initial misses are ~0, so globalopt.restore/optimise (run_ialns.settle) can start at once. None if a leg
    has no Lambert solution."""
    from .lambert import lambert_best
    from .linleg import LinLeg
    tL = float(tour['t_launch'])
    asts = [int(l['ast']) for l in tour['legs']]; tf = np.array([float(l['t_flyby']) for l in tour['legs']])
    rE, vE = eph.earth_state(tL)
    R, _ = eph.ast_states_at(np.array(asts) - 1, tf)
    r_prev = rE; v_prev = vE + np.asarray(tour['vinf'], float); t_prev = tL
    imp_t = []; imp_dv = []; vinf = None
    for k in range(len(asts)):
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
            if LIN_FALLBACK is not None and (not np.isfinite(dvj) or dvj > LIN_FALLBACK):
                # no single junction impulse does this leg (near-180 deg / long phasing geometry): seed it with the
                # linearised distributed-thrust profile instead, exactly the leg model the planner used
                tof = tf[k] - t_prev
                L = LinLeg(r_prev, v_prev, np.array([tof]), dtau=10 * DAY)
                cost, varr, _, U = L.solve(0, R[k][None], TMAX / 1000.0 * 1e-3, ub=0.9, return_U=True)
                if np.isfinite(cost[0]) and cost[0] < dvj and L.check(0, U[0], R[k]) < LIN_MISS_MAX:
                    K = L.K[0]
                    for kk in range(K):
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


def settle_tour(eph, tour, settle, lags=(1.0 * DAY, 0.25 * DAY, 0.05 * DAY), tol_km=150.0, screen_tank=2000.0):
    """from_tour + settle with a retry ladder over the junction lag (a long lag gives a large initial miss on some tours,
    a short one breaks the linearisation when a flyby shifts past its departure impulse). settle(ip) returns the max miss.
    `screen_tank` rejects a tour whose Lambert chain is degenerate before paying for the settle.
    Returns (ip, miss, lag) of the first success, or (None, best miss, None)."""
    best = np.inf
    for lag in lags:
        ip = from_tour(eph, tour, lag=lag)
        if ip is None:
            return None, np.inf, None
        # cheap screen: a degenerate Lambert leg shows up at once as an absurd junction impulse (the initial tank runs
        # to 3700-4100 kg against 780-1070 kg for a healthy tour). Settling those costs ~80 s and never converges.
        if screen_tank and float(ip.tank()) > screen_tank:
            return None, float('inf'), None
        miss = float(settle(ip))
        if miss <= tol_km:
            return ip, miss, lag
        best = min(best, miss)
    return None, best, None


def to_exact(ip, eph=None, h=1.0 * DAY, hstep=0.25 * DAY):
    """Mass-free continuous-thrust problem (globalopt) from an impulsive solution: each impulse becomes a constant
    acceleration over its dtau bin on an h sample grid (N-equivalent at the tank mass). Follow with
    globalopt.polish_ml-style restoration / optimisation / to_thrust."""
    from .globalopt import GlobalProblem
    eph = eph or ip.eph
    m0 = ip.tank()
    t_last = max(ip.tf.max(), ip.ts[-1] + ip.h / 2)
    n = int(np.ceil((min(t_last + 5 * DAY, T_MISSION - DAY) - ip.tL) / h))
    ts = ip.tL + h * (np.arange(n) + 1)
    b = np.clip(np.searchsorted(ip.ts - ip.h / 2, ts, side='right') - 1, 0, len(ip.ts) - 1)
    Ts = ip.Ts[b] / ip.h * 1e3 * m0                                   # km/s over dtau -> m/s^2 -> N at mass m0
    inside = np.abs(ts - ip.ts[b]) <= ip.h / 2
    Ts[~inside] = 0.0
    rE, vE = eph.earth_state(ip.tL)
    craft = dict(tL=ip.tL, r0=rE, v0=vE + ip.vinf, m0=m0, ts=ts, Ts=Ts, flybys=list(zip(ip.tf, ip.asts)), m_end=np.nan)
    gp = GlobalProblem(eph, craft, hstep=hstep); gp.massless = True
    return gp


def enforce_cap(ip, stages=6, log=None):
    """Bring every impulse under the continuous-thrust cap (ImpulsiveProblem.tcap) by a homotopy on the cap: clip to a
    cap lowered in `stages` steps from the current largest impulse, restoring the flybys after each clip.
    Returns the final misses."""
    from .globalopt import restore, optimise
    q = lambda s: None
    top = float(np.linalg.norm(ip.Ts, axis=1).max()) if len(ip.Ts) else 0.0
    d = np.zeros(0)
    for k in range(1, stages + 1):
        cap = top + (ip.tcap - top) * k / stages
        imp = np.linalg.norm(ip.Ts, axis=1); over = imp > cap
        if np.any(over):
            ip.Ts[over] *= (cap / imp[over])[:, None]
        d = restore(ip, tol_km=100, iters=30, rho=1e3, log=log or q)
        if len(d) and d.max() > 1e5:
            break
    return d
