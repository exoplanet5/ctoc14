"""Insert uncovered targets into existing planned tours with the Lambert-junction model of search.py.

A planned tour (search.tour_json) is a chain of Lambert arcs: launch at Earth (t_launch, v_inf), flybys of asteroids
asts[i] at times[i]; leg i ends at flyby i. Junction dv at the start of leg i (identical to search.roots / search.expand):
    dv[0] = max(0, |v1[0] - vE| - VINF_MAX)         (excess over the free launch v_inf)
    dv[i] = |v1[i] - v2[i-1]|                        (i >= 1; v1/v2 = Lambert departure/arrival velocities)
    feasible if dv <= eta * a * tof,  a = TMAX / m (m = mass at the start of the leg)
    propellant dm = m (1 - exp(-kappa dv / VE)),  kappa = 1 + dv / (2 a tof)   (leverage penalty)

Insertion of target j into gap (k, k+1), k = -1 (launch -> first flyby) .. n-1 (after the last flyby): detour
k -> j -> k+1 with the intermediate flyby time t_j on a grid; optionally the flyby time of k+1 is shifted by delta
(|delta| <= delta_max); the following leg k+1 -> k+2 is then shortened by delta and its Lambert arc is recomputed, so the
junctions at k, j, k+1 and k+2 change and the rest of the timeline is untouched. Every leg of the modified tour is
re-checked with the feasibility rule and the propellant is re-accumulated along the whole tour (later legs start with a
lower mass). All (gap, t_j, delta) candidates of one (tour, target) pair are evaluated in one vectorised pass.

Budget rule: new total propellant <= m0 - M_DRY - margin (margin default 60 kg), i.e. the extra propellant must fit the
tour's remaining propellant minus the margin; new flyby times must stay <= T_MISSION - t_end_buffer.
"""
import json, pathlib
import numpy as np
from dataclasses import dataclass, field
from .constants import VE, TMAX, T_MISSION, DAY, M_DRY, M0_MAX, VINF_MAX, cost_sc
from .lambert import lambert
from .search import UNREACHABLE

SLACK_TOL = 1e-7   # km/s tolerance on eta*a*tof - dv (the planner enforces >= 0 with the same formulas)


@dataclass
class InsParams:
    eta: float = 0.6                      # planner feasibility rule dv <= eta * a * tof
    margin: float = 60.0                  # kg kept above dry mass: fuel_new <= m0 - M_DRY - margin
    t_step: float = 1.0 * DAY             # grid step of the intermediate flyby time
    tof_min: float = 15.0 * DAY           # minimum leg duration (planner grid starts at 15 d)
    tof_max_last: float = 400.0 * DAY     # maximum leg duration for an insertion after the last flyby
    delta_max: float = 10.0 * DAY         # |shift| of the flyby time of k+1
    delta_step: float = 2.5 * DAY
    t_end_buffer: float = 5.0 * DAY       # every flyby must satisfy t <= T_MISSION - t_end_buffer
    launch_gap: bool = True               # allow the detour launch -> j -> first flyby (launch epoch fixed, v_inf re-derived)
    u_max: float = 1.0                    # accept only candidates whose modified legs use <= u_max of eta*a*tof
    try_unreachable: bool = False         # also try 131 / 144

    def deltas(self):
        if self.delta_max <= 0 or self.delta_step <= 0:
            return np.array([0.0])
        n = int(round(self.delta_max / self.delta_step))
        return np.arange(-n, n + 1) * self.delta_step


# ---------------------------------------------------------------------------------------------- planner mass model
def mass_chain(m0, dv, tof, eta):
    """Propellant along a chain of legs with the planner's rule. dv, tof: (..., L). m0: scalar or (...).
    Returns (dm, m_start, slack, m_end) with slack = eta*a*tof - dv per leg (feasible iff slack >= 0)."""
    dv = np.asarray(dv, float); tof = np.asarray(tof, float)
    L = dv.shape[-1]
    m = np.array(np.broadcast_to(np.asarray(m0, float), dv.shape[:-1]), dtype=float)
    dm = np.zeros_like(dv); ms = np.zeros_like(dv); slack = np.zeros_like(dv)
    for i in range(L):
        d = dv[..., i]; T = tof[..., i]
        with np.errstate(invalid='ignore', divide='ignore', over='ignore'):
            a = TMAX / m * 1e-3                               # km/s^2
            slack[..., i] = eta * a * T - d
            kappa = 1.0 + d / (2.0 * a * T)
            dmi = m * (1.0 - np.exp(-kappa * d / VE))
        ms[..., i] = m; dm[..., i] = dmi
        m = m - dmi
    return dm, ms, slack, m


# ---------------------------------------------------------------------------------------------- tour reconstruction
@dataclass
class TourModel:
    name: str
    t_launch: float
    m0: float
    asts: np.ndarray        # (n,) asteroid IDs (1-based)
    times: np.ndarray       # (n,) flyby times [s]
    r: np.ndarray           # (n,3) asteroid positions at the flyby times
    v1: np.ndarray          # (n,3) Lambert departure velocity of leg i
    v2: np.ndarray          # (n,3) Lambert arrival velocity of leg i (= spacecraft velocity at flyby i)
    dv: np.ndarray          # (n,) junction dv at the start of leg i
    tof: np.ndarray         # (n,) leg durations
    dm: np.ndarray          # (n,) propellant of leg i
    m_start: np.ndarray     # (n,) mass at the start of leg i
    slack: np.ndarray       # (n,) eta*a*tof - dv
    vinf: np.ndarray        # (3,) launch v_inf vector (|.| <= VINF_MAX)
    fuel: float
    eta: float
    rE: np.ndarray = field(repr=False, default=None)
    vE: np.ndarray = field(repr=False, default=None)
    extra: dict = field(default_factory=dict)   # bookkeeping (inserted targets, ...)

    @property
    def n(self):
        return len(self.asts)

    @property
    def t_end(self):
        return float(self.times[-1]) if self.n else self.t_launch

    @property
    def feasible(self):
        return bool(np.all(self.slack >= -SLACK_TOL)) and np.isfinite(self.fuel)

    def budget(self, margin):
        """Propellant still available for insertions under the budget rule."""
        return self.m0 - M_DRY - margin - self.fuel

    def utilisation(self):
        """dv / (eta a tof) per leg (1 = at the planner's limit)."""
        return self.dv / np.maximum(self.dv + self.slack, 1e-12)

    def to_json(self):
        return dict(m0=float(self.m0), t_launch=float(self.t_launch), vinf=[float(x) for x in self.vinf],
                    legs=[dict(ast=int(a), t_flyby=float(t), dv_est=float(d), tof=float(T))
                          for a, t, d, T in zip(self.asts, self.times, self.dv, self.tof)],
                    dv_total_est=float(np.sum(self.dv)), fuel_est=float(self.fuel), n_flybys=int(self.n),
                    t_end=float(self.t_end), inserted=self.extra.get('inserted', []))


def build_tour(eph, name, t_launch, m0, asts, times, eta=0.6):
    """Reconstruct the Lambert chain, junction dv's and the mass profile of a tour from (t_launch, m0, asts, times)."""
    asts = np.asarray(asts, int); times = np.asarray(times, float); n = len(asts)
    rE, vE = eph.earth_state(float(t_launch))
    r = np.array([eph.ast[a - 1].state(t)[0] for a, t in zip(asts, times)]).reshape(n, 3)
    r1 = np.vstack([rE[None, :], r[:-1]]) if n else np.zeros((0, 3))
    tof = np.diff(np.concatenate([[float(t_launch)], times]))
    v1, v2 = lambert(r1, r, tof)
    dv = np.empty(n)
    if n:
        vinf_vec = v1[0] - vE
        vinf_n = np.linalg.norm(vinf_vec)
        dv[0] = max(0.0, vinf_n - VINF_MAX) if np.isfinite(vinf_n) else np.nan
        if np.isfinite(vinf_n) and vinf_n > VINF_MAX:
            vinf_vec = vinf_vec * (VINF_MAX / vinf_n)
        dv[1:] = np.linalg.norm(v1[1:] - v2[:-1], axis=1)
    else:
        vinf_vec = np.zeros(3)
    dm, ms, slack, m_end = mass_chain(m0, dv, tof, eta)
    return TourModel(name=name, t_launch=float(t_launch), m0=float(m0), asts=asts, times=times, r=r, v1=v1, v2=v2,
                     dv=dv, tof=tof, dm=dm, m_start=ms, slack=slack, vinf=vinf_vec, fuel=float(m0 - m_end), eta=eta,
                     rE=rE, vE=vE)


def tour_from_json(eph, tour, name='tour', eta=0.6, m0=None):
    """TourModel from a search.tour_json dict (m0 may be overridden)."""
    TM = build_tour(eph, name, tour['t_launch'], tour['m0'] if m0 is None else m0,
                    [l['ast'] for l in tour['legs']], [l['t_flyby'] for l in tour['legs']], eta=eta)
    TM.extra['inserted'] = list(tour.get('inserted', []))
    return TM


def tour_from_converted(eph, tour, conv, name='tour', eta=0.6, m0=None):
    """TourModel of a CONVERTED fragment: flyby sequence and times from the conversion info (conv['legs'] with
    ast / t_flyby, conv['m0'], conv['fuel'] as produced by convert_fleet.py / lowthrust.convert_tour), launch epoch from
    the planned tour JSON. The Lambert chain through the actual flyby points is the planner's model of that trajectory;
    its propellant is compared with the actual one in TM.extra['fuel_actual']."""
    legs = conv['legs']
    TM = build_tour(eph, name, tour['t_launch'], conv.get('m0', tour['m0']) if m0 is None else m0,
                    [l['ast'] for l in legs], [l['t_flyby'] for l in legs], eta=eta)
    TM.extra['inserted'] = []
    TM.extra['fuel_actual'] = float(conv.get('fuel', np.nan))
    TM.extra['base'] = 'converted'
    return TM


def n_infeasible(TM):
    return int(np.sum(~(TM.slack >= -SLACK_TOL)))


def check_reconstruction(TM, tour):
    """Compare the reconstructed junction dv's / propellant / v_inf with the values stored in the tour JSON."""
    dv_est = np.array([l['dv_est'] for l in tour['legs']], float)
    tof_est = np.array([l['tof'] for l in tour['legs']], float)
    d_dv = np.abs(TM.dv - dv_est)
    out = dict(n=TM.n, max_ddv=float(np.nanmax(d_dv)) if TM.n else 0.0, worst_leg=int(np.nanargmax(d_dv)) if TM.n else -1,
               max_dtof=float(np.max(np.abs(TM.tof - tof_est))) if TM.n else 0.0,
               fuel_recomputed=TM.fuel, fuel_est=float(tour.get('fuel_est', np.nan)),
               dfuel=float(TM.fuel - tour.get('fuel_est', np.nan)),
               dvinf=float(np.linalg.norm(TM.vinf - np.asarray(tour['vinf'], float))) if 'vinf' in tour else np.nan,
               min_slack=float(np.min(TM.slack)) if TM.n else np.nan, feasible=TM.feasible,
               n_nan=int(np.sum(~np.isfinite(TM.dv))))
    return out


# ---------------------------------------------------------------------------------------------- insertion candidates
def evaluate_insertions(eph, TM, j, P):
    """Evaluate every (gap k, t_j, delta) insertion of asteroid j into tour TM in one vectorised pass.
    Returns a dict of flat candidate arrays, or None if there is no candidate at all."""
    n = TM.n
    if n == 0 or j in set(TM.asts.tolist()):
        return None
    body = eph.ast[j - 1]
    deltas = P.deltas(); ND = len(deltas)
    t_prev = np.concatenate([[TM.t_launch], TM.times])          # t_prev[k+1] = start time of gap k
    t_last_ok = T_MISSION - P.t_end_buffer
    ks = list(range(-1 if P.launch_gap else 0, n))
    K = []; TJ = []
    for k in ks:
        t0 = t_prev[k + 1]
        if k < n - 1:
            t1 = min(TM.times[k + 1] + deltas.max() - P.tof_min, t_last_ok)
        else:
            t1 = min(t0 + P.tof_max_last, t_last_ok)
        if t1 < t0 + P.tof_min - 1e-6:
            continue
        tj = np.arange(t0 + P.tof_min, t1 + 1e-6, P.t_step)
        K.append(np.full(len(tj), k, dtype=int)); TJ.append(tj)
    if not K:
        return None
    K = np.concatenate(K); TJ = np.concatenate(TJ); NA = len(K)
    # ---- leg A: k -> j (independent of delta)
    rj, _ = body.state(TJ)
    Kc0 = np.maximum(K, 0); isL = K < 0
    r_from = np.where(isL[:, None], TM.rE[None, :], TM.r[Kc0])
    v_arr = np.where(isL[:, None], TM.vE[None, :], TM.v2[Kc0])
    tofA = TJ - t_prev[K + 1]
    v1A, v2A = lambert(r_from, rj, tofA)
    dvA = np.linalg.norm(v1A - v_arr, axis=1)
    dvA = np.where(isL, np.maximum(0.0, dvA - VINF_MAX), dvA)
    # ---- expand over delta
    A = np.repeat(np.arange(NA), ND); D = np.tile(np.arange(ND), NA)
    Kc = K[A]; TJc = TJ[A]; DEL = deltas[D]; Nc = len(A)
    has_next = Kc < n - 1
    valid = has_next | (DEL == 0.0)
    kn = np.minimum(Kc + 1, n - 1)                                 # index of flyby k+1 (clipped)
    # shifted positions of every flyby body for every delta
    R_shift = np.empty((n, ND, 3))
    for i in range(n):
        R_shift[i] = eph.ast[TM.asts[i] - 1].state(TM.times[i] + deltas)[0]
    # ---- leg B: j -> k+1 (shifted)
    tB_end = TM.times[kn] + DEL
    tofB = tB_end - TJc
    dvB = np.full(Nc, np.nan); v2B = np.full((Nc, 3), np.nan)
    selB = has_next & valid & (tofB >= P.tof_min - 1e-6)
    if selB.any():
        v1b, v2b = lambert(rj[A[selB]], R_shift[kn[selB], D[selB]], tofB[selB])
        dvB[selB] = np.linalg.norm(v1b - v2A[A[selB]], axis=1); v2B[selB] = v2b
    # ---- leg C: k+1 (shifted) -> k+2, depends on (k+1, delta) only
    v1C = np.full((n, ND, 3), np.nan); v2C = np.full((n, ND, 3), np.nan)
    for i in range(n - 1):
        v1c, v2c = lambert(R_shift[i], np.broadcast_to(TM.r[i + 1], (ND, 3)), TM.tof[i + 1] - deltas)
        v1C[i] = v1c; v2C[i] = v2c
        z = np.where(deltas == 0.0)[0]
        if len(z):
            v1C[i, z] = TM.v1[i + 1]; v2C[i, z] = TM.v2[i + 1]
    has_C = (Kc + 2 <= n - 1) & has_next
    dvC = np.full(Nc, np.nan); tofC = np.full(Nc, np.nan)
    if has_C.any():
        s = has_C
        dvC[s] = np.linalg.norm(v1C[kn[s], D[s]] - v2B[s], axis=1)
        tofC[s] = TM.tof[np.minimum(kn[s] + 1, n - 1)] - DEL[s]
    # ---- leg D: junction at k+2 (its departure arc is unchanged, its arrival velocity v2C changed with delta)
    has_D = (Kc + 3 <= n - 1) & has_next
    dvD = np.full(Nc, np.nan)
    if has_D.any():
        s = has_D
        dvD[s] = np.linalg.norm(TM.v1[np.minimum(Kc[s] + 3, n - 1)] - v2C[kn[s], D[s]], axis=1)
    # ---- assemble the modified dv / tof sequences (Nc, n+1)
    L = n + 1
    base = np.arange(L)[None, :]
    idx = np.where(base <= Kc[:, None], base, base - 1)
    idx = np.clip(idx, 0, n - 1)
    DV = TM.dv[idx]; TOF = TM.tof[idx]
    rows = np.arange(Nc)
    DV[rows, Kc + 1] = dvA[A]; TOF[rows, Kc + 1] = tofA[A]
    s = has_next; DV[rows[s], Kc[s] + 2] = dvB[s]; TOF[rows[s], Kc[s] + 2] = tofB[s]
    s = has_C; DV[rows[s], Kc[s] + 3] = dvC[s]; TOF[rows[s], Kc[s] + 3] = tofC[s]
    s = has_D; DV[rows[s], Kc[s] + 4] = dvD[s]
    dm, ms, slack, m_end = mass_chain(TM.m0, DV, TOF, P.eta)
    fuel = TM.m0 - m_end
    with np.errstate(invalid='ignore', divide='ignore'):
        util = DV / np.maximum(DV + slack, 1e-12)                  # dv / (eta a tof)
    # utilisation of the modified legs only (positions k+1 .. k+4)
    mod = (base >= (Kc + 1)[:, None]) & (base <= (Kc + 4)[:, None])
    u_mod = np.where(mod, util, 0.0)
    u_mod = np.where(np.isfinite(u_mod), u_mod, np.inf).max(axis=1)
    # modified legs must satisfy the rule; unchanged legs must not get worse than in the base tour (they cannot: the mass
    # only decreases downstream) -- this keeps bases whose reconstruction is slightly infeasible usable
    slack_base = TM.slack[idx]
    leg_ok = np.all(np.where(mod, slack >= -SLACK_TOL, slack >= np.minimum(slack_base, 0.0) - SLACK_TOL), axis=1) & np.isfinite(fuel)
    time_ok = (TJc <= t_last_ok) & (~has_next | (tB_end <= t_last_ok))
    tof_ok = (tofA[A] >= P.tof_min - 1e-6) & (~has_next | (tofB >= P.tof_min - 1e-6)) & (~has_C | (tofC >= P.tof_min - 1e-6))
    ok_nb = leg_ok & time_ok & tof_ok & valid & (u_mod <= P.u_max + 1e-9)
    ok = ok_nb & (fuel <= TM.m0 - M_DRY - P.margin)
    extra = fuel - TM.fuel
    return dict(K=Kc, tj=TJc, delta=DEL, extra=extra, fuel=fuel, ok=ok, ok_nobudget=ok_nb, leg_ok=leg_ok,
                dvA=dvA[A], dvB=dvB, dvC=dvC, dvD=dvD, tofA=tofA[A], tofB=tofB, tofC=tofC, u_mod=u_mod, DV=DV, TOF=TOF,
                SLACK=slack)


def _summ(TM, j, res, c):
    """Summary dict of candidate c of an evaluate_insertions result."""
    k = int(res['K'][c]); n = TM.n
    return dict(target=int(j), tour=TM.name, k=k, pos=k + 1,
                after_ast=int(TM.asts[k]) if k >= 0 else 0, before_ast=int(TM.asts[k + 1]) if k + 1 < n else 0,
                t_j=float(res['tj'][c]), t_j_d=float(res['tj'][c] / DAY), delta_d=float(res['delta'][c] / DAY),
                extra_kg=float(res['extra'][c]), fuel_new=float(res['fuel'][c]),
                dvA=float(res['dvA'][c]), dvB=float(res['dvB'][c]), dvC=float(res['dvC'][c]), dvD=float(res['dvD'][c]),
                tofA_d=float(res['tofA'][c] / DAY), tofB_d=float(res['tofB'][c] / DAY), u_mod=float(res['u_mod'][c]),
                dv_seq=res['DV'][c].tolist(), tof_seq=res['TOF'][c].tolist(), slack_seq=res['SLACK'][c].tolist())


def best_insertion(eph, TM, j, P):
    """Cheapest feasible insertion of j into TM (by extra propellant) and the cheapest one ignoring the budget.
    Returns dict(ok, best, best_nobudget, n_ok, n_ok_nobudget, budget) (best entries None if nothing is feasible)."""
    res = evaluate_insertions(eph, TM, j, P)
    out = dict(target=int(j), tour=TM.name, ok=False, best=None, best_nobudget=None, n_ok=0, n_ok_nobudget=0,
               budget=float(TM.budget(P.margin)), reason='')
    if res is None:
        out['reason'] = 'already in tour' if j in set(TM.asts.tolist()) else 'no gap'
        return out
    ok = res['ok']; ok_nb = res['ok_nobudget']
    out['n_ok'] = int(ok.sum()); out['n_ok_nobudget'] = int(ok_nb.sum())
    if ok_nb.any():
        c = int(np.argmin(np.where(ok_nb, res['extra'], np.inf)))
        out['best_nobudget'] = _summ(TM, j, res, c)
    if ok.any():
        c = int(np.argmin(np.where(ok, res['extra'], np.inf)))
        out['best'] = _summ(TM, j, res, c); out['ok'] = True
    else:
        out['reason'] = 'over budget' if ok_nb.any() else 'no feasible detour'
    return out


def apply_insertion(eph, TM, cand, P, m0=None):
    """Insert cand['target'] at position cand['pos'] (after flyby k), shift flyby k+1 by delta, rebuild the tour."""
    k = cand['k']; n = TM.n
    asts = list(TM.asts); times = list(TM.times)
    if k + 1 < n:
        times[k + 1] = times[k + 1] + cand['delta_d'] * DAY
    asts.insert(k + 1, cand['target']); times.insert(k + 1, cand['t_j'])
    new = build_tour(eph, TM.name, TM.t_launch, TM.m0 if m0 is None else m0, asts, times, eta=P.eta)
    new.extra = dict(TM.extra); new.extra['inserted'] = list(TM.extra.get('inserted', [])) + [int(cand['target'])]
    return new


def min_m0_for(dv_seq, tof_seq, eta, margin, m_lo, m_hi=M0_MAX, tol=0.05, base_slack=None):
    """Smallest m0 in [m_lo, m_hi] such that the chain (dv_seq, tof_seq) is leg-feasible and its propellant fits
    m0 - M_DRY - margin. Returns None if impossible. (fuel(m0) grows slower than m0, so the fit condition is monotone;
    leg feasibility only gets worse with m0, so it is checked at the smallest fitting m0.)"""
    dv = np.asarray(dv_seq, float); tof = np.asarray(tof_seq, float)
    base_bad = ~(np.asarray(base_slack, float) >= -SLACK_TOL) if base_slack is not None else np.zeros(len(dv), bool)
    def fits(m):
        _, _, slack, m_end = mass_chain(m, dv, tof, eta)
        return (m - m_end) <= m - M_DRY - margin, bool(np.all((slack >= -SLACK_TOL) | base_bad))
    f_hi, _ = fits(m_hi)
    if not f_hi:
        return None
    lo, hi = m_lo, m_hi
    f_lo, _ = fits(lo)
    if f_lo:
        hi = lo
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if fits(mid)[0]:
            hi = mid
        else:
            lo = mid
    ok_fit, ok_leg = fits(hi)
    return float(hi) if (ok_fit and ok_leg) else None


# ---------------------------------------------------------------------------------------------- greedy driver
_EPH = None
def _init_worker(eph):
    global _EPH
    _EPH = eph

def _job(args):
    TM, j, P = args
    return best_insertion(_EPH, TM, j, P)


def insert_targets(eph, tours, targets, P=None, nproc=1, verbose=True, grow_m0=False):
    """Greedy insertion: evaluate every (tour, target), apply the globally cheapest feasible insertion (extra propellant),
    re-evaluate the modified tour, repeat. With grow_m0, a second phase enlarges the tank of non-full tours (m0 < 2000)
    when that makes the cheapest leg-feasible insertion fit, ranked by the cost increase dJ = J(m0') - J(m0) (< 1).
    tours: dict name -> TourModel (modified in place). Returns (tours, history, table) where table[(name, j)] holds the
    last evaluation for every remaining target (cheapest candidate ignoring the budget)."""
    global _EPH
    _EPH = eph
    P = P or InsParams()
    remaining = [int(j) for j in targets if P.try_unreachable or int(j) not in UNREACHABLE]
    skipped = [int(j) for j in targets if int(j) not in remaining]
    if verbose and skipped:
        print(f'skipping unreachable targets {skipped}', flush=True)
    pool = None
    if nproc > 1:
        import multiprocessing as mp
        pool = mp.get_context('fork').Pool(nproc, initializer=_init_worker, initargs=(eph,))
    table = {}; history = []
    try:
        def evaluate(names, js):
            jobs = [(tours[nm], j, P) for nm in names for j in js]
            if not jobs:
                return
            res = pool.map(_job, jobs, chunksize=1) if pool is not None else [_job(a) for a in jobs]
            for (TM, j, _), r in zip(jobs, res):
                table[(TM.name, j)] = r
        stale = list(tours)
        while True:
            evaluate(stale, remaining)
            cands = [(r['best']['extra_kg'], nm, j, r['best']) for (nm, j), r in table.items() if r['ok'] and j in remaining]
            if not cands:
                break
            extra, nm, j, c = min(cands, key=lambda x: x[0])
            old = tours[nm]
            new = apply_insertion(eph, old, c, P)
            k = c['k']; mod_ok = bool(np.all(new.slack[k + 1:k + 5] >= -SLACK_TOL))
            if abs(new.fuel - c['fuel_new']) > 1e-3 or not mod_ok or n_infeasible(new) > n_infeasible(old):
                # the rebuilt chain disagrees with the vectorised evaluation (should not happen): discard this candidate
                if verbose:
                    print(f'  !! rebuilt tour disagrees for {j} in {nm} (fuel {new.fuel:.3f} vs {c["fuel_new"]:.3f}, feasible={new.feasible}); skipped', flush=True)
                table[(nm, j)]['ok'] = False; table[(nm, j)]['reason'] = 'rebuild mismatch'
                continue
            tours[nm] = new
            rec = dict(phase='budget', **c, budget_before=float(old.budget(P.margin)), budget_after=float(new.budget(P.margin)),
                       fuel_before=float(old.fuel), n_before=int(old.n), n_after=int(new.n), m0=float(new.m0), dJ=0.0)
            for key in ('dv_seq', 'tof_seq', 'slack_seq'): rec.pop(key, None)
            history.append(rec)
            remaining.remove(j)
            if verbose:
                print(f'insert {j:4d} into {nm} after {c["after_ast"]:3d} before {c["before_ast"]:3d} at t={c["t_j_d"]:7.1f} d '
                      f'delta={c["delta_d"]:+5.1f} d: extra {c["extra_kg"]:6.2f} kg (dv {c["dvA"]:.2f}/{c["dvB"]:.2f}/{c["dvC"]:.2f} km/s, '
                      f'u={c["u_mod"]:.2f}); budget left {new.budget(P.margin):7.1f} kg', flush=True)
            table = {key: v for key, v in table.items() if key[0] != nm and key[1] != j}
            stale = [nm]
        if grow_m0:
            stale = [nm for nm in tours if tours[nm].m0 < M0_MAX - 1e-6]
            while True:
                evaluate(stale, remaining)
                cands = []
                for (nm, j), r in table.items():
                    if j not in remaining or r['best_nobudget'] is None or tours[nm].m0 >= M0_MAX - 1e-6:
                        continue
                    c = r['best_nobudget']
                    m0n = min_m0_for(c['dv_seq'], c['tof_seq'], P.eta, P.margin, tours[nm].m0, base_slack=c.get('slack_seq'))
                    if m0n is None:
                        continue
                    dJ = float(cost_sc(m0n) - cost_sc(tours[nm].m0))
                    if dJ < 1.0:
                        cands.append((dJ, nm, j, c, m0n))
                if not cands:
                    break
                dJ, nm, j, c, m0n = min(cands, key=lambda x: x[0])
                old = tours[nm]
                new = apply_insertion(eph, old, c, P, m0=m0n)
                k = c['k']; mod_ok = bool(np.all(new.slack[k + 1:k + 5] >= -SLACK_TOL))
                if not mod_ok or n_infeasible(new) > n_infeasible(old) or new.budget(P.margin) < -1e-6:
                    if verbose:
                        print(f'  !! grow-m0 rebuild infeasible for {j} in {nm}; skipped', flush=True)
                    table[(nm, j)]['best_nobudget'] = None
                    continue
                tours[nm] = new
                rec = dict(phase='grow_m0', **c, budget_before=float(old.budget(P.margin)), budget_after=float(new.budget(P.margin)),
                           fuel_before=float(old.fuel), n_before=int(old.n), n_after=int(new.n), m0=float(new.m0), m0_before=float(old.m0), dJ=dJ)
                for key in ('dv_seq', 'tof_seq', 'slack_seq'): rec.pop(key, None)
                history.append(rec)
                remaining.remove(j)
                if verbose:
                    print(f'insert {j:4d} into {nm} after {c["after_ast"]:3d} before {c["before_ast"]:3d} at t={c["t_j_d"]:7.1f} d '
                          f'delta={c["delta_d"]:+5.1f} d with m0 {old.m0:.1f} -> {m0n:.1f} kg (dJ = {dJ:.3f}); fuel {new.fuel:.1f} kg', flush=True)
                table = {key: v for key, v in table.items() if key[0] != nm and key[1] != j}
                stale = [nm]
    finally:
        if pool is not None:
            pool.close(); pool.join()
    return tours, history, table, remaining


def load_tours(eph, paths, eta=0.6, names=None):
    """Load tour JSON files into TourModels keyed by name (file stem by default)."""
    tours = {}
    for i, p in enumerate(paths):
        tour = json.load(open(p))
        nm = names[i] if names else pathlib.Path(p).stem
        tours[nm] = tour_from_json(eph, tour, name=nm, eta=eta)
    return tours
