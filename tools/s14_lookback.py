"""Stage 14, track A, stage 2 -- THE LOOK-BACK PLANNER: a twin-native forward beam, priced and extended WITHOUT any
settle in the loop.

Label   a whole impulsive-twin trajectory in run_ialns format: launch epoch + v_inf, impulses on the standard twin node
        grid tL + 10 d + 20 d k (the grid of every settled route in the project), flyby list, and the states on it
        (just before every node impulse, at every flyby).  The carried velocity is the velocity OF THAT TRAJECTORY.
Expand  L-leg look-back linearised leg.  For a label with n flybys, anchor a = n - L (the launch when n < L; v_inf is
        then free as well, |v_inf| <= 4).  ONE linearised min-L1 problem over the window [t_a, t_s] for every candidate
        (target s, epoch t_s on a TOF grid):  variables = the impulses on every window node (the last L-1 legs are
        re-optimised, the new leg gets its thrust); rows = the B-plane of each of the L-1 window flybys (epochs free,
        the along-track part gives the epoch shift) + 3 position rows at the new target.  Sensitivities are Kepler
        STMs chained along the label's own trajectory (central differences, symplectic inverse).  IRLS on the L1 norm,
        batched over all candidates of a parent (they share the matrix; only the right-hand side differs), node cap
        |T_k| <= tcap (0.43 N over 20 d at --mcap kg, the twin's own cap).  Marginal price = window L1 cost - the
        parent's L1 cost on the same window.  L = 1 is the single-leg twin-grid linear model (baseline).
Refine  only the survivors of each level: the window is re-integrated nonlinearly from the anchor and corrected by at
        most --ref-iters linearised min-fuel restoration steps (globalopt.irls_min_fuel with a proximal weight, epochs
        free, line search on the miss).  Nothing before the anchor ever moves.  A child whose window misses stay above
        --acc-tol km is dropped.
Score   t_s/T + w_fuel F/1400 - prizes, F = 1600 (1 - exp(-dv/ve)) (production convention, w_fuel 0.52).
Final   the Pareto set (best labels per depth) is settled deepest first with run_ialns.settle (the whole route, once):
        a planner label IS a twin trajectory, so no conversion is needed.

Modes
  free    beam over a target pool:  --pool own:SRC_FLEET:ROUTE | pools:FILE:ROUTE | all | list:1,2,3
  guided  calibration along a source route's own order and epochs (SRC_FLEET ROUTE): every source leg is priced from
          the planner's OWN label (the one this planner built along the same order), refined and appended.
Outputs <out>/log.txt, ckpt.pkl (restart point, every level), result.json, route_<name>.npz (settled, run_ialns
        format), planned.npz (deepest planned label, unsettled)
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, signal, pathlib, argparse, pickle, warnings
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from greedy_cover import BP, ALL, CM, _timeout
from ctoc14.search import roots as search_roots, mask_of, UNREACHABLE
from ctoc14.impulsive import ImpulsiveProblem
from ctoc14.kepler import propagate_batch
from ctoc14.lambert import lambert
from ctoc14.globalopt import irls_min_fuel, _bplane
from ctoc14.constants import DAY, AU, VE, T_MISSION, M_DRY, FUEL_MAX, VINF_MAX
warnings.filterwarnings('ignore')

DTAU = 20 * DAY            # twin node spacing
OFF = 10 * DAY             # first node 10 d after launch (tL + 10 d + 20 d k): the project's standard grid
M0P = 1600.0               # planner mass convention of the score
TOL = 150.0                # km: settled
T_END = T_MISSION - 2 * DAY


def node_times(tL, k):
    return tL + OFF + DTAU * np.asarray(k, float)


def n_before(tL, t):
    """Number of grid nodes strictly before time t."""
    return max(0, int(np.ceil((t - tL - OFF) / DTAU - 1e-12)))


def sinv(P):
    """Inverse of symplectic 6x6 matrices (Kepler STMs): [[A, B], [C, D]]^-1 = [[D^T, -B^T], [-C^T, A^T]]."""
    A = P[..., :3, :3]; B = P[..., :3, 3:]; C = P[..., 3:, :3]; D = P[..., 3:, 3:]
    out = np.empty_like(P); sw = lambda M: np.swapaxes(M, -1, -2)
    out[..., :3, :3] = sw(D); out[..., :3, 3:] = -sw(B); out[..., 3:, :3] = -sw(C); out[..., 3:, 3:] = sw(A)
    return out


_ER, _EV = 100.0, 1e-4
_PERT = np.zeros((13, 6))
for _i in range(6):
    _e = _ER if _i < 3 else _EV
    _PERT[2 * _i, _i] = _e; _PERT[2 * _i + 1, _i] = -_e
_STEP = np.array([_ER] * 3 + [_EV] * 3)


def stms(X, dt):
    """Kepler STMs Phi(t + dt, t) (G, 6, 6) from states X (G, 6) (central differences, one batched propagation) and the
    nominal end states (G, 6)."""
    X = np.atleast_2d(X); G = len(X)
    B0 = np.repeat(X, 13, axis=0) + np.tile(_PERT, (G, 1))
    rp, vp = propagate_batch(B0[:, :3], B0[:, 3:], np.repeat(np.asarray(dt, float), 13))
    O = np.concatenate([rp, vp], 1).reshape(G, 13, 6)
    S = ((O[:, 0:12:2, :] - O[:, 1:12:2, :]) / (2 * _STEP)[None, :, None]).transpose(0, 2, 1)
    return S, O[:, 12, :]


class WindowProblem(ImpulsiveProblem):
    """ImpulsiveProblem whose 'launch' is a fixed anchor state (no v_inf freedom): the look-back window."""
    def __init__(self, eph, t_a, X_a, ts, Ts, tf, asts, tcap):
        self.eph = eph; self.tL = float(t_a); self.vinf = np.zeros(3)
        self.ts = np.asarray(ts, float); self.Ts = np.asarray(Ts, float).reshape(-1, 3)
        self.tf = np.asarray(tf, float); self.asts = [int(x) for x in asts]
        self.margin = 1.5; self.tcap_fixed = tcap; self.h = DTAU
        self.rE = np.asarray(X_a[:3], float).copy(); self.vE = np.asarray(X_a[3:6], float).copy()
        self.m0 = M_DRY + 1.5

    def linearise(self, Yf, Ys, Ts=None, vinf=None):
        B, Lv = ImpulsiveProblem.linearise(self, Yf, Ys, Ts, vinf)
        return B, np.zeros_like(Lv)


# ------------------------------------------------------------------ labels
def make_label(tL, vinf, T, Xn, fa, ft, Xf, miss):
    T = np.asarray(T, float).reshape(-1, 3)
    return dict(tL=float(tL), vinf=np.asarray(vinf, float).copy(), T=T, Xn=np.asarray(Xn, float).reshape(-1, 6),
                fa=[int(x) for x in fa], ft=np.asarray(ft, float), Xf=np.asarray(Xf, float).reshape(-1, 6),
                miss=np.asarray(miss, float), dv=float(np.linalg.norm(T, axis=1).sum()) if len(T) else 0.0)


def label_tank(lab):
    return (M_DRY + 1.5) * np.exp(lab['dv'] / VE)


def label_ist(lab):
    """run_ialns ist dict of a label (the twin trajectory itself)."""
    M = len(lab['T'])
    return dict(tL=float(lab['tL']), vinf=np.array(lab['vinf'], float), ts=node_times(lab['tL'], np.arange(M)),
                Ts=np.array(lab['T'], float), tf=np.array(lab['ft'], float), asts=np.array(lab['fa'], int))


def root_label(tL, vinf, X, t1):
    """Zero-thrust twin of a launch leg (exact Lambert arc, |v_inf| <= 4)."""
    E = RI.eph(); M = n_before(tL, t1)
    ip = ImpulsiveProblem(E, tL, vinf, node_times(tL, np.arange(M)), np.zeros((M, 3)), [t1], [int(X)])
    Yf, Ys = ip.integrate(); dm, _ = ip.misses(Yf)
    return make_label(tL, vinf, np.zeros((M, 3)), Ys[:, :6], [X], [t1], Yf[:, :6], np.linalg.norm(dm, axis=1))


def score_of(dv, t, nprize, P):
    F = M0P * (1.0 - np.exp(-np.asarray(dv, float) / VE))
    return P.w_t * np.asarray(t, float) / T_MISSION + P.w_fuel * F / FUEL_MAX - nprize


def lab_score(lab, P):
    return float(score_of(lab['dv'], lab['ft'][-1], sum(P.prize[x - 1] for x in lab['fa']), P))


# ------------------------------------------------------------------ look-back pricing (batched over candidates)
class Opt:
    """Planner options (argparse namespace subset)."""
    def __init__(self, a, prize):
        self.L = a.L; self.dvcap = a.dvcap; self.drmax = a.drmax * AU; self.tcap = min(1.0, 0.43 * DTAU / (a.mcap * 1e3))
        self.tofs = np.unique(np.concatenate([np.arange(a.tof_min, 150 + 1e-9, a.tof_step1),
                                              np.arange(150 + a.tof_step2, a.tof_max + 1e-9, a.tof_step2)])) * DAY
        self.npt = a.npt; self.mc = a.mc; self.sep = a.sep * DAY
        self.w_t = 1.0; self.w_fuel = a.w_fuel; self.prize = prize
        self.eps = tuple(float(x) for x in a.eps.split(','))
        self.ref_iters = a.ref_iters; self.ref_tol = a.ref_tol; self.acc_tol = a.acc_tol; self.ref_rho = a.ref_rho
        self.opt_iters = a.opt_iters
        self.mcap = a.mcap; self.cap_frac = a.cap_frac
        self.pol_L = a.polish_L; self.pol_iters = a.polish_iters
        self.tank_max = a.tank_max


def tcap_of(lab, O):
    """Node cap for the children of a label: the twin's 0.43 N x 20 d at max(--mcap, cap_frac x planned tank)."""
    if O.cap_frac <= 0:
        return O.tcap
    return min(1.0, 0.43 * DTAU / (max(O.mcap, O.cap_frac * float(label_tank(lab))) * 1e3))


def window_setup(lab, L, tofs):
    """Sensitivities of the L-leg window of a label for every target epoch t_last + tofs.
    Returns a dict with the window node list, intermediate B-plane rows, target rows and the coast reference."""
    E = RI.eph()
    n = len(lab['fa']); a = n - L; tL = lab['tL']; M = len(lab['T'])
    t_last = float(lab['ft'][-1]); X_last = lab['Xf'][-1]
    if a < 0:
        t_a = tL; rE, vE = E.earth_state(tL); X_a = np.concatenate([rE, vE + lab['vinf']]); k0 = 0
    else:
        t_a = float(lab['ft'][a]); X_a = lab['Xf'][a]; k0 = n_before(tL, t_a)
    kw = np.arange(k0, M)                                  # window nodes before t_last
    ints = list(range(max(a, -1) + 1, n))                  # window flybys (the last one included)
    tw = node_times(tL, kw)
    # events: anchor, window nodes, window flybys (impulse before flyby at equal times)
    ev = [(t_a, -1, 0, -1)] + [(float(t), 0, 1, int(k)) for t, k in zip(tw, kw)] + [(float(lab['ft'][i]), 1, 2, i) for i in ints]
    ev.sort(key=lambda z: (z[0], z[1]))
    Xa = []
    for (t, pri, kind, idx) in ev:
        if kind == 0:
            Xa.append(X_a)
        elif kind == 1:
            x = lab['Xn'][idx].copy(); x[3:] += lab['T'][idx]; Xa.append(x)
        else:
            Xa.append(lab['Xf'][idx])
    Xa = np.array(Xa); et = np.array([z[0] for z in ev])
    Pi = np.zeros((len(ev), 6, 6)); Pi[0] = np.eye(6)
    if len(ev) > 1:
        S, _ = stms(Xa[:-1], np.diff(et))
        for g in range(len(ev) - 1):
            Pi[g + 1] = S[g] @ Pi[g]
    e_node = {z[3]: g for g, z in enumerate(ev) if z[2] == 1}
    e_fb = {z[3]: g for g, z in enumerate(ev) if z[2] == 2}
    eL = e_fb.get(n - 1, 0)                              # L = 1: the anchor is the last flyby itself
    PinvN = sinv(Pi[[e_node[int(k)] for k in kw]]) if len(kw) else np.zeros((0, 6, 6))
    W = np.einsum('ab,kbc->kac', Pi[eL], PinvN[:, :, 3:]) if len(kw) else np.zeros((0, 6, 3))   # d X(t_last)/d T_k
    # window flyby rows (B-plane)
    nI = len(ints)
    ra, va = E.ast_states_at(np.array([lab['fa'][i] - 1 for i in ints]), np.array([lab['ft'][i] for i in ints]))
    vrel = np.array([lab['Xf'][i][3:] for i in ints]).reshape(-1, 3) - va.reshape(-1, 3)
    miss_i = np.array([lab['Xf'][i][:3] for i in ints]).reshape(-1, 3) - ra.reshape(-1, 3)
    Eb = _bplane(vrel)                                                      # (nI, 2, 3)
    Bi = np.zeros((nI, len(kw), 3, 3))
    Lvi = np.zeros((nI, 3, 3))
    for q, i in enumerate(ints):
        g = e_fb[i]
        act = np.array([e_node[int(k)] < g for k in kw], bool) if len(kw) else np.zeros(0, bool)
        if act.any():
            Bi[q, act] = np.einsum('ab,kbc->kac', Pi[g][:3], PinvN[act][:, :, 3:])
        Lvi[q] = Pi[g][:3, 3:]
    # coast after the last flyby
    t_s = t_last + tofs; t_s = t_s[t_s <= T_END]
    kc = np.arange(M, n_before(tL, t_s.max()) if len(t_s) else M)
    tc = node_times(tL, kc)
    tau = np.concatenate([tc, t_s])
    C, Xnom = stms(np.repeat(X_last[None], len(tau), 0), tau - t_last) if len(tau) else (np.zeros((0, 6, 6)), np.zeros((0, 6)))
    Cn, Ct = C[:len(tc)], C[len(tc):]; r_ref = Xnom[len(tc):, :3]
    J = len(t_s); K1 = len(kw); K2 = len(kc); K = K1 + K2; m = 2 * nI + 3
    Ab = np.zeros((J, K, m, 3))
    if K1:
        Ab[:, :K1, :2 * nI, :] = np.einsum('qpa,qkac->kqpc', Eb, Bi).reshape(K1, 2 * nI, 3)[None]
        Ab[:, :K1, 2 * nI:, :] = np.einsum('jab,kbc->jkac', Ct[:, :3, :], W)
    if K2:
        D2 = np.einsum('jab,qbc->jqac', Ct[:, :3, :], sinv(Cn)[:, :, 3:])
        D2 = D2 * (tc[None, :] < t_s[:, None])[:, :, None, None]
        Ab[:, K1:, 2 * nI:, :] = D2
    Tref = lab['T'][kw] if K1 else np.zeros((0, 3))
    rhs_int = (np.einsum('qpa,qa->qp', Eb, np.einsum('qkac,kc->qa', Bi, Tref) - miss_i).reshape(-1) if nI
               else np.zeros(0))
    Dref = np.einsum('jkac,kc->ja', Ab[:, :K1, 2 * nI:, :], Tref) if K1 else np.zeros((J, 3))
    Av = None
    if a < 0:
        Av = np.zeros((J, m, 3))
        Av[:, :2 * nI, :] = np.einsum('qpa,qab->qpb', Eb, Lvi).reshape(2 * nI, 3)[None]
        Av[:, 2 * nI:, :] = np.einsum('jab,bc->jac', Ct[:, :3, :], Pi[eL][:, 3:])
    return dict(a=a, t_a=t_a, k0=k0, kw=kw, kc=kc, tc=tc, ints=ints, t_s=t_s, r_ref=r_ref, Ab=Ab, rhs_int=rhs_int,
                Dref=Dref, Av=Av, Tref=Tref, cost_ref=float(np.linalg.norm(Tref, axis=1).sum()) if K1 else 0.0,
                Eb=Eb, Bi=Bi, Lvi=Lvi, miss_i=miss_i, vrel=vrel, m=m, K1=K1, K2=K2)


def irls_batch(Ab, AA, rhs, jc, tcap, eps_sched, Av=None, w_v=1.0):
    """Batched weighted-least-norm IRLS for min sum_k |T_k|  s.t.  sum_k A_{j,k} T_k (+ A_v dv) = rhs, per candidate.
    Returns best feasible (|T_k| <= tcap) iterate: T (Nc, K, 3), cost (Nc,), feasible (Nc,), dv (Nc, 3)."""
    Nc = len(jc); K = Ab.shape[1]; m = Ab.shape[2]
    if K == 0 and Av is None:
        return np.zeros((Nc, 0, 3)), np.full(Nc, np.inf), np.zeros(Nc, bool), np.zeros((Nc, 3))
    d = np.ones((Nc, K)); capw = np.ones((Nc, K))
    bestT = np.zeros((Nc, K, 3)); bestc = np.full(Nc, np.inf); bestv = np.zeros((Nc, 3))
    AAc = AA[jc]; Abc = Ab[jc]
    AvvT = None
    if Av is not None:
        AvvT = np.einsum('jma,jna->jmn', Av, Av)[jc] / w_v; Avc = Av[jc]
    I = np.eye(m) * 1e-12
    for eps in (None,) + tuple(eps_sched):
        Dinv = 1.0 / d
        G = np.einsum('ck,ckab->cab', Dinv, AAc)
        if AvvT is not None:
            G = G + AvvT
        sc = 1.0 / np.sqrt(np.maximum(np.einsum('cii->ci', G), 1e-300))
        Gs = G * sc[:, :, None] * sc[:, None, :]
        try:
            y = np.linalg.solve(Gs + I, (rhs * sc)[..., None])[..., 0] * sc
        except np.linalg.LinAlgError:
            y = np.einsum('cij,cj->ci', np.linalg.pinv(Gs + I), rhs * sc) * sc
        T = Dinv[:, :, None] * np.einsum('ckma,cm->cka', Abc, y)
        dv = np.einsum('cma,cm->ca', Avc, y) / w_v if Av is not None else np.zeros((Nc, 3))
        Tn = np.linalg.norm(T, axis=2); cost = Tn.sum(1)
        feas = Tn.max(1) <= tcap * 1.0005 if K else np.ones(Nc, bool)
        # rank check: a window with too few nodes cannot meet all rows (the least-norm formula then returns garbage)
        res = np.einsum('ckma,cka->cm', Abc, T) - rhs
        if Av is not None:
            res = res + np.einsum('cma,ca->cm', Avc, dv)
        feas &= np.linalg.norm(res, axis=1) <= np.maximum(1000.0, 1e-3 * np.linalg.norm(rhs, axis=1))
        bet = feas & (cost < bestc)
        bestT[bet] = T[bet]; bestc[bet] = cost[bet]; bestv[bet] = dv[bet]
        over = Tn > tcap
        capw = np.where(over, capw * 4.0, np.maximum(1.0, capw * 0.8))
        e = eps if eps is not None else eps_sched[0]
        d = capw / np.maximum(Tn, e)
    return bestT, bestc, np.isfinite(bestc), bestv


def price(lab, O, cand, tofs=None, keep_all=False):
    """Look-back prices of all candidate children of a label. cand: 1-based asteroid ids allowed as the next target.
    Returns a list of child dicts (score-ranked, n_per_target / max_children pruned unless keep_all)."""
    E = RI.eph()
    tofs = O.tofs if tofs is None else tofs
    if len(cand) == 0:
        return []
    W = window_setup(lab, O.L, tofs)
    J = len(W['t_s'])
    if J == 0:
        return []
    idx = np.array(sorted(cand), int) - 1
    R = E.ast_states_at(idx[None, :], W['t_s'][:, None])[0]            # (J, nX, 3)
    dist = np.linalg.norm(R - W['r_ref'][:, None, :], axis=2)
    jj, xx = np.nonzero(dist <= O.drmax)
    if len(jj) == 0:
        return []
    Ab = W['Ab']; AA = np.einsum('jkma,jkna->jkmn', Ab, Ab)
    nI2 = 2 * len(W['ints'])
    rhs = np.zeros((len(jj), W['m']))
    rhs[:, :nI2] = W['rhs_int'][None]
    rhs[:, nI2:] = W['Dref'][jj] - W['r_ref'][jj] + R[jj, xx]
    w_v = 1.0 if (W['a'] < 0 and np.linalg.norm(lab['vinf']) < VINF_MAX - 0.05) else 1e8
    out_T = []; out_c = []; out_f = []; out_v = []
    CH = 1500
    for s0 in range(0, len(jj), CH):
        sl = slice(s0, s0 + CH)
        T, c, f, v = irls_batch(Ab, AA, rhs[sl], jj[sl], tcap_of(lab, O), O.eps, W['Av'], w_v)
        out_T.append(T); out_c.append(c); out_f.append(f); out_v.append(v)
    T = np.concatenate(out_T); c = np.concatenate(out_c); f = np.concatenate(out_f); v = np.concatenate(out_v)
    marg = c - W['cost_ref']
    ok = f & (marg <= O.dvcap)
    if W['a'] < 0 and np.any(ok):
        nv = np.linalg.norm(lab['vinf'][None] + v, axis=1)
        ok &= nv <= VINF_MAX + 0.3                           # small overshoot: the refinement clips v_inf
    dv_new = lab['dv'] + marg
    ok &= (M_DRY + 1.5) * np.exp(np.maximum(dv_new, 0) / VE) <= O.tank_max
    if not np.any(ok):
        return []
    sel = np.nonzero(ok)[0]
    npz = sum(O.prize[x - 1] for x in lab['fa'])
    X1 = idx[xx[sel]] + 1
    sc = score_of(dv_new[sel], W['t_s'][jj[sel]], npz + O.prize[X1 - 1], O)
    order = np.argsort(sc, kind='stable')
    kept = []; per = {}
    for o in order:
        X = int(X1[o]); t = float(W['t_s'][jj[sel[o]]])
        lst = per.setdefault(X, [])
        if not keep_all:
            if len(lst) >= O.npt or any(abs(t - u) < O.sep for u in lst):
                continue
        lst.append(t); kept.append(o)
        if not keep_all and len(kept) >= O.mc:
            break
    res = []
    for o in kept:
        q = sel[o]; j = int(jj[q])
        Tq = T[q]; dvi = v[q]
        # epoch shifts of the window flybys from the along-track part of the predicted relative position
        dt_int = np.zeros(len(W['ints']))
        if len(W['ints']):
            dT = Tq[:W['K1']] - W['Tref']
            dr = W['miss_i'] + np.einsum('qkac,kc->qa', W['Bi'], dT) + np.einsum('qab,b->qa', W['Lvi'], dvi if W['a'] < 0 else np.zeros(3))
            dt_int = -np.einsum('qa,qa->q', dr, W['vrel']) / np.einsum('qa,qa->q', W['vrel'], W['vrel'])
        nK = W['K1'] + int(np.sum(W['tc'] < W['t_s'][j]))
        res.append(dict(X=int(X1[o]), t=float(W['t_s'][j]), marg=float(marg[q]), score=float(sc[o]),
                        T=Tq[:nK].copy(), dvi=dvi.copy(), dt_int=dt_int, a=W['a'], k0=W['k0'],
                        coast_au=float(dist[j, xx[q]] / AU)))
    return res


# ------------------------------------------------------------------ nonlinear refinement of a survivor
def refine(lab, ch, O):
    """Nonlinear correction of the child's window; returns the child label or None."""
    E = RI.eph()
    n = len(lab['fa']); a = ch['a']; tL = lab['tL']; k0 = ch['k0']
    ints = list(range(max(a, -1) + 1, n))
    t_s = ch['t']
    # window nodes: up to two grid nodes past t_s (zero impulse) so that a later final epoch still has node states
    kend = n_before(tL, min(t_s + 2 * DTAU, T_END + DTAU))
    nodes = np.arange(k0, kend); ts = node_times(tL, nodes)
    Ts = np.zeros((len(nodes), 3)); nT = min(len(ch['T']), len(nodes)); Ts[:nT] = ch['T'][:nT]
    tf = np.append(np.array([lab['ft'][i] for i in ints]) + np.clip(ch['dt_int'], -20 * DAY, 20 * DAY), t_s)
    asts = [lab['fa'][i] for i in ints] + [ch['X']]
    t_lo = tL + 1 * DAY if a < 0 else float(lab['ft'][a]) + 1 * DAY

    def fix_tf(tf_):
        tf_ = np.minimum(np.asarray(tf_, float), T_END)
        lo = t_lo
        for q in range(len(tf_)):
            tf_[q] = max(tf_[q], lo); lo = tf_[q] + 0.5 * DAY
        return tf_
    tf = fix_tf(tf)
    if a < 0:
        vi = lab['vinf'] + ch['dvi']; nv = np.linalg.norm(vi)
        if nv > VINF_MAX - 1e-6: vi = vi * (VINF_MAX - 1e-6) / nv
        wp = ImpulsiveProblem(E, tL, vi, ts, Ts, tf, asts); free_v = True
    else:
        wp = WindowProblem(E, lab['ft'][a], lab['Xf'][a], ts, Ts, tf, asts, O.tcap); free_v = False
    tcap_r = tcap_of(lab, O)
    Yf, Ys = wp.integrate(); dm, vr = wp.misses(Yf); d = np.linalg.norm(dm, axis=1)
    d0 = float(d.max()); cst = dict(Yf=Yf, Ys=Ys, dm=dm, vr=vr, d=d)

    def step(rho, merit):
        """One linearised step (irls_min_fuel, proximal weight rho) with a line search; merit(d, Ts) decides."""
        B, Lv = wp.linearise(cst['Yf'], cst['Ys'])
        wv = (1.0 if np.linalg.norm(wp.vinf) < VINF_MAX - 0.05 else 1e8) if free_v else 1e8
        T2, dt2, dv2 = irls_min_fuel(B, Lv, cst['vr'], cst['dm'], wp.Ts, wp.vinf, tcap_r, rho=rho,
                                     eps_sched=(0.02, 0.01, 0.005), inner=2, w_v=wv)
        m_cur = merit(cst['d'], wp.Ts)
        for al in (1.0, 0.5, 0.25):
            Ts_a = wp.Ts + al * (T2 - wp.Ts)
            tf_a = fix_tf(wp.tf + al * np.clip(dt2, -20, 20) * DAY)
            vi_a = wp.vinf + al * dv2 if free_v else wp.vinf
            nv = np.linalg.norm(vi_a)
            if nv > VINF_MAX - 1e-6: vi_a = vi_a * (VINF_MAX - 1e-6) / nv
            Yf2, Ys2 = wp.integrate(Ts=Ts_a, tf=tf_a, vinf=vi_a); dm2, vr2 = wp.misses(Yf2, tf_a)
            d2 = np.linalg.norm(dm2, axis=1)
            if np.all(np.isfinite(d2)) and merit(d2, Ts_a) < m_cur:
                wp.Ts, wp.tf, wp.vinf = Ts_a, tf_a, vi_a
                cst.update(Yf=Yf2, Ys=Ys2, dm=dm2, vr=vr2, d=d2); return True
        return False

    def restore_(n_it):
        k = 0
        while cst['d'].max() > O.ref_tol and k < n_it:
            k += 1
            if not step(1e5, lambda d_, T_: float(d_.max())):         # Newton: min-norm correction of the misses
                break
        return k
    it = restore_(O.ref_iters)
    if not np.all(np.isfinite(cst['d'])) or cst['d'].max() > O.acc_tol:
        return None, dict(ok=False, miss0=d0, miss=float(cst['d'].max()) if np.all(np.isfinite(cst['d'])) else None, it=it)
    dv_lin = float(np.linalg.norm(wp.Ts, axis=1).sum())
    pen = lambda d_: float(np.maximum(0.0, d_ - O.ref_tol).sum()) / 2e4
    for q in range(O.opt_iters):                                    # local SCP polish of the window, then restore
        if not step(O.ref_rho, lambda d_, T_: float(np.linalg.norm(T_, axis=1).sum()) + pen(d_)):
            break
        it += restore_(2)
    d = cst['d']; Yf = cst['Yf']; Ys = cst['Ys']
    if not np.all(np.isfinite(d)) or d.max() > O.acc_tol:
        return None, dict(ok=False, miss0=d0, miss=float(d.max()) if np.all(np.isfinite(d)) else None, it=it, stage='opt')
    t_new = float(wp.tf[-1])
    kend2 = n_before(tL, t_new)                                  # nodes before the (final) last flyby
    nw = kend2 - k0
    if nw < 0 or nw > len(nodes):
        return None, dict(ok=False, reason='nodes')
    T_child = np.concatenate([lab['T'][:k0], wp.Ts[:nw]])
    Xn_child = np.concatenate([lab['Xn'][:k0], Ys[:nw, :6]])
    na = a + 1 if a >= 0 else 0
    child = make_label(tL, wp.vinf if free_v else lab['vinf'], T_child, Xn_child,
                       lab['fa'][:na] + asts, np.concatenate([lab['ft'][:na], wp.tf]),
                       np.concatenate([lab['Xf'][:na], Yf[:, :6]]), np.concatenate([lab['miss'][:na], d]))
    return child, dict(ok=True, miss0=d0, miss=float(d.max()), it=it, marg=ch['marg'],
                       marg_ref=child['dv'] - lab['dv'], gain_opt=round(dv_lin - float(np.linalg.norm(wp.Ts, axis=1).sum()), 4))


def polish(lab, O):
    """Longer-window SCP of a kept label: the last --polish-L legs (from flyby n - polish_L, or the launch) are
    re-optimised by --polish-iters min-fuel steps (merit = window L1 cost + miss penalty), then restored.  Returns the
    polished label (or the input label if nothing improved)."""
    n = len(lab['fa']); a = n - O.pol_L
    ch = dict(a=a, k0=0 if a < 0 else n_before(lab['tL'], lab['ft'][a]), t=float(lab['ft'][-1]), X=lab['fa'][-1],
              dt_int=np.zeros(n - 1 - max(a, -1) - 1 + 0), dvi=np.zeros(3), marg=0.0)
    k0 = ch['k0']; M = len(lab['T'])
    ch['T'] = lab['T'][k0:M]
    # refine() treats the LAST flyby as the new target and the others as window flybys: exactly a polish of the window
    base = dict(lab); base['fa'] = lab['fa'][:-1]; base['ft'] = lab['ft'][:-1]; base['Xf'] = lab['Xf'][:-1]
    base['miss'] = lab['miss'][:-1]
    O2 = Opt.__new__(Opt); O2.__dict__.update(O.__dict__); O2.opt_iters = O.pol_iters
    new, info = refine(base, ch, O2)
    if new is None or new['dv'] >= lab['dv'] - 1e-6:
        return lab, 0.0
    return new, lab['dv'] - new['dv']


# ------------------------------------------------------------------ final settle
def _cap_ok(ip):
    return (not len(ip.Ts)) or float(np.linalg.norm(ip.Ts, axis=1).max()) <= 1.02 * ip.tcap


def _feasible(ip, d):
    return len(d) == 0 or (np.all(np.isfinite(d)) and float(np.max(d)) <= TOL and _cap_ok(ip))


def settle_label(lab, iters, timeout):
    """Robust whole-route settle of a planned label (a planned label is already a twin trajectory with window misses
    <= --acc-tol km):  A = restore only (feasible fallback near the planned tank; thrust-cap homotopy if a node is
    over the twin's cap at that tank);  B = run_ialns.settle (restore + SCP + restore + cap), with one extra restore
    if it ends above 150 km.  The lightest feasible of A and B is returned."""
    from ctoc14.globalopt import restore
    from ctoc14.impulsive import enforce_cap
    QUIET = lambda *x, **k: None
    st = label_ist(lab); tic = time.time(); best = None; info = dict(ok=False)
    try:
        signal.signal(signal.SIGALRM, _timeout); signal.setitimer(signal.ITIMER_REAL, timeout)
        ipA = RI.ipr(st)
        d = restore(ipA, tol_km=100, iters=30, rho=1e3, log=QUIET)
        if not _cap_ok(ipA):
            enforce_cap(ipA, stages=6); d = restore(ipA, tol_km=100, iters=30, log=QUIET)
        if _feasible(ipA, d):
            best = ipA; info.update(A_tank=round(float(ipA.tank()), 2))
        ipB = RI.ipr(st)
        miss = float(RI.settle(ipB, iters))
        if miss > TOL:
            d = restore(ipB, tol_km=100, iters=40, rho=1e3, log=QUIET); miss = float(np.max(d)) if len(d) else 0.0
        Yf, _ = ipB.integrate(); dm, _ = ipB.misses(Yf); dB = np.linalg.norm(dm, axis=1)
        info.update(B_tank=round(float(ipB.tank()), 2), B_miss=round(float(dB.max()), 1))
        if _feasible(ipB, dB) and (best is None or ipB.tank() < best.tank()):
            best = ipB
    except Exception as e:
        info.update(err=type(e).__name__)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    info['dt_s'] = round(time.time() - tic, 1)
    if best is None:
        return None, info
    Yf, _ = best.integrate(); dm, _ = best.misses(Yf)
    info.update(ok=True, miss=round(float(np.linalg.norm(dm, axis=1).max()), 1), tank=round(float(best.tank()), 2),
                used='B' if best is not None and info.get('B_tank') == round(float(best.tank()), 2) else 'A')
    return best, info


def save_ckpt(out, obj):
    tmp = out / 'ckpt.pkl.tmp'
    with open(tmp, 'wb') as f:
        pickle.dump(obj, f)
    os.replace(tmp, out / 'ckpt.pkl')


def _pool_map(fn, jobs, nproc):
    if nproc <= 1 or len(jobs) <= 1:
        return [fn(j) for j in jobs]
    import multiprocessing as mp
    with mp.get_context('fork').Pool(min(nproc, len(jobs))) as pool:
        return pool.map(fn, jobs, chunksize=max(1, len(jobs) // (4 * nproc)))


_G = {}


def _price_job(i):
    lab = _G['beam'][i]
    rem = [x for x in _G['pool'] if x not in set(lab['fa'])]
    try:
        return price(lab, _G['O'], rem)
    except Exception as e:
        return []


def _polish_job(i):
    try:
        return polish(_G['nb'][i], _G['O'])
    except Exception as e:
        return _G['nb'][i], 0.0


def _refine_job(q):
    i, ch = _G['todo'][q]
    try:
        return refine(_G['beam'][i], ch, _G['O'])
    except Exception as e:
        return None, dict(ok=False, err=repr(e)[:80])


# ------------------------------------------------------------------ free beam
def run_free(a, pool, O, say, out):
    E = RI.eph()
    ck = out / 'ckpt.pkl'
    pset = set(pool)
    if ck.exists():
        S = pickle.load(open(ck, 'rb'))
        say(f'RESUME at level {S["level"]}: beam {len(S["beam"])}, best {len(S["best"]["fa"])} @ {label_tank(S["best"]):.1f} kg')
    else:
        P = BP(beam=a.beam, w_fuel=a.w_fuel, m_margin=40.0, dv_max=1e9, tofs=np.arange(15, 400 + 1e-9, 5) * DAY,
               vinf_cap=4.0, tof_refine=False)
        P.prize = O.prize
        g = [float(x) for x in a.grid.split(',')]; grid = np.arange(g[0], g[1] + 1e-9, g[2]) * DAY
        excl = sorted(set(ALL) - pset)
        rts = search_roots(E, grid, M0P, excl, P)
        rts.sort(key=lambda s: s.score(P)); seen = set(); b2 = []
        for s in rts:
            key = (s.seq[0][0], int(s.t_launch // (20 * DAY)))
            if key in seen: continue
            seen.add(key); b2.append(s)
        b2 = b2[:a.roots]
        beam = [root_label(s.t_launch, s.vinf, s.seq[0][0], s.seq[0][1]) for s in b2]
        S = dict(level=1, beam=beam, best=None, pareto={}, hist=[], wall=0.0, end=None, n_ref=0, n_ref_ok=0)
        for b in beam:
            S['best'] = better(S['best'], b); add_pareto(S['pareto'], b, a.pareto_k)
        say(f'roots: {len(rts)} launch legs -> {len(beam)} kept (max root miss {max(b["miss"].max() for b in beam):.1f} km)')
        save_ckpt(out, S)
    while S['end'] is None:
        tic = time.time(); beam = S['beam']
        _G.update(beam=beam, pool=pool, O=O)
        res = _pool_map(_price_job, list(range(len(beam))), a.nproc)
        t_price = time.time() - tic
        kids = []
        for i, lst in enumerate(res):
            vm = mask_of(beam[i]['fa'])
            for ch in lst:
                kids.append((ch['score'], i, vm | (1 << (ch['X'] - 1)), ch))
        nk = len(kids)
        kids.sort(key=lambda z: z[0])
        seen = set(); cand = []
        for kd in kids:
            key = (kd[2], int(kd[3]['t'] // (5 * DAY)))
            if key in seen: continue
            seen.add(key); cand.append(kd)
        nu = len(cand)
        todo = [(kd[1], kd[3]) for kd in cand[:int(np.ceil(a.beam * a.over))]]
        _G['todo'] = todo
        rr = _pool_map(_refine_job, list(range(len(todo))), a.nproc)
        t_ref = time.time() - tic - t_price
        new = []; infos = []
        for (i, ch), (lab, info) in zip(todo, rr):
            S['n_ref'] += 1; infos.append(info)
            if lab is not None:
                S['n_ref_ok'] += 1; new.append(lab)
        new.sort(key=lambda l: lab_score(l, O))
        seen = set(); nb = []
        for l in new:
            key = (mask_of(l['fa']), int(l['ft'][-1] // (5 * DAY)))
            if key in seen: continue
            seen.add(key); nb.append(l)
        nb = nb[:a.beam]
        gain = 0.0; t_pol = 0.0
        if O.pol_L > 0 and O.pol_iters > 0 and nb:
            t2 = time.time(); _G['nb'] = nb
            pr = _pool_map(_polish_job, list(range(len(nb))), a.nproc)
            nb = [l for l, g in pr]; gain = float(np.mean([g for l, g in pr])); t_pol = time.time() - t2
        depth = S['level'] + 1
        okr = [i for i in infos if i.get('ok')]
        rec = dict(depth=depth, parents=len(beam), children=nk, unique=nu, refined=len(todo), ok=len(okr), kept=len(nb),
                   t_price_s=round(t_price, 1), t_ref_s=round(t_ref, 1), t_pol_s=round(t_pol, 1), pol_gain=round(gain, 4),
                   marg_lin_vs_ref=[(round(i['marg'], 4), round(i['marg_ref'], 4)) for i in okr[:30]],
                   ref_miss0_med=float(np.median([i['miss0'] for i in okr])) if okr else None,
                   fails=[i for i in infos if not i.get('ok')][:8])
        if nb:
            rec.update(t_end_d=[round(min(l['ft'][-1] for l in nb) / DAY, 1), round(max(l['ft'][-1] for l in nb) / DAY, 1)],
                       tank=[round(min(label_tank(l) for l in nb), 1), round(max(label_tank(l) for l in nb), 1)])
            for l in nb:
                S['best'] = better(S['best'], l); add_pareto(S['pareto'], l, a.pareto_k)
        S['hist'].append(rec)
        say(f'depth {depth:2d}: {len(beam)} parents -> {nk} children ({nu} unique) -> refined {len(todo)}, ok {len(okr)}, '
            f'kept {len(nb)}' + (f'; t_end {rec["t_end_d"][0]:.0f}-{rec["t_end_d"][1]:.0f} d, tank '
                                 f'{rec["tank"][0]:.0f}-{rec["tank"][1]:.0f}' if nb else '') +
            f'; best {len(S["best"]["fa"])} @ {label_tank(S["best"]):.1f} kg (price {t_price:.1f} s, refine {t_ref:.1f} s)')
        if not nb:
            S['end'] = dict(kind='no_children' if nk == 0 else 'refine_wall', depth=S['level'],
                            time_left_d=round((T_MISSION - max(l['ft'][-1] for l in beam)) / DAY, 1))
        elif len(nb[0]['fa']) >= len(pool) or depth >= a.max_depth:
            S['end'] = dict(kind='pool_exhausted' if len(nb[0]['fa']) >= len(pool) else 'max_depth', depth=depth)
        S['beam'] = nb if nb else beam; S['level'] = depth if nb else S['level']
        S['wall'] += time.time() - tic
        save_ckpt(out, S)
    return S


def better(a_, b_):
    """Deeper wins; equal depth: lighter."""
    if a_ is None: return b_
    if b_ is None: return a_
    return a_ if (len(a_['fa']), -a_['dv']) >= (len(b_['fa']), -b_['dv']) else b_


def add_pareto(par, lab, k):
    n = len(lab['fa']); lst = par.setdefault(n, [])
    key = frozenset(lab['fa'])
    for q, l in enumerate(lst):
        if frozenset(l['fa']) == key:
            if lab['dv'] < l['dv']: lst[q] = lab
            break
    else:
        lst.append(lab)
    lst.sort(key=lambda l: l['dv']); del lst[k:]


def final_settle(S, a, say, out, name, pool):
    """Settle the Pareto labels deepest first (--settle-per per depth, all of them settled); the lightest settled route
    of the deepest depth that settles at all is the route."""
    tries = []; got = None
    depths = sorted(S['pareto'], reverse=True)
    for n in depths[:a.settle_depths]:
        for q, lab in enumerate(S['pareto'][n][:a.settle_per]):
            ip, info = settle_label(lab, a.iters, a.timeout)
            info.update(n=n, alt=q, planned_tank=round(float(label_tank(lab)), 2), planned_miss=round(float(lab['miss'].max()), 1))
            tries.append(info)
            say(f'  settle depth {n} alt {q}: planned {info["planned_tank"]:.1f} kg (max miss {info["planned_miss"]:.0f} km) -> '
                + (f'twin {info["tank"]:.1f} kg, miss {info["miss"]:.0f} km [{info.get("used")}; A {info.get("A_tank")} '
                   f'B {info.get("B_tank")}/{info.get("B_miss")} km]' if info.get('ok') else f'FAILED {info}')
                + f' ({info["dt_s"]:.0f} s)')
            if ip is not None and (got is None or ip.tank() < got[1].tank()):
                got = (lab, ip, info)
        if got is not None:
            break
    if got is not None:
        lab, ip, info = got
        np.savez(out / f'route_{name}.npz', **RI.ist(ip))
    return got, tries


# ------------------------------------------------------------------ guided calibration
def run_guided(a, src, O, say, out):
    """Along the source route's order: price the true next leg from the planner's own label (source epoch and the best
    grid epoch within +-15 d), refine, append.  If the look-back cannot add a leg, the chain goes on from a
    whole-prefix settle of the leg (s14_twinbeam.settle_child) and the leg is counted as a failure."""
    E = RI.eph()
    o = np.argsort(src['tf']); asts = [int(src['asts'][i]) for i in o]; tfs = [float(src['tf'][i]) for i in o]
    ck = out / 'ckpt.pkl'
    if ck.exists():
        S = pickle.load(open(ck, 'rb')); say(f'RESUME guided at step {S["k"]}')
    else:
        tL = float(src['tL']); rE, vE = E.earth_state(tL)
        R1 = E.ast_states_at(np.array([asts[0] - 1]), np.array([tfs[0]]))[0][0]
        v1, _ = lambert(rE[None], R1[None], np.array([tfs[0] - tL]))
        vinf = v1[0] - vE; nv = np.linalg.norm(vinf)
        lab = root_label(tL, vinf * min(1.0, (VINF_MAX - 1e-6) / nv), asts[0], tfs[0])
        S = dict(k=1, cur=lab, steps=[dict(k=0, ast=asts[0], vinf=round(float(nv), 3), miss=round(float(lab['miss'].max()), 1))],
                 wall=0.0, fails=[])
        save_ckpt(out, S)
    while S['k'] < len(asts):
        tic = time.time(); k = S['k']; lab = S['cur']; X = asts[k]; t_src = tfs[k]
        rec = dict(k=k, ast=X, t_src_d=round(t_src / DAY, 1), t_prev_d=round(float(lab['ft'][-1]) / DAY, 1))
        tof0 = t_src - float(lab['ft'][-1])
        grid = tof0 + np.arange(-15.0, 15.0 + 1e-9, 2.5) * DAY
        grid = grid[(grid >= 5 * DAY) & (float(lab['ft'][-1]) + grid <= T_END)]
        O2 = Opt.__new__(Opt); O2.__dict__.update(O.__dict__); O2.dvcap = 50.0; O2.drmax = 5.0 * AU
        chs = price(lab, O2, [X], tofs=grid, keep_all=True) if len(grid) else []
        src_ch = min(chs, key=lambda c: abs(c['t'] - t_src)) if chs else None
        best_ch = min(chs, key=lambda c: c['marg']) if chs else None
        rec['price_src'] = round(src_ch['marg'], 4) if src_ch else None
        rec['coast_au'] = round(src_ch['coast_au'], 4) if src_ch else None
        rec['price_best'] = round(best_ch['marg'], 4) if best_ch else None
        rec['best_dt_d'] = round((best_ch['t'] - t_src) / DAY, 1) if best_ch else None
        got = None; tries = []
        for tag, ch in (('best', best_ch), ('src', src_ch)):
            if ch is None or (got is not None): continue
            if tag == 'src' and best_ch is not None and abs(ch['t'] - best_ch['t']) < 0.5 * DAY: continue
            child, info = refine(lab, ch, O)
            tries.append(dict(at=tag, **{kk: (round(v, 4) if isinstance(v, float) else v) for kk, v in info.items()}))
            if child is not None:
                got = child; rec['used'] = tag; rec['marg_ref'] = round(child['dv'] - lab['dv'], 4)
        rec['refine'] = tries
        if got is None:
            from s14_twinbeam import settle_child, pack
            par = dict(st=label_ist(lab), nreg=len(lab['T']), dv=lab['dv'], t_end=float(lab['ft'][-1]),
                       r_end=lab['Xf'][-1][:3], v_end=lab['Xf'][-1][3:], asts=lab['fa'], tfs=list(lab['ft']))
            tw, info = settle_child((par, X, t_src, 100, 240.0))
            rec['fallback'] = dict(ok=info.get('ok'), tank=info.get('tank'))
            S['fails'].append(k)
            if tw is None:                      # skip the target (as s14_twinbeam guided does) and go on
                rec['used'] = 'skipped'; S['skipped'] = S.get('skipped', []) + [X]
                S['steps'].append(rec); S['k'] = k + 1; S['wall'] += time.time() - tic
                say(f'leg {k:2d} -> {X:3d}: look-back and fallback settle both failed -- skipped'); save_ckpt(out, S)
                continue
            ip = RI.ipr(tw['st']); ip = regrid(ip)
            Yf, Ys = ip.integrate(); dm, _ = ip.misses(Yf)
            got = make_label(ip.tL, ip.vinf, ip.Ts, Ys[:, :6], ip.asts, ip.tf, Yf[:, :6], np.linalg.norm(dm, axis=1))
            rec['used'] = 'fallback'; rec['marg_ref'] = round(got['dv'] - lab['dv'], 4)
        S['cur'] = got; rec['tank'] = round(float(label_tank(got)), 2); rec['miss_max'] = round(float(got['miss'].max()), 1)
        rec['wall_s'] = round(time.time() - tic, 2)
        S['steps'].append(rec); S['k'] = k + 1; S['wall'] += time.time() - tic
        say(f'leg {k:2d} -> {X:3d} (tof {tof0 / DAY:5.1f} d): L{O.L} price src {rec["price_src"]} best {rec["price_best"]} '
            f'(dt {rec["best_dt_d"]}) coast {rec["coast_au"]} AU | used {rec.get("used")} marg {rec.get("marg_ref")} | '
            f'tank {rec["tank"]:.1f} miss {rec["miss_max"]:.0f} km ({rec["wall_s"]:.2f} s)')
        save_ckpt(out, S)
    return S


def regrid(ip):
    """Put a settled twin on the standard node grid (tL + 10 d + 20 d k): impulses summed per 20-d bin (for the
    guided fallback only; s14_twinbeam children carry extra seed nodes)."""
    tL = ip.tL; t_last = float(np.max(ip.tf)); M = n_before(tL, t_last)
    tsg = node_times(tL, np.arange(M))
    b = np.clip(np.round((ip.ts - tL - OFF) / DTAU).astype(int), 0, max(M - 1, 0))
    Tg = np.zeros((M, 3))
    keep = ip.ts < t_last
    np.add.at(Tg, b[keep], ip.Ts[keep])
    o = np.argsort(ip.tf)
    ip2 = ImpulsiveProblem(RI.eph(), tL, ip.vinf, tsg, Tg, np.asarray(ip.tf)[o], [ip.asts[i] for i in o])
    RI.settle(ip2, 40)
    return ip2


# ------------------------------------------------------------------ main
def parse_pool(spec):
    """own:SRC_FLEET:ROUTE | pools:FILE:ROUTE | all | list:1,2,3 -> (sorted pool, name, source info)"""
    if spec == 'all':
        return list(ALL), 'all', {}
    kind, rest = spec.split(':', 1)
    if kind == 'own':
        src, r = rest.rsplit(':', 1); F = RI.IFleet(src)
        st = F.routes[r]['st']
        return sorted(int(x) for x in st['asts']), r, dict(src=src, route=r, src_tank=F.routes[r]['tank'], src_st=st)
    if kind == 'pools':
        f, r = rest.rsplit(':', 1)
        return sorted(int(x) for x in json.load(open(f))['pools'][r]), r, dict(pools=f, route=r)
    if kind == 'list':
        return sorted(int(x) for x in rest.split(',')), 'list', {}
    raise ValueError(spec)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['free', 'guided', 'resettle']); ap.add_argument('out')
    ap.add_argument('--pool', default='all'); ap.add_argument('--src'); ap.add_argument('--route')
    ap.add_argument('--name', default=None)
    ap.add_argument('--L', type=int, default=2, help='look-back window in legs (1 = single-leg twin-grid linear model)')
    ap.add_argument('--beam', type=int, default=100); ap.add_argument('--over', type=float, default=1.5,
                                                                      help='refinements per level = over x beam')
    ap.add_argument('--roots', type=int, default=300); ap.add_argument('--grid', default='0,800,20')
    ap.add_argument('--dvcap', type=float, default=3.0, help='max marginal price of a leg [km/s]')
    ap.add_argument('--drmax', type=float, default=1.0, help='max coast miss of a candidate [AU]')
    ap.add_argument('--mcap', type=float, default=1100.0, help='mass [kg] of the node thrust cap 0.43 N x 20 d / m')
    ap.add_argument('--tank-max', type=float, default=1400.0)
    ap.add_argument('--tof-min', type=float, default=10.0); ap.add_argument('--tof-max', type=float, default=600.0)
    ap.add_argument('--tof-step1', type=float, default=5.0); ap.add_argument('--tof-step2', type=float, default=10.0)
    ap.add_argument('--npt', type=int, default=2); ap.add_argument('--mc', type=int, default=80)
    ap.add_argument('--sep', type=float, default=15.0, help='min epoch separation [d] of two children of one target')
    ap.add_argument('--w-fuel', type=float, default=round(CM.w_fuel(480.0, M0P), 4))
    ap.add_argument('--eps', default='0.05,0.02,0.01,0.005,0.003,0.002')
    ap.add_argument('--ref-iters', type=int, default=6); ap.add_argument('--ref-tol', type=float, default=100.0)
    ap.add_argument('--acc-tol', type=float, default=1000.0); ap.add_argument('--ref-rho', type=float, default=30.0)
    ap.add_argument('--opt-iters', type=int, default=2, help='local SCP polish steps of a survivor window')
    ap.add_argument('--cap-frac', type=float, default=0.0,
                    help='> 0: node cap at max(mcap, cap_frac x planned tank) (twin-consistent heavy routes)')
    ap.add_argument('--polish-L', type=int, default=0, help='> 0: SCP polish of the kept labels over their last legs')
    ap.add_argument('--polish-iters', type=int, default=3)
    ap.add_argument('--max-depth', type=int, default=80); ap.add_argument('--nproc', type=int, default=1)
    ap.add_argument('--pareto-k', type=int, default=3); ap.add_argument('--settle-depths', type=int, default=6)
    ap.add_argument('--settle-per', type=int, default=3); ap.add_argument('--no-settle', action='store_true')
    ap.add_argument('--iters', type=int, default=100); ap.add_argument('--timeout', type=float, default=900.0)
    a = ap.parse_args()
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt')
    t0 = time.time()
    if a.mode == 'guided':
        F = RI.IFleet(a.src); src = F.routes[a.route]['st']; tank0 = F.routes[a.route]['tank']
        pool = sorted(int(x) for x in src['asts'])
        prize = np.zeros(300); prize[np.array(pool) - 1] = 1.0
        O = Opt(a, prize)
        say(f'=== s14_lookback guided route {a.route}: {len(pool)} fb @ {tank0:.1f} kg; L {a.L} mcap {a.mcap} '
            f'(tcap {O.tcap:.3f} km/s/node) eps {a.eps} ref {a.ref_iters}/{a.ref_tol}/{a.acc_tol}/rho {a.ref_rho}')
        S = run_guided(a, src, O, say, out)
        st = [r for r in S['steps'][1:] if 'price_src' in r]
        lab = S['cur']
        res = dict(mode='guided', route=a.route, L=a.L, legs=len(st), src_tank=round(tank0, 2),
                   fails=S['fails'], skipped=S.get('skipped', []), planned_n=len(lab['fa']), planned_tank=round(float(label_tank(lab)), 2),
                   planned_miss=round(float(lab['miss'].max()), 1), wall_s=round(S['wall'], 1), steps=S['steps'])
        for cap in (1.2, 2.0, 3.0):
            res[f'adm_src_{cap}'] = sum(1 for r in st if r['price_src'] is not None and r['price_src'] <= cap)
            res[f'adm_best_{cap}'] = sum(1 for r in st if r['price_best'] is not None and r['price_best'] <= cap)
        if not a.no_settle:
            ip, info = settle_label(lab, a.iters, a.timeout); res['settle'] = info
            if ip is not None:
                np.savez(out / f'route_{a.route}.npz', **RI.ist(ip))
        np.savez(out / 'planned.npz', **label_ist(lab))
        say(f'RESULT guided {a.route} L{a.L}: {len(st)} legs; admissible at the source epoch <=1.2/2.0/3.0 km/s: '
            f'{res["adm_src_1.2"]}/{res["adm_src_2.0"]}/{res["adm_src_3.0"]}; best epoch +-15 d: {res["adm_best_1.2"]}/'
            f'{res["adm_best_2.0"]}/{res["adm_best_3.0"]}; refine fails {len(S["fails"])}; planned {res["planned_n"]} fb @ '
            f'{res["planned_tank"]:.1f} kg; settled {res.get("settle")}; source {tank0:.1f} kg; wall {S["wall"]:.0f} s')
        json.dump(res, open(out / 'result.json', 'w'), indent=1, default=float)
        return
    if a.mode == 'resettle':
        res = json.load(open(out / 'result.json')); S = pickle.load(open(out / 'ckpt.pkl', 'rb'))
        pool, name, info = parse_pool(res['pool_spec']); name = a.name or res['name']
        for f in out.glob('route_*.npz'):
            f.unlink()
        say(f'=== resettle {name}: Pareto depths {sorted(S["pareto"])[-3:]}')
        t1 = time.time(); got, tries = final_settle(S, a, say, out, name, pool)
        res['settle_tries'] = tries; res['settle_wall_s'] = round(time.time() - t1, 1); res.pop('settled', None)
        if got is not None:
            lab, ip, inf = got; tg = sorted(int(x) for x in ip.asts)
            res['settled'] = dict(n=len(tg), tank=round(float(ip.tank()), 2), J_i=round(float(RI.cost(ip.tank())), 5),
                                  kg_per_fb=round((float(ip.tank()) - 600.0) / len(tg), 3), miss_km=inf['miss'],
                                  planned_tank=inf['planned_tank'], targets=tg, missing=sorted(set(pool) - set(tg)))
            say(f'RESULT {name}: settled {len(tg)}/{len(pool)} @ {ip.tank():.1f} kg ({res["settled"]["kg_per_fb"]:.2f} kg/fb)')
        res['wall_s'] = round(res.get('plan_wall_s', 0) + res['settle_wall_s'], 1)
        json.dump(res, open(out / 'result.json', 'w'), indent=1, default=float)
        return
    pool, name, info = parse_pool(a.pool)
    name = a.name or name
    prize = np.zeros(300); prize[np.array(pool) - 1] = 1.0
    O = Opt(a, prize)
    say(f'=== s14_lookback free {name}: pool {len(pool)}; L {a.L} beam {a.beam} over {a.over} roots {a.roots} dvcap {a.dvcap} '
        f'drmax {a.drmax} mcap {a.mcap} (tcap {O.tcap:.3f}) tofs {len(O.tofs)} npt {a.npt} mc {a.mc} w_fuel {a.w_fuel} '
        f'ref {a.ref_iters}/{a.ref_tol}/{a.acc_tol}/rho {a.ref_rho}' + (f'; source {info.get("src_tank"):.1f} kg' if info.get('src_tank') else ''))
    S = run_free(a, pool, O, say, out)
    b = S['best']
    np.savez(out / 'planned.npz', **label_ist(b))
    res = dict(mode='free', name=name, pool_spec=a.pool, pool_size=len(pool), L=a.L, beam=a.beam, dvcap=a.dvcap,
               drmax=a.drmax, mcap=a.mcap, planned_n=len(b['fa']), planned_tank=round(float(label_tank(b)), 2),
               planned_miss=round(float(b['miss'].max()), 1), order=b['fa'], epochs_d=[round(t / DAY, 1) for t in b['ft']],
               t_launch_d=round(b['tL'] / DAY, 1), end=S['end'], refines=S['n_ref'], refine_ok=S['n_ref_ok'],
               plan_wall_s=round(S['wall'], 1), hist=S['hist'],
               pareto={int(n): [round(float(label_tank(l)), 1) for l in v] for n, v in sorted(S['pareto'].items())})
    if info.get('src_tank'):
        res['src_tank'] = round(info['src_tank'], 2)
    json.dump(res, open(out / 'result.json', 'w'), indent=1, default=float)
    say(f'PLANNED {name}: {len(b["fa"])}/{len(pool)} @ {label_tank(b):.1f} kg (max miss {b["miss"].max():.0f} km); end {S["end"]}; '
        f'refines {S["n_ref"]} ({S["n_ref_ok"]} ok); wall {S["wall"]:.0f} s')
    if not a.no_settle:
        t1 = time.time()
        got, tries = final_settle(S, a, say, out, name, pool)
        res['settle_tries'] = tries; res['settle_wall_s'] = round(time.time() - t1, 1)
        if got is not None:
            lab, ip, inf = got
            tg = sorted(int(x) for x in ip.asts)
            res['settled'] = dict(n=len(tg), tank=round(float(ip.tank()), 2), J_i=round(float(RI.cost(ip.tank())), 5),
                                  kg_per_fb=round((float(ip.tank()) - 600.0) / len(tg), 3), miss_km=inf['miss'],
                                  planned_tank=inf['planned_tank'], targets=tg, missing=sorted(set(pool) - set(tg)))
            say(f'RESULT {name}: settled {len(tg)}/{len(pool)} @ {ip.tank():.1f} kg ({res["settled"]["kg_per_fb"]:.2f} kg/fb, '
                f'J_i {res["settled"]["J_i"]:.4f}); planned {inf["planned_tank"]:.1f} kg; total wall {time.time() - t0:.0f} s')
        else:
            say(f'RESULT {name}: nothing settled')
    res['wall_s'] = round(time.time() - t0, 1)
    json.dump(res, open(out / 'result.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
