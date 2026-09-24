"""Plan-level levers: coordinate-descent re-timing of planned tours + tank right-sizing (hybrid planner model).

The 12 planned tours of results/phasing/camp{1,2}/milp_n12/tour_sc*.json are HYBRID chains of the beam search
(search.expand): some legs are single-revolution Lambert arcs (junction dv, kappa leverage rule), the others come from the
linearised continuous-thrust model (linleg.LinLeg: cost = L1 thrust integral, no kappa, arrival velocity from Phi_vv).
insertion.build_tour (Lambert only) does NOT reproduce them (reconstruction dv of 2-55 km/s on the linear legs), so this
script re-implements the planner's chain exactly: every leg's model is identified by re-evaluating both models at the
search m0 (camp1 1000 kg, camp2 1100 kg; the pool's fuel_est was computed there, the stored m0 = 600 + 1.3 fuel_est + 20
is the tank the MILP assigned afterwards) and matching the stored dv_est.

Refinement (objective = total planned propellant at the STORED m0, the tank that is actually flown): each flyby time in
turn +-10 d at 1 d then +-1 d at 0.25 d; launch epoch +-30 d at 2.5 d with the free-v_inf rule (|v_inf| <= 4 km/s, excess
paid by thrust); passes repeated until a pass saves < 0.5 kg. Constraints per leg: utilisation <= 1 (Lambert:
dv/(eta a tof); linear: max|u|/(ub a)) -- legs already above 1 in the base may not get worse; on re-evaluated legs
dv <= 1.2 km/s (legs above the cap may not get worse) and linear coast miss |dR| <= lin_drmax; tof >= 15 d (Lambert) /
20 d (linear, <= 600 d); flybys <= T_MISSION - 5 d; t_launch >= 0. Moving flyby i re-evaluates legs i and i+1 and every
following leg up to (and including) the first Lambert leg after i+1 (a Lambert arc's arrival velocity is fixed by the
end points, a linear leg's arrival velocity depends on its start velocity); the mass profile is re-accumulated along the
whole chain. Accepted moves are re-evaluated from scratch.

Outputs: results/retime/<camp>/tour_scK.json (search.tour_json keys + per-leg 'model' + 'retime' bookkeeping),
results/retime/retime_summary.json, analysis/retime_probe.md.
Usage: exp_retime.py [--camps camp1 camp2] [--nproc 2] [--only K ...] [--check] [--out results/retime]
"""
import sys, json, time, pathlib, argparse
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from ctoc14.constants import DAY, AU, T_MISSION, VINF_MAX, M_DRY, M0_MAX, TMAX, VE, cost_sc
from ctoc14.kepler import Ephemeris
from ctoc14.lambert import lambert
from ctoc14.linleg import LinLeg

ROOT = pathlib.Path(__file__).resolve().parents[1]
ETA = 0.6; LIN_DTAU = 10.0 * DAY; LIN_UB = 0.8; LIN_DRMAX = 0.15 * AU       # planner settings of the pool campaigns
DV_MAX = 1.2                                                               # km/s cap on re-evaluated legs
TOF_MIN_LAMBERT = 15.0 * DAY; TOF_MIN_LIN = 20.0 * DAY; TOF_MAX_LIN = 600.0 * DAY
T_LAST_OK = T_MISSION - 5.0 * DAY
UTIL_TOL = 1e-7
SEARCH_M0 = {'camp1': 1000.0, 'camp2': 1100.0}
TANK_RULES = {'1.3f+20': (1.3, 20.0), '1.1f+10': (1.1, 10.0)}
LAMBERT, LINEAR = 0, 1
MODEL_NAME = {LAMBERT: 'lambert', LINEAR: 'linear'}


# ------------------------------------------------------------------------------------------- leg models
def lin_solve(L, j, R, a, ub=LIN_UB, iters=5):
    """LinLeg.solve with the planner's selection (cheapest feasible IRLS iterate) plus the max |u| of the selected
    profile; if no iterate is feasible the least-saturated iterate is returned instead of inf (utilisation > 1)."""
    K = L.K[j]; Mk = L.Prv[j, :K]; dR = R - L.rc[j]; T = len(R)
    if K < 1:
        return np.full(T, np.inf), np.broadcast_to(L.vc[j], (T, 3)).copy(), np.full(T, np.inf), np.linalg.norm(dR, axis=-1)
    MMk = np.einsum('kij,klj->kil', Mk, Mk)
    w = np.ones((T, K)); best_cost = np.full(T, np.inf); best_U = np.zeros((T, K, 3)); best_umax = np.full(T, np.inf)
    alt_cost = np.full(T, np.inf); alt_U = np.zeros((T, K, 3)); alt_umax = np.full(T, np.inf)
    for it in range(iters):
        G = np.einsum('tk,kij->tij', 1 / w, MMk)
        try:
            y = np.linalg.solve(G, dR[..., None])[..., 0]
        except np.linalg.LinAlgError:
            y = np.einsum('tij,tj->ti', np.linalg.pinv(G), dR)
        U = np.einsum('tk,kji,tj->tki', 1 / w, Mk, y)
        un = np.linalg.norm(U, axis=-1)
        cost = un.sum(1) * L.dtau; umax = un.max(1); feas = umax <= ub * a
        better = feas & (cost < best_cost)
        best_cost[better] = cost[better]; best_U[better] = U[better]; best_umax[better] = umax[better]
        alt = umax < alt_umax
        alt_cost[alt] = cost[alt]; alt_U[alt] = U[alt]; alt_umax[alt] = umax[alt]
        w = np.maximum(un, 1e-9)
    fe = np.isfinite(best_cost)
    cost = np.where(fe, best_cost, alt_cost); umax = np.where(fe, best_umax, alt_umax)
    Uo = np.where(fe[:, None, None], best_U, alt_U)
    dv_arr = np.einsum('kij,tkj->ti', L.Pvv[j, :K], Uo)
    return cost, L.vc[j] + dv_arr, umax, np.linalg.norm(dR, axis=-1)


def eval_leg(eph, r_from, v_from, m, ast, t_to, tof, typ, launch=False, vE=None):
    """One leg from state (r_from, v_from, m) to asteroid `ast` at t_to (duration tof) with the given model.
    Returns dict(r, v_arr, dv, umax, dr, vinf) -- utilisation / propellant are computed from (dv, umax, tof, m) later."""
    r_to = eph.ast[ast - 1].state(t_to)[0]
    out = dict(r=r_to, tof=float(tof), typ=typ, umax=np.nan, dr=0.0, vinf=None)
    if typ == LAMBERT:
        v1, v2 = lambert(r_from[None, :], r_to[None, :], np.array([tof]))
        v1 = v1[0]; v2 = v2[0]
        if launch:
            vinf = v1 - vE; vn = np.linalg.norm(vinf)
            out['dv'] = max(0.0, vn - VINF_MAX) if np.isfinite(vn) else np.nan
            out['vinf'] = vinf * (VINF_MAX / vn) if (np.isfinite(vn) and vn > VINF_MAX) else vinf
        else:
            out['dv'] = float(np.linalg.norm(v1 - v_from))
        out['v_arr'] = v2
    else:
        a = TMAX / m * 1e-3
        L = LinLeg(r_from, v_from, np.array([tof]), dtau=LIN_DTAU)
        cost, varr, umax, dr = lin_solve(L, 0, r_to[None, :], a)
        out['dv'] = float(cost[0]); out['v_arr'] = varr[0]; out['umax'] = float(umax[0]); out['dr'] = float(dr[0])
    return out


def mass_profile(m0, dv, tof, umax, typ):
    """Planner mass rule along the chain: Lambert legs dm = m(1-exp(-kappa dv/VE)), kappa = 1 + dv/(2 a tof), util = dv/(eta a tof);
    linear legs dm = m(1-exp(-dv/VE)), util = umax/(ub a). Returns (dm, m_start, util)."""
    n = len(dv); dm = np.zeros(n); ms = np.zeros(n); util = np.zeros(n); m = float(m0)
    for i in range(n):
        a = TMAX / m * 1e-3
        with np.errstate(invalid='ignore', divide='ignore', over='ignore'):
            if typ[i] == LAMBERT:
                util[i] = dv[i] / (ETA * a * tof[i]); kappa = 1.0 + dv[i] / (2.0 * a * tof[i])
                d = m * (1.0 - np.exp(-kappa * dv[i] / VE))
            else:
                util[i] = umax[i] / (LIN_UB * a); d = m * (1.0 - np.exp(-dv[i] / VE))
        ms[i] = m; dm[i] = d; m = m - d
    return dm, ms, util


class Chain:
    """A planned tour under the planner's hybrid leg model."""
    def __init__(self, eph, name, t_launch, m0, asts, times, typ, camp=''):
        self.eph = eph; self.name = name; self.camp = camp
        self.t_launch = float(t_launch); self.m0 = float(m0)
        self.asts = np.asarray(asts, int); self.times = np.asarray(times, float); self.typ = np.asarray(typ, int)
        self.n = len(self.asts)
        self.legs = [None] * self.n
        self.evaluate_from(0)

    @property
    def tof(self):
        return np.diff(np.concatenate([[self.t_launch], self.times]))

    def evaluate_from(self, i0, stop=None, m0=None):
        """Re-evaluate legs i0..stop (default: all) and the mass profile of the whole chain."""
        stop = self.n - 1 if stop is None else stop
        rE, vE = self.eph.earth_state(self.t_launch); self.rE = rE; self.vE = vE
        tof = self.tof
        if m0 is not None: self.m0 = float(m0)
        m = self.m0 if i0 == 0 else float(self.m_start[i0])      # legs < i0 are unchanged, so their mass profile is exact
        for i in range(i0, stop + 1):
            if i == 0:
                r_from, v_from = rE, None
            else:
                r_from, v_from = self.legs[i - 1]['r'], self.legs[i - 1]['v_arr']
            lg = eval_leg(self.eph, r_from, v_from, m, int(self.asts[i]), float(self.times[i]), float(tof[i]),
                          int(self.typ[i]), launch=(i == 0), vE=vE)
            self.legs[i] = lg
            a = TMAX / m * 1e-3
            with np.errstate(invalid='ignore', over='ignore'):
                kappa = 1.0 + lg['dv'] / (2.0 * a * tof[i]) if lg['typ'] == LAMBERT else 1.0
                m = m - m * (1.0 - np.exp(-kappa * lg['dv'] / VE))
        self.refresh_mass()

    def refresh_mass(self):
        self.dv = np.array([l['dv'] for l in self.legs]); self.umax = np.array([l['umax'] for l in self.legs])
        self.dr = np.array([l['dr'] for l in self.legs])
        self.dm, self.m_start, self.util = mass_profile(self.m0, self.dv, self.tof, self.umax, self.typ)
        self.fuel = float(np.sum(self.dm))
        self.vinf = self.legs[0]['vinf'] if self.n else np.zeros(3)

    def copy(self):
        c = Chain.__new__(Chain); c.__dict__.update(self.__dict__)
        c.times = self.times.copy(); c.legs = list(self.legs)
        return c

    def stats(self):
        fin = np.isfinite(self.dv)
        return dict(n=int(self.n), fuel=self.fuel, m0=self.m0, sum_dv=float(np.sum(self.dv[fin])), max_dv=float(np.max(self.dv[fin])),
                    n_dv_gt09=int(np.sum(self.dv > 0.9)), n_dv_gt12=int(np.sum(self.dv > 1.2 + 1e-9)), max_util=float(np.max(self.util)),
                    n_util_gt1=int(np.sum(self.util > 1 + UTIL_TOL)), n_linear=int(np.sum(self.typ == LINEAR)),
                    vinf=float(np.linalg.norm(self.vinf)), t_launch_d=self.t_launch / DAY, t_end_d=float(self.times[-1] / DAY),
                    median_leg_d=float(np.median(self.tof) / DAY), max_leg_d=float(np.max(self.tof) / DAY),
                    tanks={nm: dict(zip(('m0', 'fuel', 'max_util'), fixed_point_m0(self, *kc))) for nm, kc in TANK_RULES.items()})

    def to_json(self):
        return dict(m0=self.m0, t_launch=self.t_launch, vinf=[float(x) for x in self.vinf],
                    legs=[dict(ast=int(a), t_flyby=float(t), dv_est=float(d), tof=float(T), model=MODEL_NAME[int(y)])
                          for a, t, d, T, y in zip(self.asts, self.times, self.dv, self.tof, self.typ)],
                    dv_total_est=float(np.sum(self.dv)), fuel_est=self.fuel, n_flybys=int(self.n), t_end=float(self.times[-1]))


def fixed_point_m0(C, k, c, it=40):
    """m0 = 600 + k*fuel(m0) + c with fuel(m0) from the planner mass rule on the chain's (dv, tof, umax)."""
    m0 = C.m0
    for _ in range(it):
        dm, _, _ = mass_profile(m0, C.dv, C.tof, C.umax, C.typ)
        m0n = min(M0_MAX, M_DRY + k * float(np.sum(dm)) + c)
        if abs(m0n - m0) < 1e-6:
            m0 = m0n; break
        m0 = m0n
    dm, _, util = mass_profile(m0, C.dv, C.tof, C.umax, C.typ)
    return float(m0), float(np.sum(dm)), float(np.max(util))


# ------------------------------------------------------------------------------------------- reconstruction
def identify_models(eph, tour, m0_search, name='tour', tol=1e-4):
    """Walk the stored chain at the search m0, evaluating both leg models; pick the one matching dv_est. Returns
    (types, report)."""
    asts = [l['ast'] for l in tour['legs']]; times = [l['t_flyby'] for l in tour['legs']]
    dv_est = np.array([l['dv_est'] for l in tour['legs']], float)
    tL = float(tour['t_launch']); rE, vE = eph.earth_state(tL)
    n = len(asts); tof = np.diff(np.concatenate([[tL], times]))
    typ = np.zeros(n, int); m = m0_search; r_from, v_from = rE, None; mism = []; dvs = np.zeros(n)
    for i in range(n):
        cands = {}
        for ty in ((LAMBERT,) if i == 0 else (LAMBERT, LINEAR)):
            lg = eval_leg(eph, r_from, v_from, m, asts[i], times[i], tof[i], ty, launch=(i == 0), vE=vE)
            cands[ty] = lg
        errs = {ty: abs(lg['dv'] - dv_est[i]) if np.isfinite(lg['dv']) else np.inf for ty, lg in cands.items()}
        ty = min(errs, key=errs.get)
        if errs[ty] > tol:
            mism.append((i, float(errs[ty]), {MODEL_NAME[k]: float(v['dv']) for k, v in cands.items()}, float(dv_est[i])))
        typ[i] = ty; lg = cands[ty]; dvs[i] = lg['dv']
        a = TMAX / m * 1e-3
        if ty == LAMBERT:
            kappa = 1.0 + lg['dv'] / (2 * a * tof[i]); d = m * (1 - np.exp(-kappa * lg['dv'] / VE))
        else:
            d = m * (1 - np.exp(-lg['dv'] / VE))
        m -= d; r_from, v_from = lg['r'], lg['v_arr']
    rep = dict(name=name, n=n, n_linear=int(np.sum(typ == LINEAR)), fuel_search_m0=float(m0_search - m), fuel_est=float(tour['fuel_est']),
               max_ddv=float(np.max(np.abs(dvs - dv_est))), mismatches=mism)
    return typ, rep


def load_chain(eph, camp, k, path=None):
    path = path or (ROOT / 'results' / 'phasing' / camp / 'milp_n12' / f'tour_sc{k}.json')
    tour = json.load(open(path))
    typ, rep = identify_models(eph, tour, SEARCH_M0[camp], name=f'{camp}_sc{k}')
    C = Chain(eph, f'{camp}_sc{k}', tour['t_launch'], tour['m0'], [l['ast'] for l in tour['legs']], [l['t_flyby'] for l in tour['legs']], typ, camp=camp)
    rep['fuel_stored_m0'] = C.fuel; rep['m0_stored'] = C.m0; rep['max_util_stored_m0'] = float(np.max(C.util))
    rep['n_util_gt1_stored_m0'] = int(np.sum(C.util > 1 + UTIL_TOL)); rep['path'] = str(path)
    return C, tour, rep


# ------------------------------------------------------------------------------------------- moves
def _stop_index(C, i_first_changed_state):
    """Last leg to re-evaluate when the arrival state of leg `i_first_changed_state` changes: leg i+1 always (its junction
    or, if linear, its whole solution), then every following linear leg (their arrival velocity depends on the start
    velocity) up to and including the first Lambert leg (its arrival velocity is fixed by the end points)."""
    j = i_first_changed_state + 1
    while j < C.n and C.typ[j] == LINEAR:
        j += 1
    return min(j, C.n - 1)


def trial_chain(C, kind, i, delta):
    """Copy of C with flyby i (kind 'flyby') or the launch (kind 'launch') shifted by delta [s]; partial re-evaluation."""
    T = C.copy()
    if kind == 'flyby':
        T.times[i] += delta; i0 = i; stop = _stop_index(C, i + 1)      # legs i, i+1 change; leg i+1's arrival state changes -> i+2 (+ linear cascade up to the first Lambert leg)
    else:
        T.t_launch += delta; i0 = 0; stop = _stop_index(C, 0)
    T.evaluate_from(i0, stop)
    return T, i0, stop


def feasible(T, C, i0, stop):
    """Constraints of trial T against base C (legs i0..stop re-evaluated)."""
    if not np.isfinite(T.fuel) or not np.all(np.isfinite(T.dv)):
        return False
    tof = T.tof
    if T.t_launch < 0 or np.any(T.times > T_LAST_OK):
        return False
    tmin = np.where(T.typ == LAMBERT, TOF_MIN_LAMBERT, TOF_MIN_LIN)
    if np.any(tof < tmin - 1e-6) or np.any((T.typ == LINEAR) & (tof > TOF_MAX_LIN + 1e-6)):
        return False
    if np.any(T.util > np.maximum(1.0, C.util) + UTIL_TOL):
        return False
    mod = slice(i0, stop + 1)
    if np.any(T.dv[mod] > np.maximum(DV_MAX, C.dv[mod]) + 1e-9):
        return False
    lin = T.typ[mod] == LINEAR
    if np.any(T.dr[mod][lin] > np.maximum(LIN_DRMAX, C.dr[mod][lin]) + 1e-6):
        return False
    return True


def try_move(C, kind, i, deltas, log):
    best = None
    for d in deltas:
        if abs(d) < 1e-9:
            continue
        T, i0, stop = trial_chain(C, kind, i, d)
        if feasible(T, C, i0, stop) and T.fuel < C.fuel - 1e-9 and (best is None or T.fuel < best[0].fuel):
            best = (T, d, i0, stop)
    if best is None:
        return C, 0.0
    T, d, i0, stop = best
    fuel_partial = T.fuel
    T.evaluate_from(0)                                                    # full re-evaluation of the accepted move
    if abs(T.fuel - fuel_partial) > 1e-3:
        log.append(f'{kind} {i} delta {d/DAY:+.2f} d: partial {fuel_partial:.3f} vs full {T.fuel:.3f} kg')
    if not (feasible(T, C, 0, T.n - 1) and T.fuel < C.fuel - 1e-9):
        log.append(f'{kind} {i} delta {d/DAY:+.2f} d rejected after full re-evaluation')
        return C, 0.0
    return T, float(d)


def refine_chain(C, coarse=10.0, coarse_step=1.0, fine=1.0, fine_step=0.25, launch=30.0, launch_step=2.5,
                 min_gain=0.5, max_passes=40, verbose=False):
    tic = time.time(); log = []
    d_coarse = np.arange(-coarse, coarse + 1e-9, coarse_step) * DAY
    d_fine = np.arange(-fine, fine + 1e-9, fine_step) * DAY
    d_launch = np.arange(-launch, launch + 1e-9, launch_step) * DAY
    t0 = C.times.copy(); tL0 = C.t_launch; fuel0 = C.fuel; passes = []
    for p in range(max_passes):
        f_start = C.fuel
        C, _ = try_move(C, 'launch', -1, d_launch, log)
        for i in range(C.n):
            C, _ = try_move(C, 'flyby', i, d_coarse, log)
            C, _ = try_move(C, 'flyby', i, d_fine, log)
        gain = f_start - C.fuel
        passes.append(dict(pass_=p + 1, fuel=float(C.fuel), gain=float(gain)))
        if verbose:
            print(f'  {C.name} pass {p+1}: fuel {C.fuel:.2f} kg (gain {gain:.2f}) {time.time()-tic:.0f} s', flush=True)
        if gain < min_gain:
            break
    shifts = (C.times - t0) / DAY
    info = dict(name=C.name, n=int(C.n), fuel_before=float(fuel0), fuel_after=float(C.fuel), saving=float(fuel0 - C.fuel),
                launch_shift_d=float((C.t_launch - tL0) / DAY), max_shift_d=float(np.max(np.abs(shifts))),
                mean_abs_shift_d=float(np.mean(np.abs(shifts))), n_shifted=int(np.sum(np.abs(shifts) > 1e-6)),
                shifts_d=[float(x) for x in shifts], n_passes=len(passes), passes=passes, runtime_s=float(time.time() - tic), log=log)
    return C, info


# ------------------------------------------------------------------------------------------- driver
_EPH = None
def _init():
    global _EPH; _EPH = Ephemeris()

def _job(args):
    camp, k, opts = args
    C, tour, rep = load_chain(_EPH, camp, k)
    before = C.stats()
    Cr, info = refine_chain(C, verbose=opts.get('verbose', False), min_gain=opts['min_gain'], max_passes=opts['max_passes'])
    after = Cr.stats()
    js = Cr.to_json()
    js['inserted'] = list(tour.get('inserted', []))
    js['retime'] = dict(source=rep['path'], fuel_est_source=float(tour['fuel_est']), m0_search=SEARCH_M0[camp],
                        fuel_before_stored_m0=before['fuel'], launch_shift_d=info['launch_shift_d'], max_shift_d=info['max_shift_d'],
                        shifts_d=info['shifts_d'], t_launch_before=float(tour['t_launch']), tanks=after['tanks'])
    return dict(camp=camp, k=k, path=rep['path'], recon=rep, before=before, after=after, info=info, json=js)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--camps', nargs='+', default=['camp1', 'camp2'])
    ap.add_argument('--nproc', type=int, default=2)
    ap.add_argument('--out', default='results/retime')
    ap.add_argument('--report', default='analysis/retime_probe.md')
    ap.add_argument('--only', type=int, nargs='*', default=None, help='restrict to these sc numbers')
    ap.add_argument('--min-gain', type=float, default=0.5)
    ap.add_argument('--max-passes', type=int, default=40)
    ap.add_argument('--check', action='store_true', help='only print the reconstruction of the stored tours')
    ap.add_argument('--verbose', action='store_true')
    a = ap.parse_args()
    jobs = []
    for camp in a.camps:
        d = ROOT / 'results' / 'phasing' / camp / 'milp_n12'
        ks = sorted(int(f.stem.split('sc')[1]) for f in d.glob('tour_sc*.json') if f.stem.split('sc')[1].isdigit())
        jobs += [(camp, k, dict(min_gain=a.min_gain, max_passes=a.max_passes, verbose=a.verbose)) for k in ks if not a.only or k in a.only]
    if a.check:
        _init()
        for camp, k, _ in jobs:
            C, tour, rep = load_chain(_EPH, camp, k)
            print(f"{camp} sc{k:2d}: n={rep['n']} linear={rep['n_linear']} fuel@search m0 {rep['fuel_search_m0']:.2f} vs fuel_est {rep['fuel_est']:.2f} "
                  f"(max ddv {rep['max_ddv']:.1e}); fuel@stored m0 {C.m0:.0f}: {C.fuel:.1f} kg, max util {rep['max_util_stored_m0']:.3f} "
                  f"({rep['n_util_gt1_stored_m0']} legs > 1); mismatches {rep['mismatches'][:3]}", flush=True)
        return
    out = ROOT / a.out; out.mkdir(parents=True, exist_ok=True)
    tic = time.time()
    if a.nproc > 1:
        import multiprocessing as mp
        with mp.get_context('fork').Pool(a.nproc, initializer=_init) as pool:
            res = pool.map(_job, jobs, chunksize=1)
    else:
        _init(); res = [_job(j) for j in jobs]
    print(f'{len(res)} tours refined in {time.time() - tic:.0f} s', flush=True)
    for r in res:
        (out / r['camp']).mkdir(parents=True, exist_ok=True)
        json.dump(r['json'], open(out / r['camp'] / f'tour_sc{r["k"]}.json', 'w'), indent=1)
    summ = dict(runs=[{key: r[key] for key in ('camp', 'k', 'path', 'recon', 'before', 'after', 'info')} for r in res], fleet={})
    for camp in a.camps:
        rs = [r for r in res if r['camp'] == camp]
        if not rs: continue
        sp = ROOT / 'results' / 'phasing' / camp / 'milp_n12' / 'summary.json'
        n_miss = len(json.load(open(sp))['missed']) if sp.exists() else None
        fl = dict(n_tours=len(rs), n_miss_planned=n_miss, n_flybys=int(sum(r['before']['n'] for r in rs)),
                  fuel_est_source=float(sum(r['recon']['fuel_est'] for r in rs)),
                  fuel_before=float(sum(r['before']['fuel'] for r in rs)), fuel_after=float(sum(r['after']['fuel'] for r in rs)),
                  J_stored_m0=float(sum(cost_sc(r['before']['m0']) for r in rs)))
        for nm in TANK_RULES:
            for stage in ('before', 'after'):
                fl[f'sumJ_{nm}_{stage}'] = float(sum(cost_sc(r[stage]['tanks'][nm]['m0']) for r in rs))
                fl[f'fuel_{nm}_{stage}'] = float(sum(r[stage]['tanks'][nm]['fuel'] for r in rs))
        summ['fleet'][camp] = fl
    json.dump(summ, open(out / 'retime_summary.json', 'w'), indent=1, default=float)
    write_report(ROOT / a.report, summ, res, a)
    for r in res:
        b, f, i = r['before'], r['after'], r['info']
        print(f"{r['camp']} sc{r['k']:2d}: n={b['n']:2d} lin={b['n_linear']:2d} fuel {b['fuel']:6.1f} -> {f['fuel']:6.1f} kg (-{i['saving']:5.1f}) "
              f"sum dv {b['sum_dv']:5.2f} -> {f['sum_dv']:5.2f} max shift {i['max_shift_d']:5.2f} d launch {i['launch_shift_d']:+5.1f} d "
              f"passes {i['n_passes']} {i['runtime_s']:5.0f} s max util {b['max_util']:.3f}->{f['max_util']:.3f} log {len(i['log'])}", flush=True)
    for camp, fl in summ['fleet'].items():
        print(f"{camp}: fuel {fl['fuel_before']:.1f} -> {fl['fuel_after']:.1f} kg; sum J_i (1.3f+20) {fl['sumJ_1.3f+20_before']:.3f} -> {fl['sumJ_1.3f+20_after']:.3f}; "
              f"(1.1f+10) {fl['sumJ_1.1f+10_before']:.3f} -> {fl['sumJ_1.1f+10_after']:.3f}; planned misses {fl['n_miss_planned']}")


def write_report(path, summ, res, a):
    L = [f'# Re-timing probe: coordinate-descent flyby-time refinement and tank right-sizing of planned tours\n',
         f'Script `tools/exp_retime.py` (min gain per pass {a.min_gain} kg, max passes {a.max_passes}); source tours '
         f'`results/phasing/<camp>/milp_n12/tour_sc*.json`; refined tours `results/retime/<camp>/tour_scK.json`; raw numbers '
         f'`results/retime/retime_summary.json`. (Auto-generated tables; the discussion is appended by hand below.)\n',
         '## Per-tour results (propellant at the stored m0)\n',
         '| tour | n | linear legs | fuel_est (search m0) | fuel before | fuel after | saving | sum dv before | after | max dv before/after | legs dv>0.9 before/after | max util before/after | max shift | launch shift | passes | runtime |',
         '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
    for r in res:
        b, f, i, rc = r['before'], r['after'], r['info'], r['recon']
        L.append(f"| {r['camp']} sc{r['k']} | {b['n']} | {b['n_linear']} | {rc['fuel_est']:.1f} | {b['fuel']:.1f} | {f['fuel']:.1f} | {i['saving']:.1f} ({100*i['saving']/b['fuel']:.1f}%) | "
                 f"{b['sum_dv']:.2f} | {f['sum_dv']:.2f} | {b['max_dv']:.2f}/{f['max_dv']:.2f} | {b['n_dv_gt09']}/{f['n_dv_gt09']} | {b['max_util']:.2f}/{f['max_util']:.2f} | "
                 f"{i['max_shift_d']:.2f} d | {i['launch_shift_d']:+.1f} d | {i['n_passes']} | {i['runtime_s']:.0f} s |")
    L += ['', '## Fleet-level planned J (12 craft per camp; misses = planned misses of the MILP selection)\n',
          '| camp | flybys | planned misses | sum fuel_est (search m0) | fuel before (stored m0) | fuel after | sum J_i stored m0 | sum J_i 1.3f+20 before | after | sum J_i 1.1f+10 before | after | J(1.3f+20) after + misses | J(1.1f+10) after + misses |',
          '|---|---|---|---|---|---|---|---|---|---|---|---|---|']
    for camp, fl in summ['fleet'].items():
        nm = fl['n_miss_planned'] or 0
        L.append(f"| {camp} | {fl['n_flybys']} | {fl['n_miss_planned']} | {fl['fuel_est_source']:.1f} | {fl['fuel_before']:.1f} | {fl['fuel_after']:.1f} | {fl['J_stored_m0']:.3f} | "
                 f"{fl['sumJ_1.3f+20_before']:.3f} | {fl['sumJ_1.3f+20_after']:.3f} | {fl['sumJ_1.1f+10_before']:.3f} | {fl['sumJ_1.1f+10_after']:.3f} | "
                 f"{fl['sumJ_1.3f+20_after'] + nm:.3f} | {fl['sumJ_1.1f+10_after'] + nm:.3f} |")
    L += ['', 'Per-tour tank sizes (fixed point m0 = 600 + k fuel(m0) + c) and propellant under each rule:\n',
          '| tour | m0 stored | J_i stored | m0 1.3f+20 before | after | fuel after | J_i after | m0 1.1f+10 before | after | fuel after | max util after | J_i after |',
          '|---|---|---|---|---|---|---|---|---|---|---|---|']
    for r in res:
        b, f = r['before'], r['after']
        t13b, t13a, t11b, t11a = b['tanks']['1.3f+20'], f['tanks']['1.3f+20'], b['tanks']['1.1f+10'], f['tanks']['1.1f+10']
        L.append(f"| {r['camp']} sc{r['k']} | {b['m0']:.1f} | {cost_sc(b['m0']):.3f} | {t13b['m0']:.1f} | {t13a['m0']:.1f} | {t13a['fuel']:.1f} | {cost_sc(t13a['m0']):.3f} | "
                 f"{t11b['m0']:.1f} | {t11a['m0']:.1f} | {t11a['fuel']:.1f} | {t11a['max_util']:.2f} | {cost_sc(t11a['m0']):.3f} |")
    L.append('')
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('\n'.join(L) + '\n')


if __name__ == '__main__':
    main()
