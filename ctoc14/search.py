"""Beam search for one spacecraft's flyby sequence using a Lambert-junction low-thrust feasibility model.

Model: after flying by asteroid A at time t with velocity v (Lambert arrival velocity of the previous arc), a leg of
duration tof to asteroid B is a Lambert arc; the junction delta-v dv = |v1_lambert - v| must be deliverable by low thrust
within the leg: dv <= eta * a(m) * tof (a = Tmax/m). Propellant: dm = m (1 - exp(-kappa dv / VE)), kappa = 1 + dv/(2 a tof)
(leverage penalty of position-only targeting). Child velocity = Lambert arrival velocity.
Launch: position = Earth, v = v_Earth + v_inf, |v_inf| <= 4 km/s free; excess |v1 - vE| - 4 must be paid by thrust.
"""
import json, numpy as np
from dataclasses import dataclass, field
from .constants import VE, TMAX, T_MISSION, DAY, M_DRY, FUEL_MAX, VINF_MAX, AU
from .kepler import Ephemeris, propagate_twobody
from .lambert import lambert
from .linleg import LinLeg

UNREACHABLE = {131, 144}

@dataclass
class Params:
    tofs: np.ndarray = field(default_factory=lambda: np.arange(15, 401, 5) * DAY)
    eta: float = 0.6            # max junction dv fraction of a*tof
    w_fuel: float = 1.5         # beam score weight on fuel fraction vs elapsed-time fraction
    beam: int = 300
    n_per_target: int = 2       # best TOFs kept per target per parent
    max_children: int = 60      # per parent
    prefilter_margin: float = 1.4
    m_margin: float = 2.0       # kg kept above dry mass in the plan
    max_depth: int = 120
    abs_time: bool = True       # score on absolute mission time (window is absolute)
    w_rare: float = 0.0         # bonus weight for visiting poorly-connected ('rare') targets
    rarity: object = None       # array (300,) in [0,1]; 1 = hardest to reach; None = no weighting
    dv_max: float = 1e9         # cap on junction dv per leg (km/s) to limit conversion losses
    lin_tofs: object = None     # TOF grid (s) for the linearised continuous-thrust leg model (None = Lambert only)
    lin_ub: float = 0.8         # thrust bound fraction for the linear model (|u_k| <= lin_ub * Tmax/m)
    lin_drmax: float = 0.25 * AU  # max coast-endpoint displacement for the linear model (linearisation validity), km
    lin_dtau: float = 10 * DAY  # segment length of the linear model
    vinf_cap: float = VINF_MAX  # cap on |v_inf| of launch legs (< 4 keeps the craft near Earth's orbit ring)
    tof_refine: bool = False    # refine the flyby time of Lambert children locally (1 d then 0.25 d) before scoring
    prize: object = None        # array (300,) of per-target prizes in J units (column-generation duals); overrides rarity/w_rare,
                                # enters the beam score with weight 1 and the child pruning (None = off, old behaviour)
    w_t: float = 1.0            # weight of the (absolute) mission-time fraction in the score
    first_targets: object = None  # optional set of asteroid IDs allowed as the FIRST target (roots); None = all
    dv_max_t: object = None     # optional array (300,) of per-TARGET leg caps (km/s) replacing dv_max for legs ending there
    win_lo: object = None       # optional (300,) earliest flyby time (s, mission time) per target; -inf = free
    win_hi: object = None       # optional (300,) latest flyby time per target; +inf = free.  A windowed target is MANDATORY:
                                # a state whose time passes win_hi without having visited it is dropped (deadline prune)
    mandatory: object = None    # optional list of asteroid IDs that MUST be visited (deadline prune); None = every target
                                # with a finite win_hi.  Set it to use windows on non-mandatory targets too (segment ends).
    lookahead: bool = True      # with windows: rank children by score + Lambert fuel to the next mandatory target
    lookahead_max: float = 250 * DAY  # ... only when that target's window is within this horizon (single-rev Lambert validity)

@dataclass
class State:
    t: float; r: np.ndarray; v: np.ndarray; m: float
    visited: int                # bitmask over asteroid index (bit k = ID k+1)
    seq: tuple                  # ((ast_id, t_flyby, dv_est, tof), ...)
    fuel: float; t_launch: float; vinf: np.ndarray; m0: float
    rare: float = 0.0
    def score(self, P):
        tt = self.t if P.abs_time else (self.t - self.t_launch)
        return P.w_t * tt / T_MISSION + P.w_fuel * self.fuel / FUEL_MAX - (self.rare if P.prize is not None else P.w_rare * self.rare)
    def n(self):
        return len(self.seq)

def _cap(P, ia):
    """Leg cap for legs ending at asteroid indices ia (0-based)."""
    return P.dv_max if P.dv_max_t is None else P.dv_max_t[ia]

def _in_window(P, ia, t):
    """Mask of candidate legs (asteroid indices ia, flyby times t) that respect the per-target windows."""
    if P.win_lo is None:
        return np.ones(len(ia), bool)
    return (t >= P.win_lo[ia]) & (t <= P.win_hi[ia])

def _mandatory(P):
    if P.win_hi is None:
        return None
    if P.mandatory is not None:
        return np.array([int(i) - 1 for i in P.mandatory], dtype=int)
    return np.where(np.isfinite(P.win_hi))[0]

def _alive(P, s, mand):
    """False if state s has passed the deadline of a mandatory target it has not visited."""
    if mand is None or len(mand) == 0:
        return True
    vis = _visited_array(s.visited)
    return not np.any((P.win_hi[mand] < s.t) & ~vis[mand])

def _lookahead(eph, P, states, mand):
    """A*-style heuristic (score units): fuel fraction of the cheapest single Lambert leg from each state to its NEXT
    unvisited mandatory target (earliest deadline first), evaluated at the window start, middle and end.  States that
    can no longer reach their next waypoint cheaply sink in the ranking instead of crowding the beam until the deadline."""
    h = np.zeros(len(states))
    if mand is None or len(mand) == 0 or not states:
        return h
    hi = P.win_hi[mand]; lo = P.win_lo[mand]
    r1 = []; r2 = []; tof = []; who = []; m = []
    for i, s in enumerate(states):
        vis = _visited_array(s.visited)[mand]
        if vis.all():
            continue
        j = np.argmin(np.where(vis, np.inf, hi)); k = mand[j]
        for tt in (max(lo[j], s.t + 15 * DAY), 0.5 * (max(lo[j], s.t + 15 * DAY) + hi[j]), hi[j]):
            if tt - s.t < 10 * DAY or tt - s.t > P.lookahead_max:   # single-rev Lambert is garbage beyond ~1 rev
                continue
            r1.append(s.r); tof.append(tt - s.t); who.append(i); m.append(s.m)
            r2.append((k, tt))
    if not who:
        return h
    r1 = np.array(r1); tof = np.array(tof); who = np.array(who); m = np.array(m)
    R = eph.ast_states_at(np.array([k for k, _ in r2]), np.array([t for _, t in r2]))[0]
    v1, _ = lambert(r1, R, tof)
    dv = np.linalg.norm(v1 - np.array([states[i].v for i in who]), axis=-1)
    dv = np.where(np.isfinite(dv), dv, 50.0)
    dm = m * (1.0 - np.exp(-dv / VE))
    best = np.full(len(states), np.inf)
    np.minimum.at(best, who, dm)
    has = np.isfinite(best)
    h[has] = P.w_fuel * best[has] / FUEL_MAX
    return h

def _rare_vec(P):
    """Per-target bonus accumulated in State.rare: the prize (J units) if set, else the rarity array (weighted by w_rare in score)."""
    return P.prize if P.prize is not None else P.rarity

def mask_of(ids):
    m = 0
    for i in ids: m |= 1 << (int(i) - 1)
    return m

def _visited_array(mask, n=300):
    return np.array([(mask >> k) & 1 for k in range(n)], dtype=bool)

def _expand_linear(eph, s, P, a):
    """Candidate legs from the linearised continuous-thrust model (coast + distributed correction). Returns lists."""
    tofs = P.lin_tofs[(s.t + P.lin_tofs) <= T_MISSION]
    if len(tofs) == 0:
        return [], [], [], [], []
    L = LinLeg(s.r, s.v, tofs, dtau=P.lin_dtau)
    vis = _visited_array(s.visited)
    ia, v2, dv, tof, dm = [], [], [], [], []
    for j, tj in enumerate(tofs):
        R = eph.all_ast_states(s.t + tj)[0]
        near = (np.linalg.norm(R - L.rc[j], axis=-1) <= P.lin_drmax) & ~vis
        if not near.any():
            continue
        idx = np.where(near)[0]
        cost, varr, _ = L.solve(j, R[idx], a, ub=P.lin_ub)
        f = np.isfinite(cost) & (cost <= _cap(P, idx)) & _in_window(P, idx, s.t + tj)
        if not f.any():
            continue
        dmj = s.m * (1 - np.exp(-cost[f] / VE))
        okm = (s.m - dmj) >= M_DRY + P.m_margin
        for k, vv, c, d in zip(idx[f][okm], varr[f][okm], cost[f][okm], dmj[okm]):
            ia.append(int(k)); v2.append(vv); dv.append(float(c)); tof.append(float(tj)); dm.append(float(d))
    return ia, v2, dv, tof, dm

def expand(eph, s, P):
    """Children of state s."""
    a = TMAX / s.m * 1e-3                          # km/s^2
    tofs = P.tofs[(s.t + P.tofs) <= T_MISSION]
    if len(tofs) == 0:
        return []
    rc, _ = propagate_twobody(s.r, s.v, tofs)      # coast positions (nT,3)
    R = np.empty((len(tofs), eph.n_ast, 3))
    for i, tof in enumerate(tofs):
        R[i] = eph.all_ast_states(s.t + tof)[0]
    d = np.linalg.norm(R - rc[:, None, :], axis=-1)
    lim = P.prefilter_margin * P.eta * a * tofs ** 2
    cand = d <= lim[:, None]
    cand[:, _visited_array(s.visited)] = False
    idx = np.argwhere(cand)
    it, ia = idx[:, 0], idx[:, 1]
    v1, v2 = lambert(np.broadcast_to(s.r, (len(it), 3)), R[it, ia], tofs[it]) if len(it) else (np.zeros((0, 3)), np.zeros((0, 3)))
    dv = np.linalg.norm(v1 - s.v, axis=-1)
    tof = tofs[it]
    if P.tof_refine and len(ia):
        # refine the encounter time of every plausible candidate (junction dv within 3x the cap at the coarse grid):
        # coarse +-4 d at 1 d, then +-0.75 d at 0.25 d; keeps the best time per candidate
        pl = np.where(np.isfinite(dv) & (dv <= 3.0 * np.minimum(_cap(P, ia), P.eta * a * tof.max())))[0]
        if len(pl):
            tof = tof.astype(float).copy(); dv = dv.copy(); v2 = v2.copy()
            for step, half in ((1.0, 4.0), (0.25, 0.75)):
                dts = np.arange(-half, half + 1e-9, step) * DAY; nd = len(dts)
                tt = tof[pl][:, None] + dts[None, :]; okt = (tt >= 10 * DAY) & (s.t + tt <= T_MISSION)
                tt_f = tt.ravel(); ia_f = np.repeat(ia[pl], nd); Rf = np.zeros((len(tt_f), 3)); okf = okt.ravel()
                if okf.any():
                    Rf[okf] = eph.ast_states_at(ia_f[okf], s.t + tt_f[okf])[0]
                v1f, v2f = lambert(np.broadcast_to(s.r, (len(tt_f), 3)), Rf, np.maximum(tt_f, 1.0))
                dvf = np.linalg.norm(v1f - s.v, axis=-1); dvf[~okt.ravel() | ~np.isfinite(dvf)] = np.inf
                dvf = dvf.reshape(len(pl), nd); jb = np.argmin(dvf, axis=1); bet = dvf[np.arange(len(pl)), jb] < dv[pl]
                if bet.any():
                    q = pl[bet]; jq = jb[bet]
                    tof[q] = tt[bet, jq]; dv[q] = dvf[bet, jq]; v2[q] = v2f.reshape(len(pl), nd, 3)[bet, jq]
    feas = np.isfinite(dv) & (dv <= P.eta * a * tof) & (dv <= _cap(P, ia)) & _in_window(P, ia, s.t + tof)
    it, ia, v2, dv, tof = it[feas], ia[feas], v2[feas], dv[feas], tof[feas]
    kappa = 1 + dv / (2 * a * tof)
    dm = s.m * (1 - np.exp(-kappa * dv / VE))
    ok = (s.m - dm) >= M_DRY + P.m_margin
    it, ia, v2, dv, tof, dm = it[ok], ia[ok], v2[ok], dv[ok], tof[ok], dm[ok]
    ia, v2, dv, tof, dm = list(ia), list(v2), list(dv), list(tof), list(dm)
    if P.lin_tofs is not None:
        ia2, v22, dv2, tof2, dm2 = _expand_linear(eph, s, P, a)
        ia += ia2; v2 += v22; dv += dv2; tof += tof2; dm += dm2
    if len(ia) == 0:
        return []
    ia = np.array(ia); v2 = np.array(v2).reshape(-1, 3); dv = np.array(dv); tof = np.array(tof); dm = np.array(dm)
    R = None
    cost = P.w_t * tof / T_MISSION + P.w_fuel * dm / FUEL_MAX
    if P.prize is not None:
        cost = cost - P.prize[ia]          # prize-aware child pruning (J units): priced targets survive max_children
    # keep best n_per_target per target
    order = np.lexsort((cost, ia))
    keep = []; last = -1; cnt = 0
    for o in order:
        if ia[o] != last:
            last = ia[o]; cnt = 0
        if cnt < P.n_per_target:
            keep.append(o); cnt += 1
    keep = sorted(keep, key=lambda o: cost[o])[:P.max_children]
    children = []
    rv = _rare_vec(P)
    keep = np.array(keep, dtype=int)
    Rk = eph.ast_states_at(ia[keep], s.t + tof[keep])[0] if len(keep) else np.zeros((0, 3))
    for q, o in enumerate(keep):
        k = ia[o]
        children.append(State(t=s.t + tof[o], r=Rk[q], v=v2[o], m=s.m - dm[o], visited=s.visited | (1 << int(k)),
                              seq=s.seq + ((int(k) + 1, float(s.t + tof[o]), float(dv[o]), float(tof[o])),),
                              fuel=s.fuel + float(dm[o]), t_launch=s.t_launch, vinf=s.vinf, m0=s.m0,
                              rare=s.rare + (float(rv[k]) if rv is not None else 0.0)))
    return children

def roots(eph, t_launch_grid, m0, excluded, P, visited0=None):
    """Launch children: one per (launch date, target, tof) with |v_inf| <= 4 (+ thrust for the excess)."""
    out = []
    a0 = TMAX / m0 * 1e-3
    vis0 = mask_of(set(excluded) | UNREACHABLE) if visited0 is None else (visited0 | mask_of(UNREACHABLE))
    rv = _rare_vec(P)
    if P.first_targets is not None:
        ft_mask = np.zeros(eph.n_ast, bool); ft_mask[[int(i) - 1 for i in P.first_targets]] = True
    for tL in t_launch_grid:
        rE, vE = eph.earth_state(tL)
        tofs = P.tofs[(tL + P.tofs) <= T_MISSION]
        for tof in tofs:
            R, _ = eph.all_ast_states(tL + tof)
            v1, v2 = lambert(np.broadcast_to(rE, R.shape), R, np.full(R.shape[0], tof))
            vinf = np.linalg.norm(v1 - vE, axis=-1)
            excess = np.maximum(0.0, vinf - VINF_MAX)
            feas = np.isfinite(vinf) & (excess <= P.eta * a0 * tof) & (excess <= _cap(P, np.arange(len(excess)))) & ~_visited_array(vis0) & (vinf <= P.vinf_cap)
            if P.first_targets is not None:
                feas &= ft_mask
            feas &= _in_window(P, np.arange(len(excess)), tL + tof)
            for k in np.where(feas)[0]:
                dv = excess[k]
                kappa = 1 + dv / (2 * a0 * tof) if dv > 0 else 1.0
                dm = m0 * (1 - np.exp(-kappa * dv / VE))
                vinf_vec = v1[k] - vE
                if vinf[k] > VINF_MAX:
                    vinf_vec = vinf_vec * (VINF_MAX / vinf[k])
                out.append(State(t=tL + tof, r=R[k], v=v2[k], m=m0 - dm, visited=vis0 | (1 << int(k)),
                                 seq=((int(k) + 1, float(tL + tof), float(dv), float(tof)),), fuel=float(dm),
                                 t_launch=float(tL), vinf=vinf_vec, m0=m0,
                                 rare=(float(rv[k]) if rv is not None else 0.0)))
    return out

def _expand_batch(args):
    states, P = args
    return [c for s in states for c in expand(_EPH, s, P)]

_EPH = None
def _init_worker(eph):
    global _EPH
    _EPH = eph

def _beam_loop(eph, beam, P, n_proc=1, verbose=True, pool=None):
    """Depth-synchronous beam search from an initial list of states. Returns (best_state, final_beam)."""
    global _EPH
    _EPH = eph
    mand = _mandatory(P)
    if mand is not None:
        beam = [s for s in beam if _alive(P, s, mand)]
    beam = sorted(beam, key=lambda s: s.score(P))[:P.beam * 3]
    best = max(beam, key=lambda s: (s.n(), -s.score(P))) if beam else None
    own_pool = False
    if pool is None and n_proc > 1:
        import multiprocessing as mp
        ctx = mp.get_context('fork')
        pool = ctx.Pool(n_proc, initializer=_init_worker, initargs=(eph,)); own_pool = True
    try:
        for depth in range(1, P.max_depth):
            if not beam: break
            if pool is not None and len(beam) >= 2 * n_proc:
                chunks = [beam[i::n_proc] for i in range(n_proc)]
                res = pool.map(_expand_batch, [(c, P) for c in chunks])
                children = [c for r in res for c in r]
            else:
                children = [c for s in beam for c in expand(eph, s, P)]
            if mand is not None:
                children = [c for c in children if _alive(P, c, mand)]
            if not children:
                break
            if mand is not None and P.lookahead:
                hh = _lookahead(eph, P, children, mand)
                order = np.argsort(np.array([c.score(P) for c in children]) + hh, kind='stable')
                children = [children[i] for i in order]
            else:
                children.sort(key=lambda s: s.score(P))
            seen = set(); nb = []
            for c in children:
                key = (c.visited, int(c.t // (5 * DAY)))
                if key in seen: continue
                seen.add(key); nb.append(c)
                if len(nb) >= P.beam: break
            beam = nb
            coll = getattr(P, 'collect', None)
            if coll is not None: coll.extend(beam)
            top = beam[0]
            if top.n() > best.n() or (top.n() == best.n() and top.score(P) < best.score(P)):
                best = top
            if verbose:
                print(f'depth {depth+1:3d}: {len(children):6d} children -> beam {len(beam):4d}; best n={best.n()} '
                      f't={top.t/DAY/365.25:5.2f}yr fuel={top.fuel:6.0f}kg dv_sum={sum(x[2] for x in top.seq):5.1f} launch={top.t_launch/DAY:5.0f}d', flush=True)
    finally:
        if own_pool:
            pool.close(); pool.join()
    return best, beam

def beam_search(eph, excluded=(), m0=2000.0, t_launch_grid=None, P=None, n_proc=8, verbose=True):
    """Full search from Earth launches. Returns (best_state, final_beam)."""
    P = P or Params()
    if t_launch_grid is None:
        t_launch_grid = np.arange(0, 3 * 365.25, 20) * DAY
    beam = roots(eph, t_launch_grid, m0, excluded, P)
    if verbose: print(f'roots: {len(beam)} launch options')
    beam.sort(key=lambda s: s.score(P))
    seen = set(); b2 = []
    for s in beam:
        key = (s.seq[0][0], int(s.t_launch // (20 * DAY)))
        if key in seen: continue
        seen.add(key); b2.append(s)
    return _beam_loop(eph, b2, P, n_proc, verbose)

def beam_search_from_state(eph, t, r, v, m, visited_ids, m0, t_launch, vinf, P=None, n_proc=1, verbose=False, seq=()):
    """Continue a tour from an arbitrary spacecraft state. visited_ids: asteroid IDs that must not be targeted."""
    P = P or Params()
    s0 = State(t=float(t), r=np.asarray(r, float), v=np.asarray(v, float), m=float(m), visited=mask_of(set(visited_ids) | UNREACHABLE),
               seq=tuple(seq), fuel=float(m0 - m), t_launch=float(t_launch), vinf=np.asarray(vinf, float), m0=float(m0))
    return _beam_loop(eph, [s0], P, n_proc, verbose)

def tour_json(s):
    return dict(m0=s.m0, t_launch=s.t_launch, vinf=list(map(float, s.vinf)),
                legs=[dict(ast=a, t_flyby=t, dv_est=dv, tof=tof) for (a, t, dv, tof) in s.seq],
                dv_total_est=float(sum(x[2] for x in s.seq)), fuel_est=float(s.fuel), n_flybys=len(s.seq), t_end=float(s.t))

def save_tour(s, path):
    with open(path, 'w') as f:
        json.dump(tour_json(s), f, indent=1)
