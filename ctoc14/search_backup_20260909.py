"""Beam search for one spacecraft's flyby sequence using a Lambert-junction low-thrust feasibility model.

Model: after flying by asteroid A at time t with velocity v (Lambert arrival velocity of the previous arc), a leg of
duration tof to asteroid B is a Lambert arc; the junction delta-v dv = |v1_lambert - v| must be deliverable by low thrust
within the leg: dv <= eta * a(m) * tof (a = Tmax/m). Propellant: dm = m (1 - exp(-kappa dv / VE)), kappa = 1 + dv/(2 a tof)
(leverage penalty of position-only targeting). Child velocity = Lambert arrival velocity.
Launch: position = Earth, v = v_Earth + v_inf, |v_inf| <= 4 km/s free; excess |v1 - vE| - 4 must be paid by thrust.
"""
import json, numpy as np
from dataclasses import dataclass, field
from .constants import VE, TMAX, T_MISSION, DAY, M_DRY, FUEL_MAX, VINF_MAX
from .kepler import Ephemeris, propagate_twobody
from .lambert import lambert

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

@dataclass
class State:
    t: float; r: np.ndarray; v: np.ndarray; m: float
    visited: int                # bitmask over asteroid index (bit k = ID k+1)
    seq: tuple                  # ((ast_id, t_flyby, dv_est, tof), ...)
    fuel: float; t_launch: float; vinf: np.ndarray; m0: float
    rare: float = 0.0
    def score(self, P):
        tt = self.t if P.abs_time else (self.t - self.t_launch)
        return tt / T_MISSION + P.w_fuel * self.fuel / FUEL_MAX - P.w_rare * self.rare
    def n(self):
        return len(self.seq)

def mask_of(ids):
    m = 0
    for i in ids: m |= 1 << (int(i) - 1)
    return m

def _visited_array(mask, n=300):
    return np.array([(mask >> k) & 1 for k in range(n)], dtype=bool)

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
    if len(idx) == 0:
        return []
    it, ia = idx[:, 0], idx[:, 1]
    v1, v2 = lambert(np.broadcast_to(s.r, (len(it), 3)), R[it, ia], tofs[it])
    dv = np.linalg.norm(v1 - s.v, axis=-1)
    tof = tofs[it]
    feas = np.isfinite(dv) & (dv <= P.eta * a * tof) & (dv <= P.dv_max)
    if not feas.any():
        return []
    it, ia, v2, dv, tof = it[feas], ia[feas], v2[feas], dv[feas], tof[feas]
    kappa = 1 + dv / (2 * a * tof)
    dm = s.m * (1 - np.exp(-kappa * dv / VE))
    ok = (s.m - dm) >= M_DRY + P.m_margin
    it, ia, v2, dv, tof, dm = it[ok], ia[ok], v2[ok], dv[ok], tof[ok], dm[ok]
    if len(ia) == 0:
        return []
    cost = tof / T_MISSION + P.w_fuel * dm / FUEL_MAX
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
    for o in keep:
        k = ia[o]
        children.append(State(t=s.t + tof[o], r=R[it[o], k], v=v2[o], m=s.m - dm[o], visited=s.visited | (1 << int(k)),
                              seq=s.seq + ((int(k) + 1, float(s.t + tof[o]), float(dv[o]), float(tof[o])),),
                              fuel=s.fuel + float(dm[o]), t_launch=s.t_launch, vinf=s.vinf, m0=s.m0,
                              rare=s.rare + (float(P.rarity[k]) if P.rarity is not None else 0.0)))
    return children

def roots(eph, t_launch_grid, m0, excluded, P):
    """Launch children: one per (launch date, target, tof) with |v_inf| <= 4 (+ thrust for the excess)."""
    out = []
    a0 = TMAX / m0 * 1e-3
    vis0 = mask_of(set(excluded) | UNREACHABLE)
    for tL in t_launch_grid:
        rE, vE = eph.earth_state(tL)
        tofs = P.tofs[(tL + P.tofs) <= T_MISSION]
        for tof in tofs:
            R, _ = eph.all_ast_states(tL + tof)
            v1, v2 = lambert(np.broadcast_to(rE, R.shape), R, np.full(R.shape[0], tof))
            vinf = np.linalg.norm(v1 - vE, axis=-1)
            excess = np.maximum(0.0, vinf - VINF_MAX)
            feas = np.isfinite(vinf) & (excess <= P.eta * a0 * tof) & (excess <= P.dv_max) & ~_visited_array(vis0)
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
                                 rare=(float(P.rarity[k]) if P.rarity is not None else 0.0)))
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
            if not children:
                break
            children.sort(key=lambda s: s.score(P))
            seen = set(); nb = []
            for c in children:
                key = (c.visited, int(c.t // (5 * DAY)))
                if key in seen: continue
                seen.add(key); nb.append(c)
                if len(nb) >= P.beam: break
            beam = nb
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
