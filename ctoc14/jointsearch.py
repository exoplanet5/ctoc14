"""Joint multi-spacecraft beam search ("swarm" planner).

A joint state holds K spacecraft states (search.State) sharing ONE visited mask. At every step the active spacecraft with
the earliest current time is expanded (search.expand with the shared mask); an unlaunched spacecraft is expanded through
search.roots over its launch grid. Beam pruning happens among joint states with the same total flyby count, on
score = w_fuel * sum(fuel_i)/1400 + sum(elapsed_i)/T_MISSION.  Spacecraft that cannot expand are retired.
This assigns every target to the spacecraft that catches it cheapest at that moment, instead of the greedy sequential
allocation that leaves the hardest targets to the last spacecraft."""
import numpy as np
from dataclasses import dataclass
from .constants import DAY, T_MISSION, FUEL_MAX
from .search import State, Params, expand, roots, mask_of, UNREACHABLE, tour_json
from .kepler import propagate_twobody
from .constants import M_DRY
WAIT = 120 * DAY   # coast step when a craft has no feasible leg right now (patience)

@dataclass
class Joint:
    craft: tuple            # tuple of State or None (unlaunched)
    active: tuple           # bool per craft
    visited: int            # shared mask
    launch_t: tuple         # earliest launch time per unlaunched craft
    parent: int = 0         # id of the parent joint state (diversity cap)
    def n(self): return sum(s.n() for s in self.craft if s is not None)
    def fuel(self): return sum(s.fuel for s in self.craft if s is not None)
    def score(self, P):
        el = sum((s.t - s.t_launch) for s in self.craft if s is not None) / T_MISSION
        rare = sum(s.rare for s in self.craft if s is not None)
        return P.w_fuel * self.fuel() / FUEL_MAX + el - P.w_rare * rare
    def next_index(self):
        """Index of the active craft with the earliest time (unlaunched craft count with their launch time)."""
        best, bi = None, -1
        for i, (s, act) in enumerate(zip(self.craft, self.active)):
            if not act: continue
            t = s.t if s is not None else self.launch_t[i]
            if best is None or t < best: best, bi = t, i
        return bi

def _replace(j, i, s_new, visited_new):
    craft = list(j.craft); craft[i] = s_new
    return Joint(craft=tuple(craft), active=j.active, visited=visited_new, launch_t=j.launch_t)

def _retire(j, i):
    act = list(j.active); act[i] = False
    return Joint(craft=j.craft, active=tuple(act), visited=j.visited, launch_t=j.launch_t)

def expand_joint(eph, j, P, m0s, launch_grid, excluded_mask):
    """Children of joint state j (exactly one more flyby each), or [] if terminal. Retires craft that cannot move."""
    while True:
        i = j.next_index()
        if i < 0:
            return []
        s = j.craft[i]
        if s is None:
            pool_i = _CTX['roots'][m0s[i]]
            vis = j.visited | excluded_mask
            kids = [c for c in pool_i if c.t_launch >= j.launch_t[i] and not (vis & (1 << (c.seq[0][0] - 1)))]
            if len(kids) > P.max_children:      # keep the cheapest launch legs per target (fuel, then arrival time), then the best overall
                kids.sort(key=lambda c: (c.fuel, c.t)); per = {}; sel = []
                for c in kids:
                    a0 = c.seq[0][0]
                    if per.get(a0, 0) < P.n_per_target: per[a0] = per.get(a0, 0) + 1; sel.append(c)
                kids = sel[:P.max_children]
        else:
            s2 = State(t=s.t, r=s.r, v=s.v, m=s.m, visited=j.visited | excluded_mask, seq=s.seq, fuel=s.fuel,
                       t_launch=s.t_launch, vinf=s.vinf, m0=s.m0, rare=s.rare)
            kids = expand(eph, s2, P)
        if kids:
            out = []
            for c in kids:
                new_id = c.seq[-1][0]
                out.append(_replace(j, i, c, j.visited | (1 << (new_id - 1))))
            return out
        # no feasible leg now: a launched craft with fuel and time left simply coasts (waits) and tries again later
        wait = getattr(P, 'wait', WAIT)
        if s is not None and s.m > M_DRY + P.m_margin + 5.0 and s.t + wait + 15 * DAY < T_MISSION:
            r2, v2 = propagate_twobody(s.r, s.v, np.array([wait]))
            sw = State(t=s.t + wait, r=r2[0], v=v2[0], m=s.m, visited=s.visited, seq=s.seq, fuel=s.fuel,
                       t_launch=s.t_launch, vinf=s.vinf, m0=s.m0, rare=s.rare)
            j = _replace(j, i, sw, j.visited)
            continue
        j = _retire(j, i)

_CTX = {}
def _init(eph, P, m0s, grid, exm, root_pool=None):
    _CTX.update(eph=eph, P=P, m0s=m0s, grid=grid, exm=exm)
    if root_pool is None:
        root_pool = {}
        for m0 in sorted(set(m0s)):
            rp = roots(eph, grid, m0, (), P); rp.sort(key=lambda c: (c.fuel, c.t)); root_pool[m0] = rp
    _CTX['roots'] = root_pool
def _batch(states):
    c = _CTX; out = []
    for pid, j in enumerate(states):
        for k in expand_joint(c['eph'], j, c['P'], c['m0s'], c['grid'], c['exm']):
            k.parent = id(j); out.append(k)
    return out

def joint_beam_search(eph, m0s, launch_grid, P=None, excluded=(), n_proc=8, verbose=True, launch_stagger=0.0, max_children_per_state=None):
    """m0s: list of initial masses (one per spacecraft). launch_grid: array of allowed launch times (s).
    Returns the best terminal/last Joint state (max flybys, then min score)."""
    P = P or Params(); K = len(m0s)
    exm = mask_of(set(excluded) | UNREACHABLE)
    j0 = Joint(craft=(None,) * K, active=(True,) * K, visited=0, launch_t=tuple(i * launch_stagger for i in range(K)))
    beam = [j0]; best = j0
    pool = None
    _init(eph, P, m0s, launch_grid, exm)
    if verbose: print(f'root pool: {[len(v) for v in _CTX["roots"].values()]} launch legs', flush=True)
    if n_proc > 1:
        import multiprocessing as mp
        ctx = mp.get_context('fork'); pool = ctx.Pool(n_proc, initializer=_init, initargs=(eph, P, m0s, launch_grid, exm, _CTX['roots']))
    try:
        for depth in range(1, 400):
            if pool is not None and len(beam) >= 2 * n_proc:
                chunks = [beam[i::n_proc] for i in range(n_proc)]
                res = pool.map(_batch, chunks); children = [c for r in res for c in r]
            else:
                children = _batch(beam)
            if not children:
                break
            children.sort(key=lambda j: j.score(P))
            seen = set(); nb = []; per_parent = {}; cap = getattr(P, 'parent_cap', 0)
            for c in children:
                key = (c.visited, tuple(int((s.t if s is not None else 0) // (10 * DAY)) for s in c.craft))
                if key in seen: continue
                pid = getattr(c, 'parent', None)
                if cap and per_parent.get(pid, 0) >= cap: continue
                seen.add(key); nb.append(c); per_parent[pid] = per_parent.get(pid, 0) + 1
                if len(nb) >= P.beam: break
            beam = nb; top = beam[0]
            coll = getattr(P, 'collect_joint', None)
            if coll is not None: coll.extend(beam)
            if top.n() > best.n() or (top.n() == best.n() and top.score(P) < best.score(P)): best = top
            if verbose and (depth % 5 == 0 or depth < 5):
                act = sum(top.active); launched = sum(s is not None for s in top.craft)
                print(f'depth {depth:3d}: {len(children):6d} children -> beam {len(beam):4d}; n={top.n()} fuel={top.fuel():.0f} kg '
                      f'launched {launched}/{K} active {act}; per craft {[s.n() if s else 0 for s in top.craft]}', flush=True)
    finally:
        if pool is not None: pool.close(); pool.join()
    return best

def joint_tours(j):
    return [tour_json(s) for s in j.craft if s is not None and s.n() > 0]
