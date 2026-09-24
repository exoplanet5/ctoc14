"""Column generation for the fleet (docs/design_colgen.md).

Master problem (RMP): set cover with misses and a fleet-size cap over planned routes (columns)
    min  sum_r c_r y_r + sum_a w_a   s.t.  sum_{r covers a} y_r + w_a >= 1  (a in rows),  sum_r y_r <= N,  0 <= y, w <= 1.
Duals pi_a in [0, 1] (value of covering target a), mu >= 0 (value of one more craft). Reduced cost of a route:
    rc_r = c_r + mu - sum_{a in r} pi_a.
Pricing: the single-craft beam search (search.beam_search) with the duals as per-target prizes in J units
(Params.prize), fuel weighted by the marginal J per kg of the cost model, time weighted by w_t; every beam state of
every depth is a candidate column, its reduced cost is exact from its own fuel and target set.
Integer fleet: dive = fix the route with the best one-step look-ahead LP value, exclude its targets (zero overlap),
re-price the residual problem with one craft fewer, repeat.

Column cost: c_r = cost_sc(m0) with the tank fixed point m0 = 600 + R + s * F(m0), F(m0) = fuel_est * m0 / m0_search
(planned propellant scales with the tank mass for a fixed Delta-v profile): m0 = (600 + R) / (1 - s * fuel_est / m0_search).
"""
import json, time, numpy as np
from dataclasses import dataclass, field
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from scipy.sparse import csr_matrix
from .constants import cost_sc, M_DRY, M0_MAX, DAY, AU, T_MISSION, FUEL_MAX
from .search import Params, beam_search, tour_json, UNREACHABLE

TARGETS = np.array([t for t in range(1, 301) if t not in UNREACHABLE])      # 298 coverable IDs


# ------------------------------------------------------------------------------------------------ cost model
@dataclass
class CostModel:
    s: float = 1.15          # actual/planned propellant factor (conversion + tail), calibrated median ~1.03, p75 ~1.13
    reserve: float = 20.0    # kg: launch-leg allowance + reserve above dry mass

    def tank(self, fuel_est, m0_search):
        k = self.s * np.asarray(fuel_est, float) / np.asarray(m0_search, float)
        with np.errstate(divide='ignore'):
            m0 = (M_DRY + self.reserve) / (1.0 - k)
        return np.where((k < 0.99) & (m0 <= M0_MAX), m0, np.inf)

    def cost(self, fuel_est, m0_search):
        m0 = self.tank(fuel_est, m0_search)
        return np.where(np.isfinite(m0), cost_sc(np.minimum(m0, M0_MAX)), 10.0)

    def w_fuel(self, fuel_typ=360.0, m0_search=1000.0):
        """Beam fuel weight = FUEL_MAX * dc/dF at a typical route (J per full-tank fraction of planned propellant)."""
        h = 1.0
        return float(FUEL_MAX * (self.cost(fuel_typ + h, m0_search) - self.cost(fuel_typ - h, m0_search)) / (2 * h))


# ------------------------------------------------------------------------------------------------ column store
class ColumnStore:
    """All known routes as target-index lists + cost; pool routes are referenced by (file, byte offset), generated
    routes keep their tour JSON inline. Exact-duplicate target sets keep the cheapest route."""
    def __init__(self, cm):
        self.cm = cm
        self.tidx = []; self.cost = []; self.fuel = []; self.m0s = []; self.n = []; self.t_launch = []
        self.src = []            # (file, offset) or ('inline', tour_dict)
        self.tag = []            # origin label
        self.key = {}            # packed target mask (bytes) -> column id
        self._csr = None

    def __len__(self):
        return len(self.cost)

    def add(self, targets, fuel_est, m0_search, t_launch, src, tag):
        """Returns (column id, is_new_or_improved)."""
        tset = sorted(set(int(t) for t in targets if int(t) not in UNREACHABLE))
        if len(tset) == 0:
            return -1, False
        c = float(self.cm.cost(fuel_est, m0_search))
        mk = np.zeros(304, np.uint8); mk[np.array(tset) - 1] = 1
        kb = np.packbits(mk).tobytes()                  # 38-byte key (a frozenset would cost ~1.5 kB per route)
        j = self.key.get(kb)
        if j is not None:
            if c < self.cost[j] - 1e-9:
                self.cost[j] = c; self.fuel[j] = float(fuel_est); self.m0s[j] = float(m0_search)
                self.t_launch[j] = float(t_launch); self.src[j] = src; self.tag[j] = tag
                return j, True
            return j, False
        j = len(self.cost); self.key[kb] = j
        self.tidx.append(np.array([t - 1 for t in tset], np.int16)); self.cost.append(c); self.fuel.append(float(fuel_est))
        self.m0s.append(float(m0_search)); self.n.append(len(tset)); self.t_launch.append(float(t_launch))
        self.src.append(src); self.tag.append(tag); self._csr = None
        return j, True

    def recost(self):
        self.cost = [float(self.cm.cost(f, m)) for f, m in zip(self.fuel, self.m0s)]

    def csr(self):
        """(ncol x 300) 0/1 matrix, cached until the next add."""
        if self._csr is None or self._csr.shape[0] != len(self.cost):
            lens = np.array([len(t) for t in self.tidx]); indptr = np.concatenate([[0], np.cumsum(lens)])
            ind = np.concatenate(self.tidx).astype(np.int32) if len(self.tidx) else np.zeros(0, np.int32)
            self._csr = csr_matrix((np.ones(len(ind), np.float64), ind, indptr), shape=(len(self.cost), 300))
        return self._csr

    def costs(self):
        return np.asarray(self.cost)

    def tour(self, j):
        s = self.src[j]
        if s[0] == 'inline':
            return dict(s[1])
        f, off = s
        with open(f, 'rb') as fh:
            fh.seek(off); return json.loads(fh.readline())

    def targets(self, j):
        return [int(i) + 1 for i in self.tidx[j]]

    def load_jsonl(self, path, tag=None, min_targets=5):
        n0 = len(self); tag = tag or str(path); spath = str(path)
        with open(path, 'rb') as fh:
            off = 0
            for line in fh:
                if line.strip():
                    t = json.loads(line)
                    ids = [l['ast'] for l in t['legs']]
                    if len(set(ids)) >= min_targets:
                        self.add(ids, t['fuel_est'], t.get('m0_search', t['m0']), t['t_launch'],
                                 (spath, off), tag)
                off += len(line)
        return len(self) - n0

    def load_columns_file(self, path):
        """Generated columns (inline tours, see column_record)."""
        n0 = len(self)
        for line in open(path):
            if line.strip():
                try:
                    t = json.loads(line)
                except json.JSONDecodeError:          # partial last line of a file another process is appending to
                    continue
                self.add([l['ast'] for l in t['legs']], t['fuel_est'], t['m0_search'], t['t_launch'], ('inline', t), t.get('tag', 'gen'))
        return len(self) - n0


# ------------------------------------------------------------------------------------------------ master problem
@dataclass
class LPResult:
    value: float            # planned J of the LP incl. fixed routes, 131/144 and nothing else
    y: np.ndarray           # (n_active,) route values
    cols: np.ndarray        # active column ids
    pi: np.ndarray          # (300,) duals, 0 outside rows
    mu: float
    w: np.ndarray           # (nR,) miss values
    rows: np.ndarray        # target indices (0-based) in the rows
    status: int = 0


class RMP:
    """Restricted master problem over a working set of columns; pricing of the whole store by the LP duals."""
    def __init__(self, store, craft_cost=0.0):
        self.store = store
        self.active = np.zeros(0, int)
        self.craft_cost = craft_cost      # fixed cost per used column (penalty homotopy instead of a hard craft cap)

    def seed(self, rows, per_target=30):
        """Initial working set: per target in rows, the `per_target` cheapest-per-target columns containing it."""
        C = self.store.csr(); cost = self.store.costs(); n = np.asarray(C.sum(axis=1)).ravel()
        cpt = cost / np.maximum(n, 1)
        CT = C.T.tocsr(); sel = set(self.active.tolist())
        for a in rows:
            idx = CT.indices[CT.indptr[a]:CT.indptr[a + 1]]
            if len(idx):
                sel.update(idx[np.argsort(cpt[idx])[:per_target]].tolist())
        self.active = np.array(sorted(sel), int)

    def _lp(self, rows, N, cols):
        C = self.store.csr()[cols][:, rows]           # (n, nR)
        n, nR = C.shape
        cost = self.store.costs()[cols] + self.craft_cost
        A_cov = -C.T.tocsr()                          # -sum y - w <= -1
        from scipy.sparse import hstack, vstack, identity
        A1 = hstack([A_cov, -identity(nR, format='csr')])
        A2 = csr_matrix((np.ones(n), (np.zeros(n, int), np.arange(n))), shape=(1, n + nR))
        A = vstack([A1, A2]).tocsr()
        b = np.concatenate([-np.ones(nR), [N]])
        res = linprog(np.concatenate([cost, np.ones(nR)]), A_ub=A, b_ub=b, bounds=(0, 1), method='highs')
        if res.status != 0:
            raise RuntimeError(f'LP failed: {res.message}')
        marg = -res.ineqlin.marginals
        pi = np.zeros(300); pi[rows] = marg[:nR]; mu = float(marg[nR])
        return res.fun, res.x[:n], res.x[n:], pi, mu

    def solve(self, rows, N, fixed_cost=0.0, max_iter=30, add_per_iter=2000, eps=1e-6, verbose=None):
        """LP over the working set, repeatedly adding the most negative reduced-cost columns of the WHOLE store until
        none is left (the LP is then exact over the store). rows: 0-based target indices still to cover."""
        rows = np.asarray(rows, int)
        if len(self.active) == 0:
            self.seed(rows)
        C = self.store.csr(); cost = self.store.costs()
        for it in range(max_iter):
            fun, y, w, pi, mu = self._lp(rows, N, self.active)
            rc = cost + self.craft_cost - C @ pi + mu
            rc[self.active] = np.inf
            neg = np.where(rc < -eps)[0]
            if verbose:
                verbose(f'    master it {it}: LP {fun + fixed_cost + 2:.4f} sum y {y.sum():.2f} mu {mu:.3f} #pi>0.99 {(pi > 0.99).sum()} '
                        f'active {len(self.active)} negative {len(neg)}')
            if len(neg) == 0:
                break
            add = neg[np.argsort(rc[neg])[:add_per_iter]]
            self.active = np.concatenate([self.active, add])
        return LPResult(value=fun + fixed_cost + 2.0 - self.craft_cost * float(y.sum()), y=y, cols=self.active.copy(), pi=pi, mu=mu, w=w, rows=rows)

    def lp_value_fixed(self, rows, N, fixed_cost):
        """Quick LP value (working set only) — for look-ahead."""
        rows = np.asarray(rows, int)
        if len(rows) == 0:
            return fixed_cost + 2.0
        fun, y, w, pi, mu = self._lp(rows, max(N, 0), self.active)
        return fun + fixed_cost + 2.0

    def mip(self, rows, N, time_limit=120.0, fixed_cost=0.0, cols=None):
        """Integer set cover over the working set (or given cols). Returns (J, selected column ids, missed rows)."""
        rows = np.asarray(rows, int); cols = self.active if cols is None else np.asarray(cols, int)
        C = self.store.csr()[cols][:, rows]; n, nR = C.shape; cost = self.store.costs()[cols]
        from scipy.sparse import hstack, vstack, identity
        A1 = hstack([-C.T.tocsr(), -identity(nR, format='csr')])
        A2 = csr_matrix((np.ones(n), (np.zeros(n, int), np.arange(n))), shape=(1, n + nR))
        A = vstack([A1, A2]).tocsr()
        lo = np.full(nR + 1, -np.inf); hi = np.concatenate([-np.ones(nR), [N]])
        res = milp(np.concatenate([cost, np.ones(nR)]), constraints=LinearConstraint(A, lo, hi),
                   integrality=np.concatenate([np.ones(n), np.zeros(nR)]), bounds=Bounds(0, 1),
                   options=dict(time_limit=time_limit, mip_rel_gap=1e-3))
        if res.x is None:
            return np.inf, [], rows.tolist()
        y = np.round(res.x[:n]).astype(int); sel = cols[y > 0]
        cov = np.zeros(300, bool)
        for j in sel: cov[self.store.tidx[j]] = True
        missed = [int(a) for a in rows if not cov[a]]
        J = float(self.store.costs()[sel].sum()) + len(missed) + fixed_cost + 2.0
        return J, sel.tolist(), missed


# ------------------------------------------------------------------------------------------------ pricing
@dataclass
class PricingSpec:
    tag: str = 'A'
    launch: tuple = (0.0, 600.0, 20.0)      # days: start, stop, step of the launch grid
    beam: int = 300
    m0: float = 1000.0                       # search mass (acceleration + planned propellant scale)
    w_t: float = 0.9
    w_fuel: float = None                     # None = cost model derivative
    dv_max: float = 1.2
    lin_tofmax: float = 600.0
    lin_drmax: float = 0.15                  # AU
    vinf_cap: float = 2.0
    m_margin: float = 40.0
    tof_refine: bool = True
    perturb: float = 0.0                     # relative Gaussian noise on the prizes
    first_top: int = 0                       # restrict first targets to the `first_top` highest-prize targets
    nproc: int = 5
    seed: int = 0
    max_cols: int = 400
    jaccard: float = 0.9
    hard_cap: float = 0.0                    # > 0: legs ending at targets with prize >= hard_thresh may use this cap (km/s)
    hard_thresh: float = 0.5


def column_record(s, m0_search, tag, rc=None, rnd=None):
    t = tour_json(s)
    t['m0_search'] = float(m0_search); t['tag'] = tag
    t['r_end'] = [float(x) for x in s.r]; t['v_end'] = [float(x) for x in s.v]; t['m_end'] = float(s.m)
    if rc is not None: t['rc'] = float(rc)
    if rnd is not None: t['round'] = rnd
    return t


def price(eph, store, pi_prize, pi_true, mu, excluded, spec, log=print, rnd=None, eps=0.02):
    """One pricing beam. Returns (list of added column ids, stats dict)."""
    tic = time.time()
    rng = np.random.default_rng(spec.seed)
    prize = np.array(pi_prize, float)
    if spec.perturb > 0:
        prize = np.clip(prize * (1 + spec.perturb * rng.standard_normal(300)), 0, 1)
    ex = set(int(x) for x in excluded) | UNREACHABLE
    for t in ex: prize[t - 1] = 0.0
    wf = spec.w_fuel if spec.w_fuel is not None else store.cm.w_fuel(360.0 * spec.m0 / 1000.0, spec.m0)
    first = None
    if spec.first_top > 0:
        order = [int(i) + 1 for i in np.argsort(-prize) if prize[i] > 0][:spec.first_top]
        first = set(order) if order else None
    P = Params(beam=spec.beam, w_fuel=wf, w_t=spec.w_t, m_margin=spec.m_margin, dv_max=spec.dv_max,
               lin_tofs=np.arange(20, spec.lin_tofmax + 1e-9, 10) * DAY, lin_drmax=spec.lin_drmax * AU, vinf_cap=spec.vinf_cap,
               tof_refine=spec.tof_refine, prize=prize, first_targets=first)
    if spec.hard_cap > 0:
        P.dv_max_t = np.where(prize >= spec.hard_thresh, max(spec.hard_cap, spec.dv_max), spec.dv_max)
    P.collect = []
    grid = np.arange(*spec.launch) * DAY
    best, beam = beam_search(eph, excluded=sorted(ex - UNREACHABLE), m0=spec.m0, t_launch_grid=grid, P=P, n_proc=spec.nproc, verbose=False)
    states = P.collect + list(beam or [])
    # reduced cost of every state w.r.t. the TRUE duals
    best_by_set = {}
    for s in states:
        ids = [x[0] for x in s.seq]
        tset = frozenset(ids)
        if len(tset) < 5: continue
        c = float(store.cm.cost(s.fuel, spec.m0))
        rc = c + mu - float(sum(pi_true[i - 1] for i in tset))
        o = best_by_set.get(tset)
        if o is None or rc < o[0]:
            best_by_set[tset] = (rc, s)
    cands = sorted(best_by_set.values(), key=lambda x: x[0])
    rc_min = cands[0][0] if cands else np.inf
    neg = [(rc, s) for rc, s in cands if rc < -eps]
    # greedy diversity filter (Jaccard among the selected)
    sel = []; sel_masks = []
    for rc, s in neg:
        m = np.zeros(300, bool); m[[x[0] - 1 for x in s.seq]] = True
        if sel_masks:
            M = np.array(sel_masks); inter = (M & m).sum(1); uni = (M | m).sum(1)
            if np.max(inter / uni) > spec.jaccard: continue
        sel.append((rc, s)); sel_masks.append(m)
        if len(sel) >= spec.max_cols: break
    added = []
    for rc, s in sel:
        rec = column_record(s, spec.m0, spec.tag, rc=rc, rnd=rnd)
        j, new = store.add([x[0] for x in s.seq], s.fuel, spec.m0, s.t_launch, ('inline', rec), spec.tag)
        if new: added.append(j)
    bn = max((s.n() for s in states), default=0)
    st = dict(tag=spec.tag, states=len(states), sets=len(cands), negative=len(neg), added=len(added), rc_min=float(rc_min),
              best_n=int(bn), w_fuel=float(wf), runtime=time.time() - tic)
    log(f'  pricing {spec.tag}: {len(states)} states, {len(cands)} sets, {len(neg)} with rc<-{eps}, added {len(added)}; '
        f'rc_min {rc_min:.3f}, deepest {bn}, w_fuel {wf:.2f} ({time.time()-tic:.0f} s)')
    return added, st


# ------------------------------------------------------------------------------------------------ fleet helpers
def fleet_summary(store, cols):
    cov = set()
    for j in cols: cov |= set(store.targets(j))
    sumJ = float(sum(store.cost[j] for j in cols))
    missed = sorted(set(TARGETS.tolist()) - cov)
    return dict(n=len(cols), covered=len(cov), sumJ=sumJ, J=sumJ + len(missed) + 2, missed=missed)


def write_fleet(store, cols, outdir, extra=None):
    import pathlib
    out = pathlib.Path(outdir); out.mkdir(parents=True, exist_ok=True)
    for f in out.glob('tour_sc*.json'): f.unlink()
    order = sorted(cols, key=lambda j: -store.n[j])
    for k, j in enumerate(order, 1):
        t = store.tour(j); t = dict(t)
        t['m0_search'] = float(store.m0s[j]); t['m0'] = float(min(M0_MAX, store.cm.tank(store.fuel[j], store.m0s[j])))
        t['cost_J'] = float(store.cost[j]); t['col'] = int(j); t['tag'] = store.tag[j]
        json.dump(t, open(out / f'tour_sc{k}.json', 'w'), indent=1)
    summ = fleet_summary(store, cols); summ.update(extra or {})
    json.dump(summ, open(out / 'summary.json', 'w'), indent=1)
    return summ
