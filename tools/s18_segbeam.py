"""Stage 18 / WP1: the per-craft SEGMENT BEAM -- one slice of proposals for one craft (or one root beam).

What it is.  The lintwin twin beam of tools/s14_twinbeam.py (state = a SETTLED impulsive twin; children = every
remaining target at its single-leg encounter epochs, priced by the linearised whole-prefix twin; kept children
settled) restricted to a time horizon: candidate encounter epochs in (T, T + H + L].  Every settled state whose last
flyby is <= T + H is a COLUMN (the committed part); states beyond T + H are LOOKAHEAD and only give their column
ancestor its `potential` (the extra flybys a descendant reached in (T + H, T + H + L]).  Lookahead flybys are never
committed.

Two kinds of job (JOB schema below):
  craft  the beam starts from the craft's settled, regridded prefix twin (route npz); its columns belong to that craft;
  root   the beam starts from launch legs (ctoc14.search.roots over launch epochs in the job's `grid`, |v_inf| <= 4, the
         exact zero-thrust Lambert arc as the twin: s14_twinbeam.root_state) whose first flyby is <= T + H + L; its
         columns have craft = None and are interchangeable among the unlaunched craft.  The driver splits the launch
         window (T, T + H] into n_root contiguous sub-grids (one job each) so that the root front is diverse by
         launch date by construction; inside a job the per-level selection is round-robin over `root_bin`-day launch
         bins (diverse_select).

Reused from s14_twinbeam (imported, never edited): params, root_state, search_state, seeds, child_problem,
settle_child, lin_price, twin_score, pack, save_ckpt.  epoch_candidates and _price_parent are COPIED here with the
horizon (t_lo, t_hi] and a per-target prize (the beam score is s15_twintail.pscore: t_end/T + w_fuel F/1400 - sum of
prizes; prize = 1 + edf when the driver passes `prize`, else 1).

Integrator changes (09-24, measured on results/s18/g0 and g0b): (1) BRANCH DIVERSITY `div_cap` (JOB_DEF): at most
div_cap settle tries / kept states per branch = (launch bin, first target) for a root beam, first NEW target for a craft
beam, so that every job offers several mutually compatible columns to the auction (g0 r0: 4-fb best -> 6-fb best,
16 vs 12 disjoint targets); (2) a STRANDED prefix (t_end far before the window) gets a LinLeg tof grid reaching 600 d
into the window (epoch_candidates_h), so an idle craft can always propose again; (3) the tie-break RNG is stored in the
checkpoint (bit-for-bit resume).

THE WORKER IS SINGLE-THREADED: run_job never opens a Pool (TB._pool_map(fn, jobs, 1)); the driver's fork Pool runs
one job per process.  Everything the job needs is in the JOB dict (JSON) and on disk; the result is RESULT (JSON) plus
the column twins col_<k>.npz in the job's out dir.  run_job is idempotent (result.json exists -> returned as is) and
checkpoints the beam every level (ckpt.pkl, s14_twinbeam.save_ckpt) so a killed slice resumes.

JOB (slices/s<kk>/jobs/<job_id>.json; all times SECONDS; targets 1-based)
  job_id      "s03_c2" (craft 2, slice 3) | "s03_r1" (root sub-grid 1)
  kind        "craft" | "root"
  craft       int | None
  route       str | None   path (relative to run_dir) of the prefix twin (kind craft)
  tank_prev   float        honest tank of the prefix (root: 601.5)
  counts_prev {"pin","med","fil"} of the prefix (root: zeros)   [diagnostic; quotas are applied by the auction]
  remaining   [int]        targets the beam may take (the fleet's `remaining`, minus nothing: the auction resolves
                           conflicts between craft)
  T, H, L     float s
  grid        [t_lo, t_hi, step]  s (kind root): launch epochs tL in (t_lo, t_hi], step; the driver gives each root
                           job one contiguous part of (T, T + H]
  roots       int          kind root: max root states kept (s14: 60), after (first target, 20 d launch bin) dedupe
  prize       {"<id>": float} | None   per-target prize in the beam score (default 1.0 for every remaining target)
  params      {B, tries, k_epochs, lin_cap, res_max, iters, timeout, max_levels, wall_budget, root_bin (DAYS),
               w_t, w_fuel (None = s14 production), dv, drmax, mc, price_cap}
  seed        int          (ordering ties only; the beam is otherwise deterministic)
  run_dir     str          absolute run directory
  out         str          job directory relative to run_dir ("slices/s03/c2"): log.txt, ckpt.pkl, col_<k>.npz, result.json

RESULT (out/result.json)
  job_id, kind, craft, ok (bool), error (None | str), columns [COLUMN ...] (s18_common.COLUMN_SCHEMA; npz paths
  relative to run_dir), stats {levels, n_states, n_columns, settles, settle_ok, price_s, settle_s, wall_s, end
  {kind: no_children | settle_wall | horizon | max_levels | wall_budget | pool_exhausted | lookahead, depth}, hist [per level:
  {depth, parents, children, unique, tried, settled, kept, in_H, wall_s}]}

`max_levels` counts EXPANSION ROUNDS: the start states are depth 1, after `max_levels` rounds the deepest state has
depth max_levels + 1 (e.g. --levels 3 from a 4-flyby prefix reaches 7 flybys).

usage: s18_segbeam.py run JOB.json            (single process; prints the RESULT summary)
       s18_segbeam.py test-craft OUT [--route results/s16/best/fleet/route_r1.npz --T 540 --H 540 --L 360 --B 3 --tries 4 --levels 3]
       s18_segbeam.py test-root  OUT [--T 0 --H 540 --L 360 --step 30 --roots 10 --B 3 --tries 4 --levels 3]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pickle, pathlib, argparse, warnings, shutil, traceback
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import numpy as np
import run_ialns as RI
import s14_twinbeam as TB
import s18_common as C
from ctoc14.search import mask_of, UNREACHABLE
from ctoc14.search import roots as search_roots
from ctoc14.lambert import lambert
from ctoc14.linleg import LinLeg
from ctoc14.constants import DAY, VE, T_MISSION, FUEL_MAX, AU, TMAX
warnings.filterwarnings('ignore')

JOB_DEF = dict(B=8, tries=16, k_epochs=3, lin_cap=3.0, res_max=1.5e8, iters=100, timeout=240.0, max_levels=14,
               wall_budget=2400.0, root_bin=60.0, w_t=None, w_fuel=None, dv=1.2, drmax=0.15, mc=60, price_cap=250,
               div_cap=3, look_stop=2, settle_cap=8, settle_ok_min=0.3)
# Fixer 09-24 (reviewer finding 6, MEASURED on results/s18/g0b: 54 % of a deep job's wall went to levels with in_H = 0):
#   look_stop    the beam ends after this many CONSECUTIVE levels whose settled children all lie past T + H (t_end is
#                monotone along a branch, so no further column can appear; the potential of every column is then known
#                to depth + look_stop).  wall_budget (2400 s) stays as the safety only.
#   settle_cap / settle_ok_min   when the previous level settled fewer than settle_ok_min of its tries, this level tries
#                at most settle_cap children (g0b s02_c0: 83 settles, 26 ok, 625 s of settle).
# Reviewer finding 4: a craft job may carry job['quota'] = {'fil': head-room, 'med': head-room} (None = unlimited);
#   when a parent's NEW targets already fill a class's head-room, that class is excluded from its candidates
#   (price_parent's `remaining`), so the beam never proposes what the auction must drop for the quota.
# div_cap: BRANCH DIVERSITY (integrator fix 09-24 01:35, measured on results/s18/g0 slice 0: all 30 columns of root job
# r0 launched on day 0 and shared target 139/174, so the auction found only 2 disjoint offers among 51 columns).  The
# branch key of a state is (launch bin, first target) for a root beam and the first NEW target for a craft beam; at most
# div_cap settle tries and at most div_cap kept states per branch, filled in score order, the rest by score (0 = off).
# The auction can adopt only one column per craft, so columns that share the first new target are near-duplicates
# for it; the cap makes every job offer several mutually compatible columns.
HORIZON_MARGIN = 15 * DAY       # a state with t_end > T + H + L - 15 d cannot get a child inside the horizon
EDGE_LO = 5 * DAY               # candidate epochs start at T + 5 d: the settle moves flyby epochs by a few days, and a
                                # settled child whose new flyby left (T, T + H + L] is rejected (would fail check_column)
PRICE_SLOW_S = 180.0            # a level whose pricing exceeded this -> k_epochs 2 for parents deeper than 25
DEEP_PARENT = 25
DEDUP_BIN = 5 * DAY             # (visited set, 5 d epoch bin) dedupe of children and survivors (s14 / s15)
ROOT_DEDUP_BIN = 20 * DAY       # (first target, 20 d launch bin) dedupe of launch legs (s14 / s15)


# ------------------------------------------------------------------ pricing with a horizon
def epoch_candidates_h(par, remaining, t_lo, t_hi, K=3, sep=30 * DAY):
    """Horizon-limited copy of s14_twinbeam.epoch_candidates: the K lowest local minima (>= sep apart) of the single-leg
    cost curves from the twin end state (par['t_end'], r_end, v_end) -- single-rev Lambert junction over 15-400 d / 5 d
    and LinLeg L1 cost over 20-600 d / 10 d -- but only epochs t with t_lo < t <= t_hi (and t <= T_MISSION - 2 d).
    Returns [(ast, t [s], single-leg cost [km/s])].  `remaining` are 1-based ids.
    The tof grids are masked BEFORE the Lambert / LinLeg batches (this is what makes a slice cheap).  The grids keep
    ONE point beyond each horizon edge so that a cut curve does not fake a local minimum at the edge (the settle would
    then slide the flyby out of the horizon: MEASURED 2 of 4 drifts below T in the first acceptance run); the minima
    returned are all inside (t_lo, t_hi]."""
    E = RI.eph(); t0 = float(par['t_end']); r0 = np.asarray(par['r_end'], float); v0 = np.asarray(par['v_end'], float)
    idx = np.array(sorted(set(int(x) for x in remaining)), int) - 1
    if len(idx) == 0:
        return []
    t_end_m = T_MISSION - 2 * DAY; t_max = min(float(t_hi), t_end_m); t_lo = float(t_lo)
    curves = {int(i) + 1: [] for i in idx}
    tl = np.arange(15, 400 + 1e-9, 5) * DAY; tl = tl[(t0 + tl > t_lo - 5 * DAY) & (t0 + tl <= min(t_max + 5 * DAY, t_end_m))]
    # A STRANDED prefix (t_end far before the window, e.g. a craft whose column was rejected or that got no column for a
    # slice) needs a long coast: the LinLeg grid then reaches 600 d INTO the window instead of 600 d past t_end
    # (integrator fix 09-24: results/s18/g0b s02_c0 had 0 candidates with t_end 444 d and the window at 1085 d).
    # The single-rev Lambert grid is not extended (multi-rev arcs are not modelled); the L1 LinLeg cost covers it.
    tn_max = max(600 * DAY, (t_lo - t0) + 600 * DAY)
    tn = np.arange(20, tn_max / DAY + 1e-9, 10) * DAY; tn = tn[(t0 + tn > t_lo - 10 * DAY) & (t0 + tn <= min(t_max + 10 * DAY, t_end_m))]
    nX = len(idx)
    if len(tl):
        R = E.ast_states_at(idx[None, :], (t0 + tl)[:, None])[0]
        v1, _ = lambert(np.broadcast_to(r0, (len(tl) * nX, 3)), R.reshape(-1, 3), np.repeat(tl, nX))
        dv = np.linalg.norm(v1 - v0, axis=-1).reshape(len(tl), nX)
        dv = np.where(np.isfinite(dv), dv, np.inf)
        for q in range(nX):
            curves[int(idx[q]) + 1].append((tl, dv[:, q]))
    if len(tn):
        L = LinLeg(r0, v0, tn, dtau=10 * DAY); acc = 50.0 * TMAX / TB.M0P * 1e-3       # loose bound: epochs only
        Cc = np.full((len(tn), nX), np.inf)
        for j in range(len(tn)):
            R = E.ast_states_at(idx, np.full(nX, t0 + tn[j]))[0]
            c, _, _ = L.solve(j, R, acc, ub=0.8); Cc[j] = c
        for q in range(nX):
            curves[int(idx[q]) + 1].append((tn, Cc[:, q]))
    out = []
    for X, cl in curves.items():
        mins = []
        for tt, cc in cl:
            for i in range(len(tt)):
                if not np.isfinite(cc[i]):
                    continue
                lo = cc[i - 1] if i > 0 else np.inf; hi = cc[i + 1] if i + 1 < len(tt) else np.inf
                if cc[i] <= lo and cc[i] <= hi and t_lo < t0 + tt[i] <= t_max:
                    mins.append((float(cc[i]), float(t0 + tt[i])))
        mins.sort(); keep = []
        for c, t in mins:
            if all(abs(t - u) >= sep for _, u in keep):
                keep.append((c, t))
            if len(keep) >= K:
                break
        out += [(X, t, c) for c, t in keep]
    return out


def prize_of(par, prize):
    """Sum of the prizes of the targets a twin state has flown."""
    return float(sum(float(prize[int(a) - 1]) for a in par['asts']))


def price_parent_full(job):
    """price_parent with the pricing statistics: (candidates, info); info = {n_epochs (candidates inside the horizon),
    n_priced (after the price_cap cut), n_ok (finite, <= lin_cap, residual <= res_max), secs}."""
    par, remaining, t_lo, t_hi, K, cap, wt, wf, res_max, prize = job[:10]
    price_cap = int(job[10]) if len(job) > 10 and job[10] else 0
    tic = time.time()
    flown = set(int(a) for a in par['asts'])
    rem = [int(x) for x in remaining if int(x) not in flown]
    ec = epoch_candidates_h(par, rem, t_lo, t_hi, K)
    n_ep = len(ec)
    if price_cap and n_ep > price_cap:
        ec = sorted(ec, key=lambda z: (z[2], z[1], z[0]))[:price_cap]
    ppar = prize_of(par, prize); out = []
    for X, t, c1 in ec:
        try:
            dvm, lin, cm = TB.lin_price(par, X, t)
        except Exception:
            continue
        if not np.isfinite(dvm) or dvm > cap or not (lin['res_km'] <= res_max):
            continue
        F = TB.M0P * (1.0 - np.exp(-(par['dv'] + max(dvm, 0.0)) / VE))
        sc = wt * t / T_MISSION + wf * F / FUEL_MAX - (ppar + float(prize[int(X) - 1]))
        out.append((sc, int(X), float(t), float(dvm), lin, float(cm), float(c1)))
    return out, dict(n_epochs=n_ep, n_priced=len(ec), n_ok=len(out), secs=round(time.time() - tic, 2))


def price_parent(job):
    """All lintwin candidates of one parent inside the horizon, with prizes:
    job = (par, remaining, t_lo, t_hi, K, lin_cap, w_t, w_fuel, res_max, prize (len-300 array)[, price_cap])
    -> [(score, ast, t, marginal dv, lin step, coast miss AU, single-leg cost)] where
    score = w_t t/T_MISSION + w_fuel F/1400 - (prize sum of parent + prize[ast]), F = M0P (1 - exp(-(dv_par + dvm)/ve)).
    Candidates with non-finite dvm, dvm > lin_cap or lin residual > res_max are dropped (s14 _price_parent); with a
    price_cap only the `price_cap` cheapest single-leg candidates are priced by the linearised twin."""
    return price_parent_full(job)[0]


def pscore(ts_, w_t, w_fuel, prize):
    """s15_twintail.pscore: the beam score of a settled twin state with prizes (lower is better)."""
    F = TB.M0P * (1.0 - np.exp(-float(ts_['dv']) / VE))
    return float(w_t) * float(ts_['t_end']) / T_MISSION + float(w_fuel) * F / FUEL_MAX - prize_of(ts_, prize)


def beam_weights(prm):
    """(P, w_t, w_fuel): the s14 production planner parameters and the score weights (None = production)."""
    P = TB.params(float(prm.get('dv', JOB_DEF['dv'])), float(prm.get('drmax', JOB_DEF['drmax'])), int(prm.get('mc', JOB_DEF['mc'])))
    w_t = float(P.w_t if prm.get('w_t') is None else prm['w_t'])
    w_fuel = float(P.w_fuel if prm.get('w_fuel') is None else prm['w_fuel'])
    P.w_t = w_t; P.w_fuel = w_fuel
    return P, w_t, w_fuel


# ------------------------------------------------------------------ start states
def prefix_state(route_npz):
    """Twin STATE dict (s14 pack format) of a committed prefix: C.load_twin -> C.twin_of -> C.pack_state (nreg = the
    regular nodes present, see s18_common.nreg_of).  The prefix is regridded and settled by the driver, so no settle here;
    assert miss <= 150 km (else raise: the driver's commit was dishonest)."""
    st = C.load_twin(route_npz); ip = C.twin_of(st); miss = C.miss_of(ip)
    if not miss <= C.TOL_MISS:
        raise ValueError(f'prefix {route_npz}: miss {miss:.1f} km > {C.TOL_MISS} (not a settled twin)')
    return C.pack_state(ip, miss, 'prefix')


def root_states(grid, remaining, t_first_max, n_keep, w_t=None, w_fuel=None, prize=None, dv=None, drmax=None, mc=None):
    """Launch legs as settled twins.  ctoc14.search.roots(E, launch epochs, M0P, excluded = everything not in
    `remaining`, P = TB.params(...) with P.prize) -> keep states with first flyby <= t_first_max, sort by score, dedupe on
    (first target, 20 d launch bin), keep n_keep, and turn each into TB.root_state(tL, vinf, ast, t1) (exact zero-thrust
    Lambert twin, miss ~0).  Returns [twin state] with how = 'root'.
    grid = (t_lo, t_hi, step) in SECONDS: launch epochs t_lo + step, ..., <= t_hi."""
    E = RI.eph()
    P, w_t, w_fuel = beam_weights(dict(w_t=w_t, w_fuel=w_fuel, dv=JOB_DEF['dv'] if dv is None else dv,
                                       drmax=JOB_DEF['drmax'] if drmax is None else drmax, mc=JOB_DEF['mc'] if mc is None else mc))
    rem = set(int(x) for x in remaining)
    if prize is None:
        prize = np.zeros(300)
        for t in rem:
            prize[t - 1] = 1.0
    P.prize = np.asarray(prize, float)
    t_lo, t_hi, step = [float(x) for x in grid]
    launches = np.arange(t_lo + step, t_hi + 1e-6, step)
    launches = launches[(launches > t_lo) & (launches <= t_hi) & (launches < T_MISSION - 20 * DAY)]
    if len(launches) == 0 or not rem:
        return []
    excl = sorted(set(C.TARGETS) - rem)
    rts = search_roots(E, launches, TB.M0P, excl, P)
    rts = [s for s in rts if s.t <= float(t_first_max)]
    rts.sort(key=lambda s: (s.score(P), s.t, s.t_launch))
    seen = set(); keep = []
    for s in rts:
        key = (int(s.seq[0][0]), int(s.t_launch // ROOT_DEDUP_BIN))
        if key in seen:
            continue
        seen.add(key); keep.append(s)
        if len(keep) >= int(n_keep):
            break
    out = []
    for s in keep:
        ts_ = TB.root_state(float(s.t_launch), np.asarray(s.vinf, float), int(s.seq[0][0]), float(s.seq[0][1]))
        if ts_['miss'] <= C.TOL_MISS:
            out.append(ts_)
    return out


def diverse_select(states, B, key, score):
    """Keep at most B states, round-robin over the groups key(state) in order of `score` (lower better) inside each
    group: the best of every group first, then the second best of every group, ...  Used for the root front (key = launch
    epoch // root_bin) so that one launch date cannot fill the beam.  Groups are visited in order of their best score."""
    groups = {}
    for s in states:
        groups.setdefault(key(s), []).append(s)
    for g in groups.values():
        g.sort(key=score)
    order = sorted(groups, key=lambda k: score(groups[k][0]))
    out = []; rank = 0
    while len(out) < int(B):
        got = False
        for k in order:
            if rank < len(groups[k]):
                out.append(groups[k][rank]); got = True
                if len(out) >= int(B):
                    break
        if not got:
            break
        rank += 1
    return out


def capped_select(items, n, key, score, cap):
    """The best `n` items by `score` (lower better) with at most `cap` per group key(item): pass 1 takes items in score
    order while their group has fewer than cap, pass 2 fills up to n from the rest in score order.  cap <= 0: plain
    best n.  Used for the settle order of the children and for the kept beam (branch diversity, see JOB_DEF)."""
    order = sorted(items, key=score)
    if int(cap) <= 0 or n <= 0:
        return order[:int(n)]
    cnt = {}; first = []; rest = []
    for it in order:
        k = key(it)
        if cnt.get(k, 0) < int(cap):
            cnt[k] = cnt.get(k, 0) + 1; first.append(it)
        else:
            rest.append(it)
    out = first[:int(n)]
    if len(out) < n:
        out += rest[:int(n) - len(out)]
    return out


# ------------------------------------------------------------------ the beam
def _fmt_range(vals, f):
    return f'{f(min(vals))}-{f(max(vals))}' if vals else '-'


def beam(start, remaining, T, H, L, prm, prize, say, out, seed=0, kind='craft', quota=None, cls=None):
    """Depth-synchronous lintwin beam from the start states over `remaining`, epochs in (T, T + H + L].

    Per level: price every beam state (price_parent, single process), sort children by score, dedupe on (visited set,
    5 d bin of the epoch), settle in rank order (TB.settle_child with the lin step as first seed) until B succeed or
    `tries` are spent; the survivors are re-scored by pscore, deduped, and the next beam is the best B (craft) or
    diverse_select over launch bins (root).  EVERY settled state is recorded in `states` with id, parent id, level,
    t_end, tank, miss, asts, tfs, and its ist dict; states with t_end <= T + H are column candidates.
    Stops when no child, when the settle wall is hit, when every beam state has t_end > T + H + L - 15 d (horizon),
    at max_levels, at wall_budget (the beam then still returns what it has), or when `remaining` is exhausted.
    Checkpoint every level: out/ckpt.pkl = dict(level, beam, states, hist, end, wall, settles, settle_ok); resumes.
    Returns dict(states=[...], hist=[...], end={...}, settles, settle_ok, wall, price_s, settle_s, levels).
    NOTE the state list can hold ~tries x levels twins (each ist dict ~100 kB); fine for B <= 8, levels <= 14."""
    out = pathlib.Path(out); out.mkdir(parents=True, exist_ok=True)
    prm = {k: prm.get(k, JOB_DEF[k]) for k in JOB_DEF}
    P, w_t, w_fuel = beam_weights(prm)
    prize = np.asarray(prize, float)
    B = int(prm['B']); tries = int(prm['tries']); K0 = int(prm['k_epochs']); max_levels = int(prm['max_levels'])
    wall_budget = float(prm['wall_budget']); root_bin = float(prm['root_bin']) * DAY; price_cap = int(prm['price_cap'] or 0)
    T = float(T); H = float(H); L = float(L); t_hi = T + H + L
    rem = sorted(set(int(x) for x in remaining)); rem_set = set(rem)
    excl_mask = mask_of((set(C.TARGETS) - rem_set) | set(UNREACHABLE))
    start_asts = set(int(a) for s in start for a in s['asts']) if kind == 'craft' else set()
    rng = np.random.RandomState(int(seed))
    key_launch = lambda s: int(float(s['st']['tL']) // root_bin)
    sc_of = lambda s: pscore(s, w_t, w_fuel, prize)
    div_cap = int(prm.get('div_cap') or 0)
    look_stop = int(prm.get('look_stop') or 0); settle_cap = int(prm.get('settle_cap') or 0)
    settle_ok_min = float(prm.get('settle_ok_min') or 0.0)
    quota = {k: int(v) for k, v in (quota or {}).items() if v is not None}
    cls = cls or {}

    def rem_for(s):
        """The parent's admissible targets: `rem` minus every class whose head-room the parent's new targets fill."""
        if not quota:
            return rem
        new = [int(a) for a in s['asts'] if int(a) not in start_asts]
        full = {k for k, q in quota.items() if sum(1 for a in new if cls.get(a) == k) >= q}
        return [t for t in rem if cls.get(t) not in full] if full else rem

    def branch_of(s):
        """(launch bin, first target) of a root state; the first NEW target of a craft state (None for the prefix)."""
        new = [int(a) for a in s['asts'] if int(a) not in start_asts]
        if kind == 'root':
            return (key_launch(s), new[0] if new else None)
        return new[0] if new else None

    def kid_branch(kd):
        b = branch_of(kd[2])
        return b if (kind == 'root' or b is not None) else int(kd[3])
    ck = out / 'ckpt.pkl'
    if ck.exists():
        with open(ck, 'rb') as f:
            S = pickle.load(f)
        say(f'RESUME level {S["level"]}: beam {len(S["beam"])}, states {len(S["states"])}, wall so far {S["wall"]:.0f} s')
        if S.get('rng') is not None:
            rng.set_state(S['rng'])              # tie-break stream continues where the checkpoint left it (bit-for-bit resume)
    else:
        states = []
        for i, s in enumerate(start):
            s = dict(s); s.update(id=i, parent=None, level=1, score=sc_of(s)); states.append(s)
        S = dict(level=1, beam=list(states), states=states, hist=[], end=None, wall=0.0, settles=0, settle_ok=0,
                 price_s=0.0, settle_s=0.0, next_id=len(states), last_price_s=0.0)
        if not states:
            S['end'] = dict(kind='no_children', depth=0)
        elif not rem_set:
            S['end'] = dict(kind='pool_exhausted', depth=1)
        say(f'start: {len(states)} {kind} state(s), remaining {len(rem)}, horizon ({C.s2d(T):.0f}, {C.s2d(T + H):.0f}] + L '
            f'{C.s2d(L):.0f} d, B {B} tries {tries} k_epochs {K0} price_cap {price_cap} max_levels {max_levels} wall {wall_budget:.0f} s')
        TB.save_ckpt(out, S)
    while S['end'] is None:
        if S['wall'] >= wall_budget:
            S['end'] = dict(kind='wall_budget', depth=S['level']); break
        if S['level'] - 1 >= max_levels:
            S['end'] = dict(kind='max_levels', depth=S['level']); break
        tic = time.time(); bm = S['beam']; depth = S['level'] + 1
        parents = [s for s in bm if s['t_end'] <= t_hi - HORIZON_MARGIN]
        if not parents:
            S['end'] = dict(kind='horizon', depth=S['level']); break
        # ---- price (single process); k_epochs 2 for deep parents when the last level's pricing was slow
        slow = S.get('last_price_s', 0.0) > PRICE_SLOW_S
        jobs = []
        for s in parents:
            K = 2 if (slow and len(s['asts']) > DEEP_PARENT) else K0
            jobs.append((s, rem_for(s), T + EDGE_LO, t_hi, K, float(prm['lin_cap']), w_t, w_fuel, float(prm['res_max']), prize, price_cap))
        if slow:
            say(f'  pricing was slow last level ({S["last_price_s"]:.0f} s): k_epochs 2 for parents deeper than {DEEP_PARENT}')
        res = TB._pool_map(price_parent_full, jobs, 1)
        kids = []; n_ep = 0; n_pr = 0
        for s, (lst, info) in zip(parents, res):
            n_ep += info['n_epochs']; n_pr += info['n_priced']
            vm = excl_mask | mask_of(s['asts'])
            for sc, X, t, dvm, lin, cm, c1 in lst:
                kids.append((sc, rng.random_sample(), s, X, t, vm | (1 << (X - 1)), lin, dvm))
        t_price = time.time() - tic; S['price_s'] += t_price; S['last_price_s'] = t_price
        if n_pr < n_ep:
            say(f'  price_cap: {n_pr} of {n_ep} horizon candidates priced ({len(parents)} parents)')
        kids.sort(key=lambda x: (x[0], x[1]))
        seen = set(); cand = []
        for kd in kids:
            key = (kd[5], int(kd[4] // DEDUP_BIN))
            if key in seen:
                continue
            seen.add(key); cand.append(kd)
        if div_cap > 0:                          # branch diversity of the settle order (JOB_DEF div_cap)
            cand = capped_select(cand, len(cand), kid_branch, lambda x: (x[0], x[1]), div_cap)
        # ---- settle in rank order until B succeed or tries spent (one at a time: the worker is single-threaded)
        new = []; tried = 0; nk = len(kids); nu = len(cand); t_settle = 0.0; fails = []; n_drift = 0
        tries_lvl = tries
        if settle_cap > 0 and S['hist']:
            h0 = S['hist'][-1]
            if h0.get('tried', 0) > 0 and h0.get('settled', 0) < settle_ok_min * h0['tried']:
                tries_lvl = min(tries, settle_cap)
                say(f'  settle cap {tries_lvl}: last level settled {h0.get("settled", 0)}/{h0["tried"]}')
        while cand and len(new) < B and tried < tries_lvl:
            sc, tb, par, X, t, vm, lin, dvm = cand.pop(0); tried += 1
            tw, info = TB._pool_map(TB.settle_child, [(par, X, t, int(prm['iters']), float(prm['timeout']), lin)], 1)[0]
            S['settles'] += 1; t_settle += float(info.get('dt_s', 0.0))
            if tw is None:
                fails.append((X, round(C.s2d(t)), info.get('tries')))
                continue
            # the settle moves every flyby epoch: a new flyby that left (T, T + H + L] would fail check_column
            new_eps = [float(e) for a, e in zip(tw['asts'], tw['tfs']) if int(a) not in start_asts]
            if any(not (T < e <= t_hi) for e in new_eps):
                fails.append((X, round(C.s2d(t)), 'drift', [round(C.s2d(e), 1) for e in new_eps if not (T < e <= t_hi)]))
                n_drift += 1
                continue
            S['settle_ok'] += 1
            tw = dict(tw); tw.update(id=S['next_id'], parent=par['id'], level=depth, score=sc_of(tw), dv_lin=round(dvm, 4))
            S['next_id'] += 1; new.append(tw); S['states'].append(tw)
        S['settle_s'] += t_settle
        new.sort(key=lambda s: s['score'])
        seen = set(); nb = []
        for s in new:
            key = (mask_of(s['asts']), int(s['t_end'] // DEDUP_BIN))
            if key in seen:
                continue
            seen.add(key); nb.append(s)
        if div_cap > 0:
            nb = capped_select(nb, B, branch_of, lambda s: s['score'], div_cap)
        else:
            nb = diverse_select(nb, B, key_launch, lambda s: s['score']) if kind == 'root' else nb[:B]
        n_branch = len(set(branch_of(s) for s in nb))
        in_H = sum(1 for s in new if s['t_end'] <= T + H)
        rec = dict(depth=depth, parents=len(parents), children=nk, unique=nu, tried=tried, settled=len(new), kept=len(nb),
                   in_H=in_H, wall_s=round(time.time() - tic, 1), price_s=round(t_price, 1), settle_s=round(t_settle, 1),
                   n_epochs=n_ep, n_priced=n_pr, drift=n_drift, branches=n_branch, fails=fails[:8])
        S['hist'].append(rec)
        say(f'depth {depth:2d}: {len(parents)}/{len(bm)} parents -> {nk} children ({nu} unique) -> tried {tried}, settled '
            f'{len(new)} ({in_H} in H' + (f', {n_drift} drifted out' if n_drift else '') + f'), kept {len(nb)} ({n_branch} branches)'
            + (f'; t_end {_fmt_range([s["t_end"] for s in nb], lambda v: f"{C.s2d(v):.0f}")} d, tank '
               f'{_fmt_range([s["tank"] for s in nb], lambda v: f"{v:.0f}")}' if nb else '')
            + f' ({rec["wall_s"]:.0f} s: price {t_price:.0f}, settle {t_settle:.0f})')
        S['n_noH'] = (S.get('n_noH', 0) + 1) if (in_H == 0 and len(new) > 0) else 0
        if not nb:
            S['end'] = dict(kind='no_children' if nk == 0 else 'settle_wall', depth=S['level'], tried=tried)
        elif any(rem_set <= set(s['asts']) for s in nb):
            S['end'] = dict(kind='pool_exhausted', depth=depth)
        elif look_stop > 0 and S['n_noH'] >= look_stop:
            S['end'] = dict(kind='lookahead', depth=depth, n_noH=S['n_noH'])
        S['beam'] = nb if nb else bm; S['level'] = depth if nb else S['level']
        S['wall'] += time.time() - tic
        S['rng'] = rng.get_state()
        TB.save_ckpt(out, S)
    say(f'beam end {S["end"]}: {len(S["states"])} settled states, {S["settles"]} settles ({S["settle_ok"]} ok), '
        f'wall {S["wall"]:.0f} s (price {S["price_s"]:.0f}, settle {S["settle_s"]:.0f})')
    return dict(states=S['states'], hist=S['hist'], end=S['end'], settles=S['settles'], settle_ok=S['settle_ok'],
                wall=S['wall'], price_s=S['price_s'], settle_s=S['settle_s'], levels=S['level'] - 1)


def columns_of(res, T, H, L, craft, tank_prev, W, run_dir, out_rel, job_id, prefix_asts=()):
    """Turn the beam's settled states into COLUMNs (s18_common.new_column / COLUMN_SCHEMA):
    - a state is a column iff its last flyby <= T + H and it has >= 1 new target (the prefix itself is not a column: the
      auction's "choose nothing" is the empty choice);
    - potential(col) = max over descendants d of col (parent chain) of n(d) - n(col) counting only flybys in
      (T + H, T + H + L]; 0 when no descendant went past T + H;
    - dedupe columns on the set of new targets (keep the lightest tank; the potential is the max over the duplicates,
      they are the same committed proposal);
    - write each column's twin to run_dir/out_rel/col_<k>.npz (C.save_twin) and set npz = "<out_rel>/col_<k>.npz";
    - counts = C.count_classes(W, new targets); epochs from the twin in flyby order; check_column(c, T, H) must pass.
    Returns the list sorted by (-n_new, tank)."""
    T = float(T); H = float(H); L = float(L); pre = set(int(a) for a in prefix_asts)
    states = res['states']; by_id = {s['id']: s for s in states}
    kids = {}
    for s in states:
        if s['parent'] is not None:
            kids.setdefault(s['parent'], []).append(s['id'])

    def look(s):
        return sum(1 for t in s['tfs'] if T + H < float(t) <= T + H + L)

    pot = {}

    def potential(sid):
        if sid in pot:
            return pot[sid]
        best = look(by_id[sid])
        for c in kids.get(sid, []):
            best = max(best, potential(c))
        pot[sid] = best
        return best
    best = {}
    for s in states:
        if float(s['t_end']) > T + H + 1e-6:
            continue
        new = [int(a) for a in s['asts'] if int(a) not in pre]
        if not new:
            continue
        key = frozenset(new); p = potential(s['id'])
        if key not in best:
            best[key] = (s, p)
        else:
            s0, p0 = best[key]
            keep = s if (float(s['tank']), float(s['t_end'])) < (float(s0['tank']), float(s0['t_end'])) else s0
            best[key] = (keep, max(p, p0))
    cols = []
    od = pathlib.Path(run_dir) / out_rel
    ordered = sorted(best.values(), key=lambda z: (-len(z[0]['asts']), float(z[0]['tank']), float(z[0]['t_end'])))
    for k, (s, p) in enumerate(ordered):
        st = s['st']; o = C.order_of(st)
        all_t = [int(st['asts'][i]) for i in o]; all_e = [float(st['tf'][i]) for i in o]
        new = [(a, e) for a, e in zip(all_t, all_e) if a not in pre]
        rel = f'{out_rel}/col_{k}.npz'
        C.save_twin(od / f'col_{k}.npz', st)
        c = C.new_column(id=f'{job_id}_k{k}', craft=craft, targets=[a for a, _ in new], epochs=[e for _, e in new],
                         all_targets=all_t, npz=rel, tank=float(s['tank']), tank_prev=float(tank_prev), miss=float(s['miss']),
                         t_launch=float(st['tL']), t_end=float(s['t_end']), potential=int(p),
                         counts=C.count_classes(W, [a for a, _ in new]), depth=int(s['level']), score=float(s['score']),
                         job_id=job_id)
        C.check_column(c, T, H)
        cols.append(c)
    cols.sort(key=lambda c: (-c['n_new'], c['tank']))
    return cols


def run_job(job):
    """Execute one JOB (dict) in THIS process, single-threaded; write out/result.json; return the RESULT dict.
    Idempotent: an existing result.json is returned unchanged.  Any exception is caught and reported as ok False with
    the traceback text in `error` (the driver then treats the craft as "no proposal this slice").
    Steps: set_threads; RI.eph(); W = C.load_windows(); start = [prefix_state(route)] (craft) or root_states(...) (root);
    prize vector from job['prize'] (default 1.0 on remaining); res = beam(...); columns = columns_of(...); RESULT."""
    rd = pathlib.Path(job['run_dir']); od = rd / job['out']; od.mkdir(parents=True, exist_ok=True)
    rf = od / 'result.json'; retried = False
    if rf.exists():
        R0 = C.jload(rf)
        if R0.get('ok') or R0.get('retried'):
            R0['cached'] = True
            return R0
        # fixer 09-24 (reviewer finding 5): a cached FAILURE is re-run once (a corrupt checkpoint is discarded first)
        retried = True
        ck = od / 'ckpt.pkl'
        if ck.exists():
            try:
                with open(ck, 'rb') as f:
                    pickle.load(f)
            except Exception:
                ck.unlink()
        try:
            rf.rename(od / 'result_failed.json')
        except Exception:
            pass
    C.set_threads()
    say = C.logger(od / 'log.txt'); tic = time.time()
    kind = job['kind']; craft = job.get('craft'); job_id = job['job_id']
    if retried:
        say(f'job {job_id}: cached FAILURE -> re-running once (result_failed.json kept)')
    R = dict(job_id=job_id, kind=kind, craft=craft, ok=False, error=None, columns=[], stats={}, retried=retried)
    try:
        RI.eph(); W = C.load_windows()
        prm = {k: (job.get('params') or {}).get(k, JOB_DEF[k]) for k in JOB_DEF}
        T = float(job['T']); H = float(job['H']); L = float(job['L'])
        remaining = sorted(set(int(x) for x in job['remaining']) - set(UNREACHABLE))
        prize = np.zeros(300); pz = job.get('prize') or {}
        for t in remaining:
            prize[t - 1] = float(pz.get(str(t), pz.get(t, 1.0)))
        say(f'job {job_id} ({kind}, craft {craft}): T {C.s2d(T):.0f} d H {C.s2d(H):.0f} L {C.s2d(L):.0f}; remaining {len(remaining)}; '
            f'prize >1 on {int(np.sum(prize > 1.0))}, tank_prev {float(job.get("tank_prev", C.TANK_DRY)):.1f}')
        P, w_t, w_fuel = beam_weights(prm)
        if kind == 'craft':
            route = job['route']; rp = pathlib.Path(route)
            start = [prefix_state(rp if rp.is_absolute() else rd / rp)]
            prefix_asts = list(start[0]['asts'])
            say(f'prefix {route}: {len(prefix_asts)} flybys, t_end {C.s2d(start[0]["t_end"]):.1f} d, tank {start[0]["tank"]:.1f} kg '
                f'(planner view), miss {start[0]["miss"]:.1f} km, nreg {start[0]["nreg"]}')
        else:
            tic_r = time.time()
            start = root_states(job['grid'], remaining, T + H + L, int(job.get('roots', 60)), w_t, w_fuel, prize,
                                prm['dv'], prm['drmax'], prm['mc'])
            prefix_asts = []
            g = job['grid']
            say(f'roots: grid ({C.s2d(g[0]):.0f}, {C.s2d(g[1]):.0f}] step {C.s2d(g[2]):.0f} d -> {len(start)} launch legs kept '
                f'({time.time() - tic_r:.1f} s)' + (f'; launch {_fmt_range([float(s["st"]["tL"]) for s in start], lambda v: f"{C.s2d(v):.0f}")} d,'
                                                  f' first flyby {_fmt_range([s["t_end"] for s in start], lambda v: f"{C.s2d(v):.0f}")} d' if start else ''))
        quota = job.get('quota') or None
        cls = {t: C.class_of(W, t) for t in remaining} if quota else None
        if quota:
            say(f'quota head-room {quota} (targets of a full class are excluded from a parent that filled it)')
        res = beam(start, remaining, T, H, L, prm, prize, say, od, seed=int(job.get('seed', 0)), kind=kind, quota=quota, cls=cls)
        cols = columns_of(res, T, H, L, craft, float(job.get('tank_prev', C.TANK_DRY)), W, rd, job['out'], job_id, prefix_asts)
        R['columns'] = cols; R['ok'] = True
        R['stats'] = dict(levels=res['levels'], n_states=len(res['states']), n_columns=len(cols), settles=res['settles'],
                          settle_ok=res['settle_ok'], price_s=round(res['price_s'], 1), settle_s=round(res['settle_s'], 1),
                          wall_s=round(time.time() - tic, 1), end=res['end'], hist=res['hist'])
        say(f'RESULT {job_id}: {len(cols)} columns (n_new {_fmt_range([c["n_new"] for c in cols], str)}, potential '
            f'{_fmt_range([c["potential"] for c in cols], str)}), {len(res["states"])} states, end {res["end"]["kind"]}, '
            f'wall {time.time() - tic:.0f} s')
    except Exception:
        R['error'] = traceback.format_exc()[-2000:]
        R['stats'] = dict(wall_s=round(time.time() - tic, 1))
        say(f'ERROR {job_id}: {R["error"].strip().splitlines()[-1]}')
    C.jdump(R, rf)
    return R


def make_job(kind, run_dir, slice_k, T, H, L, remaining, params, craft=None, route=None, tank_prev=C.TANK_DRY,
             counts_prev=None, grid=None, roots=60, prize=None, seed=0, idx=0, quota=None):
    """Build a JOB dict (schema above) for the driver; job_id = f"s{slice_k:02d}_c{craft}" or f"s{slice_k:02d}_r{idx}";
    out = f"slices/s{slice_k:02d}/{c<craft>|r<idx>}".  params = {k: params.get(k, JOB_DEF[k]) for k in JOB_DEF}."""
    params = params or {}
    tag = f'c{int(craft)}' if kind == 'craft' else f'r{int(idx)}'
    job_id = f's{int(slice_k):02d}_{tag}'
    if kind == 'craft':
        assert craft is not None and route is not None, 'a craft job needs craft and route'
    else:
        assert grid is not None and len(grid) == 3, 'a root job needs grid = [t_lo, t_hi, step] (s)'
    return dict(job_id=job_id, kind=kind, craft=(int(craft) if kind == 'craft' else None),
                route=(str(route) if kind == 'craft' else None), tank_prev=float(tank_prev if kind == 'craft' else C.TANK_DRY),
                counts_prev=dict(counts_prev or dict(pin=0, med=0, fil=0)), remaining=[int(x) for x in remaining],
                T=float(T), H=float(H), L=float(L), grid=([float(x) for x in grid] if grid is not None else None),
                roots=int(roots), prize=({str(k): float(v) for k, v in prize.items()} if prize else None),
                params={k: params.get(k, JOB_DEF[k]) for k in JOB_DEF}, seed=int(seed),
                run_dir=str(pathlib.Path(run_dir).resolve()), out=f'slices/s{int(slice_k):02d}/{tag}',
                quota=({k: int(v) for k, v in quota.items()} if quota else None))


# ------------------------------------------------------------------ tests (WP1 acceptance)
def _column_table(cols, T, H):
    L = ['  k  id                n_new  depth  pot   tank    dtank   miss  t_launch  t_end   classes         targets @ epochs [d]']
    for c in cols:
        eps = ' '.join(f'{a}@{C.s2d(e):.0f}' for a, e in zip(c['targets'], c['epochs']))
        L.append(f'  {c["id"]:>18s}  {c["n_new"]:5d}  {c["depth"]:5d}  {c["potential"]:3d}  {c["tank"]:6.1f}  {c["dtank"]:6.1f}  '
                 f'{c["miss"]:5.1f}  {C.s2d(c["t_launch"]):8.1f}  {C.s2d(c["t_end"]):6.1f}  '
                 f'p{c["counts"]["pin"]}/m{c["counts"]["med"]}/f{c["counts"]["fil"]}  {eps}')
    return '\n'.join(L)


def _check_columns(R, rd, T, H, kind, craft):
    """Common assertions of the two tests: schema, twin reloads with miss <= 150, epochs in (T, T+H], dedupe."""
    cols = R['columns']
    assert R['ok'], f'job failed: {R.get("error")}'
    assert len(cols) >= 1, 'no column'
    seen = set()
    for c in cols:
        assert c['craft'] == craft and c['kind'] == kind, f'column {c["id"]}: craft {c["craft"]} kind {c["kind"]}'
        assert c['n_new'] >= 1, f'column {c["id"]}: n_new 0'
        C.check_column(c, T, H)
        st = C.load_twin(pathlib.Path(rd) / c['npz']); ip = C.twin_of(st); miss = C.miss_of(ip)
        assert miss <= C.TOL_MISS, f'column {c["id"]}: reloaded miss {miss:.1f} km'
        assert abs(float(ip.tank()) - c['tank']) < 0.5, f'column {c["id"]}: tank {c["tank"]} vs reload {ip.tank()}'
        key = frozenset(c['targets']); assert key not in seen, f'column {c["id"]}: duplicate target set'
        seen.add(key)
        assert c['potential'] >= 0
    return cols


def test_craft(out, route, T_d, H_d, L_d, B, tries, levels, w_t=None, w_fuel=None):
    """Acceptance test 1 (< 8 min on 1 core, MEASURE and print wall): truncate `route` at T_d (C.truncate_twin), save it
    as OUT/craft_0.npz, build a craft JOB over remaining = TARGETS - prefix targets, run_job, then assert:
    >= 1 column with n_new >= 1; every column's npz reloads with miss <= 150 km (C.miss_of); every epoch in (T, T+H];
    potential >= 0 and potential > 0 for at least one column when levels >= 3; columns deduped by target set;
    result.json idempotent (second run_job call returns without computing).  Prints the column table."""
    C.set_threads(); tic0 = time.time()
    rd = C.run_dir(out); rd.mkdir(parents=True, exist_ok=True)
    RI.eph(); W = C.load_windows()
    T = C.d2s(T_d); H = C.d2s(H_d); L = C.d2s(L_d)
    rp = pathlib.Path(route); rp = rp if rp.is_absolute() else ROOT / rp
    st = C.load_twin(rp)
    tic = time.time(); R0, nk = C.truncate_twin(st, T)
    t_tr = time.time() - tic
    assert R0['ok'], f'truncation not honest: {R0}'
    C.save_twin(rd / 'craft_0.npz', R0['st'])
    pre = [int(a) for a in R0['st']['asts']]
    print(f'[{C.now()}] prefix: {rp.name} truncated at {T_d:.0f} d -> {nk} flybys, honest tank {R0["tank"]:.1f} kg, miss '
          f'{R0["miss"]:.1f} km, honest {R0["honest"]["ok"]} (MEASURED {t_tr:.1f} s)', flush=True)
    remaining = [t for t in C.TARGETS if t not in set(pre)]
    job = make_job('craft', rd, 0, T, H, L, remaining, dict(B=B, tries=tries, max_levels=levels, w_t=w_t, w_fuel=w_fuel),
                   craft=0, route='craft_0.npz',
                   tank_prev=R0['tank'], counts_prev=C.count_classes(W, pre), seed=0)
    jd = rd / job['out']
    if jd.exists():
        shutil.rmtree(jd)                        # the test measures a fresh run
    C.jdump(job, rd / 'slices/s00/jobs' / f'{job["job_id"]}.json')
    tic = time.time(); R = run_job(job); wall = time.time() - tic
    tic = time.time(); R2 = run_job(job); w2 = time.time() - tic
    cols = _check_columns(R, rd, T, H, 'craft', 0)
    assert R2 == C.jload(jd / 'result.json') and w2 < 1.0, f'second run_job not the cached result ({w2:.2f} s)'
    for c in cols:
        assert set(c['all_targets']) == set(pre) | set(c['targets']), f'column {c["id"]}: all_targets != prefix + new'
        assert not (set(c['targets']) & set(pre))
    npos = sum(1 for c in cols if c['potential'] > 0)
    print(f'[{C.now()}] test-craft: {len(cols)} columns, n_new {min(c["n_new"] for c in cols)}-{max(c["n_new"] for c in cols)}, '
          f'{npos} with potential > 0; stats {json.dumps({k: v for k, v in R["stats"].items() if k != "hist"}, default=float)}')
    print(_column_table(cols, T, H))
    for h in R['stats']['hist']:
        print('  level', {k: v for k, v in h.items() if k != 'fails'})
    if levels >= 3:
        if npos == 0:
            # the potential machinery needs a beam state past T + H; MEASURED 09-24: the production score (w_t 1.0) advances
            # 80-120 d per level, so 3 levels from the 464 d prefix reach ~835 d < 1080 d; --levels 8 passes in 130 s
            far = max(c['t_end'] for c in cols)
            print(f'[{C.now()}] test-craft FAIL on the potential criterion: no beam state reached T+H = {C.s2d(T + H):.0f} d '
                  f'after {levels} levels (deepest column ends {C.s2d(far):.0f} d, {(C.s2d(far) - C.s2d(R0["st"]["tf"].max())) / max(1, levels):.0f} '
                  f'd per level); rerun with more --levels (8 passes) or a smaller --w-t', flush=True)
        assert npos >= 1, 'no column with potential > 0'
    print(f'[{C.now()}] test-craft PASS: run_job MEASURED {wall:.1f} s (cached call {w2:.2f} s); total {time.time() - tic0:.1f} s; '
          f'result {jd / "result.json"}', flush=True)
    return R


def test_root(out, T_d, H_d, L_d, step_d, roots, B, tries, levels, w_t=None, w_fuel=None):
    """Acceptance test 2 (< 6 min on 1 core): root JOB over grid (T, T+H] step step_d, remaining = TARGETS; assert >= 1
    column, craft None, t_launch in (T, T+H], first epoch <= T+H, and >= 2 distinct root_bin launch bins among the
    columns when the beam kept >= 2 states at level 1."""
    C.set_threads(); tic0 = time.time()
    rd = C.run_dir(out); rd.mkdir(parents=True, exist_ok=True)
    RI.eph(); C.load_windows()
    T = C.d2s(T_d); H = C.d2s(H_d); L = C.d2s(L_d)
    job = make_job('root', rd, 0, T, H, L, list(C.TARGETS), dict(B=B, tries=tries, max_levels=levels, w_t=w_t, w_fuel=w_fuel),
                   grid=[T, T + H, C.d2s(step_d)], roots=roots, seed=0, idx=0)
    jd = rd / job['out']
    if jd.exists():
        shutil.rmtree(jd)
    C.jdump(job, rd / 'slices/s00/jobs' / f'{job["job_id"]}.json')
    tic = time.time(); R = run_job(job); wall = time.time() - tic
    tic = time.time(); R2 = run_job(job); w2 = time.time() - tic
    cols = _check_columns(R, rd, T, H, 'root', None)
    assert R2 == C.jload(jd / 'result.json') and w2 < 1.0, f'second run_job not the cached result ({w2:.2f} s)'
    root_bin = C.d2s(job['params']['root_bin'])
    bins = set()
    for c in cols:
        assert T < c['t_launch'] <= T + H, f'column {c["id"]}: t_launch {C.s2d(c["t_launch"]):.1f} d outside (T, T+H]'
        assert c['epochs'][0] <= T + H and c['epochs'][0] > c['t_launch']
        assert set(c['all_targets']) == set(c['targets'])
        bins.add(int(c['t_launch'] // root_bin))
    n_start = R['stats']['hist'][0]['parents'] if R['stats']['hist'] else 0
    npos = sum(1 for c in cols if c['potential'] > 0)
    print(f'[{C.now()}] test-root: {len(cols)} columns, n_new {min(c["n_new"] for c in cols)}-{max(c["n_new"] for c in cols)}, '
          f'{npos} with potential > 0, launch bins ({root_bin / DAY:.0f} d) {sorted(bins)}; '
          f'stats {json.dumps({k: v for k, v in R["stats"].items() if k != "hist"}, default=float)}')
    print(_column_table(cols, T, H))
    for h in R['stats']['hist']:
        print('  level', {k: v for k, v in h.items() if k != 'fails'})
    if n_start >= 2:
        assert len(bins) >= 2, f'only {len(bins)} launch bin(s) among the columns'
    print(f'[{C.now()}] test-root PASS: run_job MEASURED {wall:.1f} s (cached call {w2:.2f} s); total {time.time() - tic0:.1f} s; '
          f'result {jd / "result.json"}', flush=True)
    return R


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('cmd', choices=['run', 'test-craft', 'test-root']); ap.add_argument('arg')
    ap.add_argument('--route', default='results/s16/best/fleet/route_r1.npz')
    ap.add_argument('--T', type=float, default=540.0); ap.add_argument('--H', type=float, default=540.0)
    ap.add_argument('--L', type=float, default=360.0); ap.add_argument('--B', type=int, default=3)
    ap.add_argument('--tries', type=int, default=4); ap.add_argument('--levels', type=int, default=3)
    ap.add_argument('--step', type=float, default=30.0); ap.add_argument('--roots', type=int, default=10)
    ap.add_argument('--w-t', type=float, default=None, help='beam time weight (default None = s14 production 1.0)')
    ap.add_argument('--w-fuel', type=float, default=None, help='beam fuel weight (default None = s14 production 0.52)')
    a = ap.parse_args()
    if a.cmd == 'run':
        r = run_job(C.jload(a.arg))
        print(json.dumps(dict(job_id=r['job_id'], ok=r['ok'], n_columns=len(r['columns']), stats=r.get('stats'), error=r.get('error')),
                         default=float))
    elif a.cmd == 'test-craft':
        test_craft(a.arg, a.route, a.T, a.H, a.L, a.B, a.tries, a.levels, a.w_t, a.w_fuel)
    else:
        test_root(a.arg, a.T, a.H, a.L, a.step, a.roots, a.B, a.tries, a.levels, a.w_t, a.w_fuel)


if __name__ == '__main__':
    main()
