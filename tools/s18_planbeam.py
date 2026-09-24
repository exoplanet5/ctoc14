"""Stage 18 / PLANBEAM: a drop-in PROPOSER for the rolling-horizon fleet auction that PLANS LONG and COMMITS SHORT.

Why.  The segment beam (tools/s18_segbeam.py: twin beam, width 8, every child settled, horizon T + H + L) offers per
craft and 540 d slice 1 target for 1-20 kg, 2 for 20-45 kg, 3 for 50-85 kg, 4 for 70-220 kg (MEASURED on
results/s18/g1_base/slices/s02,s03/columns.json), while s16a's deep routes fly 4.2 flybys per 540 d at ~8 kg each: a
short-horizon, narrow, settle-every-child beam cannot set up cheap sequences.  Here the PRODUCTION planner
(ctoc14.search: Lambert 15-400 d with tof refinement + LinLeg 20-600 d, m0 1600, prizes) runs a WIDE beam over a long
horizon (T .. min(T_MISSION, T + H + L + plan_extra)), and the columns are TRUNCATIONS of the planned tours at T + H:
  committed  = the new flybys with epoch <= T + H (>= 1), settled into an honest-enough twin (miss <= 150 km) by chaining
               s14_twinbeam.lin_price / settle_child from the craft's prefix twin (root jobs: s14_twinbeam.root_state of
               the planner's launch leg, then the same chain);
  potential  = the number of planned flybys in (T + H, T + H + L] (the segbeam definition; flybys planned beyond
               T + H + L shape the plan but are not counted).  A TRUNCATED candidate (k < the in-window planned flybys)
               keeps that definition in the emitted column (pot_mode 'lookahead': the auction's gamma calibration is
               segbeam's), while the proposer RANKS candidates by the continuation count pot_cont = planned flybys in
               (t_last_committed, T + H + L] (so shallow truncations of good plans rank above dead ends); stats.plan
               carries plan_cont per column; pot_mode 'continuation' emits pot_cont as the column potential.
               pot_mode 'full' (09-24, additive) emits pot_full = the planned flybys in the WHOLE continuation
               (t_last_committed, plan end] (plan end = min(T_MISSION, T + H + L + plan_extra_d)), i.e. the depth the
               plan reaches beyond the committed part, and in that mode the dedupe / frontier ranking use pot_full
               (the "best continuation" is the full one); pot_cap (default None) caps the emitted potential.
               Settable per run by S18_PLANBEAM='{"pot_mode": "full", "pot_cap": 8}'.  stats.plan carries plan_full.
Every column also carries `plan_next` (09-24, additive): the plan's flybys beyond the committed part, in plan order
(a settle-truncated column's plan_next starts with the legs the settle cut off); the driver's --claims option turns
them into the craft's claims for the next slices (s18_rhfa).  stats.plan[col].plan_next mirrors it.
Columns are deduplicated by their set of new targets (best continuation, then cheapest committed plan fuel), and
max_cols are kept: per depth 1..max_depth_col the best-continuation and the (per_depth - 1) cheapest columns, then the
rest by (continuation, fuel) round-robin over the first new target (so the auction gets mutually compatible offers);
when more than max_cols were taken, the cut keeps every depth's rank-0 column first, then rank 1, ..., deepest first
within a rank (the deep columns are the ones the fleet needs; a 1-target column is never preferred over a 5-target one).
WALL: one job_wall budget (300 s) is shared by the plan (min(plan_wall, job_wall - elapsed - settle_reserve)), the
settles and the regrids (each SIGALRM capped at the remaining budget; columns beyond it are dropped, logged); a
prefix whose regrids fail regrid_fail_max times in a row (a stranded twin, the c0 pattern) stops settling early.

Contract: EXACTLY the s18_segbeam JOB / RESULT contract (run_job(job) -> RESULT dict, out/result.json cache and
idempotence, a cached failure re-run once, ok False on error, single-threaded, column twins col_<k>.npz under the job's
out dir, s18_common.check_column passes for every column).  The driver runs it with `s18_rhfa.py run ... --proposer
planbeam` (s18_rhfa.resolve_proposer).  Extra JOB keys are optional: job['planbeam'] = {PB_DEF overrides}; the
environment variable S18_PLANBEAM (a JSON dict) overrides PB_DEF for a whole run (the driver's workers inherit it).
Checkpoint: out/plan.pkl (the planned states) so that a job killed while settling does not re-plan.

RESULT.stats: levels (plan depth = max new flybys planned), n_states (planned states), n_columns, settles, settle_ok,
regrids, regrid_ok, price_s (lin_price wall), settle_s, plan_s, wall_s, end {kind: plan_done | plan_wall | no_children |
pool_exhausted, depth}, budget {job_wall, plan_cap, settle_left_s (budget left when the settles started), hit (a
budget or streak stop: None | 'settle_wall' | 'job_wall' | 'regrid_streak'), n_unsettled}, hist [per plan depth {depth,
kept (beam states at that depth), in_H (with last flyby <= T + H)}], plan {column id -> {plan_n, plan_fuel_kg,
plan_t_end_d, commit_fuel_kg, plan_cont, cut}} (planner estimates, labelled as such), frontier [per n_new: {n_new,
dtank_min, potential_max}].

usage: s18_planbeam.py run JOB.json
       s18_planbeam.py test-slice OUT SRC_RUN K --crafts 0,1,4,7 [--nproc 2] [--width 150] [--plan-wall 240]
           rebuilds the slice-K craft jobs of SRC_RUN from SRC_RUN/slices/sKK/jobs/*.json (routes copied into OUT
           first), runs run_job for the listed crafts (nproc <= 2), and prints the frontier table (cheapest dtank per
           n_new, its potential, wall) next to segbeam's from SRC_RUN/slices/sKK/columns.json.
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
import s18_segbeam as SB
from ctoc14 import search as SR
from ctoc14.search import mask_of, UNREACHABLE, State
from ctoc14.constants import DAY, VE, T_MISSION, FUEL_MAX, AU, VINF_MAX
warnings.filterwarnings('ignore')

PB_DEF = dict(width=150,            # planner beam width (production isogen: 100; s16a routes: 100-150)
              npt=6, mc=150,        # children per target per parent / per parent (production isogen member H)
              dv=1.2, drmax=0.15,   # production leg caps (km/s, AU)
              w_t=3.0, w_fuel=None,   # None = production (w_t 1.0, w_fuel CM.w_fuel(480, 1600) = 0.52).  w_t 3: cadence
                                      # matters inside a 540 d commit window (MEASURED 09-24 c4: the 4-flyby column
                                      # 67 -> 30 kg, c7 unchanged, plan wall +20 s; results/s18/planbeam_test/var)
              job_wall=300.0,       # s: TOTAL budget of a job (plan + settles + regrids); every SIGALRM below is capped
                                    # at what is left of it (the driver's LOST deadline is 3000 s; the slice wall is the
                                    # max over the craft jobs, so one stranded prefix must not take 700 s)
              plan_wall=240.0,      # s: SIGALRM cap on the planning beam (the states collected so far are used); the
                                    # effective cap is min(plan_wall, job_wall - elapsed - settle_reserve)
              settle_reserve=90.0,  # s: budget kept back from the plan for the settles + regrids
              plan_extra_d=720.0,   # the plan horizon is min(T_MISSION, T + H + L + plan_extra_d)
              max_depth=40,         # planner expansion rounds (the horizon stops it first)
              max_cols=16, max_depth_col=8,   # columns kept per job; depths covered by the frontier
              per_depth=3,          # per depth: the best-continuation candidate + the (per_depth - 1) cheapest
              pot_mode='lookahead', # column potential: 'lookahead' = planned flybys in (T + H, T + H + L] (segbeam's
                                    # definition, the auction's gamma calibration) | 'continuation' = planned flybys in
                                    # (t_last_committed, T + H + L] (counts the in-window flybys a truncation dropped)
                                    # | 'full' (09-24) = planned flybys in (t_last_committed, plan end]: the whole
                                    # continuation of the plan.  The proposer's own ranking uses the continuation
                                    # count (pot_cont), in mode 'full' the full count.
              pot_cap=None,         # cap on the emitted column potential (None = no cap); e.g. 8 with pot_mode 'full'
              regrid_cols=1,        # 1: C.regrid_settle every settled column and emit the HONEST twin (a column whose
                                    # regrid fails is dropped: the commit would reject it anyway; MEASURED 09-24 c0: 6 of
                                    # 8 planner-seeded twins on a stranded prefix failed the regrid law)
              regrid_fail_max=4,    # consecutive regrid failures (no success in between) after which only candidates
                                    # that EXTEND a regrid-ok prefix are still settled (a stranded twin: c0 09-24 had
                                    # 6 of 7 fresh first legs fail the law at 5-40 s each; its 2-target column extends
                                    # the one that passed)
              iters=100, timeout=240.0,       # settle_child iterations / SIGALRM per settle (capped at the budget left)
              settle_wall=480.0,    # s: settle budget of a job (the job_wall budget is the binding one by default)
              lin_seed=1,           # 1: s14 lin_price step as the first settle seed (lintwin; no irregular nodes)
              # FIRST LEVEL of a craft plan (MEASURED 09-24 on g1_base slice 3: the production caps give 0-29 children
              # from a twin end state, c0 / c5 stranded at 1363 / 1130 d give 0 / 1): the first legs come from the
              # relaxed caps dv1 / drmax1 (the settle re-optimises the whole prefix, so a myopic 1.5 km/s first leg is
              # not what the twin pays) AND from s18_segbeam.price_parent_full (the whole-prefix linearised twin price,
              # the segbeam's own level-1 candidates; the child state is one integration of the linear step).
              first_leg='both',     # 'relaxed' | 'lintwin' | 'both' | 'prod'
              dv1=2.5, drmax1=0.30, k_epochs=3, lin_cap=3.0, res_max=1.5e8, price_cap=250,
              root_dedup_d=20.0,    # root legs: (first target, launch bin) dedupe, days
              roots_max=450)        # root legs kept for the planner front (search.roots gives ~13000 per 54 dates);
                                    # search._beam_loop keeps at most 3 * width of the initial front, so the effective
                                    # value is min(roots_max, 3 * width): raise width for a wider root front
LOOK_KEYS = ('dv', 'drmax', 'w_t', 'w_fuel', 'iters', 'timeout')   # JOB params (segbeam keys) that also apply here


# ------------------------------------------------------------------ parameters
def pb_params(job):
    """PB_DEF <- job['params'] (shared keys) <- env S18_PLANBEAM (JSON) <- job['planbeam']."""
    prm = dict(PB_DEF)
    jp = job.get('params') or {}
    for k in LOOK_KEYS:
        if k in jp and jp[k] is not None:
            prm[k] = jp[k]
    env = os.environ.get('S18_PLANBEAM')
    if env:
        try:
            prm.update({k: v for k, v in json.loads(env).items() if k in PB_DEF})
        except Exception:
            pass
    prm.update({k: v for k, v in (job.get('planbeam') or {}).items() if k in PB_DEF})
    return prm


def planner_params(prm, prize, t_plan, t_lo=-np.inf):
    """ctoc14.search.Params of the production planner (s14_twinbeam.params / s15_isogen.params) with the beam width,
    the prize vector and a PLAN HORIZON: every target gets win_hi = t_plan with mandatory = [] (no deadline prune, no
    lookahead heuristic), so expand() only makes legs that end by t_plan."""
    P = TB.params(float(prm['dv']), float(prm['drmax']), int(prm['mc']))
    P.beam = int(prm['width']); P.n_per_target = int(prm['npt']); P.max_depth = int(prm['max_depth'])
    if prm.get('w_t') is not None:
        P.w_t = float(prm['w_t'])
    if prm.get('w_fuel') is not None:
        P.w_fuel = float(prm['w_fuel'])
    P.prize = np.asarray(prize, float)
    P.win_lo = np.full(300, float(t_lo)); P.win_hi = np.full(300, float(t_plan)); P.mandatory = []
    P.lookahead = False
    return P


def first_level(start, job, prm, P, T, H, L, t_plan, excl_mask, say):
    """Depth-1 children of a craft plan (see PB_DEF first_leg): search.State list, deduped on (visited, 5 d)."""
    import copy
    E = RI.eph(); mode = str(prm['first_leg']); s0 = TB.search_state(start, excl_mask, P)
    s0 = State(t=s0.t, r=s0.r, v=s0.v, m=s0.m, visited=s0.visited, seq=s0.seq, fuel=s0.fuel, t_launch=s0.t_launch,
               vinf=s0.vinf, m0=s0.m0, rare=0.0)
    kids = []; n_rel = n_lin = n_prod = 0; tic = time.time()
    if mode == 'prod':
        kids += SR.expand(E, s0, P); n_prod = len(kids)
    if mode in ('relaxed', 'both'):
        Pr = copy.copy(P); Pr.dv_max = float(prm['dv1']); Pr.lin_drmax = float(prm['drmax1']) * AU; Pr.collect = None
        k = SR.expand(E, s0, Pr); n_rel = len(k); kids += k
    t_rel = time.time() - tic
    if mode in ('lintwin', 'both'):
        rem = [int(t) for t in job['remaining']]
        lst, info = SB.price_parent_full((start, rem, T + SB.EDGE_LO, min(float(t_plan), T + H + L), int(prm['k_epochs']),
                                          float(prm['lin_cap']), P.w_t, P.w_fuel, float(prm['res_max']), P.prize, int(prm['price_cap'])))
        for sc, X, t, dvm, lin, cm, c1 in lst:
            try:
                ip, _ = TB.child_problem(start, int(X), float(t), None)
                tf2 = np.minimum(ip.tf + np.asarray(lin['dt']) * DAY, T_MISSION - 2 * DAY)
                vi = ip.vinf + np.asarray(lin['dvinf']); nv = np.linalg.norm(vi)
                if nv > VINF_MAX - 1e-6:
                    vi = vi * (VINF_MAX - 1e-6) / nv
                Yf2, _ = ip.integrate(Ts=np.asarray(lin['T'], float), tf=tf2, vinf=vi)
                j = [i for i, a in enumerate(ip.asts) if int(a) == int(X)][-1]
                t_new = float(tf2[j])
                if not (T < t_new <= float(t_plan)):
                    continue
                dv = max(float(dvm), 0.0); m = s0.m * np.exp(-dv / VE)
                kids.append(State(t=t_new, r=np.asarray(Yf2[j, :3], float), v=np.asarray(Yf2[j, 3:6], float), m=float(m),
                                  visited=s0.visited | (1 << (int(X) - 1)), seq=s0.seq + ((int(X), t_new, dv, t_new - s0.t),),
                                  fuel=float(s0.fuel + s0.m - m), t_launch=s0.t_launch, vinf=s0.vinf, m0=s0.m0,
                                  rare=s0.rare + float(P.prize[int(X) - 1])))
                n_lin += 1
            except Exception:
                continue
    kids.sort(key=lambda s: s.score(P))
    seen = set(); out = []
    for c in kids:
        key = (c.visited, int(c.t // (5 * DAY)))
        if key in seen:
            continue
        seen.add(key); out.append(c)
    say(f'first level ({mode}): {n_prod} production + {n_rel} relaxed (dv {prm["dv1"]}, {prm["drmax1"]} AU; {t_rel:.1f} s) + '
        f'{n_lin} lintwin children -> {len(out)} unique ({time.time() - tic:.1f} s)')
    return s0, out


# ------------------------------------------------------------------ planning
class _Collect(list):
    """P.collect: every beam (all depths) is appended; the depth-1 start states are added by plan()."""
    def extend(self, states):
        list.extend(self, states)


def plan(kind, start, job, prm, P, T, H, L, say, plan_cap=None):
    """The long planning beam.  Returns dict(states=[search.State ...] (every beam state of every depth, start states
    included), end={kind, depth}, plan_s, n_pre, t_plan, plan_cap).
    craft: one start State from the prefix twin (TB.search_state: planner mass M0P exp(-dv_twin / ve));
    root:  search.roots over the job's launch grid (T, T + H] (launches >= T), first flyby <= T + H, sorted by score,
           deduped on (first target, root_dedup_d launch bin), at most min(roots_max, 3 * width) kept (the cap of
           search._beam_loop on its initial front).
    plan_cap: SIGALRM cap on the beam loop [s] (None: prm['plan_wall']); the caller derives it from the job budget."""
    E = RI.eph(); tic = time.time()
    coll = _Collect(); P.collect = coll
    plan_cap = float(prm['plan_wall'] if plan_cap is None else plan_cap)
    if kind == 'craft':
        ts_ = start
        excl_mask = mask_of((set(C.TARGETS) - set(int(t) for t in job['remaining'])) | set(UNREACHABLE))
        s0, beam0 = first_level(start, job, prm, P, T, H, L, float(P.win_hi[0]), excl_mask, say)
        n_pre = len(s0.seq)
        say(f'plan start: prefix {n_pre} flybys, t_end {C.s2d(s0.t):.1f} d, planner mass {s0.m:.1f} kg (dv {ts_["dv"]:.3f} km/s); '
            f'{len(beam0)} first-level states')
    else:
        t_lo, t_hi, step = [float(x) for x in job['grid']]
        launches = np.arange(t_lo + step, t_hi + 1e-6, step)
        launches = launches[(launches > t_lo) & (launches <= t_hi) & (launches >= T - 1e-6) & (launches < T_MISSION - 20 * DAY)]
        excl = sorted(set(C.TARGETS) - set(int(t) for t in job['remaining']))
        rts = SR.roots(E, launches, TB.M0P, excl, P) if len(launches) else []
        rts = [s for s in rts if s.t <= T + H]
        rts.sort(key=lambda s: (s.score(P), s.t, s.t_launch))
        seen = set(); beam0 = []; bin_s = float(prm['root_dedup_d']) * DAY
        n_keep = min(int(prm['roots_max']), 3 * int(P.beam))       # _beam_loop keeps beam[:3 * P.beam] of the front
        for s in rts:
            key = (int(s.seq[0][0]), int(s.t_launch // bin_s))
            if key in seen:
                continue
            seen.add(key); beam0.append(s)
            if len(beam0) >= n_keep:
                break
        n_pre = 0
        say(f'plan start: {len(launches)} launch dates ({C.s2d(launches.min()) if len(launches) else 0:.0f}-'
            f'{C.s2d(launches.max()) if len(launches) else 0:.0f} d) -> {len(rts)} launch legs with first flyby <= T+H, '
            f'{len(beam0)} kept (cap min(roots_max {prm["roots_max"]}, 3 x width) = {n_keep}; {time.time() - tic:.1f} s)')
    coll.extend(beam0)
    end = dict(kind='plan_done', depth=0)
    if beam0:
        try:
            C.with_timeout(SR._beam_loop, max(1.0, plan_cap), E, beam0, P, 1, False)
        except C.Timeout:
            end = dict(kind='plan_wall', depth=0)
    else:
        end = dict(kind='no_children', depth=0)
    states = list(coll)
    depth = max((len(s.seq) - n_pre for s in states), default=0); end['depth'] = int(depth)
    return dict(states=states, end=end, plan_s=time.time() - tic, n_pre=n_pre, t_plan=float(P.win_hi[0]), plan_cap=plan_cap)


def _leg_key(a, t):
    return (int(a), int(round(float(t) / DAY * 10)))


def _launch_key(s):
    return int(round(float(s.t_launch) / DAY * 10))


def _rk(c):
    """The candidate's ranking count: pot_cont (the default modes) or pot_full (pot_mode 'full'); candidates() sets
    'rank', hand-made candidates (unit tests) fall back to pot_cont."""
    return c.get('rank', c['pot_cont'])


def candidates(res, T, H, L, quota, cls, prm, say):
    """Column candidates from the planned states: for every state and every truncation k of its NEW legs (1 <= k <= the
    number of new legs with epoch <= T + H, cut where a class quota head-room would be exceeded), one candidate keyed by
    the committed target set: dict(legs=[(ast, t)], pot (planned flybys in (T + H, T + H + L]), pot_cont (planned flybys
    in (t_last_committed, T + H + L]: pot plus the in-window flybys the truncation dropped), pot_full (planned flybys in
    (t_last_committed, plan end]: the whole continuation, = plan_n - k), rank (pot_cont, or pot_full when
    prm['pot_mode'] == 'full'), fuel_c (planner kg of the
    committed part: the fuel of the beam state with exactly those legs AND this launch, else an estimate from the leg
    dv), plan_n, plan_fuel, plan_t_end, score, first, t_launch, vinf (the launch of the plan that won the dedupe: a root
    column is settled from ITS launch, not from another plan's with the same first leg)).
    Dedupe: best (rank, -fuel_c)."""
    n_pre = res['n_pre']; states = res['states']
    full_mode = str((prm or {}).get('pot_mode', 'lookahead')) == 'full'
    fuel_of = {}
    for s in states:
        fuel_of[(_launch_key(s),) + tuple(_leg_key(a, t) for a, t, _, _ in s.seq[n_pre:])] = float(s.fuel)
    quota = {k: int(v) for k, v in (quota or {}).items() if v is not None}
    best = {}; n_cut = 0
    for s in states:
        new = list(s.seq[n_pre:])
        if not new:
            continue
        com = [l for l in new if l[1] <= T + H + 1e-6]      # seq is in flyby order: a time prefix
        if quota:
            cnt = dict(pin=0, med=0, fil=0); k_ok = 0
            for a, t, _, _ in com:
                c = cls.get(int(a), 'pin'); cnt[c] += 1
                if c in quota and cnt[c] > quota[c]:
                    n_cut += 1; break
                k_ok += 1
            com = com[:k_ok]
        if not com:
            continue
        pot = sum(1 for a, t, _, _ in new if T + H < t <= T + H + L)
        lk = _launch_key(s)
        for k in range(1, len(com) + 1):
            legs = com[:k]; key = frozenset(int(a) for a, _, _, _ in legs)
            fk = fuel_of.get((lk,) + tuple(_leg_key(a, t) for a, t, _, _ in legs))
            if fk is None:
                # the ancestor with exactly these legs was not kept: the state's fuel minus the uncommitted legs' share
                # (rocket equation on the leg dv estimates from the end mass; the leverage factor kappa is ignored)
                dv_unc = sum(float(l[2]) for l in new[k:]); m_end = s.m0 - s.fuel
                fk = float(s.fuel - m_end * (np.exp(dv_unc / VE) - 1.0))
            t_last = float(legs[-1][1])
            pot_cont = sum(1 for a, t, _, _ in new if t_last < t <= T + H + L)
            pot_full = sum(1 for a, t, _, _ in new if t > t_last)          # the whole planned continuation
            cand = dict(legs=[(int(a), float(t)) for a, t, _, _ in legs], pot=int(pot), pot_cont=int(pot_cont),
                        pot_full=int(pot_full), rank=int(pot_full if full_mode else pot_cont), fuel_c=float(fk),
                        plan_n=len(s.seq) - n_pre, plan_fuel=float(s.fuel), plan_t_end=float(s.t),
                        score=float(s.score_val) if hasattr(s, 'score_val') else 0.0, first=int(legs[0][0]), n_new=k,
                        t_launch=float(s.t_launch), vinf=np.asarray(s.vinf, float),
                        plan_next=[int(a) for a, _, _, _ in new[k:]])     # the plan's continuation, in plan order (09-24)
            o = best.get(key)
            if o is None or (cand['rank'], -cand['fuel_c']) > (o['rank'], -o['fuel_c']):
                best[key] = cand
    if n_cut:
        say(f'  quota head-room {quota}: {n_cut} plans cut at the quota')
    return list(best.values())


def frontier_select(cands, prm):
    """max_cols candidates spanning depths 1..max_depth_col: per depth the best continuation (pot_cont) and the
    (per_depth - 1) cheapest committed fuel, then the rest by (continuation desc, fuel asc) round-robin over the first
    new target.  When the per-depth pass alone exceeds max_cols (more than max_cols / per_depth depths in the window),
    the cut is by (rank within its depth, deepest first, fuel): every depth keeps its rank-0 column before any depth
    keeps a rank-1 one, and within a rank the DEEPER columns survive (the old cut sorted by (n_new, fuel) and dropped
    the deepest columns: reviewer finding, 09-24)."""
    max_cols = int(prm['max_cols']); dmax = int(prm['max_depth_col'])
    by_d = {}
    for c in cands:
        if c['n_new'] <= dmax:
            by_d.setdefault(c['n_new'], []).append(c)
    chosen = []; ids = set(); rank_of = {}

    def take(c, rank):
        k = frozenset(a for a, _ in c['legs'])
        if k not in ids:
            ids.add(k); chosen.append(c); rank_of[k] = rank
            return True
        return False
    n_per = max(1, int(prm.get('per_depth', 2)))
    for d in sorted(by_d):
        lst = by_d[d]
        take(max(lst, key=lambda c: (_rk(c), -c['fuel_c'])), 0)
        r = 1
        for c in sorted(lst, key=lambda c: (c['fuel_c'], -_rk(c))):
            if r >= n_per:
                break
            if take(c, r):
                r += 1
    rest = [c for c in cands if frozenset(a for a, _ in c['legs']) not in ids and c['n_new'] <= dmax]
    rest.sort(key=lambda c: (-_rk(c), c['fuel_c']))
    groups = {}
    for c in rest:
        groups.setdefault(c['first'], []).append(c)
    order = sorted(groups, key=lambda g: (-_rk(groups[g][0]), groups[g][0]['fuel_c']))
    rank = 0
    while len(chosen) < max_cols:
        got = False
        for g in order:
            if rank < len(groups[g]):
                take(groups[g][rank], n_per + rank); got = True
                if len(chosen) >= max_cols:
                    break
        if not got:
            break
        rank += 1
    if len(chosen) > max_cols:
        chosen.sort(key=lambda c: (rank_of[frozenset(a for a, _ in c['legs'])], -c['n_new'], c['fuel_c']))
        chosen = chosen[:max_cols]
    chosen.sort(key=lambda c: (c['n_new'], c['fuel_c']))      # readability / settle order (prefix cache: shallow first)
    return chosen


# ------------------------------------------------------------------ settling the committed prefixes
def root_twin(tL, vinf, ast, t1, iters, timeout):
    """s14_twinbeam.root_state (exact zero-thrust Lambert arc when |v_inf| <= 4); a leg that needs thrust for the
    v_inf excess (miss > 150 km) is settled once."""
    rs = TB.root_state(float(tL), np.asarray(vinf, float), int(ast), float(t1))
    if rs['miss'] <= C.TOL_MISS:
        return rs, dict(ok=True, how='root')
    ip = RI.ipr(rs['st'])
    try:
        miss = float(C.with_timeout(RI.settle, float(timeout), ip, int(iters)))
    except Exception as e:
        return None, dict(ok=False, how='root_settle', err=repr(e)[:80])
    if miss <= C.TOL_MISS:
        return TB.pack(ip, rs['nreg'], miss, 'root_settled'), dict(ok=True, how='root_settled')
    return None, dict(ok=False, how='root_settled', miss=miss)


MIN_LEFT_S = 3.0        # s: a settle / regrid is not started with less budget than this


def settle_columns(chosen, kind, start, T, H, prm, say, stats, deadline=None):
    """Chain-settle every chosen candidate's committed legs from the start twin (craft: the prefix; root: the launch leg
    twin of the candidate's own planner launch).  A cache over leg prefixes shares the work between columns (settle
    AND regrid results).  A leg that fails to settle truncates the column to the settled prefix (dropped if that
    duplicates another column or is empty).  The settle moves the epochs: a new flyby that left (T, T + H] rejects the
    column (drift).  Every settle / regrid SIGALRM is capped at the budget left (min(settle_wall, deadline - now));
    columns beyond the budget are dropped (logged, stats.budget.hit).  After regrid_fail_max consecutive regrid
    failures (a stranded prefix; stats.budget.hit 'regrid_streak') only candidates extending a regrid-ok prefix are
    still settled (fresh first legs are skipped).  The emitted twin of a regridded
    column is re-packed from the regridded problem (t_end / r_end / v_end / tfs consistent with its st).
    Returns [(cand, twin state)]."""
    iters = int(prm['iters']); timeout = float(prm['timeout']); lin_seed = int(prm['lin_seed'])
    wall = float(prm['settle_wall']); tic0 = time.time()
    deadline = float('inf') if deadline is None else float(deadline)
    do_regrid = int(prm.get('regrid_cols', 0)); streak_max = int(prm.get('regrid_fail_max', 0) or 0)
    cache = {}                                    # tuple of leg keys -> twin state | None (failed)
    rg_cache = {}                                 # tuple of leg keys -> regridded twin state | None (failed)
    out = []; have = set(); n_drop = 0; streak = 0; n_left = 0
    stats.setdefault('regrids', 0); stats.setdefault('regrid_ok', 0)

    def left():
        return min(wall - (time.time() - tic0), deadline - time.time())
    for ic, c in enumerate(chosen):
        lft = left()
        if lft < MIN_LEFT_S:
            n_left = len(chosen) - ic
            stats['budget_hit'] = 'settle_wall' if wall - (time.time() - tic0) < deadline - time.time() else 'job_wall'
            say(f'  budget reached ({stats["budget_hit"]}: settle {time.time() - tic0:.0f} s of {wall:.0f}, '
                f'{max(0.0, deadline - time.time()):.0f} s of the job left): {n_left} candidates not settled')
            break
        legs = c['legs']; keys = [_leg_key(a, t) for a, t in legs]
        # the longest cached prefix
        j = 0
        for i in range(len(keys), 0, -1):
            if tuple(keys[:i]) in cache:
                j = i; break
        par = cache[tuple(keys[:j])] if j > 0 else None
        if j > 0 and par is None:
            n_drop += 1; continue                 # a failed prefix
        if streak_max and streak >= streak_max and not any(rg_cache.get(tuple(keys[:i])) is not None for i in range(1, len(keys))):
            # a stranded prefix (regrid_fail_max failures in a row): only candidates that EXTEND a regrid-ok prefix
            # are still tried (chained from a good twin they are cheap and pass; c0 09-24: its 2-target column)
            n_left += 1
            if stats.get('budget_hit') != 'regrid_streak':
                stats['budget_hit'] = 'regrid_streak'
                say(f'  {streak} consecutive regrid failures (regrid_fail_max {streak_max}): only extensions of regrid-ok prefixes are settled from here')
            continue
        ok_upto = j
        for i in range(j, len(legs)):
            a, t = legs[i]
            lft = left()
            if lft < MIN_LEFT_S:
                break                             # the column is cut at the settled prefix; the outer loop then stops
            to = min(timeout, lft)
            if i == 0 and kind == 'root':
                tw, info = root_twin(c['t_launch'], c['vinf'], a, t, iters, to)
                stats['settles'] += 1; stats['settle_ok'] += int(tw is not None)
                if tw is None:
                    say(f'    root leg {int(a):3d} @ {C.s2d(t):6.0f} d (launch {C.s2d(c["t_launch"]):.0f} d) FAIL {info}')
            elif i == 0 and kind == 'craft':
                tw, info = _settle_one(start, a, t, iters, to, lin_seed, stats)
            else:
                tw, info = _settle_one(par, a, t, iters, to, lin_seed, stats)
            if tw is not None:
                leg_set = set(x for x, _ in legs)
                new_eps = [float(e) for a_, e in zip(tw['asts'], tw['tfs']) if int(a_) in leg_set]
                if any(not (T < e <= T + H + 1e-6) for e in new_eps):
                    say(f'  drift: leg {a}@{C.s2d(t):.0f} d settled outside (T, T+H] -> column cut at {i}')
                    tw = None
            if tw is None and left() < MIN_LEFT_S:
                break                             # a settle cut by the budget is not a verdict on the leg: not cached
            cache[tuple(keys[:i + 1])] = tw
            if tw is None:
                break
            par = tw; ok_upto = i + 1
        if ok_upto == 0:
            n_drop += 1; continue
        key = frozenset(a for a, _ in legs[:ok_upto])
        if key in have:
            n_drop += 1; continue
        cc = dict(c); cc['legs'] = legs[:ok_upto]; cc['n_new'] = ok_upto; cc['cut'] = ok_upto < len(legs)
        # a truncated column's continuation starts with the legs the settle cut off (plan order kept)
        cc['plan_next'] = [int(a) for a, _ in legs[ok_upto:]] + [int(a) for a in c.get('plan_next', [])]
        tw_out = par
        if do_regrid:
            rk = tuple(keys[:ok_upto])
            if rk in rg_cache:
                tw_out = rg_cache[rk]
                if tw_out is None:
                    n_drop += 1; continue         # this prefix already failed the regrid law
            else:
                lft = left()
                if lft < MIN_LEFT_S:
                    n_left = len(chosen) - ic; stats['budget_hit'] = 'job_wall'
                    say(f'  budget reached before the regrid of column {ic}: {n_left} candidates not emitted'); break
                r = C.regrid_settle(C.twin_of(par['st']), timeout=min(timeout, lft))
                stats['regrids'] += 1
                if r['ok']:
                    tw_out = C.pack_state(r['ip'], float(r['miss']), par['how'] + '+regrid')
                    stats['regrid_ok'] += 1; streak = 0
                else:
                    say(f'  regrid FAIL ({ok_upto} fb: miss {r["miss"]:.0f} km, {r["error"]}, {r["secs"]} s) -> column dropped')
                    tw_out = None
                    if r['error'] != 'timeout' or lft >= timeout:
                        streak += 1               # a regrid cut by the budget is not counted as a law failure
                rg_cache[rk] = tw_out
                if tw_out is None:
                    n_drop += 1; continue
        have.add(key)
        out.append((cc, tw_out))
    stats['settle_s'] += time.time() - tic0; stats['n_unsettled'] = n_left
    if n_drop:
        say(f'  {n_drop} candidates dropped (settle failure / regrid failure / duplicate after truncation)')
    return out


def _settle_one(par, a, t, iters, timeout, lin_seed, stats):
    lin = None
    if lin_seed:
        tic = time.time()
        try:
            dvm, lin, cm = TB.lin_price(par, int(a), float(t))
            if not np.isfinite(dvm) or lin['res_km'] > 1.5e8:
                lin = None
        except Exception:
            lin = None
        stats['price_s'] += time.time() - tic
    tic = time.time()
    tw, info = TB.settle_child((par, int(a), float(t), iters, timeout, lin))
    stats['settles'] += 1
    if tw is None and lin is not None:
        tw, info = TB.settle_child((par, int(a), float(t), iters, timeout, None))
        stats['settles'] += 1
    if tw is not None:
        stats['settle_ok'] += 1
    say = stats.get('say')
    if say is not None:
        say(f'    settle {len(par["asts"]) + 1:2d}: {int(a):3d} @ {C.s2d(t):6.0f} d -> ' + (
            f'ok {tw["how"]} tank {tw["tank"]:.1f} (+{tw["tank"] - par["tank"]:.1f}) miss {tw["miss"]:.0f} km' if tw is not None
            else f'FAIL {info.get("tries")}') + f' ({time.time() - tic:.1f} s)')
    return tw, info


# ------------------------------------------------------------------ columns
def column_potential(c, pot_mode='lookahead', pot_cap=None):
    """The potential a candidate emits: pot ('lookahead', segbeam's definition; also any unknown mode), pot_cont
    ('continuation') or pot_full ('full': the whole planned continuation), then min(., pot_cap) when pot_cap is set."""
    mode = str(pot_mode)
    if mode == 'continuation':
        p = c.get('pot_cont', c['pot'])
    elif mode == 'full':
        p = c.get('pot_full', c.get('pot_cont', c['pot']))
    else:
        p = c['pot']
    p = int(p)
    return min(p, int(pot_cap)) if pot_cap is not None else p


def columns_of(settled, T, H, L, craft, tank_prev, W, run_dir, out_rel, job_id, prefix_asts, n_pre, pot_mode='lookahead',
               pot_cap=None):
    """COLUMN dicts (s18_common.new_column) of the settled candidates; twins to run_dir/out_rel/col_<k>.npz.
    pot_mode: the column potential is the candidate's pot ('lookahead', segbeam's definition), pot_cont ('continuation')
    or pot_full ('full'); pot_cap caps it (column_potential)."""
    pre = set(int(a) for a in prefix_asts); od = pathlib.Path(run_dir) / out_rel
    cols = []; plan = {}
    ordered = sorted(settled, key=lambda z: (-z[0]['n_new'], z[1]['tank']))
    for k, (c, tw) in enumerate(ordered):
        pot_look = int(c['pot'])                          # the lookahead count itself (stats.plan.plan_lookahead)
        c = dict(c, pot=column_potential(c, pot_mode, pot_cap))
        st = tw['st']; o = C.order_of(st)
        all_t = [int(st['asts'][i]) for i in o]; all_e = [float(st['tf'][i]) for i in o]
        new = [(a, e) for a, e in zip(all_t, all_e) if a not in pre]
        rel = f'{out_rel}/col_{k}.npz'
        C.save_twin(od / f'col_{k}.npz', st)
        cid = f'{job_id}_k{k}'
        col = C.new_column(id=cid, craft=craft, targets=[a for a, _ in new], epochs=[e for _, e in new], all_targets=all_t,
                           npz=rel, tank=float(tw['tank']), tank_prev=float(tank_prev), miss=float(tw['miss']),
                           t_launch=float(st['tL']), t_end=float(tw['t_end']), potential=int(c['pot']),
                           counts=C.count_classes(W, [a for a, _ in new]), depth=int(n_pre + len(new)) if craft is not None else int(len(new)),
                           score=float(c.get('score', 0.0)), job_id=job_id)
        # plan_next (09-24, additive): the plan's flybys beyond the committed part, in plan order (the driver's
        # --claims turns them into the craft's claims for the next slices); never a committed or prefix target
        com = set(col['targets']) | pre
        col['plan_next'] = [int(a) for a in c.get('plan_next', []) if int(a) not in com]
        C.check_column(col, T, H)
        cols.append(col)
        plan[cid] = dict(plan_n=int(c['plan_n']), plan_fuel_kg=round(float(c['plan_fuel']), 1), plan_t_end_d=round(C.s2d(c['plan_t_end']), 1),
                         commit_fuel_kg=round(float(c['fuel_c']), 1), cut=bool(c.get('cut', False)),
                         plan_cont=int(c.get('pot_cont', pot_look)), plan_lookahead=pot_look,
                         plan_full=int(c.get('pot_full', c.get('pot_cont', pot_look))), pot_mode=str(pot_mode),
                         potential=int(c['pot']),
                         t_launch_d=round(C.s2d(c['t_launch']), 1) if c.get('t_launch') is not None else None,
                         plan_next=list(col['plan_next']))
    cols.sort(key=lambda c: (-c['n_new'], c['tank']))
    return cols, plan


def frontier_of(cols):
    out = []
    for n in sorted(set(c['n_new'] for c in cols)):
        lst = [c for c in cols if c['n_new'] == n]
        cheap = min(lst, key=lambda c: c['dtank'])
        out.append(dict(n_new=n, dtank_min=round(cheap['dtank'], 1), potential_of_cheapest=cheap['potential'],
                        potential_max=max(c['potential'] for c in lst), n_cols=len(lst)))
    return out


# ------------------------------------------------------------------ the job
def run_job(job):
    """Execute one JOB (s18_segbeam schema) in THIS process, single-threaded; write out/result.json; return RESULT.
    Idempotent (result.json returned as is; a cached failure re-run once, s18_segbeam semantics); exceptions -> ok False."""
    rd = pathlib.Path(job['run_dir']); od = rd / job['out']; od.mkdir(parents=True, exist_ok=True)
    rf = od / 'result.json'; retried = False
    if rf.exists():
        R0 = C.jload(rf)
        if R0.get('ok') or R0.get('retried'):
            R0['cached'] = True
            return R0
        retried = True
        pk = od / 'plan.pkl'
        if pk.exists():
            try:
                with open(pk, 'rb') as f:
                    pickle.load(f)
            except Exception:
                pk.unlink()
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
    stats = dict(settles=0, settle_ok=0, price_s=0.0, settle_s=0.0, say=say)
    try:
        RI.eph(); W = C.load_windows()
        prm = pb_params(job)
        T = float(job['T']); H = float(job['H']); L = float(job['L'])
        remaining = sorted(set(int(x) for x in job['remaining']) - set(UNREACHABLE))
        quota = {k: int(v) for k, v in (job.get('quota') or {}).items() if v is not None} or None
        cls = {t: C.class_of(W, t) for t in remaining}
        if quota:
            full = [k for k, q in quota.items() if q <= 0]
            if full:
                remaining = [t for t in remaining if cls[t] not in full]
                say(f'quota: classes at quota {full} excluded from the pool ({len(remaining)} targets left)')
        job = dict(job, remaining=remaining)
        prize = np.zeros(300); pz = job.get('prize') or {}
        for t in remaining:
            prize[t - 1] = float(pz.get(str(t), pz.get(t, 1.0)))
        t_plan = min(T_MISSION - 2 * DAY, T + H + L + C.d2s(prm['plan_extra_d']))
        say(f'job {job_id} ({kind}, craft {craft}, planbeam): T {C.s2d(T):.0f} d H {C.s2d(H):.0f} L {C.s2d(L):.0f}; plan horizon '
            f'{C.s2d(t_plan):.0f} d; remaining {len(remaining)}; prize >1 on {int(np.sum(prize > 1.0))}; width {prm["width"]} npt '
            f'{prm["npt"]} mc {prm["mc"]} dv {prm["dv"]} drmax {prm["drmax"]} plan_wall {prm["plan_wall"]:.0f} s; quota {quota}')
        P = planner_params(prm, prize, t_plan, T + SB.EDGE_LO)     # legs end in [T + 5 d, t_plan]
        job_wall = float(prm['job_wall']); deadline = tic + job_wall
        if kind == 'craft':
            route = job['route']; rp = pathlib.Path(route)
            start = SB.prefix_state(rp if rp.is_absolute() else rd / rp)
            prefix_asts = list(start['asts'])
            say(f'prefix {route}: {len(prefix_asts)} flybys, t_end {C.s2d(start["t_end"]):.1f} d, tank {start["tank"]:.1f} kg '
                f'(planner view), miss {start["miss"]:.1f} km, nreg {start["nreg"]}')
        else:
            start = None; prefix_asts = []
        # ---- plan (checkpointed)
        pk = od / 'plan.pkl'
        res = None
        if pk.exists():
            try:
                with open(pk, 'rb') as f:
                    res = pickle.load(f)
                say(f'RESUME: plan.pkl with {len(res["states"])} states (plan {res["plan_s"]:.0f} s)')
            except Exception:
                res = None
        if res is None:
            plan_cap = min(float(prm['plan_wall']), deadline - time.time() - float(prm['settle_reserve']))
            say(f'budget: job_wall {job_wall:.0f} s, {deadline - time.time():.0f} s left before the plan -> plan cap '
                f'{max(1.0, plan_cap):.0f} s (plan_wall {prm["plan_wall"]:.0f}, settle_reserve {prm["settle_reserve"]:.0f})')
            res = plan(kind, start, job, prm, P, T, H, L, say, plan_cap=max(1.0, plan_cap))
            tmp = od / 'plan.pkl.tmp'
            with open(tmp, 'wb') as f:
                pickle.dump(res, f)
            os.replace(tmp, pk)
        states = res['states']; n_pre = res['n_pre']
        for s in states:
            s.score_val = s.score(P)
        hist = []
        for d in sorted(set(len(s.seq) - n_pre for s in states)):
            lst = [s for s in states if len(s.seq) - n_pre == d]
            hist.append(dict(depth=int(d), kept=len(lst), in_H=sum(1 for s in lst if s.t <= T + H + 1e-6),
                             t_end_d=[round(C.s2d(min(s.t for s in lst)), 0), round(C.s2d(max(s.t for s in lst)), 0)],
                             fuel_kg=[round(min(s.fuel for s in lst), 0), round(max(s.fuel for s in lst), 0)]))
        deep = max(states, key=lambda s: (len(s.seq), -s.fuel)) if states else None
        say(f'plan: {len(states)} states, depth {res["end"]["depth"]} new flybys, end {res["end"]["kind"]} ({res["plan_s"]:.0f} s)'
            + (f'; deepest {len(deep.seq) - n_pre} new to {C.s2d(deep.t):.0f} d at planner fuel {deep.fuel:.0f} kg' if deep else ''))
        for h in hist:
            say(f'  depth {h["depth"]:2d}: {h["kept"]:4d} states ({h["in_H"]} in H), t_end {h["t_end_d"][0]:.0f}-{h["t_end_d"][1]:.0f} d, '
                f'fuel {h["fuel_kg"][0]:.0f}-{h["fuel_kg"][1]:.0f} kg')
        # ---- candidates -> frontier -> settle
        cands = candidates(res, T, H, L, quota, cls, prm, say)
        say(f'candidates: {len(cands)} distinct committed target sets (n_new {SB._fmt_range([c["n_new"] for c in cands], str)}, '
            f'potential {SB._fmt_range([c["pot"] for c in cands], str)})')
        chosen = frontier_select(cands, prm)
        # every candidate carries the launch of the plan that won its dedupe (root jobs settle from it); the planner's
        # launch grid is the job's, so the settled column launches inside [T, T + H] and the grid
        say(f'frontier: {len(chosen)} candidates chosen (n fb / committed planner kg / continuation (lookahead) [full]; '
            f'pot_mode {prm.get("pot_mode", "lookahead")}, pot_cap {prm.get("pot_cap")}): '
            + ', '.join(f'{c["n_new"]}fb/{c["fuel_c"]:.0f}kg/p{c["pot_cont"]}({c["pot"]})[{c.get("pot_full", "-")}]' for c in chosen))
        settle_left = deadline - time.time()
        say(f'budget: {settle_left:.0f} s of {job_wall:.0f} left for the settles + regrids (settle_wall {prm["settle_wall"]:.0f})')
        settled = settle_columns(chosen, kind, start, T, H, prm, say, stats, deadline=deadline)
        cols, plan_info = columns_of(settled, T, H, L, craft, float(job.get('tank_prev', C.TANK_DRY)), W, rd, job['out'], job_id,
                                     prefix_asts, n_pre, pot_mode=str(prm.get('pot_mode', 'lookahead')), pot_cap=prm.get('pot_cap'))
        R['columns'] = cols; R['ok'] = True
        fr = frontier_of(cols)
        R['stats'] = dict(levels=int(res['end']['depth']), n_states=len(states), n_columns=len(cols), settles=stats['settles'],
                          settle_ok=stats['settle_ok'], regrids=stats.get('regrids', 0), regrid_ok=stats.get('regrid_ok', 0),
                          honest=bool(int(prm.get('regrid_cols', 0))),
                          price_s=round(stats['price_s'], 1), settle_s=round(stats['settle_s'], 1),
                          plan_s=round(res['plan_s'], 1), wall_s=round(time.time() - tic, 1), end=res['end'],
                          budget=dict(job_wall=job_wall, plan_cap=round(float(res.get('plan_cap', prm['plan_wall'])), 1),
                                      settle_left_s=round(settle_left, 1), hit=stats.get('budget_hit'),
                                      n_unsettled=int(stats.get('n_unsettled', 0))),
                          hist=hist, plan=plan_info, frontier=fr, n_candidates=len(cands), params=prm)
        stats.pop('say', None)
        say('RESULT ' + job_id + f': {len(cols)} columns; frontier ' + ' | '.join(
            f'{f["n_new"]}fb {f["dtank_min"]:.0f} kg p{f["potential_of_cheapest"]}/{f["potential_max"]}' for f in fr)
            + f'; settles {stats["settles"]} ({stats["settle_ok"]} ok), regrids {stats.get("regrids", 0)} ({stats.get("regrid_ok", 0)} ok); '
            f'wall {time.time() - tic:.0f} s of {job_wall:.0f} (plan {res["plan_s"]:.0f}, '
            f'price {stats["price_s"]:.0f}, settle {stats["settle_s"]:.0f}; budget hit {stats.get("budget_hit")}, '
            f'{stats.get("n_unsettled", 0)} unsettled)')
    except Exception:
        R['error'] = traceback.format_exc()[-2000:]
        R['stats'] = dict(wall_s=round(time.time() - tic, 1))
        say(f'ERROR {job_id}: {R["error"].strip().splitlines()[-1]}')
    C.jdump(R, rf)
    return R


def make_job(*a, **kw):
    """The JOB schema is s18_segbeam's (the driver builds jobs with s18_segbeam.make_job)."""
    return SB.make_job(*a, **kw)


# ------------------------------------------------------------------ test-slice
def _seg_frontier(cols, craft):
    out = {}
    for c in cols:
        if c.get('craft') != craft:
            continue
        n = c['n_new']; o = out.get(n)
        if o is None or c['dtank'] < o['dtank_min']:
            out[n] = dict(dtank_min=c['dtank'], potential_of_cheapest=c['potential'], potential_max=c['potential'], n_cols=0)
        out[n]['potential_max'] = max(out[n]['potential_max'], c['potential']); out[n]['n_cols'] += 1
    return out


def _run_one(job):
    return run_job(job)


def _regrid_check(rd, cols, tank_prev, timeout=120.0):
    """Honest view of the columns (the commit's regrid law): {column id: (ok, honest dtank, secs)}."""
    out = {}
    for c in cols:
        try:
            st = C.load_twin(pathlib.Path(rd) / c['npz']); r = C.regrid_settle(C.twin_of(st), timeout=timeout)
            out[c['id']] = dict(ok=bool(r['ok']), dtank=round(float(r['tank']) - float(tank_prev), 1), miss=round(float(r['miss']), 1),
                                secs=r['secs'], planner_dtank=round(float(c['dtank']), 1))
        except Exception as e:
            out[c['id']] = dict(ok=False, error=repr(e)[:80])
    return out


def test_slice(out, src_run, k, crafts, nproc=2, overrides=None, regrid=False):
    """Rebuild the slice-k craft jobs of src_run in OUT (routes copied), run run_job for `crafts` on <= nproc workers,
    print the frontier table next to segbeam's (src_run/slices/sKK/columns.json).  Writes OUT/frontier.json.
    regrid: also C.regrid_settle every column twin (the commit's honest view) and print honest vs planner dtank."""
    C.set_threads(); tic0 = time.time()
    src = C.run_dir(src_run); rd = C.run_dir(out); rd.mkdir(parents=True, exist_ok=True)
    sd_src = C.slice_dir(src, k); sd = C.slice_dir(rd, k)
    RI.eph(); C.load_windows()
    jobs = []
    for i in crafts:
        jf = sd_src / 'jobs' / f's{k:02d}_c{i}.json'
        job = C.jload(jf)
        route = job['route']
        if route:
            dst = rd / route; dst.parent.mkdir(parents=True, exist_ok=True)
            if not dst.exists():
                shutil.copy2(src / route, dst)
        job['run_dir'] = str(rd); job['out'] = f'slices/s{k:02d}/c{i}'
        if overrides:
            job['planbeam'] = dict(overrides)
        C.jdump(job, sd / 'jobs' / f'{job["job_id"]}.json'); jobs.append(job)
    T = float(jobs[0]['T']); H = float(jobs[0]['H']); L = float(jobs[0]['L'])
    print(f'[{C.now()}] test-slice {k}: T {C.s2d(T):.0f} d H {C.s2d(H):.0f} L {C.s2d(L):.0f}; crafts {crafts}; nproc {nproc}', flush=True)
    if nproc > 1 and len(jobs) > 1:
        import multiprocessing as mp
        ctx = mp.get_context('fork')
        with ctx.Pool(min(int(nproc), len(jobs)), maxtasksperchild=1) as pool:
            results = pool.map(_run_one, jobs, chunksize=1)
    else:
        results = [run_job(j) for j in jobs]
    seg = C.jload(sd_src / 'columns.json') if (sd_src / 'columns.json').exists() else []
    seg_res = {}
    for i in crafts:
        rf = sd_src / f'c{i}' / 'result.json'
        if rf.exists():
            seg_res[i] = C.jload(rf)
    table = []
    lines = [f'frontier (slice {k}, T {C.s2d(T):.0f} d, commit (T, T+H] = ({C.s2d(T):.0f}, {C.s2d(T + H):.0f}] d): cheapest dtank '
             f'[kg] per n_new, potential (of the cheapest / max); planbeam vs segbeam',
             f'  craft n_new | planbeam dtank  pot(cheap/max) ncol | segbeam dtank  pot(cheap/max) ncol']
    for i, R in zip(crafts, results):
        pf = {f['n_new']: f for f in (R.get('stats') or {}).get('frontier', [])}
        sf = _seg_frontier(seg, i)
        ws = (seg_res.get(i, {}).get('stats') or {}).get('wall_s')
        row = dict(craft=i, ok=R['ok'], error=(R.get('error') or '')[-200:] or None, wall_s=(R.get('stats') or {}).get('wall_s'),
                   plan_s=(R.get('stats') or {}).get('plan_s'), n_columns=len(R.get('columns') or []),
                   planbeam=[pf[n] for n in sorted(pf)], segbeam={n: sf[n] for n in sorted(sf)}, seg_wall_s=ws)
        table.append(row)
        for n in sorted(set(pf) | set(sf)):
            a = pf.get(n); b = sf.get(n)
            lines.append(f'  c{i}   {n:3d}   | ' + (f'{a["dtank_min"]:8.1f}   {a["potential_of_cheapest"]}/{a["potential_max"]}   {a["n_cols"]:3d}' if a else '       -             -    -')
                         + '   | ' + (f'{b["dtank_min"]:8.1f}   {b["potential_of_cheapest"]}/{b["potential_max"]}   {b["n_cols"]:3d}' if b else '       -             -    -'))
        lines.append(f'  c{i} wall: planbeam {row["wall_s"]} s (plan {row["plan_s"]} s, {row["n_columns"]} columns, ok {R["ok"]}'
                     + (f', error {row["error"]}' if row['error'] else '') + f'); segbeam {ws} s')
        if regrid and R.get('columns'):
            job = [j for j in jobs if j['craft'] == i][0]
            rg = _regrid_check(rd, R['columns'], job['tank_prev']); row['regrid'] = rg
            n_ok = sum(1 for v in rg.values() if v.get('ok'))
            lines.append(f'  c{i} regrid (honest view of the {len(rg)} column twins): {n_ok} ok; planner -> honest dtank: '
                         + ', '.join(f'{v["planner_dtank"]}->{v["dtank"]}' + ('' if v['ok'] else '(FAIL)') for v in rg.values() if 'dtank' in v))
    txt = '\n'.join(lines)
    print(txt, flush=True)
    C.jdump(dict(slice=k, T_d=C.s2d(T), H_d=C.s2d(H), L_d=C.s2d(L), crafts=crafts, rows=table, total_wall_s=round(time.time() - tic0, 1),
                 built=C.now()), rd / 'frontier.json')
    with open(rd / 'frontier.txt', 'w') as f:
        f.write(txt + '\n')
    print(f'[{C.now()}] test-slice done: total {time.time() - tic0:.0f} s; {rd / "frontier.json"}', flush=True)
    return table


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['run', 'test-slice']); ap.add_argument('arg'); ap.add_argument('src', nargs='?')
    ap.add_argument('k', nargs='?', type=int)
    ap.add_argument('--crafts', default='0'); ap.add_argument('--nproc', type=int, default=2)
    ap.add_argument('--width', type=int, default=None); ap.add_argument('--plan-wall', type=float, default=None)
    ap.add_argument('--max-cols', type=int, default=None); ap.add_argument('--plan-extra', type=float, default=None)
    ap.add_argument('--w-t', type=float, default=None); ap.add_argument('--w-fuel', type=float, default=None)
    ap.add_argument('--first-leg', default=None, choices=[None, 'relaxed', 'lintwin', 'both', 'prod'])
    ap.add_argument('--job-wall', type=float, default=None); ap.add_argument('--pot-mode', default=None, choices=[None, 'lookahead', 'continuation', 'full'])
    ap.add_argument('--pot-cap', type=int, default=None, help='cap on the emitted column potential (PB_DEF pot_cap)')
    ap.add_argument('--regrid', action='store_true', help='test-slice: regrid_settle every column (honest dtank)')
    a = ap.parse_args()
    if a.cmd == 'run':
        r = run_job(C.jload(a.arg))
        print(json.dumps(dict(job_id=r['job_id'], ok=r['ok'], n_columns=len(r['columns']),
                              stats={k: v for k, v in (r.get('stats') or {}).items() if k not in ('hist', 'plan')}, error=r.get('error')), default=float))
    else:
        ov = {k: v for k, v in dict(width=a.width, plan_wall=a.plan_wall, max_cols=a.max_cols, plan_extra_d=a.plan_extra,
                                    w_t=a.w_t, w_fuel=a.w_fuel, first_leg=a.first_leg, job_wall=a.job_wall, pot_mode=a.pot_mode,
                                    pot_cap=a.pot_cap).items() if v is not None}
        test_slice(a.arg, a.src, int(a.k), [int(x) for x in a.crafts.split(',')], a.nproc, ov, a.regrid)


if __name__ == '__main__':
    main()
