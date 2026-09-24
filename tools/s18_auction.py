"""Stage 18 / WP2: the AUCTION -- one exact assignment of columns to craft per slice (scipy.optimize.milp, HiGHS).

Given the columns of a slice (s18_common.COLUMN_SCHEMA; from every craft job and the root jobs) and the run STATE
(s18_common.STATE_SCHEMA), choose at most one column per launched craft and at most n_unlaunched root columns so that
no target is taken twice, maximising
    sum_c x_c [ sum_{t in new(c)} (1 + edf(t)) + gamma * potential_c - lam * dtank_c / 100 ]
with  edf(t) = s18_common.edf(W, t, T + H, beta)   (earliest deadline first: beta / beta/2 / beta/4 / 0 for 0 / 1 / 2 /
more windows after T + H), gamma = 0.3, lam = 0.15 (1 target ~ 670 kg), beta = 1.5 by default.

Pre-filter (dropped columns are returned with a reason, never silently):
  quota     cumulative class counts of the craft after adoption: fil <= qf (12), med <= qm (18); pins unlimited; a root
            column starts from zero counts; NOT applied in the last quota_free_slices slices (C.quota_active);
  fuel      column tank (planner twin) <= m0max (1070 kg = the fleet bar / N): the honest tank at commit is re-checked;
  ramp      column tank <= C.ramp_cap(T + H_eff, m0max) (fixer 09-24, finding 3: a craft may not spend its 15-yr budget
            in year 4; reason 'ramp');
  stale     a new target no longer in state['remaining'] (cannot happen inside one slice, but the resume path may
            re-run an auction after a partial commit) -> drop;
  horizon   epochs outside (T, T + H] or any other check_column failure -> drop (the reason carries the message);
  empty     n_new = 0 -> drop (choosing nothing is the free alternative);
  craft     a craft column of an unlaunched / unknown craft, or a root column when no craft is unlaunched -> drop
            (reasons 'unlaunched_craft', 'unknown_craft', 'no_unlaunched').
MILP:  x_c in {0,1};  per launched craft i: sum_{c: craft = i} x_c <= 1;  root: sum_{c: craft None} x_c <= n_unlaunched;
       per target: sum_{c: t in new(c)} x_c <= 1;  fleet fuel (fixer 09-24): sum_c x_c max(dtank_c, 0) <=
       C.slice_fuel_budget(H_eff) - fuel_used (fuel_used = what the slice's earlier commit rounds already spent).
       lam = S['lam_dual'] when the state carries it (the driver raises it while the budget binds), else params lam.
       H = S['H_eff_d'] (the slice's effective horizon, days) when present, else params H.
       LAUNCH CAP (09-24, additive, default off): with params max_launch = K > 0 the root row's rhs is
       min(K - launches_used, n_unlaunched) (launches_used = root columns already accepted by earlier commit rounds of
       the slice); a root column offered when that cap is 0 is dropped with reason 'max_launch'.  LAUNCH STAGGER
       (params launch_stagger = D days > 0): the k-th launch of the slice by epoch (k = 0, 1, ..; it goes to the k-th
       unlaunched craft, see the root assignment) must satisfy t_launch >= T + (k + launches_used) D, written as the
       cumulative rows  sum_{root c: t_launch_c < T + (k + used) D} x_c <= k  for k = 1 .. cap - 1 (rows 'stagger').
       Both are off (K = 0 / None, D = 0) unless the driver's --max-launch / --launch-stagger set them.
       PRIZE BONUS (09-24, additive, default off): params prize_bonus = {target: bonus} (the driver's --prize-json,
       frozen in state['params']) is added to the objective coefficient of every column per new target it contains
       (value_of: 1 + edf + bonus); absent / {} = the original objective.
       Solved with scipy milp (time_limit tlim, default 60 s; the problem
       is < 2000 columns x 300 rows, seconds).  A time-limited but feasible incumbent is accepted and flagged.
Root assignment: chosen root columns go to the unlaunched craft in index order (lowest index first; roots ordered by
launch epoch).

CHOICE (slices/s<kk>/auction.json)
  {chosen: {"<craft>": column id}, columns: {"<craft>": COLUMN}, objective, status ("optimal"|"time_limit"|"infeasible"|
   "empty"|"trivial"), n_columns, n_kept, dropped: [{id, reason}], n_vars, n_rows, solve_s, values: {"<column id>": value}}
  'trivial' = no constraint row was needed (every kept column independent): x_c = 1 iff value > 0, no solver call.

SELFTEST (hand-made instance, no twins; beta = 0, gamma = 0.3, lam = 0.15, qf = 10, qm = 16, m0max = 1120)
  N = 3, craft 0 and 1 launched, craft 2 unlaunched; remaining = TARGETS; T = 0, H = 540 d.
  columns (new targets; tank = 601.5 + 12 kg per target unless stated):
    c0a  craft 0  {5, 6}                value 2 - 0.036
    c0b  craft 0  {5, 7, 8} potential 2 value 3 + 0.6 - 0.054   (shares 5 with c0a, 7 with c1a)
    c1a  craft 1  {7}                   value 1 - 0.018
    r1   root     {9, 10}               value 2 - 0.036
    r2   root     {9}                   value 1 - 0.018          (conflict with r1 on 9)
    c1q  craft 1  11 filler-class targets            -> dropped: quota (fil 11 > qf 10)
    c0f  craft 0  {11, 12} tank 1300                 -> dropped: fuel (1300 > 1120)
  optimum: c0b + r1 = 5.51 (the alternative c0a + c1a + r1 = 4.91); chosen == {'0': 'c0b', '2': 'r1'} (r1 goes to the
  unlaunched craft 2), dropped == {c1q: quota, c0f: fuel}, craft 1 gets nothing (c1a conflicts on 7).

usage (tests): s18_auction.py selftest          (hand-made columns with a known optimum; < 5 s)
               s18_auction.py choose RUN_DIR SLICE   (re-run the auction of a slice from its columns.json; prints the choice)
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import csr_matrix
import s18_common as C

STATUS = {0: 'optimal', 1: 'time_limit', 2: 'infeasible', 3: 'unbounded', 4: 'other'}


def _p(prm, k):
    """Parameter with the C.DEF fallback."""
    return prm.get(k, C.DEF[k]) if prm is not None else C.DEF[k]


def edf_bonus(W, ast, T, H, beta):
    """s18_common.edf(W, ast, T + H, beta): bonus of taking `ast` in the slice that commits (T, T + H]."""
    return C.edf(W, int(ast), float(T) + float(H), float(beta))


def H_of(S, prm):
    """The slice's effective horizon [s]: S['H_eff_d'] (set by the driver, finding 2) else params H."""
    return C.d2s(S['H_eff_d']) if S.get('H_eff_d') is not None else C.d2s(_p(prm, 'H'))


def lam_of(S, prm):
    """The fuel price: S['lam_dual'] (the driver's dual of the slice budget) else params lam."""
    return float(S['lam_dual']) if S.get('lam_dual') is not None else float(_p(prm, 'lam'))


def launch_cap(prm, n_unl, launches_used=0):
    """Root columns the auction may still choose in this slice: n_unlaunched, capped at max_launch - launches_used when
    prm['max_launch'] > 0 (the driver's --max-launch; 0 / None / absent = no cap, the original behaviour)."""
    K = prm.get('max_launch') if prm is not None else None
    n = int(n_unl)
    if K is None or int(K) <= 0:
        return n
    return max(0, min(n, int(K) - int(launches_used)))


def launch_stagger(prm):
    """The launch stagger D [s] (prm['launch_stagger'] in DAYS; 0 / None / absent = off)."""
    D = prm.get('launch_stagger') if prm is not None else None
    return C.d2s(D) if D is not None and float(D) > 0 else 0.0


def prize_bonus_of(prm):
    """{target id (str): bonus} of the driver's --prize-json (prm['prize_bonus'], frozen in state['params']; 09-24,
    additive): added to the auction prize of those targets in every slice.  Absent / None / {} = off."""
    b = prm.get('prize_bonus') if prm is not None else None
    return {str(k): float(v) for k, v in b.items()} if b else {}


def value_of(col, S, W, prm):
    """Objective coefficient of a column: sum over new targets of (1 + edf + prize_bonus) + gamma * potential - lam *
    dtank / 100.  prm = S['params'] (beta, gamma; lam = lam_of(S, prm); H = H_of(S, prm); prize_bonus = the
    --prize-json table, empty by default)."""
    T = float(S['T']); H = H_of(S, prm)
    beta = float(_p(prm, 'beta')); gamma = float(_p(prm, 'gamma')); lam = lam_of(S, prm)
    v = sum(1.0 + edf_bonus(W, t, T, H, beta) for t in col['targets'])
    bonus = prize_bonus_of(prm)
    if bonus:
        v += sum(bonus.get(str(int(t)), 0.0) for t in col['targets'])
    v += gamma * float(col.get('potential', 0) or 0)
    v -= lam * float(col.get('dtank', 0.0) or 0.0) / 100.0
    return float(v)


def filter_columns(columns, S, prm, W, launches_used=0):
    """-> (kept [COLUMN], dropped [{id, reason}]) applying quota / fuel / stale / horizon / empty (module docstring).
    Quota uses S['craft'][i]['counts'] + col['counts'] (root: col['counts']).  launches_used: root columns already
    accepted in this slice (retry rounds); with a launch cap of 0 left every root column is dropped ('max_launch')."""
    T = float(S['T']); H = H_of(S, prm)
    qf = int(_p(prm, 'qf')); qm = int(_p(prm, 'qm')); m0max = float(_p(prm, 'm0max'))
    q_on = C.quota_active(T, H, prm); rcap = C.ramp_cap(T + H, m0max, prm=prm)
    remaining = set(int(t) for t in S['remaining'])
    craft = {int(c['i']): c for c in S['craft']}
    n_unl = sum(1 for c in S['craft'] if not c['launched'])
    n_cap = launch_cap(prm, n_unl, launches_used)
    kept, dropped = [], []
    for col in columns:
        cid = col.get('id')
        # schema + horizon
        try:
            C.check_column(col, T, H)
        except AssertionError as e:
            dropped.append(dict(id=cid, reason='horizon', detail=str(e)[:80])); continue
        if col['n_new'] == 0:
            dropped.append(dict(id=cid, reason='empty')); continue
        ci = col.get('craft')
        if ci is None:
            if n_unl == 0:
                dropped.append(dict(id=cid, reason='no_unlaunched')); continue
            if n_cap <= 0:
                dropped.append(dict(id=cid, reason='max_launch', detail=f'launch cap {prm.get("max_launch")} used {int(launches_used)}')); continue
            base = dict(pin=0, med=0, fil=0)
        else:
            ci = int(ci)
            if ci not in craft:
                dropped.append(dict(id=cid, reason='unknown_craft')); continue
            if not craft[ci]['launched']:
                dropped.append(dict(id=cid, reason='unlaunched_craft')); continue
            base = craft[ci]['counts']
        tot = C.add_counts(base, col['counts'])
        if q_on and (tot['fil'] > qf or tot['med'] > qm):
            dropped.append(dict(id=cid, reason='quota', detail=f"fil {tot['fil']}/{qf} med {tot['med']}/{qm}")); continue
        if float(col['tank']) > m0max:
            dropped.append(dict(id=cid, reason='fuel', detail=f"tank {float(col['tank']):.1f} > {m0max:.0f}")); continue
        if float(col['tank']) > rcap:
            dropped.append(dict(id=cid, reason='ramp', detail=f"tank {float(col['tank']):.1f} > ramp {rcap:.0f}")); continue
        stale = [int(t) for t in col['targets'] if int(t) not in remaining]
        if stale:
            dropped.append(dict(id=cid, reason='stale', detail=f'{stale[:6]}')); continue
        kept.append(col)
    return kept, dropped


def build_milp(kept, S, W, prm, fuel_used=0.0, launches_used=0):
    """-> dict(c (objective, to MAXIMISE), A (csr), lb, ub, meta) for scipy milp (which minimises: pass -c).
    Rows: one per launched craft that has columns, one root row (if root columns and n_unlaunched > 0, rhs
    launch_cap(prm, n_unlaunched, launches_used) = n_unlaunched unless --max-launch caps it; if n_unlaunched == 0 root
    columns must have been dropped by filter_columns), one per target that appears in >= 2 kept columns (targets in
    exactly one column need no row), ONE fleet fuel row (fixer 09-24, finding 3): sum_c x_c max(dtank_c, 0) <=
    C.slice_fuel_budget(H_eff, prm) - fuel_used (meta['fuel_budget']), and, with --launch-stagger D > 0, the
    cumulative launch rows of the module docstring (kinds ('stagger', k))."""
    n = len(kept)
    c = np.array([value_of(col, S, W, prm) for col in kept], float)
    rows, cols, ub, kinds = [], [], [], []
    r = 0
    # craft rows
    by_craft = {}
    for j, col in enumerate(kept):
        by_craft.setdefault(col.get('craft'), []).append(j)
    n_unl = sum(1 for cc in S['craft'] if not cc['launched'])
    n_cap = launch_cap(prm, n_unl, launches_used)
    for ci in sorted(k for k in by_craft if k is not None):
        for j in by_craft[ci]:
            rows.append(r); cols.append(j)
        ub.append(1.0); kinds.append(('craft', int(ci))); r += 1
    if None in by_craft:
        for j in by_craft[None]:
            rows.append(r); cols.append(j)
        ub.append(float(n_cap)); kinds.append(('root', None)); r += 1
    # target rows (only targets in >= 2 columns)
    by_t = {}
    for j, col in enumerate(kept):
        for t in col['targets']:
            by_t.setdefault(int(t), []).append(j)
    for t in sorted(by_t):
        js = by_t[t]
        if len(js) < 2:
            continue
        for j in js:
            rows.append(r); cols.append(j)
        ub.append(1.0); kinds.append(('target', t)); r += 1
    vals = [1.0] * len(rows)
    # fleet fuel row (knapsack): coefficients max(dtank, 0) in kg
    budget = C.slice_fuel_budget(H_of(S, prm), prm) - float(fuel_used)
    dt = [max(0.0, float(col.get('dtank', 0.0) or 0.0)) for col in kept]
    fuel_row = None
    if any(d > 0 for d in dt):
        for j, d in enumerate(dt):
            if d > 0:
                rows.append(r); cols.append(j); vals.append(d)
        ub.append(max(0.0, budget)); kinds.append(('fuel', None)); fuel_row = r; r += 1
    # launch stagger rows (--launch-stagger D > 0): at most k root columns launch before T + (k + used) D, k = 1 .. cap-1
    Ds = launch_stagger(prm); n_stagger = 0
    if Ds > 0 and None in by_craft and n_cap >= 2:
        T = float(S['T']); used = int(launches_used)
        for k in range(1, n_cap):
            thr = T + (k + used) * Ds
            js = [j for j in by_craft[None] if float(kept[j].get('t_launch') if kept[j].get('t_launch') is not None else T) < thr - 1e-6]
            if len(js) <= k:
                continue                                  # the row can never bind
            for j in js:
                rows.append(r); cols.append(j); vals.append(1.0)
            ub.append(float(k)); kinds.append(('stagger', k)); r += 1; n_stagger += 1
    A = csr_matrix((np.array(vals, float), (rows, cols)), shape=(r, n)) if r else csr_matrix((0, n))
    return dict(c=c, A=A, lb=np.zeros(r), ub=np.array(ub, float), meta=dict(n_vars=n, n_rows=r, kinds=kinds,
                n_unlaunched=n_unl, n_target_rows=sum(1 for k in kinds if k[0] == 'target'),
                fuel_budget=round(budget, 1), fuel_row=fuel_row, launch_cap=n_cap, n_stagger_rows=n_stagger))


def _assign(kept, x, S):
    """chosen / columns maps from a 0/1 vector: craft columns to their craft; root columns to the unlaunched craft in
    index order (roots sorted by launch epoch)."""
    chosen, columns = {}, {}
    unl = [int(c['i']) for c in S['craft'] if not c['launched']]
    roots = []
    for j, col in enumerate(kept):
        if x[j] < 0.5:
            continue
        if col.get('craft') is None:
            roots.append(col)
        else:
            k = str(int(col['craft']))
            assert k not in chosen, f'two columns for craft {k}'
            chosen[k] = col['id']; columns[k] = col
    roots.sort(key=lambda c: (float(c['t_launch'] if c.get('t_launch') is not None else 0.0), c['id']))
    assert len(roots) <= len(unl), 'more root columns than unlaunched craft'
    for col, i in zip(roots, unl):
        chosen[str(i)] = col['id']; columns[str(i)] = col
    return chosen, columns


def choose(columns, S, W=None, prm=None, tlim=60.0, say=None, fuel_used=0.0, launches_used=0):
    """The slice auction (module docstring).  Returns the CHOICE dict.  Never raises on an empty column list (status
    'empty', chosen {}).  W defaults to C.load_windows(), prm to S['params'].  fuel_used [kg]: propellant already
    committed in this slice by earlier rounds (the retry auction of the driver's commit); it reduces the fleet row.
    launches_used: root columns already accepted by those rounds (counts against --max-launch / --launch-stagger).
    The CHOICE also carries lam, H_eff_d, fuel_budget, fuel_chosen (sum of the chosen columns' dtank), n_ramp
    (columns dropped by the ramp) so that the driver can move the dual, and launch_cap (root columns allowed)."""
    tic = time.time()
    W = W if W is not None else C.load_windows()
    prm = prm if prm is not None else S['params']
    say = say or (lambda s: None)
    columns = list(columns or [])
    kept, dropped = filter_columns(columns, S, prm, W, launches_used=launches_used)
    H = H_of(S, prm)
    n_unl = sum(1 for c in S['craft'] if not c['launched'])
    n_cap = launch_cap(prm, n_unl, launches_used)
    out = dict(chosen={}, columns={}, objective=0.0, status='empty', n_columns=len(columns), n_kept=len(kept),
               dropped=dropped, n_vars=len(kept), n_rows=0, solve_s=0.0, values={}, lam=lam_of(S, prm),
               H_eff_d=round(C.s2d(H), 3), fuel_budget=round(C.slice_fuel_budget(H, prm) - float(fuel_used), 1),
               fuel_chosen=0.0, n_ramp=sum(1 for d in dropped if d['reason'] == 'ramp'),
               quota_active=bool(C.quota_active(float(S['T']), H, prm)),
               launch_cap=int(n_cap), max_launch=prm.get('max_launch'), launch_stagger_d=prm.get('launch_stagger'),
               launches_used=int(launches_used))
    if say:
        from collections import Counter
        say(f'auction: {len(columns)} columns, kept {len(kept)}, dropped {dict(Counter(d["reason"] for d in dropped))}; '
            f'lam {out["lam"]:.2f} H_eff {out["H_eff_d"]:.0f} d fuel budget {out["fuel_budget"]:.0f} kg quota {"on" if out["quota_active"] else "OFF"}'
            + (f' launch cap {n_cap} of {n_unl} unlaunched (max_launch {prm.get("max_launch")}, used {int(launches_used)})'
               if n_cap != n_unl else '')
            + (f' stagger {prm.get("launch_stagger")} d' if launch_stagger(prm) > 0 else '')
            + (f' prize bonus on {len(prize_bonus_of(prm))} targets' if prize_bonus_of(prm) else ''))
    if not kept:
        out['solve_s'] = round(time.time() - tic, 3)
        return out
    M = build_milp(kept, S, W, prm, fuel_used=fuel_used, launches_used=launches_used)
    c, A = M['c'], M['A']
    out['values'] = {col['id']: round(float(v), 6) for col, v in zip(kept, c)}
    out['n_rows'] = int(M['meta']['n_rows'])
    if M['meta']['n_rows'] == 0:
        x = (c > 0).astype(float); status = 'trivial'
    else:
        res = milp(-c, constraints=LinearConstraint(A, M['lb'], M['ub']), integrality=np.ones(len(c)),
                   bounds=Bounds(np.zeros(len(c)), np.ones(len(c))), options=dict(time_limit=float(tlim)))
        status = STATUS.get(int(res.status), 'other')
        if res.x is None:
            x = np.zeros(len(c)); status = status if status != 'optimal' else 'infeasible'
        else:
            x = np.round(np.asarray(res.x, float))
        # never accept a vector that breaks a row (safety; rounding of an incumbent)
        if M['meta']['n_rows'] and np.any(A @ x > M['ub'] + 1e-6):
            say('auction: incumbent violates a row -> discarded'); x = np.zeros(len(c)); status = 'infeasible'
    chosen, cols = _assign(kept, x, S)
    fuel_chosen = float(sum(max(0.0, float(col.get('dtank', 0.0) or 0.0)) for col in cols.values()))
    out.update(chosen=chosen, columns=cols, objective=float(c @ x), status=status, solve_s=round(time.time() - tic, 3),
               fuel_chosen=round(fuel_chosen, 1))
    say(f'auction: status {status} objective {out["objective"]:.3f} vars {len(c)} rows {out["n_rows"]} '
        f'chosen {len(chosen)} fuel {fuel_chosen:.0f}/{out["fuel_budget"]:.0f} kg ({out["solve_s"]} s)')
    return out


def selftest():
    """The hand-made instance of the module docstring (beta = 0 to kill the edf ties); asserts
    chosen == {'0': 'c0b', '2': 'r1'} and dropped == {c1q: quota, c0f: fuel}."""
    tic = time.time()
    W = C.load_windows()
    prm = dict(C.DEF); prm.update(N=3, H=540.0, beta=0.0, gamma=0.3, lam=0.15, qf=10, qm=16, m0max=1120.0,
               ramp_slack=600.0, fuel_bar=1e5)          # the instance predates the ramp / fleet row: keep them slack
    S = C.new_state('selftest', prm, 3)
    for i in (0, 1):
        S['craft'][i].update(launched=True, route=f'craft_{i}.npz', t_launch=0.0, t_end=0.0, tank=650.0)
    C.check_state(S)
    T = 0.0; H = C.d2s(540); d = C.DAY

    def col(cid, craft, targets, tank=None, potential=0, t_launch=10 * d):
        ep = [T + (k + 1) * 40 * d for k in range(len(targets))]     # 11 x 40 d = 440 d < H
        tank = C.TANK_DRY + 12.0 * len(targets) if tank is None else tank
        return C.new_column(id=cid, craft=craft, targets=list(targets), epochs=ep, all_targets=list(targets), npz=None,
                            tank=tank, tank_prev=C.TANK_DRY, miss=1.0, t_launch=t_launch, t_end=ep[-1],
                            potential=potential, counts=C.count_classes(W, targets), job_id='selftest')
    fillers = [t for t in C.TARGETS if C.class_of(W, t) == 'fil' and t > 12][:11]
    assert len(fillers) == 11
    cols = [col('c0a', 0, [5, 6]), col('c0b', 0, [5, 7, 8], potential=2), col('c1a', 1, [7]),
            col('r1', None, [9, 10], t_launch=20 * d), col('r2', None, [9], t_launch=15 * d),
            col('c1q', 1, fillers), col('c0f', 0, [11, 12], tank=1300.0)]
    ch = choose(cols, S, W, prm, tlim=10.0, say=print)
    print(json.dumps({k: v for k, v in ch.items() if k != 'columns'}, indent=1, default=float))
    dropped = {dd['id']: dd['reason'] for dd in ch['dropped']}
    assert dropped == {'c1q': 'quota', 'c0f': 'fuel'}, dropped
    assert ch['chosen'] == {'0': 'c0b', '2': 'r1'}, ch['chosen']
    assert ch['status'] == 'optimal', ch['status']
    assert abs(ch['objective'] - 5.51) < 1e-6, ch['objective']
    assert ch['columns']['2']['craft'] is None and ch['columns']['0']['id'] == 'c0b'
    # empty input never raises
    e = choose([], S, W, prm, say=print)
    assert e['status'] == 'empty' and e['chosen'] == {}
    # everything conflicting except one -> one column; n_unlaunched 0 -> roots dropped
    S2 = json.loads(json.dumps(S)); S2['craft'][2].update(launched=True, route='craft_2.npz')
    ch2 = choose(cols, S2, W, prm, say=print)
    assert {dd['id']: dd['reason'] for dd in ch2['dropped']}['r1'] == 'no_unlaunched'
    assert ch2['chosen'] == {'0': 'c0b'}, ch2['chosen']
    # --max-launch / --launch-stagger (09-24): craft 1 and 2 unlaunched; roots r1 {9,10} @ 20 d, r3 {13,14,15} @ 15 d,
    # r4 {16} @ 200 d (independent targets).  No cap: r1 + r3 (2 launches = n_unlaunched; r4 loses to the root row);
    # max_launch 1: r3 alone (the best root); max_launch 2 + stagger 100 d: the 2nd launch by epoch must be >= 100 d ->
    # r3 + r4 (r1 + r3 both launch before 100 d); launches_used 1 with max_launch 2: one root left -> r3.
    S3 = json.loads(json.dumps(S)); S3['craft'][1].update(launched=False, route=None, tank=C.TANK_DRY); C.check_state(S3)
    r4 = col('r4', None, [16], t_launch=200 * d); r4['epochs'] = [260 * d]; r4['t_end'] = 260 * d   # flyby after the launch
    cols3 = [col('c0a', 0, [5, 6]), col('r1', None, [9, 10], t_launch=20 * d), col('r3', None, [13, 14, 15], t_launch=15 * d), r4]
    ch3 = choose(cols3, S3, W, prm, say=print)
    assert ch3['chosen'] == {'0': 'c0a', '1': 'r3', '2': 'r1'}, ch3['chosen']      # roots by launch epoch: r3 (15 d) -> craft 1
    assert ch3['launch_cap'] == 2 and 'stagger' not in {k[0] for k in build_milp(cols3, S3, W, prm)['meta']['kinds']}
    p1 = dict(prm, max_launch=1)
    ch4 = choose(cols3, S3, W, p1, say=print)
    assert ch4['chosen'] == {'0': 'c0a', '1': 'r3'} and ch4['launch_cap'] == 1, ch4['chosen']
    p2 = dict(prm, max_launch=2, launch_stagger=100.0)
    M2 = build_milp([c for c in cols3], S3, W, p2)
    assert ('stagger', 1) in M2['meta']['kinds'] and M2['meta']['n_stagger_rows'] == 1, M2['meta']['kinds']
    ch5 = choose(cols3, S3, W, p2, say=print)
    assert ch5['chosen'] == {'0': 'c0a', '1': 'r3', '2': 'r4'}, ch5['chosen']
    ch6 = choose(cols3, S3, W, p2, say=print, launches_used=1)
    assert ch6['chosen'] == {'0': 'c0a', '1': 'r3'} and ch6['launch_cap'] == 1, ch6['chosen']
    ch7 = choose(cols3, S3, W, p1, say=print, launches_used=1)
    assert ch7['chosen'] == {'0': 'c0a'} and {dd['id']: dd['reason'] for dd in ch7['dropped']} == {'r1': 'max_launch', 'r3': 'max_launch', 'r4': 'max_launch'}, ch7
    # default parameters (no max_launch / launch_stagger keys at all) reproduce the uncapped choice
    ch8 = choose(cols3, S3, W, {k: v for k, v in prm.items() if k not in ('max_launch', 'launch_stagger')}, say=print)
    assert ch8['chosen'] == ch3['chosen'] and ch8['objective'] == ch3['objective']
    # --prize-json (09-24): prize_bonus {6: 3.0} on the original instance: c0a {5, 6} becomes 2 - 0.036 + 3 = 4.964, so
    # c0a + c1a + r1 = 7.91 beats c0b + r1 = 5.51; the values of the other columns are unchanged; an empty / absent
    # table gives the original choice bit-identically.
    pb = dict(prm, prize_bonus={'6': 3.0})
    ch9 = choose(cols, S, W, pb, say=print)
    assert ch9['chosen'] == {'0': 'c0a', '1': 'c1a', '2': 'r1'}, ch9['chosen']
    assert abs(ch9['objective'] - 7.91) < 1e-6, ch9['objective']
    assert abs(ch9['values']['c0a'] - 4.964) < 1e-6 and ch9['values']['c0b'] == ch['values']['c0b'], ch9['values']
    ch10 = choose(cols, S, W, dict(prm, prize_bonus={}), say=print)
    assert ch10['chosen'] == ch['chosen'] and ch10['objective'] == ch['objective'] and ch10['values'] == ch['values']
    assert prize_bonus_of(prm) == {} and prize_bonus_of(dict(prm, prize_bonus=None)) == {}
    print(f'selftest wall {time.time() - tic:.2f} s')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('cmd', choices=['selftest', 'choose']); ap.add_argument('run', nargs='?')
    ap.add_argument('slice', nargs='?', type=int); ap.add_argument('--tlim', type=float, default=60.0)
    a = ap.parse_args()
    if a.cmd == 'selftest':
        selftest(); print('selftest OK')
    else:
        rd = C.run_dir(a.run); S = C.load_state(rd)
        cols = C.jload(C.slice_dir(rd, a.slice) / 'columns.json')
        ch = choose(cols, S, tlim=a.tlim, say=print)
        print(json.dumps({k: v for k, v in ch.items() if k != 'columns'}, indent=1, default=float))


if __name__ == '__main__':
    main()
