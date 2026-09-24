"""Stage 18 / WP3: FIXTURE proposers -- columns without the segment beam, so that the auction and the driver (WP2) can
be built and tested before the beam (WP1) lands.  Same JOB / RESULT contract as s18_segbeam.run_job (the driver calls
`run_job(job)` of the module named by --proposer), so the driver never knows which proposer it runs.

fake     (proposer fixture_fake; use with --commit none)  columns with NO twins: for a craft job, n_cols random
         subsets (1..6 targets) of `remaining` with epochs spread in (T, T+H], tank = tank_prev + 12 kg per target (+-
         noise), potential = Poisson(2), counts from the class table, npz = None; for a root job, the same with
         t_launch in the job's grid and tank_prev = 601.5.  Deterministic in job['seed'] + slice + craft.  Two craft
         jobs share a fraction of targets on purpose (conflicts for the auction): every job of a slice draws half of
         its targets from one slice-level "hot" set of 12 targets (seeded by seed + slice only).  Every job also emits
         one column that violates the filler quota (11 fillers) and one that violates the fuel cap (tank 1300 kg) so
         the auction's filter is exercised; both pass C.check_column (the schema has no quota / fuel).  With the
         environment variable S18_FIXTURE_PLAN_NEXT=K (K > 0; default off) every fake column also carries a fake
         `plan_next` (the K smallest remaining targets outside the column, shifted by one on odd slices) so that the
         driver's --claims bookkeeping (claims, release, ttl, conflicts, resume) can be exercised without twins.
replay   (proposer fixture_replay; --commit regrid)  columns with REAL twins cut from an existing fleet
         (--fixture-fleet, default results/s16/best/fleet; craft i <-> route r<i+1>): the route truncated at T + H
         (every flyby <= T + H; cut_twin = prefix + new flybys, regrid + settle) is the main column, the truncation at
         the flyby before that the second one; potential = the route's flybys in (T + H, T + H + L].  Only the route's
         flybys in (T, T + H] are "new"; targets already taken by another craft (not in `remaining`) are excluded by
         truncating earlier.  A craft whose committed prefix skipped some of the source route's flybys (the auction
         took the shorter column, or a quota dropped a slice) keeps its prefix: the skipped flybys are REMOVED from the
         replayed twin before the settle (cut_twin keeps exactly prefix + new).  Root jobs: the routes of the
         still-unlaunched craft indices (read from run_dir/state.json when it exists, else every route not yet adopted),
         truncated the same way, when their launch lies in the job's grid (s16a launches on days 0-88, r9 on day 944:
         at H = 540 the root jobs of slice 0 propose r1..r8 and r9 appears in slice 1 -- or not at all when it has no
         flyby before T + H; an EMPTY root result is a feature the driver must survive).  Real regrid_settle at commit:
         the G0 smoke test then exercises the whole regrid law with realistic twins in ~3-5 min.

usage: s18_fixture.py fake JOB.json | replay JOB.json     (prints the RESULT summary)
       s18_fixture.py selftest [--fleet results/s16/best/fleet]   (one fake craft job, one fake root job, one replay
                                                                  craft job for r1 at T = 0, H = 540: < 3 min on 1 core)
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse, shutil, traceback
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import numpy as np
import run_ialns as RI
import s18_common as C
from ctoc14.impulsive import ImpulsiveProblem
from ctoc14.constants import DAY

DEFAULT_FLEET = ROOT / 'results/s16/best/fleet'
KG_PER_TARGET = 12.0        # fake tank growth per new target
N_HOT = 12                  # size of the slice-level shared target set (fake conflicts)
QUOTA_VIOLATOR_FIL = 11     # fillers in the deliberate quota violator (qf default 10)
FUEL_VIOLATOR_TANK = 1300.0 # tank of the deliberate fuel violator (m0max default 1120)

_FLEET_CACHE = {}


# ------------------------------------------------------------------ job helpers
def _slice_of(job_id):
    """'s03_c2' -> 3."""
    return int(str(job_id).split('_')[0][1:])


def _idx_of(job_id):
    """'s03_r1' -> 1 (root sub-window index); 's03_c2' -> 2."""
    return int(str(job_id).split('_')[1][1:])


def _state_hint(job):
    """The slice's state when present: the driver's snapshot slices/s<kk>/in/state_in.json (next to the job's out dir),
    else run_dir/state.json.  The fixture uses it for the craft's prefix targets / launch epoch (fake) and for the
    unlaunched indices (replay root).  None when absent (selftest)."""
    rd = pathlib.Path(job['run_dir'])
    for p in (rd / pathlib.Path(job['out']).parent / 'in' / 'state_in.json', rd / 'state.json'):
        try:
            if p.exists():
                return C.jload(p)
        except Exception:
            pass
    return None


def _prefix_of(job, S=None):
    """(targets, epochs, t_launch) of the craft's committed prefix: from the route npz when given and loadable, else
    from the state hint (a --commit none run writes placeholder twins), else empty."""
    rd = pathlib.Path(job['run_dir'])
    if job.get('route'):
        try:
            st = C.load_twin(rd / job['route']); o = C.order_of(st)
            if len(o):
                return ([int(st['asts'][i]) for i in o], [float(st['tf'][i]) for i in o], float(st['tL']))
        except Exception:
            pass
    if S is not None and job.get('craft') is not None:
        c = S['craft'][int(job['craft'])]
        if c['launched']:
            return (list(map(int, c['targets'])), list(map(float, c['epochs'])), float(c['t_launch']))
    return [], [], None


def _stats(n_columns, settles=0, settle_ok=0, settle_s=0.0, wall_s=0.0, depth=0, note=''):
    return dict(levels=1, n_states=int(n_columns), n_columns=int(n_columns), settles=int(settles), settle_ok=int(settle_ok),
                price_s=0.0, settle_s=round(float(settle_s), 2), wall_s=round(float(wall_s), 2),
                end=dict(kind='fixture', depth=int(depth)), hist=[], note=note)


def _result(job, columns, stats, ok=True, error=None):
    return dict(job_id=job['job_id'], kind=job['kind'], craft=job.get('craft'), ok=bool(ok), error=error,
                columns=list(columns), stats=stats)


# ------------------------------------------------------------------ fake
def _rng_of(job):
    """Deterministic generators: (per job: seed + slice + craft/idx, per slice: seed + slice only -> the shared hot set)."""
    seed = int(job.get('seed', 0) or 0); k = _slice_of(job['job_id'])
    who = int(job['craft']) if job.get('craft') is not None else 100 + _idx_of(job['job_id'])
    return np.random.default_rng([seed, k, who]), np.random.default_rng([seed, k, 999])


def _fake_epochs(rng, n, lo_d, hi_d):
    """n distinct integer days in [lo_d, hi_d] (ascending) -> seconds."""
    days = np.arange(int(lo_d), int(hi_d) + 1)
    if len(days) < n:
        return None
    pick = np.sort(rng.choice(days, size=n, replace=False))
    return [float(d) * DAY for d in pick]


def fake_columns(job, W, rng, n_cols=4):
    """[COLUMN] for one JOB without twins (module docstring).  Every column passes C.check_column(c, T, H) except the
    deliberate quota / fuel violators, which pass the schema but are dropped by s18_auction.filter_columns."""
    rng_job, rng_slice = (rng, np.random.default_rng([int(job.get('seed', 0) or 0), _slice_of(job['job_id']), 999])) \
        if not isinstance(rng, tuple) else rng
    T = float(job['T']); H = float(job['H']); kind = job['kind']; craft = job.get('craft')
    rem = sorted(int(t) for t in job['remaining'])
    if not rem:
        return []
    S = _state_hint(job)
    pre_t, pre_e, t_launch = _prefix_of(job, S) if kind == 'craft' else ([], [], None)
    tank_prev = float(job.get('tank_prev', C.TANK_DRY))
    hot = [int(x) for x in rng_slice.choice(rem, size=min(N_HOT, len(rem)), replace=False)]
    hot_set = set(hot); cold = [t for t in rem if t not in hot_set]
    H_d = int(H // DAY)
    cols = []

    def make(k, targets, tank=None, t_launch_col=None, lo_d=1, potential=None, note=''):
        n = len(targets)
        ep = _fake_epochs(rng_job, n, lo_d, H_d)
        if ep is None:
            return None
        ep = [T + e for e in ep]
        tk = tank if tank is not None else tank_prev + KG_PER_TARGET * n + float(rng_job.normal(0.0, 2.0))
        pot = int(rng_job.poisson(2.0)) if potential is None else int(potential)
        c = C.new_column(id=f'{job["job_id"]}_k{k}', craft=craft, targets=[int(t) for t in targets], epochs=ep,
                         all_targets=[int(t) for t in pre_t] + [int(t) for t in targets], npz=None, tank=float(tk),
                         tank_prev=tank_prev, miss=0.0, t_launch=float(t_launch_col), t_end=float(ep[-1]), potential=pot,
                         counts=C.count_classes(W, targets), depth=1, score=-float(n), job_id=job['job_id'])
        if note:
            c['note'] = note
        return c

    def draw(n):
        n_hot = int(rng_job.integers(1, n + 1)) if hot else 0
        n_hot = min(n_hot, len(hot)); n_cold = min(n - n_hot, len(cold))
        pick = [int(x) for x in rng_job.choice(hot, size=n_hot, replace=False)] if n_hot else []
        pick += [int(x) for x in rng_job.choice(cold, size=n_cold, replace=False)] if n_cold else []
        rng_job.shuffle(pick)
        return pick

    if kind == 'craft':
        t_launch_col = t_launch if t_launch is not None else T
        lo_d = 1
    else:
        t_lo, t_hi, step = (float(x) for x in job['grid'])
        # launch epochs t_lo + step, ..., <= t_hi, and at least 40 d of room before T + H for a first flyby
        launches = np.arange(t_lo + step, min(t_hi, T + H - 40 * DAY) + 1e-6, step)
        if len(launches) == 0:
            return []
        t_launch_col = float(rng_job.choice(launches))
        lo_d = int(np.ceil((t_launch_col - T) / DAY)) + 20
        tank_prev = C.TANK_DRY
    k = 0
    for _ in range(int(n_cols)):
        n = int(rng_job.integers(1, 7))
        c = make(k, draw(n), t_launch_col=t_launch_col, lo_d=lo_d)
        if c is not None:
            cols.append(c); k += 1
    # deliberate violators (schema-valid; dropped by the auction's pre-filter)
    fil = [t for t in rem if C.class_of(W, t) == 'fil']
    if len(fil) >= QUOTA_VIOLATOR_FIL:
        tq = [int(x) for x in rng_job.choice(fil, size=QUOTA_VIOLATOR_FIL, replace=False)]
        c = make(k, tq, t_launch_col=t_launch_col, lo_d=lo_d, note='quota violator (11 fillers)')
        if c is not None:
            cols.append(c); k += 1
    c = make(k, draw(int(rng_job.integers(1, 4))), tank=FUEL_VIOLATOR_TANK, t_launch_col=t_launch_col, lo_d=lo_d,
             note='fuel violator (tank 1300)')
    if c is not None:
        cols.append(c); k += 1
    # fake plan_next (09-24; env S18_FIXTURE_PLAN_NEXT = K > 0, default off): every column claims the K smallest
    # remaining targets outside its own targets, shifted by one on odd slices, so the driver's --claims bookkeeping
    # sees persisting claims (ttl expiry), one release + one new claim per slice, and claim conflicts between craft
    # (every job of a slice names the same targets).  Deterministic; the columns themselves are unchanged.
    K = fake_plan_next_k()
    if K > 0:
        off = _slice_of(job['job_id']) % 2
        for c in cols:
            own = set(c['targets'])
            c['plan_next'] = [t for t in rem if t not in own][off:off + K]
    for c in cols:
        C.check_column(c, T, H)
    return cols


def fake_plan_next_k():
    """K of the fake plan_next (module docstring; env S18_FIXTURE_PLAN_NEXT, 0 / unset = off)."""
    try:
        return max(0, int(os.environ.get('S18_FIXTURE_PLAN_NEXT', '0') or 0))
    except ValueError:
        return 0


# ------------------------------------------------------------------ replay
def load_fixture_fleet(fleet_dir):
    """{index (0-based): {name, st, targets, epochs, tL}} of a route_*.npz directory; r<i+1> -> index i when the names
    are r<n>, else sorted order.  Cached per process."""
    fleet_dir = pathlib.Path(fleet_dir)
    if not fleet_dir.is_absolute():
        fleet_dir = ROOT / fleet_dir
    key = str(fleet_dir)
    if key in _FLEET_CACHE:
        return _FLEET_CACHE[key]
    files = sorted(fleet_dir.glob('route_*.npz'), key=lambda f: (len(f.stem), f.stem))
    if not files:
        raise FileNotFoundError(f'no route_*.npz in {fleet_dir}')
    out = {}
    for j, f in enumerate(files):
        name = f.stem[6:]
        idx = int(name[1:]) - 1 if (name[:1] == 'r' and name[1:].isdigit()) else j
        st = C.load_twin(f); o = C.order_of(st)
        out[idx] = dict(name=name, file=str(f), st=st, targets=[int(st['asts'][i]) for i in o],
                        epochs=[float(st['tf'][i]) for i in o], tL=float(st['tL']))
    _FLEET_CACHE[key] = out
    return out


def cut_twin(st, keep, iters=100, timeout=300.0):
    """Generalised C.truncate_twin: keep the flybys with the given indices (into st['tf'] / st['asts']; any subset), the
    impulse nodes up to the last kept flyby, then THE REGRID LAW (C.regrid_settle).  Returns the regrid_settle result
    dict (ip, st, tank, miss, ok, honest, secs, error).  Removing a skipped flyby is a real re-settle (miss <= 150 km
    decides ok)."""
    keep = sorted(int(i) for i in keep)
    if not keep:
        raise ValueError('no flyby to keep')
    ip = C.twin_of(st); tf = np.asarray(st['tf'], float)
    tk = float(tf[keep].max()); m = ip.ts <= tk
    ip2 = ImpulsiveProblem(RI.eph(), ip.tL, ip.vinf, ip.ts[m], ip.Ts[m], tf[keep], [int(st['asts'][i]) for i in keep])
    return C.regrid_settle(ip2, iters=iters, timeout=timeout)


def _match_route(fleet, craft, pre_t, T):
    """The fixture route a craft replays: r<craft+1> when its flybys up to T contain the prefix targets (in order), else
    any route with that property (the auction may have handed a root column to another index); None when nothing fits."""
    order = ([craft] if craft in fleet else []) + [i for i in sorted(fleet) if i != craft]
    for i in order:
        r = fleet[i]
        upto = [t for t, e in zip(r['targets'], r['epochs']) if e <= T + 1e-6]
        if not pre_t:
            if craft is None or i == craft:
                return i
            continue
        # every prefix target is one of the route's flybys up to T, in the route's order
        pos = [upto.index(t) if t in upto else -1 for t in pre_t]
        if all(p >= 0 for p in pos) and pos == sorted(pos):
            return i
    return None


def _replay_route(job, W, r, pre_t, remaining, out, k0, say, craft, tank_prev):
    """The (up to two) columns of one source route for this job: keep = prefix flybys + new flybys in (T, T+H] that are
    still remaining (cut at the first taken one); the second column drops the last new flyby.  Returns (columns, stats)."""
    T = float(job['T']); H = float(job['H']); L = float(job['L']); rd = pathlib.Path(job['run_dir'])
    rem = set(int(t) for t in remaining); pre = set(pre_t)
    tg, ep = r['targets'], r['epochs']
    keep_pre = [i for i, (t, e) in enumerate(zip(tg, ep)) if e <= T + 1e-6 and t in pre]
    new = []
    for i, (t, e) in enumerate(zip(tg, ep)):
        if e <= T + 1e-6 or t in pre:
            continue
        if e > T + H + 1e-6:
            break
        if t not in rem:
            say(f'  {r["name"]}: target {t} at {e / DAY:.0f} d already taken -> truncate before it')
            break
        new.append(i)
    pot = sum(1 for t, e in zip(tg, ep) if T + H < e <= T + H + L and t in rem)
    cols = []; settles = 0; ok_n = 0; ssec = 0.0
    if not new:
        say(f'  {r["name"]}: no new flyby in ({T / DAY:.0f}, {(T + H) / DAY:.0f}] d -> no column')
        return cols, (settles, ok_n, ssec)
    variants = [new] + ([new[:-1]] if len(new) >= 2 else [])
    for v in variants:
        keep = keep_pre + v
        tic = time.time()
        try:
            res = cut_twin(r['st'], keep)
        except Exception as e:
            say(f'  {r["name"]}: cut_twin failed ({repr(e)[:80]})'); settles += 1; ssec += time.time() - tic
            continue
        settles += 1; ssec += res['secs']
        if not res['ok']:
            say(f'  {r["name"]}: settle not ok (miss {res["miss"]:.0f} km, err {res["error"]}) -> column dropped')
            continue
        st2 = res['st']; o = C.order_of(st2)
        all_t = [int(st2['asts'][i]) for i in o]; all_e = [float(st2['tf'][i]) for i in o]
        new_t = [tg[i] for i in v]; new_set = set(new_t)
        new_e = [e for t, e in zip(all_t, all_e) if t in new_set]
        if min(new_e) <= T + 1e-6 or max(new_e) > T + H + 1e-6:
            say(f'  {r["name"]}: a settled epoch moved outside (T, T+H] -> column dropped'); continue
        if [t for t in all_t if t in new_set] != new_t:
            say(f'  {r["name"]}: settled flyby order changed -> column dropped'); continue
        k = k0 + len(cols)
        npz_rel = f'{job["out"]}/col_{k}.npz'
        C.save_twin(rd / npz_rel, st2)
        c = C.new_column(id=f'{job["job_id"]}_k{k}', craft=craft, targets=new_t, epochs=new_e, all_targets=all_t,
                         npz=npz_rel, tank=float(res['tank']), tank_prev=float(tank_prev), miss=float(res['miss']),
                         t_launch=float(st2['tL']), t_end=float(max(new_e)), potential=int(pot),
                         counts=C.count_classes(W, new_t), depth=len(all_t), score=-float(len(new_t)), job_id=job['job_id'])
        c['source'] = r['name']; c['honest_planner_twin'] = bool(res['honest'].get('ok', False))
        C.check_column(c, T, H)
        cols.append(c); ok_n += 1
        say(f'  {r["name"]}: column {c["id"]} new {new_t} tank {c["tank"]:.1f} (+{c["dtank"]:.1f}) miss {c["miss"]:.0f} km '
            f'potential {pot} ({res["secs"]:.1f} s)')
    return cols, (settles, ok_n, ssec)


def replay_columns(job, W, fleet_dir, say):
    """[COLUMN] with real twins from `fleet_dir` (module docstring): truncations of the craft's source route at T + H
    and at the previous flyby; twins saved to run_dir/out/col_<k>.npz; miss <= 150 asserted (cut_twin re-settles).
    Returns [] when the route has no new flyby in (T, T + H] (e.g. r9 before day 944)."""
    fleet = load_fixture_fleet(fleet_dir)
    T = float(job['T']); H = float(job['H']); kind = job['kind']
    remaining = [int(t) for t in job['remaining']]
    S = _state_hint(job)
    cols = []; settles = ok_n = 0; ssec = 0.0
    if kind == 'craft':
        craft = int(job['craft'])
        pre_t, pre_e, _ = _prefix_of(job, S)
        i = _match_route(fleet, craft, pre_t, T)
        if i is None:
            say(f'  craft {craft}: no fixture route matches the prefix {pre_t} -> no column')
            return cols, _stats(0, note='no matching route')
        say(f'  craft {craft}: replaying {fleet[i]["name"]} (prefix {len(pre_t)} flybys)')
        cols, (settles, ok_n, ssec) = _replay_route(job, W, fleet[i], pre_t, remaining, job['out'], 0, say, craft,
                                                   float(job.get('tank_prev', C.TANK_DRY)))
    else:
        t_lo, t_hi, _ = (float(x) for x in job['grid'])
        if job.get('unlaunched') is not None:              # the driver names the unlaunched indices in the job
            idx = [int(i) for i in job['unlaunched']]; hinted = True
        elif S is not None:
            idx = [int(c['i']) for c in S['craft'] if not c['launched']]; hinted = True
        else:
            idx = sorted(fleet); hinted = False
        rem = set(remaining)
        # launch inside the job's sub-window; the slice's first sub-window (t_lo = T) also takes a launch AT T (s16a's
        # r2 / r4 / r6 launch on day 0, which the strict (T, T + H] would never replay)
        lo_ok = (lambda tL: tL >= t_lo - 1e-6) if abs(t_lo - T) < 1.0 else (lambda tL: tL > t_lo + 1e-6)
        for i in idx:
            r = fleet.get(i)
            if r is None:
                continue
            if not (lo_ok(r['tL']) and r['tL'] <= t_hi + 1e-6) or r['tL'] > T + H:
                continue
            if not hinted and r['targets'][0] not in rem:
                continue
            say(f'  root: replaying {r["name"]} (launch {r["tL"] / DAY:.0f} d in ({t_lo / DAY:.0f}, {t_hi / DAY:.0f}])')
            cs, (s, o, sec) = _replay_route(job, W, r, [], remaining, job['out'], len(cols), say, None, C.TANK_DRY)
            cols += cs; settles += s; ok_n += o; ssec += sec
        if not cols:
            say(f'  root: no route launches in ({t_lo / DAY:.0f}, {t_hi / DAY:.0f}] d with a flyby before {(T + H) / DAY:.0f} d -> empty result')
    return cols, _stats(len(cols), settles, ok_n, ssec, depth=max([c['depth'] for c in cols], default=0))


# ------------------------------------------------------------------ the contract
def run_job(job):
    """The proposer entry point (driver contract): job['proposer'] in ('fixture_fake', 'fixture_replay') decides;
    writes out/result.json (idempotent) and returns the RESULT dict (s18_segbeam docstring)."""
    C.set_threads()
    rd = pathlib.Path(job['run_dir']); out = rd / job['out']; out.mkdir(parents=True, exist_ok=True)
    rf = out / 'result.json'
    if rf.exists():
        return C.jload(rf)
    say = C.logger(out / 'log.txt'); tic = time.time()
    prop = job.get('proposer', 'fixture_fake')
    say(f'fixture {prop} job {job["job_id"]} kind {job["kind"]} T {float(job["T"]) / DAY:.0f} d H {float(job["H"]) / DAY:.0f} d '
        f'remaining {len(job["remaining"])}')
    try:
        W = C.load_windows()
        if prop == 'fixture_fake':
            cols = fake_columns(job, W, _rng_of(job))
            st = _stats(len(cols), note='fake: no twins')
        elif prop == 'fixture_replay':
            fleet = job.get('fixture_fleet')
            if not fleet:
                S = _state_hint(job)
                fleet = (S or {}).get('params', {}).get('fixture_fleet') or str(DEFAULT_FLEET)
            RI.eph()
            cols, st = replay_columns(job, W, fleet, say)
        else:
            raise ValueError(f'unknown fixture proposer {prop!r}')
        for c in cols:
            C.check_column(c, float(job['T']), float(job['H']))
        st['wall_s'] = round(time.time() - tic, 2)
        res = _result(job, cols, st)
    except Exception:
        err = traceback.format_exc()[-2000:]
        say('ERROR ' + err.strip().splitlines()[-1])
        res = _result(job, [], _stats(0, wall_s=time.time() - tic, note='exception'), ok=False, error=err)
    say(f'done: ok {res["ok"]} columns {len(res["columns"])} wall {res["stats"]["wall_s"]:.1f} s')
    C.jdump(res, rf)
    return res


def make_job(job_id, kind, run_dir, T_d, H_d, L_d, remaining, craft=None, route=None, tank_prev=C.TANK_DRY, grid_d=None,
             proposer='fixture_fake', seed=0, fixture_fleet=None):
    """A JOB dict of the s18_segbeam schema for the tests (days in, seconds inside)."""
    out = f'slices/s{_slice_of(job_id):02d}/' + (f'c{craft}' if kind == 'craft' else f'r{_idx_of(job_id)}')
    job = dict(job_id=job_id, kind=kind, craft=craft, route=route, tank_prev=float(tank_prev),
               counts_prev=dict(pin=0, med=0, fil=0), remaining=[int(t) for t in remaining], T=C.d2s(T_d), H=C.d2s(H_d),
               L=C.d2s(L_d), grid=[C.d2s(x) for x in grid_d] if grid_d else None, roots=60, prize=None,
               params=dict(B=8, tries=16), seed=int(seed), run_dir=str(run_dir), out=out, proposer=proposer)
    if fixture_fleet:
        job['fixture_fleet'] = str(fixture_fleet)
    return job


# ------------------------------------------------------------------ selftest
def _table(cols, T, H):
    lines = [f'  {"id":16s} {"craft":>5s} {"n":>2s} {"targets":40s} {"epochs_d":28s} {"tank":>7s} {"dtank":>6s} {"pot":>3s} {"npz":6s}']
    for c in cols:
        lines.append(f'  {c["id"]:16s} {str(c["craft"]):>5s} {c["n_new"]:2d} {str(c["targets"])[:40]:40s} '
                     f'{",".join(f"{e / DAY:.0f}" for e in c["epochs"])[:28]:28s} {c["tank"]:7.1f} {c["dtank"]:6.1f} '
                     f'{c["potential"]:3d} {"yes" if c["npz"] else "none":6s}' + (f'  {c["note"]}' if c.get('note') else ''))
    return '\n'.join(lines)


def selftest(fleet=DEFAULT_FLEET):
    """Builds the three jobs of the module docstring in a scratch run dir under results/s18/_fixture_selftest, runs
    them, asserts the schema of every column (C.check_column), that fake columns contain the violators, that replay
    columns' twins reload with miss <= 150 km and epochs in (T, T+H], and prints the column tables + wall time."""
    tic = time.time()
    rd = C.S18 / '_fixture_selftest'
    if rd.exists():
        shutil.rmtree(rd)
    rd.mkdir(parents=True)
    say = C.logger(rd / 'log.txt')
    W = C.load_windows(); RI.eph()
    T, H, L = 0.0, 540.0, 360.0
    Ts, Hs = C.d2s(T), C.d2s(H)

    # 1. fake craft job (+ a second craft of the same slice: the conflicts)
    jA = make_job('s00_c0', 'craft', rd, T, H, L, C.TARGETS, craft=0, proposer='fixture_fake')
    jA2 = make_job('s00_c1', 'craft', rd, T, H, L, C.TARGETS, craft=1, proposer='fixture_fake')
    rA = run_job(jA); rA2 = run_job(jA2)
    assert rA['ok'] and rA2['ok'], (rA.get('error'), rA2.get('error'))
    for r in (rA, rA2):
        assert len(r['columns']) >= 3
        for c in r['columns']:
            C.check_column(c, Ts, Hs); assert c['npz'] is None and c['kind'] == 'craft' and c['craft'] == r['craft']
        assert any(c['counts']['fil'] == QUOTA_VIOLATOR_FIL for c in r['columns']), 'fake: quota violator missing'
        assert any(abs(c['tank'] - FUEL_VIOLATOR_TANK) < 1e-9 for c in r['columns']), 'fake: fuel violator missing'
    uA = set(t for c in rA['columns'] for t in c['targets']); uA2 = set(t for c in rA2['columns'] for t in c['targets'])
    assert uA & uA2, 'fake: craft 0 and craft 1 share no target (no conflicts)'
    rA_again = run_job(jA); assert rA_again == rA, 'fake: result.json not idempotent'
    rA_det = fake_columns(jA, W, _rng_of(jA))
    assert [c['targets'] for c in rA_det] == [c['targets'] for c in rA['columns']], 'fake: not deterministic'
    say(f'fake craft job s00_c0: {len(rA["columns"])} columns, conflicts with c1 on {sorted(uA & uA2)}')
    print(_table(rA['columns'], Ts, Hs))

    # 2. fake root job
    jB = make_job('s00_r0', 'root', rd, T, H, L, C.TARGETS, grid_d=(0.0, 540.0, 10.0), proposer='fixture_fake')
    rB = run_job(jB); assert rB['ok'], rB.get('error')
    assert len(rB['columns']) >= 3
    for c in rB['columns']:
        C.check_column(c, Ts, Hs); assert c['craft'] is None and c['kind'] == 'root' and c['npz'] is None
        assert 0 < c['t_launch'] <= Hs and c['epochs'][0] > c['t_launch']
    assert any(c['counts']['fil'] == QUOTA_VIOLATOR_FIL for c in rB['columns'])
    assert any(abs(c['tank'] - FUEL_VIOLATOR_TANK) < 1e-9 for c in rB['columns'])
    say(f'fake root job s00_r0: {len(rB["columns"])} columns')
    print(_table(rB['columns'], Ts, Hs))

    # 3. replay craft job for r1 at T = 0 (no prefix: the craft is "assigned" r1)
    jC = make_job('s00_c0', 'craft', rd / 'replay', T, H, L, C.TARGETS, craft=0, proposer='fixture_replay', fixture_fleet=fleet)
    (rd / 'replay').mkdir(exist_ok=True)
    t0 = time.time(); rC = run_job(jC); wC = time.time() - t0
    assert rC['ok'], rC.get('error')
    assert len(rC['columns']) == 2, f'expected 2 replay columns (r1 at 540 d and one flyby earlier), got {len(rC["columns"])}'
    for c in rC['columns']:
        C.check_column(c, Ts, Hs); assert c['npz'] is not None and c['craft'] == 0
        st = C.load_twin(rd / 'replay' / c['npz']); miss = C.miss_of(C.twin_of(st))
        assert miss <= C.TOL_MISS, f'replay twin {c["npz"]} miss {miss}'
        assert all(Ts < e <= Hs for e in c['epochs']); assert c['all_targets'] == c['targets']
        assert abs(c['miss'] - miss) < 5.0
    assert rC['columns'][0]['n_new'] == rC['columns'][1]['n_new'] + 1
    t0 = time.time(); rC2 = run_job(jC); assert rC2 == rC and time.time() - t0 < 1.0, 'replay: not idempotent'
    say(f'replay craft job s00_c0 (r1, T 0, H 540): {len(rC["columns"])} columns in {wC:.1f} s '
        f'(main: {rC["columns"][0]["n_new"]} new, tank {rC["columns"][0]["tank"]:.1f}, miss {rC["columns"][0]["miss"]:.0f} km)')
    print(_table(rC['columns'], Ts, Hs))

    # 3b. replay craft job at T = 540 with the committed prefix (the driver's next slice)
    main = rC['columns'][0]
    shutil.copy(rd / 'replay' / main['npz'], rd / 'replay' / 'craft_0.npz')
    remaining = [t for t in C.TARGETS if t not in set(main['all_targets'])]
    jD = make_job('s01_c0', 'craft', rd / 'replay', 540.0, H, L, remaining, craft=0, route='craft_0.npz',
                  tank_prev=main['tank'], proposer='fixture_replay', fixture_fleet=fleet)
    t0 = time.time(); rD = run_job(jD); wD = time.time() - t0
    assert rD['ok'] and len(rD['columns']) >= 1, rD.get('error')
    for c in rD['columns']:
        C.check_column(c, C.d2s(540.0), Hs)
        st = C.load_twin(rd / 'replay' / c['npz']); miss = C.miss_of(C.twin_of(st))
        assert miss <= C.TOL_MISS and c['all_targets'][:len(main['all_targets'])] == main['all_targets']
        assert all(C.d2s(540.0) < e <= C.d2s(1080.0) for e in c['epochs'])
        assert not (set(c['targets']) & set(main['all_targets']))
    say(f'replay craft job s01_c0 (r1 prefix 4 flybys, T 540): {len(rD["columns"])} columns in {wD:.1f} s '
        f'(main: {rD["columns"][0]["n_new"]} new, tank {rD["columns"][0]["tank"]:.1f})')
    print(_table(rD['columns'], C.d2s(540.0), Hs))

    # 3c. replay root jobs: slice 0 over (0, 540] (r1..r8 launch on days 0-88) and slice 1 over (540, 1080] (r9: day 944)
    jE = make_job('s00_r0', 'root', rd / 'replay', T, H, L, C.TARGETS, grid_d=(0.0, 540.0, 10.0), proposer='fixture_replay', fixture_fleet=fleet)
    t0 = time.time(); rE = run_job(jE); wE = time.time() - t0
    assert rE['ok'] and len(rE['columns']) >= 2, rE.get('error')
    for c in rE['columns']:
        C.check_column(c, Ts, Hs); assert c['craft'] is None and 0 <= c['t_launch'] <= Hs
        st = C.load_twin(rd / 'replay' / c['npz']); assert C.miss_of(C.twin_of(st)) <= C.TOL_MISS
    srcs = sorted(set(c['source'] for c in rE['columns']))
    say(f'replay root job s00_r0 (0, 540]: {len(rE["columns"])} columns from {srcs} in {wE:.1f} s')
    jF = make_job('s01_r0', 'root', rd / 'replay', 540.0, H, L, C.TARGETS, grid_d=(540.0, 1080.0, 10.0), proposer='fixture_replay', fixture_fleet=fleet)
    t0 = time.time(); rF = run_job(jF); wF = time.time() - t0
    assert rF['ok'], rF.get('error')
    for c in rF['columns']:
        C.check_column(c, C.d2s(540.0), Hs); assert c['craft'] is None
    say(f'replay root job s01_r0 (540, 1080]: {len(rF["columns"])} columns '
        f'({", ".join(c["source"] + ":" + str(c["targets"]) for c in rF["columns"]) or "empty result: r9 has no flyby before 1080 d"}) in {wF:.1f} s')
    say(f'selftest wall {time.time() - tic:.1f} s (MEASURED); scratch dir {rd}')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('cmd', choices=['fake', 'replay', 'selftest']); ap.add_argument('job', nargs='?')
    ap.add_argument('--fleet', default=str(DEFAULT_FLEET))
    a = ap.parse_args()
    if a.cmd == 'selftest':
        selftest(a.fleet); print('selftest OK')
    else:
        job = C.jload(a.job); job['proposer'] = 'fixture_' + a.cmd
        if a.cmd == 'replay' and a.fleet:
            job.setdefault('fixture_fleet', a.fleet)
        r = run_job(job)
        print(json.dumps(dict(job_id=r['job_id'], ok=r['ok'], n_columns=len(r['columns']), error=r.get('error')), default=float))


if __name__ == '__main__':
    main()
