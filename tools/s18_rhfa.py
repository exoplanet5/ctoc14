"""Stage 18 / WP2: the RHFA DRIVER -- rolling-horizon fleet auction (docs/stage18_rhfa.md section 1).

Loop (T_0 = 0; H, L in DAYS on the CLI, seconds inside):
  1. PROPOSE   one JOB per launched craft (its prefix twin, remaining) + n_root root JOBs (n_root = max(1,
               nproc - n_launched) if any craft is unlaunched, each over a contiguous part of the launch window
               (T, T + H]); all jobs run on ONE fork Pool(nproc, maxtasksperchild=1) with <proposer>.run_job
               (imap_unordered); each worker is single-threaded.  RI.eph() is loaded in the parent before the pool is
               created (copy-on-write).  The union of the RESULTs' columns is written to slices/s<kk>/columns.json.
               Job results are cached on disk (result.json), so a killed slice re-runs only its unfinished jobs.
  2. AUCTION   s18_auction.choose -> slices/s<kk>/auction.json.
  3. COMMIT    for every chosen column, in the same pool (w_commit): load its twin, C.regrid_settle (THE REGRID LAW),
               accept iff ok and honest tank <= m0max + m0slack; a craft whose column is REJECTED is re-auctioned over
               the slice's remaining columns (COMMIT_RETRIES rounds; integrator fix 09-24, g0b slice 1); on accept: save craft_<i>.npz (regridded), update the
               craft (launched, targets, epochs, counts, tank, t_launch, t_end, history), remaining, taken; on reject:
               log, the craft keeps its previous state, the targets stay in remaining.  Write order: craft_<i>.npz,
               then state.json (T <- T + H, slice += 1, phase 'propose'), then slices/s<kk>/commit.json (the marker
               LAST, so a crash anywhere inside the commit re-runs it idempotently from the cached auction.json).
  4. STOP      when nothing usable is left (C.slice_H: the LAST slice absorbs the mission tail, H_eff = T_END - T when
               that is <= H + 120 d; fixer 09-24, finding 2), when remaining is empty, or when --stop-T (DAYS) is given
               and T >= it (gate G0b: --stop-T 1461).  Then finalize: fleet/route_<i>.npz for every launched craft (the
               regridded twins), fleet/fleet.json (RI.IFleet.save keys + honest flags + covered + misses + sumJi + fuel).

Fixer 09-24 (reviewer findings 1-6, all inside tools/s18_*.py):
  (1) PoolBox: jobs are dispatched with apply_async and polled under a deadline (wall_job + job_grace for the beam,
      regrid_timeout + 120 s for the commit); a lost / raised job is logged, counted as ok False, and the pool is
      re-created (a crashed worker no longer hangs the slice for ever).
  (2) S['H_eff_d'] per slice (start_slice; stored in jobs/index.json and the log record).
  (3) fuel pacing: m0max 1070 (= the 3764 kg bar / 8), C.ramp_cap at T + H_eff in the auction filter AND at the commit
      (+ m0slack 10 kg regrid tolerance), one fleet fuel row per auction (C.slice_fuel_budget), lam as the dual of that
      budget (S['lam_dual'], update_lam: x1.5 when the row binds or the ramp drops a column, /1.25 when < half is used).
  (4) quotas 12 / 18, lifted in the last quota_free_slices (3) slices (C.quota_active); the beam is quota-aware: classes
      at the quota are withheld from a craft job's remaining, the head-room travels as job['quota'].
  (5) a job's cached FAILURE is re-run once by s18_segbeam.run_job (retried flag) and logged here.
  (6) the beam stops after look_stop (2) consecutive levels without a column candidate (in_H = 0); wall_job 2400 is a
      safety only; settle tries are capped at settle_cap (8) after a level with < settle_ok_min (0.3) success.

Resume: `run OUT ...` with an existing OUT/state.json continues at (slice, phase); the phase files decide what is
skipped: columns.json -> propose skipped, auction.json -> auction skipped; jobs/<id>.json are frozen once written and a
job's result.json is returned by the proposer without recomputing.  Every slice starts by snapshotting the state and
the launched prefixes into slices/s<kk>/in/ (state_in.json, craft_<i>.npz; the jobs point at the snapshot), and a
resume REWINDS to the earliest slice k <= state.slice whose commit.json is missing (state := its state_in.json): the
phase files are the truth, state.json is derived.  Deleting slices/s02/{auction,commit}.json and re-running therefore
re-auctions and re-commits slice 2 from its cached columns.json and reproduces the same state.
All CLI parameters are stored in state['params'] on creation and NOT overridden on resume (a changed CLI value is
logged and ignored), except --nproc, --stop-T and --tlim which are runtime knobs (so a run stopped by --stop-T is
continued by running it again without / with a later --stop-T).

Log line per slice (OUT/log.txt): `[slice s T=1080 d] launched 6/8  covered 71  remaining 227  fuel 812 kg  sumJi 7.1
chosen: c0:5(+31) c1:4(+20) r->c6:3(+9) ...  rejected: c3 (miss 412 km)  (wall 1830 s)`.

Proposers (--proposer): segbeam (default; s18_segbeam.run_job) | planbeam (s18_planbeam.run_job) | fixture_fake | fixture_replay (s18_fixture.run_job;
WP3).  While s18_segbeam.make_job / s18_fixture.run_job are still skeletons (NotImplementedError) the driver falls back
to make_job_local and to the stub proposers below (stub_fake_run_job / stub_replay_run_job: same JOB / RESULT contract,
labelled 'stub' in the log) so that WP2 is testable before WP1 / WP3 land.  --commit regrid (default) | none
(fixture_fake only: adopt the column's tank unchanged, a placeholder twin is written, every tank flagged honest=False;
NEVER for a real run).
SIGTERM / Ctrl-C: the pool is terminated, nothing extra is saved (every phase already saved its file); exit code 130.

usage: s18_rhfa.py run OUT [--N 8 --H 540 --L 360 --B 8 --tries 16 --beta 1.5 --gamma 0.3 --lam 1.0 --qf 12 --qm 18
                            --m0max 1070 --m0slack 10 --nproc 8 --seed 0 --stop-T DAYS --proposer segbeam --commit regrid
                            --fuel-bar 3764 --ramp-alpha 0.9 --ramp-slack 60 --fleet-factor 1.3 --lam-max 6
                            --quota-free-slices 3 --look-stop 2 --settle-cap 8 --settle-ok-min 0.3 --job-grace 600
                            --k-epochs 3 --lin-cap 3.0 --res-max 1.5e8 --iters 100 --timeout 240 --grid-step 10 --roots 60
                            --max-levels 14 --wall-job 2400 --root-bin 60 --beam-edf 1 --tlim 60 --fixture-fleet DIR
                            --max-launch K --launch-stagger D --prize-json PATH --claims --claim-bonus 0.5 --claim-ttl 2]
       --max-launch K (09-24, default 0 = no cap): at most K root columns (launches) are chosen per slice (the auction's
                      root row: sum root x_c <= min(K, n_unlaunched); the commit's retry rounds count the launches the
                      slice already accepted).  --launch-stagger D (days, default 0 = off): the k-th launch of a slice
                      by epoch (it goes to the k-th unlaunched craft) must lie at >= T + k D (auction rows, see
                      s18_auction).  Both are frozen in state['params'] like every other run parameter.
       --prize-json PATH (09-24, default off): a JSON {target_id: bonus}; the bonus is added to the auction prize of
                      those targets in EVERY slice (s18_auction.value_of: 1 + edf + bonus) and to the proposer jobs'
                      prize dict (job['prize'] = 1 + edf + bonus, so the beam / planner steers toward them too).  The
                      table is frozen into state['params'] as prize_bonus (the path as prize_json); a resume of a run
                      created without it logs and ignores the flag.
       --claims (09-24, default off): plan continuations become CLAIMS.  After each slice's commit every accepted
                      column's `plan_next` (planbeam: the plan's flybys beyond T + H, in plan order; segbeam / fixtures
                      emit none) is stored as the claims of that craft (state['claims'] = {"<craft>": {"<target>": slice
                      claimed}}).  In the next slice's jobs (make_jobs) the other craft (and the root jobs) get
                      `remaining` MINUS the targets claimed by others, the claimant keeps them and gets them prized
                      +claim_bonus (--claim-bonus, 0.5).  Claims are rewritten from the chosen columns every slice
                      (update_claims: a claim not in the claimant's new plan_next, or whose craft got no column, is
                      released; a claim older than --claim-ttl (2) slices is released even when still planned; a target
                      planned by two craft goes to the lower index).  Logged per slice (claimed targets per craft),
                      persisted in state.json (and in every slice's state_in.json snapshot, so resume / rewind see the
                      same claims).  Frozen in state['params'] like every other run parameter.
       s18_rhfa.py status OUT            (one line from state.json; exit 0)
       s18_rhfa.py finalize OUT          (re-write fleet/ from state.json; idempotent)
       s18_rhfa.py misses-prize RUN_DIR OUT.json --bonus 1.0 [--prev PREV.json --decay 0.5]
                      writes {miss: bonus} for the misses of a FINISHED run (RUN_DIR/fleet/fleet.json), merged with a
                      previous prize file (prev values x decay, then + bonus for the current misses; entries that decay
                      below 1e-6 are dropped).  The miss-feedback iteration: run k+1 with --prize-json = the misses of
                      run k:  misses-prize results/s18/<run_k> results/s18/prizes_k.json --bonus 1.0 --prev
                      results/s18/prizes_k-1.json --decay 0.5;  run results/s18/<run_k+1> ... --prize-json
                      results/s18/prizes_k.json.
Long runs: nohup nice -n 5 ~/.venvs/astro313/bin/python tools/s18_rhfa.py run results/s18/<run> ... > results/s18/<run>/nohup.out 2>&1 &
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, math, shutil, signal, pathlib, argparse, traceback, multiprocessing as mp
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import numpy as np
import run_ialns as RI
import s18_common as C
import s18_auction as AU
from ctoc14.constants import DAY, T_MISSION

RUNTIME_KNOBS = ('nproc', 'stop_T', 'tlim', 'job_grace')
COMMIT_RETRIES = 3          # re-auction rounds for craft whose column failed the regrid law (commit)
POLL_S = 2.0                # pool polling period (PoolBox.map)
LAM_UP, LAM_DOWN = 1.5, 1.25   # dual step of the fuel price (update_lam)


class PoolBox:
    """The fork Pool with a LOST-job guard (fixer 09-24, reviewer finding 1 -- FATAL: a worker killed by the OOM killer
    or a segfault never delivers its task and imap_unordered waits forever; MEASURED by results/s18/review/pool_crash.py).
    map(fn, jobs, timeout_each, say, label, on_result) dispatches every job with apply_async and polls the AsyncResults;
    a job that has not returned by its deadline (timeout_each x waves after dispatch, waves = ceil(n / nproc)) or whose
    worker raised is logged 'LOST' / 'RAISED' and returned as None; after a lost job the pool is terminated and re-created
    so that no stale worker survives.  Results are handed to on_result in completion order."""

    def __init__(self, nproc):
        self.nproc = int(nproc); self.pool = None; self.create()

    def create(self):
        self.pool = mp.get_context('fork').Pool(self.nproc, initializer=_worker_init, maxtasksperchild=1)

    def recreate(self, say):
        _shutdown_pool(self.pool, say, graceful=False, timeout=10.0)
        self.create(); say('pool re-created')

    def map(self, fn, jobs, timeout_each, say, label='job', on_result=None, key=None):
        """-> list of results aligned with jobs (None for a lost / raised job)."""
        jobs = list(jobs); n = len(jobs)
        if n == 0:
            return []
        key = key or (lambda j, i: str(i))
        waves = int(math.ceil(n / max(1, self.nproc)))
        deadline = time.time() + float(timeout_each) * waves
        ars = [self.pool.apply_async(fn, (j,)) for j in jobs]
        out = [None] * n; pending = set(range(n)); lost = []
        while pending:
            for i in sorted(pending):
                ar = ars[i]
                if ar.ready():
                    pending.discard(i)
                    try:
                        out[i] = ar.get(timeout=0)
                    except Exception as e:
                        say(f'  {label} {key(jobs[i], i)} RAISED in the worker: {repr(e)[:200]}'); lost.append(i); continue
                    if on_result is not None:
                        on_result(out[i])
            if not pending:
                break
            if time.time() > deadline:
                for i in sorted(pending):
                    say(f'  {label} {key(jobs[i], i)} LOST: no result {time.time() - (deadline - timeout_each * waves):.0f} s '
                        f'after dispatch (deadline {timeout_each:.0f} s x {waves} wave(s)) -> treated as failed')
                lost += sorted(pending); pending = set()
                break
            time.sleep(POLL_S)
        if lost:
            say(f'{len(lost)} {label}(s) lost -> terminating and re-creating the pool')
            self.recreate(say)
        return out
DEFAULT_FLEET = str(ROOT / 'results/s16/best/fleet')
# copy of s18_segbeam.JOB_DEF (used only while that module is a skeleton / not importable)
JOB_DEF_LOCAL = dict(B=8, tries=16, k_epochs=3, lin_cap=3.0, res_max=1.5e8, iters=100, timeout=240.0, max_levels=14,
                     wall_budget=2400.0, root_bin=60.0, w_t=None, w_fuel=None, dv=1.2, drmax=0.15, mc=60)
_MAIN_PID = os.getpid()
_FALLBACK_NOTED = set()


class Interrupted(Exception):
    pass


# ------------------------------------------------------------------ CLI / params
def parse_args(argv=None):
    """argparse over C.DEF (every key a --flag; underscores -> dashes; H, L, grid_step, root_bin, stop_T in DAYS) plus
    cmd, out, --tlim (auction time limit, s), --fixture-fleet (fixture_replay source, default results/s16/best/fleet).
    Returns the Namespace; params_of(a) turns it into the state['params'] dict."""
    ap = argparse.ArgumentParser(description='stage 18 RHFA driver')
    ap.add_argument('cmd', choices=['run', 'status', 'finalize', 'misses-prize']); ap.add_argument('out')
    ap.add_argument('out2', nargs='?', default=None, help='misses-prize: the OUT.json to write')
    for k, v in C.DEF.items():
        flag = '--' + k.replace('_', '-')
        if v is None:
            ap.add_argument(flag, type=float, default=None, dest=k)
        elif isinstance(v, bool):
            ap.add_argument(flag, type=int, default=int(v), dest=k)
        elif isinstance(v, int):
            ap.add_argument(flag, type=int, default=v, dest=k)
        elif isinstance(v, float):
            ap.add_argument(flag, type=float, default=v, dest=k)
        else:
            ap.add_argument(flag, type=str, default=v, dest=k)
    ap.add_argument('--tlim', type=float, default=60.0)
    ap.add_argument('--fixture-fleet', default=DEFAULT_FLEET, dest='fixture_fleet')
    # launch pacing (09-24; additive, both off by default): see the module docstring and s18_auction.launch_cap
    ap.add_argument('--max-launch', type=int, default=0, dest='max_launch',
                    help='at most K launches (root columns) per slice; 0 = no cap (default)')
    ap.add_argument('--launch-stagger', type=float, default=0.0, dest='launch_stagger',
                    help='days: the k-th launch of a slice (by epoch) must be >= T + k D; 0 = off (default)')
    # prize feedback (09-24; additive, off by default): a JSON {target: bonus} added to the auction / proposer prizes
    ap.add_argument('--prize-json', default=None, dest='prize_json',
                    help='JSON {target_id: bonus} added to the auction prize and the proposer prize of those targets in every slice')
    # claims (09-24; additive, off by default): plan continuations reserve targets for the claimant (module docstring)
    ap.add_argument('--claims', action='store_true', dest='claims',
                    help='store every chosen column\'s plan_next as the craft\'s claims for the next slices (default off)')
    ap.add_argument('--claim-bonus', type=float, default=0.5, dest='claim_bonus',
                    help='prize bonus of a claimed target in the claimant\'s own proposer job (default 0.5)')
    ap.add_argument('--claim-ttl', type=int, default=2, dest='claim_ttl',
                    help='a claim older than this many slices is released (default 2)')
    # misses-prize helper arguments
    ap.add_argument('--bonus', type=float, default=1.0, help='misses-prize: bonus per miss of RUN_DIR')
    ap.add_argument('--prev', default=None, help='misses-prize: previous prize JSON to merge (values x --decay)')
    ap.add_argument('--decay', type=float, default=0.5, help='misses-prize: factor applied to the --prev values')
    return ap.parse_args(argv)


LAUNCH_KEYS = ('max_launch', 'launch_stagger')     # params added 09-24 (absent in older state.json files = off)
ADDED_KEYS = LAUNCH_KEYS + ('prize_json', 'claims')   # on/off params added after the first runs (absent in old states = off)


def params_of(a):
    """state['params'] dict from the Namespace: C.DEF keys (values as given, days stay days) + tlim + fixture_fleet
    + max_launch (int, 0 = off) + launch_stagger (days, 0 = off) + prize_json (path | None; the table itself is
    loaded into prm['prize_bonus'] by run() when the run is created) + claims (int 0/1) + claim_bonus + claim_ttl."""
    prm = {k: getattr(a, k) for k in C.DEF}
    prm['tlim'] = float(a.tlim); prm['fixture_fleet'] = str(a.fixture_fleet)
    prm['max_launch'] = int(getattr(a, 'max_launch', 0) or 0)
    prm['launch_stagger'] = float(getattr(a, 'launch_stagger', 0.0) or 0.0)
    prm['prize_json'] = str(a.prize_json) if getattr(a, 'prize_json', None) else None
    prm['claims'] = int(bool(getattr(a, 'claims', False)))
    prm['claim_bonus'] = float(getattr(a, 'claim_bonus', 0.5))
    prm['claim_ttl'] = int(getattr(a, 'claim_ttl', 2))
    return prm


def claims_on(prm):
    """True when the run was created with --claims (prm['claims']; absent in older states = off)."""
    return bool(int(prm.get('claims', 0) or 0)) if prm is not None else False


def claim_owner(S):
    """{target (int): craft (int)} of the claims in force (state['claims'] = {"<craft>": {"<target>": slice})."""
    out = {}
    for ci, d in (S.get('claims') or {}).items():
        for t in d:
            out[int(t)] = int(ci)
    return out


def update_claims(S, prm, accepted, k, say):
    """Rewrite state['claims'] after the commit of slice k from the accepted columns' plan_next (module docstring):
    per accepted craft, in index order, every plan_next target that is still remaining becomes / stays its claim
    (`since` = the slice it was first claimed, kept across rewrites); a target already claimed by a lower craft in this
    rewrite is a conflict (skipped); a claim with k - since >= claim_ttl is released (expired); every old claim not
    re-made is released.  Logs the claimed targets per craft.  Returns {craft (str): n_claims}."""
    old = S.get('claims') or {}; ttl = int(prm.get('claim_ttl', 2)); rem = set(int(t) for t in S['remaining'])
    new = {}; owner = {}; n_new = n_keep = n_exp = n_conf = 0
    for a in sorted(accepted, key=lambda a: int(a['craft'])):
        i = str(int(a['craft'])); oc = old.get(i) or {}; d = {}
        for t in (a.get('plan_next') or []):
            t = int(t)
            if t not in rem or str(t) in d:
                continue
            if t in owner:
                n_conf += 1; continue
            since = int(oc.get(str(t), k))
            if k - since >= ttl:
                n_exp += 1; continue
            d[str(t)] = since; owner[t] = i
            if since == k:
                n_new += 1
            else:
                n_keep += 1
        if d:
            new[i] = d
    n_old = sum(len(v) for v in old.values()); n_rel = max(0, n_old - n_keep - n_exp)
    S['claims'] = new
    cnt = {i: len(d) for i, d in new.items()}
    say(f'[slice {k}] claims: ' + (' '.join(f'c{i}:{n}' for i, n in sorted(cnt.items(), key=lambda kv: int(kv[0]))) or '-')
        + f'  ({sum(cnt.values())} targets: {n_keep} kept, {n_new} new; {n_rel} released, {n_exp} expired (ttl {ttl}), {n_conf} conflicts)')
    return cnt


def load_prize_json(path):
    """{target id (str): bonus (float)} from a --prize-json file; unreachable / unknown ids are dropped (logged by the
    caller through the returned dict size); None / '' -> {}."""
    if not path:
        return {}
    raw = C.jload(path)
    ok = set(C.TARGETS); out = {}
    for k, v in raw.items():
        try:
            t = int(k); b = float(v)
        except (TypeError, ValueError):
            continue
        if t in ok and b != 0.0:
            out[str(t)] = b
    return out


def misses_prize(run_dir, out, bonus=1.0, prev=None, decay=0.5):
    """Write OUT.json = {miss: bonus} for the misses of a finished run (RUN_DIR/fleet/fleet.json), merged with a
    previous prize file: prev values x decay, then + bonus for every current miss; entries below 1e-6 are dropped.
    Returns the table."""
    rd = C.run_dir(run_dir); F = C.jload(rd / 'fleet' / 'fleet.json')
    misses = [int(t) for t in F.get('misses', [])]
    tab = {}
    if prev:
        for k, v in C.jload(prev).items():
            val = float(v) * float(decay)
            if abs(val) >= 1e-6:
                tab[str(int(k))] = val
    for t in misses:
        tab[str(t)] = tab.get(str(t), 0.0) + float(bonus)
    tab = {k: round(v, 6) for k, v in sorted(tab.items(), key=lambda kv: int(kv[0])) if abs(v) >= 1e-6}
    C.jdump(tab, out)
    print(f'misses-prize: {rd.name}: {len(misses)} misses (covered {F.get("covered")}), bonus {bonus}'
          + (f', prev {prev} x {decay}' if prev else '') + f' -> {len(tab)} entries, sum {sum(tab.values()):.3f} -> {out}')
    return tab


def init_run(out, prm):
    """Create OUT (mkdir), state.json (C.new_state), copy windows.json into OUT (frozen table for the run), log header.
    Returns S.  Refuses to overwrite an existing state.json (that is a resume)."""
    rd = pathlib.Path(out); rd.mkdir(parents=True, exist_ok=True)
    if (rd / 'state.json').exists():
        raise RuntimeError(f'{rd}/state.json exists: use run to resume')
    if pathlib.Path(C.WINDOWS_FILE).exists():
        shutil.copy(C.WINDOWS_FILE, rd / 'windows.json')
    else:
        C.build_windows(out=rd / 'windows.json')
    S = C.new_state(rd.name, prm, int(prm['N']))
    C.save_state(S, rd)
    return S


# ------------------------------------------------------------------ jobs
def prizes_for(S, W, prm):
    """{target: 1 + edf(W, t, T + H, beta) [+ prize_bonus]} over remaining when prm['beam_edf'] else None (the beam's
    own prizes; the auction always applies edf regardless).  With a --prize-json table (prm['prize_bonus']) the bonus
    is added per target (and the dict is built from 1.0 even when beam_edf is off, so the bonus still reaches the
    proposer)."""
    bonus = AU.prize_bonus_of(prm)
    if not int(prm.get('beam_edf', 1)) and not bonus:
        return None
    TH = float(S['T']) + H_eff(S, prm); beta = float(prm['beta']); use_edf = bool(int(prm.get('beam_edf', 1)))
    pz = {int(t): 1.0 + (C.edf(W, int(t), TH, beta) if use_edf else 0.0) for t in S['remaining']}
    for t, b in bonus.items():
        if int(t) in pz:
            pz[int(t)] += float(b)
    return pz


def H_eff(S, prm):
    """The current slice's effective horizon [s] (finding 2): S['H_eff_d'] when set by start_slice, else C.slice_H."""
    if S.get('H_eff_d') is not None:
        return C.d2s(S['H_eff_d'])
    return C.slice_H(float(S['T']), C.d2s(prm['H']))


def start_slice(S, rd, prm, say):
    """Fix the slice's effective horizon (S['H_eff_d']) and the fuel dual (S['lam_dual'], initialised to prm['lam'])
    before anything of the slice is snapshotted; saved to state.json so a rewind sees the same values."""
    changed = False
    Hd = round(C.s2d(C.slice_H(float(S['T']), C.d2s(prm['H']))), 4)
    if S.get('H_eff_d') != Hd:
        S['H_eff_d'] = Hd; changed = True
        if abs(Hd - float(prm['H'])) > 1e-6:
            say(f'[slice {S["slice"]}] effective horizon {Hd:.1f} d (the last slice absorbs the mission tail; H {prm["H"]})')
    if S.get('lam_dual') is None:
        S['lam_dual'] = float(prm['lam']); changed = True
    if changed:
        C.save_state(S, rd)
    return S


def quota_headroom(c, prm):
    """{'fil': n, 'med': n} the craft may still take under the quotas (root / unlaunched: the full quotas)."""
    cnt = c.get('counts') or dict(pin=0, med=0, fil=0)
    return dict(fil=max(0, int(prm['qf']) - int(cnt.get('fil', 0))), med=max(0, int(prm['qm']) - int(cnt.get('med', 0))))


def update_lam(S, prm, choice, say):
    """The fuel price as the dual of the slice budget (finding 3): after the auction, raise lam_dual by LAM_UP when the
    fleet row is binding (chosen fuel >= 0.95 x budget) or a column was dropped by the ramp; lower it by LAM_DOWN
    (never below prm['lam']) when less than half the budget was used and nothing was dropped by the ramp.  Applies to
    the NEXT slice (the current auction is already decided and cached)."""
    lam0 = float(S.get('lam_dual', prm['lam'])); lam_max = float(prm.get('lam_max', C.DEF['lam_max']))
    budget = float(choice.get('fuel_budget', 0.0) or 0.0); used = float(choice.get('fuel_chosen', 0.0) or 0.0)
    n_ramp = int(choice.get('n_ramp', 0) or 0); n_kept = int(choice.get('n_kept', 0) or 0)
    f_ramp = n_ramp / max(1, n_kept + n_ramp)          # share of the admissible columns the ramp removed
    if budget > 0 and (used >= 0.95 * budget or f_ramp >= 0.2):
        lam = min(lam_max, lam0 * LAM_UP)
    elif budget > 0 and used < 0.5 * budget and f_ramp < 0.05:
        lam = max(float(prm['lam']), lam0 / LAM_DOWN)
    else:
        lam = lam0
    if abs(lam - lam0) > 1e-9:
        say(f'[slice {S["slice"]}] fuel dual: lam {lam0:.3f} -> {lam:.3f} (chosen {used:.0f} / budget {budget:.0f} kg, ramp drops {n_ramp} of {n_kept + n_ramp})')
    S['lam_dual'] = float(lam)
    return S


def _job_def():
    try:
        import s18_segbeam as SB
        return dict(SB.JOB_DEF)
    except Exception:
        return dict(JOB_DEF_LOCAL)


def make_job_local(kind, run_dir, slice_k, T, H, L, remaining, params, craft=None, route=None, tank_prev=C.TANK_DRY,
                   counts_prev=None, grid=None, roots=60, prize=None, seed=0, idx=0, quota=None):
    """Local copy of s18_segbeam.make_job (the JOB schema of its docstring); used while that one is a skeleton."""
    JD = _job_def()
    job_id = f's{int(slice_k):02d}_c{int(craft)}' if kind == 'craft' else f's{int(slice_k):02d}_r{int(idx)}'
    out = f'slices/s{int(slice_k):02d}/' + (f'c{int(craft)}' if kind == 'craft' else f'r{int(idx)}')
    return dict(job_id=job_id, kind=kind, craft=(int(craft) if craft is not None else None), route=route,
                tank_prev=float(tank_prev), counts_prev=dict(counts_prev or dict(pin=0, med=0, fil=0)),
                remaining=[int(t) for t in remaining], T=float(T), H=float(H), L=float(L),
                grid=([float(g) for g in grid] if grid is not None else None), roots=int(roots),
                prize=({str(int(k)): float(v) for k, v in prize.items()} if prize else None),
                params={k: params.get(k, JD[k]) for k in JD}, seed=int(seed), run_dir=str(run_dir), out=out,
                quota=(dict(quota) if quota else None))


def _make_job(say, **kw):
    """s18_segbeam.make_job when implemented, else make_job_local (logged once)."""
    try:
        import s18_segbeam as SB
        return SB.make_job(**kw)
    except (ImportError, NotImplementedError, AttributeError, SyntaxError) as e:
        if 'make_job' not in _FALLBACK_NOTED:
            _FALLBACK_NOTED.add('make_job'); say(f'note: s18_segbeam.make_job unavailable ({type(e).__name__}) -> make_job_local')
        return make_job_local(**kw)


def make_jobs(S, rd, prm, W, say):
    """The slice's JOB list (s18_segbeam.make_job): one craft job per launched craft; if unlaunched craft exist, n_root =
    max(1, prm['nproc'] - n_launched) root jobs over contiguous sub-grids of (T, T + H] with step grid_step (DAYS ->
    s; slice 0's first sub-grid also contains the epoch 0, a launch at mission start); root `roots` per job =
    ceil(prm['roots'] / n_root) but >= 10.  Writes each job to slices/s<kk>/jobs/<job_id>.json
    (+ jobs/index.json last); an existing index freezes the job set (resume).  Snapshots state_in.json and the launched
    prefixes into slices/s<kk>/in/ first.  prize = prizes_for(...)."""
    k = int(S['slice']); sd = C.slice_dir(rd, k); jd = sd / 'jobs'; ind = jd / 'index.json'
    if ind.exists():
        return [C.jload(jd / f'{jid}.json') for jid in C.jload(ind)['jobs']]
    # snapshot (the rewind target and the frozen prefixes of this slice)
    sin = sd / 'in'; sin.mkdir(parents=True, exist_ok=True)
    C.jdump(S, sin / 'state_in.json')
    T = float(S['T']); H = H_eff(S, prm); L = C.d2s(prm['L']); rem = [int(t) for t in S['remaining']]
    prize = prizes_for(S, W, prm)
    params = dict(prm); params['wall_budget'] = float(prm.get('wall_job', JOB_DEF_LOCAL['wall_budget']))
    extra = dict(proposer=prm['proposer'], fixture_fleet=prm.get('fixture_fleet', DEFAULT_FLEET), slice=k, N=int(S['N']),
                 unlaunched=[int(c['i']) for c in S['craft'] if not c['launched']], H_eff_d=C.s2d(H))
    q_on = C.quota_active(T, H, prm)
    # claims (09-24, --claims): targets claimed by another craft are withheld from a job's remaining; the claimant's
    # own claims are prized +claim_bonus in its job
    cl_on = claims_on(prm); owner = claim_owner(S) if cl_on else {}; cbonus = float(prm.get('claim_bonus', 0.5))
    if cl_on:
        cnt = {i: len(d) for i, d in (S.get('claims') or {}).items()}
        say(f'[slice {k}] claims in force: ' + (' '.join(f'c{i}:{n}' for i, n in sorted(cnt.items(), key=lambda kv: int(kv[0]))) or '-')
            + f'  ({len(owner)} targets; claim_bonus {cbonus}, ttl {prm.get("claim_ttl", 2)})')
    jobs = []
    for c in S['craft']:
        if not c['launched']:
            continue
        i = int(c['i']); snap = sin / f'craft_{i}.npz'
        shutil.copy(C.craft_npz(rd, i), snap)
        # quota-aware beam (finding 4): classes already at the quota are removed from the job's remaining, the head-room
        # of the others is passed as job['quota']; no quota in the last quota_free_slices slices
        quota = quota_headroom(c, prm) if q_on else None
        rem_i = rem; prize_i = prize; own = []
        if owner:
            others = {t for t, ci in owner.items() if ci != i}
            own = sorted(t for t, ci in owner.items() if ci == i and t in rem)
            rem_i = [t for t in rem if t not in others]
            if own:
                prize_i = dict(prize) if prize else {int(t): 1.0 for t in rem}
                for t in own:
                    prize_i[int(t)] = float(prize_i.get(int(t), 1.0)) + cbonus
            say(f'  c{i}: claims: {len(own)} own (prized +{cbonus}), {len(rem) - len(rem_i)} claimed by others withheld')
        if quota:
            full = {cl for cl, n in quota.items() if n <= 0}
            if full:
                n0 = len(rem_i); rem_i = [t for t in rem_i if C.class_of(W, t) not in full]
                say(f'  c{i}: classes at quota {sorted(full)} -> {n0 - len(rem_i)} targets withheld from its beam')
            quota = {cl: n for cl, n in quota.items() if n > 0}
        j = _make_job(say, kind='craft', run_dir=str(rd), slice_k=k, T=T, H=H, L=L, remaining=rem_i, params=params, craft=i,
                      route=str(snap.relative_to(rd)), tank_prev=float(c['tank']), counts_prev=dict(c['counts']),
                      prize=prize_i, seed=int(prm['seed']), quota=quota)
        j.update(extra)
        if cl_on:
            j['claims'] = [int(t) for t in own]
        jobs.append(j)
    n_launched = len(jobs); n_unl = len(extra['unlaunched'])
    if owner:
        rem = [t for t in rem if t not in owner]                 # root jobs: nothing that a launched craft has claimed
    if n_unl > 0:
        step = C.d2s(prm['grid_step']); M = int(math.floor(H / step + 1e-9))
        n_root = max(1, min(max(1, int(prm['nproc']) - n_launched), max(1, M)))
        roots = max(10, int(math.ceil(float(prm['roots']) / n_root)))
        bounds = [int(round(M * j / n_root)) for j in range(n_root + 1)]
        idx = 0
        for j in range(n_root):
            i0, i1 = bounds[j], bounds[j + 1]
            if i1 <= i0:
                continue
            # launch epochs T + step (i0 + 1) .. T + step i1; the very first grid of the run also holds the epoch T = 0
            # (a launch at mission start is legitimate: s16a's r2 launches at day 0)
            t_lo = T + step * i0 - (step if (k == 0 and i0 == 0 and T <= 0.0) else 0.0)
            grid = [t_lo, T + step * i1, step]
            jb = _make_job(say, kind='root', run_dir=str(rd), slice_k=k, T=T, H=H, L=L, remaining=rem, params=params,
                           grid=grid, roots=roots, prize=prize, seed=int(prm['seed']), idx=idx,
                           quota=(dict(fil=int(prm['qf']), med=int(prm['qm'])) if q_on else None))
            jb.update(extra); jb['idx'] = idx; jobs.append(jb); idx += 1
    jd.mkdir(parents=True, exist_ok=True)
    for j in jobs:
        C.jdump(j, jd / f"{j['job_id']}.json")
    index = dict(jobs=[j['job_id'] for j in jobs], T_d=C.s2d(T), H_eff_d=C.s2d(H), n_launched=n_launched, n_unlaunched=n_unl,
                 quota_active=bool(q_on), lam_dual=S.get('lam_dual'))
    if cl_on:
        index['claims'] = {i: sorted(int(t) for t in d) for i, d in (S.get('claims') or {}).items()}
    C.jdump(index, ind)
    return jobs


# ------------------------------------------------------------------ proposers (+ stubs while WP1 / WP3 are skeletons)
def _result_cached(job):
    rf = pathlib.Path(job['run_dir']) / job['out'] / 'result.json'
    return C.jload(rf) if rf.exists() else None


def _result(job, columns, ok=True, error=None, tic=None, note='stub'):
    r = dict(job_id=job['job_id'], kind=job['kind'], craft=job.get('craft'), ok=bool(ok), error=error, columns=columns,
             stats=dict(levels=1, n_states=len(columns), n_columns=len(columns), settles=0, settle_ok=0, price_s=0.0,
                        settle_s=0.0, wall_s=round(time.time() - tic, 2) if tic else 0.0, end=dict(kind=note, depth=1),
                        hist=[]))
    out = pathlib.Path(job['run_dir']) / job['out']; out.mkdir(parents=True, exist_ok=True)
    C.jdump(r, out / 'result.json')
    return r


def _windows_of(job):
    p = pathlib.Path(job['run_dir']) / 'windows.json'
    return C.load_windows(p if p.exists() else C.WINDOWS_FILE)


def stub_fake_run_job(job):
    """STUB fake proposer (used only while s18_fixture.run_job is a skeleton; fake test only): columns without twins,
    deterministic in (seed, slice, craft | root idx); ~half of every column's targets come from a shared 'hot' pool so
    that craft conflict; one quota violator (11 fillers) and one fuel violator (tank 1300) per job; npz = None."""
    r = _result_cached(job)
    if r is not None:
        return r
    tic = time.time()
    try:
        W = _windows_of(job); rd = pathlib.Path(job['run_dir'])
        k = int(job.get('slice', 0)); craft = job.get('craft'); idx = int(job.get('idx', 0))
        rng = np.random.default_rng([int(job.get('seed', 0)), k, (int(craft) if craft is not None else 50 + idx)])
        rem = [int(t) for t in job['remaining']]; T = float(job['T']); H = float(job['H'])
        hot, cold = rem[:24], rem[24:]
        prefix, tL = [], T
        if job.get('route') and (rd / job['route']).exists():
            st = C.load_twin(rd / job['route']); prefix = [int(a) for a in st['asts'][C.order_of(st)]]; tL = float(st['tL'])
        tank_prev = float(job.get('tank_prev', C.TANK_DRY)); cols = []

        def mk(targets, tank):
            n = len(cols)
            if craft is None:
                lo, hi, _ = job['grid']; t_launch = float(rng.uniform(lo, hi)); base = t_launch
            else:
                t_launch, base = tL, T
            ep = base + (T + H - base) * np.sort(rng.uniform(0.05, 1.0, len(targets)))
            for q in range(1, len(ep)):
                ep[q] = max(ep[q], ep[q - 1] + 1.0)
            cols.append(C.new_column(id=f"{job['job_id']}_k{n}", craft=craft, targets=list(targets), epochs=ep.tolist(),
                                     all_targets=prefix + list(targets), npz=None, tank=float(tank), tank_prev=tank_prev,
                                     miss=float(rng.uniform(0, 50)), t_launch=t_launch, t_end=float(ep[-1]),
                                     potential=int(rng.poisson(2)), counts=C.count_classes(W, targets),
                                     depth=len(prefix) + len(targets), score=0.0, job_id=job['job_id']))
        for _ in range(4):
            m = int(rng.integers(1, 7)); nh = min(len(hot), (m + 1) // 2); nc = min(len(cold), m - nh)
            tg = list(rng.choice(hot, nh, replace=False)) + list(rng.choice(cold, nc, replace=False)) if (nh + nc) else []
            if tg:
                mk([int(t) for t in rng.permutation(tg)], tank_prev + 12.0 * len(tg) + rng.normal(0, 2))
        fil = [t for t in rem if C.class_of(W, t) == 'fil'][:11]
        if len(fil) == 11:
            mk(fil, tank_prev + 12.0 * 11)                              # quota violator (fil 11 > qf 10)
        if len(cold) >= 2:
            mk([int(t) for t in rng.choice(cold, 2, replace=False)], 1300.0)   # fuel violator
        return _result(job, cols, tic=tic, note='stub_fake')
    except Exception:
        return _result(job, [], ok=False, error=traceback.format_exc()[-800:], tic=tic, note='stub_fake')


def stub_replay_run_job(job):
    """STUB replay proposer (used only while s18_fixture.run_job is a skeleton): real twins = C.truncate_twin of the
    fleet routes (--fixture-fleet, r1..rN) at the last new flyby in (T, T+H] and at the flyby before; a craft job replays
    the route that contains its prefix; a root job replays every route not yet started (first target still remaining)
    whose launch is inside the job's grid; potential = the route's flybys in (T+H, T+H+L]."""
    r = _result_cached(job)
    if r is not None:
        return r
    tic = time.time()
    try:
        C.set_threads(); RI.eph(); W = _windows_of(job); rd = pathlib.Path(job['run_dir'])
        fleet = pathlib.Path(job.get('fixture_fleet') or DEFAULT_FLEET); N = int(job.get('N', 8))
        T = float(job['T']); H = float(job['H']); L = float(job['L']); rem = set(int(t) for t in job['remaining'])
        craft = job.get('craft'); routes = sorted(fleet.glob('route_r*.npz'), key=lambda p: int(p.stem[7:]))[:N]
        src = []
        for p in routes:
            st = C.load_twin(p); o = C.order_of(st); tg = [int(st['asts'][i]) for i in o]
            if craft is not None:
                pre = C.load_twin(rd / job['route'])
                if set(int(a) for a in pre['asts']) <= set(tg):
                    src.append(st); break
            else:
                lo, hi, _ = job['grid']
                if tg[0] in rem and lo < float(st['tL']) <= hi:
                    src.append(st)
        cols = []
        for st in src:
            o = C.order_of(st); ep = [float(st['tf'][i]) for i in o]; tg = [int(st['asts'][i]) for i in o]
            cuts = []
            for e, t in zip(ep, tg):
                if e <= T:
                    continue
                if e > T + H or t not in rem:
                    break
                cuts.append(e)
            if not cuts:
                continue
            potential = sum(1 for e in ep if T + H < e <= T + H + L)
            for t_cut in ([cuts[-1]] + ([cuts[-2]] if len(cuts) >= 2 else [])):
                res, n_kept = C.truncate_twin(st, t_cut)
                if not res['ok']:
                    continue
                s2 = res['st']; o2 = C.order_of(s2); ep2 = [float(s2['tf'][i]) for i in o2]; tg2 = [int(s2['asts'][i]) for i in o2]
                new = [(t, e) for t, e in zip(tg2, ep2) if e > T]
                n = len(cols); npz = f"{job['out']}/col_{n}.npz"
                c = C.new_column(id=f"{job['job_id']}_k{n}", craft=craft, targets=[t for t, _ in new], epochs=[e for _, e in new],
                                 all_targets=tg2, npz=npz, tank=float(res['tank']), tank_prev=float(job.get('tank_prev', C.TANK_DRY)),
                                 miss=float(res['miss']), t_launch=float(s2['tL']), t_end=ep2[-1], potential=potential,
                                 counts=C.count_classes(W, [t for t, _ in new]), depth=n_kept, score=0.0, job_id=job['job_id'])
                C.check_column(c, T, H); C.save_twin(rd / npz, s2); cols.append(c)
        return _result(job, cols, tic=tic, note='stub_replay')
    except Exception:
        return _result(job, [], ok=False, error=traceback.format_exc()[-800:], tic=tic, note='stub_replay')


def fixture_run_job(job):
    """s18_fixture.run_job (WP3) with the stub fallback while it is a skeleton (NotImplementedError, raised or caught
    into the RESULT's error)."""
    stub = stub_fake_run_job if job['proposer'] == 'fixture_fake' else stub_replay_run_job
    try:
        import s18_fixture as FX
        r = FX.run_job(job)
        if r.get('ok') or 'NotImplementedError' not in str(r.get('error', '')):
            return r
    except (NotImplementedError, ImportError, AttributeError, SyntaxError):
        pass
    return stub(job)


def resolve_proposer(name):
    """The callable mapped over the JOBs (importing its module in the parent BEFORE the fork Pool is created)."""
    if name == 'segbeam':
        import s18_segbeam as SB
        return SB.run_job
    if name == 'planbeam':                                  # plan long, commit short (tools/s18_planbeam.py; same contract)
        import s18_planbeam as PB
        return PB.run_job
    if name in ('fixture_fake', 'fixture_replay'):
        try:
            import s18_fixture  # noqa: F401  (imported here so the forked workers inherit it)
        except Exception:
            pass
        return fixture_run_job
    raise ValueError(f'unknown proposer {name}')


# ------------------------------------------------------------------ phases
def propose(S, rd, prm, pool, say):
    """Phase 1.  If slices/s<kk>/columns.json exists return it.  Else make_jobs, dispatch through `pool` with the
    proposer (segbeam: s18_segbeam.run_job; fixture_*: s18_fixture.run_job), collect RESULTs (ok False -> logged, no
    columns), concatenate columns, C.check_column each (bad ones logged and dropped), write columns.json, return list.
    Logs per job: id, kind, n_columns, levels, settles/ok, wall."""
    k = int(S['slice']); sd = C.slice_dir(rd, k); cf = sd / 'columns.json'
    if cf.exists():
        cols = C.jload(cf); say(f'[slice {k}] propose: cached columns.json ({len(cols)} columns)'); return cols
    W = C.load_windows(rd / 'windows.json')
    jobs = make_jobs(S, rd, prm, W, say)
    T = float(S['T']); H = H_eff(S, prm)
    say(f'[slice {k} T={C.s2d(T):.0f} d] propose: {len(jobs)} jobs ({sum(1 for j in jobs if j["kind"] == "craft")} craft, '
        f'{sum(1 for j in jobs if j["kind"] == "root")} root) on {prm["nproc"]} workers, proposer {prm["proposer"]}, '
        f'H_eff {C.s2d(H):.0f} d, lam {S.get("lam_dual", prm["lam"])}')
    fn = resolve_proposer(prm['proposer']); columns = []; n_bad = 0; tic = time.time()

    def on_result(res):
        nonlocal n_bad
        st = res.get('stats') or {}
        if not res.get('ok'):
            say(f'  job {res.get("job_id")} {"cached FAILURE (already retried once, not retried again)" if res.get("cached") else "FAILED"}: '
                f'{str(res.get("error"))[-200:]}'); return
        good = []
        for c in res.get('columns', []):
            try:
                C.check_column(c, T, H); good.append(c)
            except AssertionError as e:
                n_bad += 1; say(f'  job {res["job_id"]}: column {c.get("id")} dropped ({e})')
        columns.extend(good)
        say(f'  job {res["job_id"]} {res["kind"]:5s} columns {len(good):3d}  levels {st.get("levels", "-")}  settles '
            f'{st.get("settle_ok", "-")}/{st.get("settles", "-")}  end {(st.get("end") or {}).get("kind", "-")}  wall {st.get("wall_s", 0):.1f} s'
            + ('  (cached)' if res.get('cached') else '') + ('  (retried)' if res.get('retried') else ''))
    if jobs:
        # LOST guard (finding 1): a job may take wall_job (beam) + one level + the settle timeout; job_grace covers that
        timeout = float(prm.get('wall_job', JOB_DEF_LOCAL['wall_budget'])) + float(prm.get('job_grace', C.DEF['job_grace']))
        pool.map(fn, jobs, timeout, say, label='job', on_result=on_result, key=lambda j, i: j['job_id'])
    say(f'[slice {k}] propose: {len(columns)} columns ({n_bad} bad dropped) in {time.time() - tic:.1f} s')
    C.jdump(columns, cf)
    return columns


def auction(S, rd, prm, columns, W, say):
    """Phase 2.  If slices/s<kk>/auction.json exists return it.  Else s18_auction.choose(columns, S, W, prm, tlim) ->
    write and return the CHOICE."""
    k = int(S['slice']); af = C.slice_dir(rd, k) / 'auction.json'
    if af.exists():
        ch = C.jload(af); say(f'[slice {k}] auction: cached ({len(ch["chosen"])} chosen, status {ch["status"]})'); return ch
    ch = AU.choose(columns, S, W, prm, tlim=float(prm.get('tlim', 60.0)), say=lambda s: say(f'[slice {k}] ' + s))
    C.jdump(ch, af)
    return ch


def _lost_commit(i, col):
    return dict(craft=int(i), column_id=col['id'], ok=False, tank=float(col['tank']), miss=float(col.get('miss', 0.0)),
                secs=0.0, error='LOST (worker died or exceeded the commit deadline)', st=None, honest=dict(ok=False))


def _placeholder_twin(col, prev_epochs, all_targets):
    """A syntactically valid ist dict for --commit none (NO trajectory; flagged honest=False everywhere)."""
    ep = list(prev_epochs) + [float(e) for e in col['epochs']]
    return dict(tL=float(col['t_launch'] if col.get('t_launch') is not None else (ep[0] - 30 * DAY if ep else 0.0)),
                vinf=np.zeros(3), ts=np.zeros(0), Ts=np.zeros((0, 3)), tf=np.asarray(ep, float), asts=np.asarray(all_targets, int))


def w_commit(job):
    """Pool worker: job = (craft index, column dict, run_dir, regrid_iters, regrid_timeout, mode) ->
    dict(craft, column_id, ok, tank, miss, secs, error, st (ist dict, regridded) | None, honest (dict)).
    mode 'regrid': C.regrid_settle on C.twin_of(C.load_twin(run_dir / column['npz'])); mode 'none' (fixture only):
    ok True, tank = column['tank'], st = the column's twin if any, honest = {'ok': False, 'note': 'commit none'}."""
    i, col, rd, iters, timeout, mode = job
    tic = time.time(); rd = pathlib.Path(rd)
    out = dict(craft=int(i), column_id=col['id'], ok=False, tank=float(col['tank']), miss=float(col.get('miss', 0.0)),
               secs=0.0, error=None, st=None, honest=dict(ok=False))
    try:
        C.set_threads()
        if mode == 'none':
            st = C.load_twin(rd / col['npz']) if col.get('npz') else None
            out.update(ok=True, st=st, honest=dict(ok=False, note='commit none'))
        else:
            if not col.get('npz'):
                raise ValueError('column has no twin (npz None) and commit mode is regrid')
            ip = C.twin_of(C.load_twin(rd / col['npz']))
            res = C.regrid_settle(ip, iters=int(iters), timeout=float(timeout))
            out.update(ok=bool(res['ok']), tank=float(res['tank']), miss=float(res['miss']), st=res['st'], honest=res['honest'],
                       error=res['error'], tank_in=res['tank_in'], over_in=res['over_in'], over_out=res['over_out'])
    except Exception:
        out['error'] = traceback.format_exc()[-600:]
    out['secs'] = round(time.time() - tic, 2)
    return out


def _fleet_numbers(S):
    tanks = [c['tank'] for c in S['craft'] if c['launched']]
    return dict(launched=len(tanks), covered=len(S['taken']), fuel=C.fuel_of(tanks), sumJi=C.sum_Ji(tanks))


def commit(S, rd, prm, choice, W, pool, say, tic_slice=None):
    """Phase 3 (module docstring).  Accept iff res['ok'] and res['tank'] <= m0max + m0slack (mode none: always).
    Root columns are assigned to unlaunched craft in index order (the CHOICE already carries the craft).  Writes
    craft_<i>.npz (C.save_twin) BEFORE state.json, and slices/s<kk>/commit.json after it (the marker; a crash before it
    re-runs the commit idempotently).  Then advances T, slice, phase and saves the state; appends the slice's log
    record.  Returns S."""
    k = int(S['slice']); sd = C.slice_dir(rd, k); cf = sd / 'commit.json'; tic = time.time()
    if cf.exists():
        say(f'[slice {k}] commit: commit.json exists but the state was not advanced -> re-running (idempotent)')
    mode = prm['commit']; H = H_eff(S, prm)
    # acceptance cap (finding 3): the ramp at T + H_eff (never above m0max) + m0slack, the honest-vs-planner tolerance
    cap = C.ramp_cap(float(S['T']) + H, float(prm['m0max']), prm=prm) + float(prm['m0slack'])
    items = sorted(((int(i), col) for i, col in choice['columns'].items()), key=lambda x: x[0])
    accepted, rejected, chosen_txt, results, rounds = [], [], [], {}, []
    # RETRY (integrator fix 09-24; results/s18/g0b slice 1: c0's column failed the regrid law, miss 1e5 km, and the craft
    # idled for a whole slice): a craft whose column is rejected is re-auctioned over the slice's remaining columns
    # (its own ids minus the rejected ones, root columns while it is unlaunched; targets accepted meanwhile are stale
    # for the auction because S['remaining'] is already updated), up to COMMIT_RETRIES rounds.
    all_cols = C.jload(sd / 'columns.json') if (sd / 'columns.json').exists() else []
    tried_ids = set()
    for rnd in range(1 + COMMIT_RETRIES):
        if rnd > 0:
            rej_craft = {int(r['craft']) for r in rounds[-1]}
            if not rej_craft or not all_cols:
                break
            unl = {int(c['i']) for c in S['craft'] if not c['launched']}
            sub = [c for c in all_cols if c['id'] not in tried_ids and
                   ((c.get('craft') is not None and int(c['craft']) in rej_craft) or (c.get('craft') is None and rej_craft & unl))]
            if not sub:
                say(f'[slice {k}] commit retry {rnd}: no columns left for craft {sorted(rej_craft)}'); break
            fuel_used = float(sum(max(0.0, a['dtank']) for a in accepted))
            launches_used = sum(1 for a in accepted if a['root'])      # --max-launch counts the slice's accepted launches
            ch2 = AU.choose(sub, S, W, prm, tlim=float(prm.get('tlim', 60.0)), say=lambda s: say(f'[slice {k}] retry {rnd} ' + s),
                            fuel_used=fuel_used, launches_used=launches_used)
            items = sorted(((int(i), col) for i, col in ch2['columns'].items()), key=lambda x: x[0])
            if not items:
                break
            choice['retries'] = choice.get('retries', []) + [dict(round=rnd, chosen=dict(ch2['chosen']), status=ch2['status'], n_sub=len(sub))]
        tried_ids |= {col['id'] for _, col in items}
        acc_r, rej_r = _commit_round(S, rd, prm, items, W, pool, say, k, mode, cap, results)
        accepted += acc_r; rejected += rej_r; rounds.append(rej_r)
        chosen_txt += [f'{"r->" if a["root"] else ""}c{a["craft"]}:{a["n_new"]}({a["dtank"]:+.0f})' + (f'[retry {rnd}]' if rnd else '')
                       for a in acc_r]
        if not rej_r:
            break
    rejected = [r for r in rejected if r['craft'] not in {a['craft'] for a in accepted}] + \
               [dict(r, note='recovered by a retry') for r in rejected if r['craft'] in {a['craft'] for a in accepted}]
    T_slice = float(S['T'])
    fn = _fleet_numbers(S); wall = round(time.time() - (tic_slice or tic), 1)
    rec = dict(slice=k, T_d=C.s2d(T_slice), H_eff_d=C.s2d(H), n_jobs=len(list((sd / 'jobs').glob('s*.json'))), n_columns=choice['n_columns'],
               n_kept=choice['n_kept'], chosen={str(a['craft']): a['column_id'] for a in accepted}, rejected=rejected, launched=fn['launched'],
               covered=fn['covered'], remaining=len(S['remaining']), fuel=round(fn['fuel'], 1), sumJi=round(fn['sumJi'], 4),
               wall_s=wall, auction_status=choice['status'], retries=choice.get('retries', []), lam=choice.get('lam'),
               fuel_budget=choice.get('fuel_budget'), fuel_slice=round(float(sum(a['dtank'] for a in accepted)), 1), cap=round(cap, 1))
    if claims_on(prm):
        rec['claims'] = update_claims(S, prm, accepted, k, say)       # --claims: rewrite the claims from the chosen columns
    S['log'].append(rec)
    update_lam(S, prm, choice, say)
    S['T'] = T_slice + H; S['slice'] = k + 1; S['phase'] = 'propose'; S['H_eff_d'] = None
    C.save_state(S, rd)
    C.jdump(dict(slice=k, T_d=C.s2d(T_slice), mode=mode, accepted=accepted, rejected=rejected, results={
        str(i): {kk: v for kk, v in r.items() if kk != 'st'} for i, r in results.items()}, record=rec), cf)
    say(f'[slice {k} T={C.s2d(T_slice):.0f} d] launched {fn["launched"]}/{S["N"]}  covered {fn["covered"]}  remaining '
        f'{len(S["remaining"])}  fuel {fn["fuel"]:.0f} kg  sumJi {fn["sumJi"]:.4f}  chosen: {" ".join(chosen_txt) or "-"}  '
        f'rejected: {" ".join(f"c{r['craft']} ({r['reason']})" for r in rejected if not r.get("note")) or "-"}'
        f'{"  recovered by retry: " + " ".join(f"c{r['craft']}" for r in rejected if r.get("note")) if any(r.get("note") for r in rejected) else ""}'
        f'  (wall {wall:.0f} s)')
    return S


def _commit_round(S, rd, prm, items, W, pool, say, k, mode, cap, results):
    """One commit round over items = [(craft, column)]: w_commit in the pool, then accept / reject each (module
    docstring rules); accepted craft states, remaining and taken are updated in S.  Returns (accepted, rejected)."""
    jobs = [(i, col, str(rd), prm['regrid_iters'], prm['regrid_timeout'], mode) for i, col in items]
    res_r = {}
    if jobs:
        # LOST guard (finding 1): regrid_settle is bounded by regrid_timeout (SIGALRM) + the regrid itself
        outs = pool.map(w_commit, jobs, float(prm['regrid_timeout']) + 120.0, say, label='commit', key=lambda j, i: f'c{j[0]}')
        for (i, col), res in zip(items, outs):
            res_r[i] = res if res is not None else _lost_commit(i, col)
    results.update({f'{i}:{r["column_id"]}': r for i, r in res_r.items()})
    accepted, rejected = [], []
    for i, col in items:
        res = res_r[i]; c = S['craft'][i]; prev_t = list(c['targets']); prev_e = list(c['epochs']); new = [int(t) for t in col['targets']]
        reason = None
        if not res['ok']:
            reason = f"regrid_settle failed ({res.get('error') or 'miss ' + format(res['miss'], '.0f') + ' km / cap ' + format(res['honest'].get('max_window_cap', -1), '.3f')})"
        elif mode != 'none' and res['tank'] > cap:
            reason = f"honest tank {res['tank']:.1f} > {cap:.0f}"
        elif any(t not in set(S['remaining']) for t in new):
            reason = 'stale target'
        if reason is None:
            if mode == 'none':
                st = res['st'] if res['st'] is not None else _placeholder_twin(col, prev_e, prev_t + new)
                tg, ep = prev_t + new, prev_e + [float(e) for e in col['epochs']]
            else:
                st = res['st']; o = C.order_of(st); tg = [int(st['asts'][j]) for j in o]; ep = [float(st['tf'][j]) for j in o]
                if set(tg) != set(prev_t) | set(new) or len(tg) != len(prev_t) + len(new):
                    reason = f'twin targets != prefix + new ({len(tg)} vs {len(prev_t)}+{len(new)})'
            if reason is None and not all(a < b for a, b in zip(ep, ep[1:])):
                reason = 'epochs not increasing'
        if reason is not None:
            rejected.append(dict(craft=i, column_id=col['id'], reason=reason, secs=res['secs']))
            say(f'[slice {k}] commit: c{i} column {col["id"]} REJECTED: {reason}  (targets stay in remaining)'); continue
        tank_prev = float(c['tank']); tank = float(res['tank']) if mode != 'none' else float(col['tank'])
        C.save_twin(C.craft_npz(rd, i), st)
        was_root = not c['launched']
        c.update(launched=True, route=f'craft_{i}.npz', targets=tg, epochs=ep, counts=C.count_classes(W, tg), tank=tank,
                 t_launch=float(st['tL']), t_end=float(max(ep)), honest=bool(mode != 'none' and res['honest'].get('ok', False)))
        c['history'].append(dict(slice=k, column_id=col['id'], n_new=len(new), tank=tank, dtank=round(tank - tank_prev, 3),
                                 secs=res['secs'], miss=round(float(res['miss']), 1)))
        rem = set(S['remaining'])
        for t in new:
            rem.discard(t); S['taken'][str(t)] = i
        S['remaining'] = sorted(rem)
        acc = dict(craft=i, column_id=col['id'], n_new=len(new), tank=tank, dtank=round(tank - tank_prev, 3),
                   miss=round(float(res['miss']), 1), secs=res['secs'], honest=c['honest'], root=was_root)
        if claims_on(prm):
            acc['plan_next'] = [int(t) for t in (col.get('plan_next') or [])]     # --claims: the column's continuation
        accepted.append(acc)
    return accepted, rejected


def should_stop(S, prm):
    """(stop: bool, reason dict) per the STOP rule of the module docstring."""
    T = float(S['T']); Td = C.s2d(T)
    if not S['remaining']:
        return True, dict(kind='covered', T_d=Td)
    # finding 2: the last slice absorbs the tail (C.slice_H), so the loop ends only when nothing usable is left
    if C.slice_H(T, C.d2s(prm['H'])) < C.H_TAIL_MIN - 1e-6:
        return True, dict(kind='horizon', T_d=Td)
    st = prm.get('stop_T')
    if st is not None and T >= C.d2s(st) - 1e-6:
        return True, dict(kind='stop_T', T_d=Td, stop_T_d=float(st))
    return False, None


def finalize(S, rd, prm, say):
    """fleet/route_<i>.npz (copy of craft_<i>.npz for launched craft) + fleet/fleet.json:
    {J (raw: sumJi + misses + 2), n, covered, misses [ids], sumJi, fuel, honest (all routes pass C.honest_check),
     routes: {"<i>": {tank, J_i, flybys, t_launch_d, t_end_d, counts, cadence_d, honest, check}}, params, stop, note}.
    Idempotent.  Never writes outside OUT."""
    rd = pathlib.Path(rd); fd = rd / 'fleet'; fd.mkdir(parents=True, exist_ok=True)
    mode = prm.get('commit', 'regrid'); routes = {}; tanks = []; all_ok = True
    for c in S['craft']:
        if not c['launched']:
            continue
        i = int(c['i']); src = C.craft_npz(rd, i); dst = fd / f'route_{i}.npz'
        shutil.copy(src, dst)
        st = C.load_twin(src); n = len(c['targets'])
        if mode == 'none':
            chk = dict(ok=False, note='commit none: placeholder twin, tank NOT honest')
        else:
            chk = C.honest_check(C.twin_of(st))
        ok = bool(chk.get('ok', False)); all_ok &= ok; tanks.append(float(c['tank']))
        span_d = C.s2d(float(c['t_end']) - float(c['t_launch'])) if c['t_end'] is not None else 0.0
        routes[str(i)] = dict(tank=float(c['tank']), J_i=C.cost(c['tank']), flybys=n, t_launch_d=C.s2d(c['t_launch']),
                              t_end_d=C.s2d(c['t_end']), counts=dict(c['counts']), cadence_d=round(span_d / n, 1) if n else None,
                              honest=ok, check=chk, targets=list(c['targets']))
    covered = len(S['taken']); misses = sorted(int(t) for t in S['remaining'])
    F = dict(J=C.fleet_J(tanks, covered), n=len(routes), covered=covered, misses=misses, sumJi=C.sum_Ji(tanks),
             fuel=C.fuel_of(tanks), honest=bool(all_ok and mode != 'none'), routes=routes, params=dict(S['params']),
             stop=S.get('stop'), note=f'stage 18 RHFA run {S["run"]}: commit mode {mode}; tanks after regrid + settle'
             if mode != 'none' else f'stage 18 RHFA run {S["run"]}: commit mode none (FIXTURE): tanks NOT honest')
    C.jdump(F, fd / 'fleet.json')
    say(f'finalize: fleet/ {len(routes)} routes  covered {covered}  misses {len(misses)}  fuel {F["fuel"]:.0f} kg  '
        f'sumJi {F["sumJi"]:.4f}  J {F["J"]:.4f}  honest {F["honest"]}')
    return F


# ------------------------------------------------------------------ the loop
def resume_point(S, rd, say):
    """The phase files are the truth: rewind to the earliest slice k < state.slice whose commit.json is missing (its
    in/state_in.json); returns the state to continue from (runtime knobs are re-applied by the caller)."""
    for k in range(int(S['slice'])):
        sd = C.slice_dir(rd, k)
        if not (sd / 'commit.json').exists():
            sin = sd / 'in' / 'state_in.json'
            if sin.exists():
                S2 = C.jload(sin); C.check_state(S2); S2['params'] = S['params']
                say(f'resume: slice {k} has no commit.json -> rewind to its state_in (T={C.s2d(S2["T"]):.0f} d)')
                return S2
            say(f'resume: slice {k} has no commit.json and no state_in.json -> cannot rewind, continuing from state.json')
    return S


def _install_sigterm():
    def h(sig, frm):
        if os.getpid() != _MAIN_PID:                      # a forked worker: default action
            signal.signal(signal.SIGTERM, signal.SIG_DFL); os.kill(os.getpid(), signal.SIGTERM); return
        raise Interrupted('SIGTERM')
    signal.signal(signal.SIGTERM, h)


def _worker_init():
    """Pool initializer, run first in EVERY fresh worker.  With maxtasksperchild=1 the replacement workers are forked
    by the pool's handler thread while the parent's main thread may be inside a print (holding the stdout buffer
    lock): a child that inherits that lock blocks forever on its first log line and can then not even be terminated
    (measured 09-24 01:10: worker stuck in PyThread_acquire_lock_timed, parent stuck in Pool.join).  So: fresh
    stdout / stderr objects on duplicated fds (the inherited objects are kept referenced, never flushed by us; a fork
    worker exits through os._exit), default SIGTERM, BLAS threads = 1."""
    C.set_threads()
    try:
        signal.signal(signal.SIGTERM, signal.SIG_DFL)
    except Exception:
        pass
    for name, fd in (('stdout', 1), ('stderr', 2)):
        try:
            setattr(sys, '_inherited_' + name, getattr(sys, name))
            setattr(sys, name, os.fdopen(os.dup(fd), 'w', buffering=1))
        except Exception:
            pass


def _shutdown_pool(pool, say, graceful=True, timeout=30.0):
    """close()/terminate() + join() under a watchdog: a worker stuck in a C lock ignores SIGTERM and would hang the
    join forever, so after `timeout` s the remaining workers get SIGKILL.  Returns True when the pool shut down
    cleanly (the caller may then exit normally; otherwise it should os._exit)."""
    import threading

    def go():
        try:
            if graceful:
                pool.close()
            else:
                pool.terminate()
            pool.join()
        except Exception:
            pass
    th = threading.Thread(target=go, daemon=True); th.start(); th.join(timeout)
    if not th.is_alive():
        return True
    left = [p for p in list(getattr(pool, '_pool', [])) if p.is_alive()]
    say(f'pool shutdown did not finish in {timeout:.0f} s -> SIGKILL to {len(left)} worker(s)')
    for p in left:
        try:
            p.kill()
        except Exception:
            pass
    return False


def run(a):
    """The loop.  Creates or resumes the run; opens the fork Pool once (after RI.eph()); phases with on-disk caching;
    handles KeyboardInterrupt / SIGTERM by saving nothing extra (every phase already saved its file).  Sets
    S['done'], S['stop'], finalize at the end."""
    rd = C.run_dir(a.out); rd.mkdir(parents=True, exist_ok=True); say = C.logger(rd / 'log.txt')
    cli = params_of(a)
    if (rd / 'state.json').exists():
        S = C.load_state(rd); prm = S['params']
        for k, v in cli.items():
            if k in RUNTIME_KNOBS:
                if prm.get(k) != v:
                    say(f'resume: runtime knob {k}: {prm.get(k)} -> {v}')
                prm[k] = v
            elif k in prm and prm[k] != v and k != 'fixture_fleet':
                say(f'resume: CLI --{k.replace("_", "-")} {v} ignored (frozen {prm[k]})')
            elif k not in prm and k in ADDED_KEYS and v:
                say(f'resume: CLI --{k.replace("_", "-")} {v} ignored (the run was created without this option: off)')
        S['params'] = prm
        S = resume_point(S, rd, say); S['params'] = prm
        say(f'resume {rd.name}: slice {S["slice"]} phase {S["phase"]} T={C.s2d(S["T"]):.0f} d done={S["done"]}')
    else:
        prm = cli
        # --prize-json (09-24): the table is frozen into the params at creation (the file is read once, here)
        prm['prize_bonus'] = load_prize_json(prm.get('prize_json'))
        S = init_run(rd, prm)
        if claims_on(prm):
            S['claims'] = {}; C.save_state(S, rd)
        say(f'new run {rd.name}: N {prm["N"]} H {prm["H"]} d L {prm["L"]} d proposer {prm["proposer"]} commit {prm["commit"]} '
            f'nproc {prm["nproc"]} stop_T {prm["stop_T"]}'
            + (f' max_launch {prm["max_launch"]}' if prm.get('max_launch') else '')
            + (f' launch_stagger {prm["launch_stagger"]} d' if prm.get('launch_stagger') else '')
            + (f' prize_json {prm["prize_json"]} ({len(prm["prize_bonus"])} targets, bonus sum {sum(prm["prize_bonus"].values()):.2f})'
               if prm.get('prize_json') else '')
            + (f' claims on (bonus {prm["claim_bonus"]}, ttl {prm["claim_ttl"]})' if claims_on(prm) else ''))
    if prm['commit'] == 'none' and prm['proposer'] != 'fixture_fake':
        raise SystemExit('--commit none is allowed only with --proposer fixture_fake')
    W = C.load_windows(rd / 'windows.json')
    RI.eph()                                              # loaded once in the parent, inherited by the fork
    resolve_proposer(prm['proposer'])                     # import the proposer module before the fork
    pool = PoolBox(int(prm['nproc']))
    _install_sigterm(); code = 0; clean = True
    try:
        while True:
            stop, reason = should_stop(S, prm)
            if stop:
                S['done'] = True; S['stop'] = reason; C.save_state(S, rd)
                say(f'stop: {reason["kind"]} at T={reason["T_d"]:.0f} d (slice {S["slice"]})'); break
            if S.get('done'):
                S['done'] = False; S['stop'] = None; C.save_state(S, rd); say('resume: stop rule no longer holds -> continuing')
            tic = time.time()
            S = start_slice(S, rd, prm, say)
            columns = propose(S, rd, prm, pool, say)
            choice = auction(S, rd, prm, columns, W, say)
            S = commit(S, rd, prm, choice, W, pool, say, tic_slice=tic)
        finalize(S, rd, prm, say)
        say(status(rd))
    except (KeyboardInterrupt, Interrupted) as e:
        say(f'interrupted ({type(e).__name__}): slice {S["slice"]} -- phase files on disk are the resume point'); code = 130
    finally:
        clean = _shutdown_pool(pool.pool, say, graceful=(code == 0), timeout=30.0 if code == 0 else 10.0)
    if code or not clean:
        sys.stdout.flush(); sys.stderr.flush()
        os._exit(code)                                    # a killed worker must not be joined again at interpreter exit


def status(rd):
    """One line: run, slice, phase, T_d, launched, covered, remaining, fuel, sumJi, done."""
    rd = pathlib.Path(rd); p = rd / 'state.json'
    if not p.exists():
        return f'{rd.name}: no state.json'
    S = C.jload(p); fn = _fleet_numbers(S)
    return (f'{S["run"]}: slice {S["slice"]} phase {S["phase"]} T={C.s2d(S["T"]):.0f} d  launched {fn["launched"]}/{S["N"]}  '
            f'covered {fn["covered"]}  remaining {len(S["remaining"])}  fuel {fn["fuel"]:.0f} kg  sumJi {fn["sumJi"]:.4f}  '
            f'done {S["done"]}' + (f' ({S["stop"]["kind"]})' if S.get('stop') else ''))


def main(argv=None):
    a = parse_args(argv)
    if a.cmd == 'run':
        run(a)
    elif a.cmd == 'status':
        print(status(C.run_dir(a.out)))
    elif a.cmd == 'misses-prize':
        if not a.out2:
            raise SystemExit('usage: s18_rhfa.py misses-prize RUN_DIR OUT.json [--bonus 1.0 --prev PREV.json --decay 0.5]')
        misses_prize(a.out, a.out2, bonus=a.bonus, prev=a.prev, decay=a.decay)
    else:
        rd = C.run_dir(a.out); S = C.load_state(rd); finalize(S, rd, S['params'], C.logger(rd / 'log.txt'))


if __name__ == '__main__':
    main()
