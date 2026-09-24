# Stage 18: RHFA, the rolling-horizon fleet auction (2026-09-24 01:00 CST)

A from-launch, concurrent construction of an N-craft fleet. All craft are advanced together through the 15-year
window in time slices; inside each slice every craft proposes segment plans built by a single-craft twin beam, and
one exact assignment (the "auction") picks one plan per craft so that no target is taken twice.

## 0. Why this, and what it must beat
Facts from stages 1-17 (docs/stage17_brief.md, docs/stage17/PLAN.md):
- The leaders' J < 12 = 7-8 routes ALL at full cadence (~130 d per flyby for 15 yr) and mutually disjoint. Their
  total delta-v equals ours; they distribute it evenly. Our 4 deep routes have that cadence; our 5 tails (starved
  of "fillers") fly at 175-210 d.
- Every method that reuses existing tracks tops out at 285-286 targets with 8 craft (law 9). The leftovers are
  "pins": targets with <= 4 node windows in 15 yr.
- Sequential from-launch construction starves late routes (fills 228-249) because early routes hoard fillers.
- A route's depth is bought with fillers (>= 10 windows); pins cost nothing extra inside a deep route.
- Depth needs a candidate pool ~2.5x the route: partition-first (assign, then plan) collapses depth.
- Fuel does not bind at N = 8: the 09-26 bar allows 3764 kg in total (mean m0 1070); J < 12 needs about 2800 kg.

RHFA is designed against exactly these: no route is "last" (all advance together), no pre-assignment (each beam
sees the whole remaining pool), fillers are rationed by a per-craft quota, and pins are taken when their last
windows arrive (earliest-deadline prizes), by whichever craft can take them cheapest (the auction).

Targets: N = 8 covering 298 at <= 3764 kg (raw < 13.609, submit 09-26 12:00). Stretch: <= 3000 kg (raw ~12.9).

## 1. Architecture

### State (per run, on disk, resumable): `results/s18/<run>/state.json` + `craft_<i>.npz`
- `T` current commit epoch (s); `slice` index; `params`.
- `craft[i]`: `launched` (bool), `route` (npz, run_ialns ist format = settled impulsive twin: tL, vinf, ts, Ts, tf,
  asts), `targets` (list), `counts` {pin, med, fil}, `tank` (kg, honest = after regrid + settle).
- `remaining`: targets not yet taken (start: all 298; 131/144 excluded).
- `log`: per slice, the columns offered and chosen, coverage, fuel.

### Slice loop  (T_0 = 0; horizon H; lookahead L; both in days; defaults H = 540, L = 360)
1. **Propose.** In parallel (one process per craft, 8 cores):
   - launched craft i: a twin beam from its settled prefix over `remaining`, with candidate encounter epochs limited
     to (T, T + H + L]. Reuse `tools/s14_twinbeam.py`: `epoch_candidates`, `lin_price`, `child_problem`,
     `settle_child`, `search_state`, `twin_score`. Beam width B (8), tries per level (16), levels until no child.
     Children are settled twins (max miss <= 150 km).
   - unlaunched craft: ONE shared root beam (`s14_twinbeam.root_state` / `ctoc14.search.roots`, launch grid
     10 d over (T, T + H], |v_inf| <= 4) whose front is kept DIVERSE (all settled states per depth, binned by
     launch date); its columns are interchangeable among the unlaunched craft.
   - Output per beam: every settled state = (targets, epochs, twin file, tank). A **column** is a state whose last
     flyby epoch <= T + H (the committed part); its `potential` = max number of extra flybys any descendant reaches
     in (T + H, T + H + L]. The committed twin is what gets adopted; lookahead flybys are never committed.
2. **Auction** (exact, scipy `milp`, seconds):
   - x_c in {0,1} per column; per launched craft sum x_c <= 1; root columns: sum x_c <= number of unlaunched craft;
     per target sum_{c contains t} x_c <= 1;
   - cumulative class quotas per craft: fillers <= q_f (10), medium <= q_m (16); pins unlimited; a column that
     would exceed a quota is dropped before the MILP;
   - fuel cap per craft: committed tank <= m0_max (1120 kg default);
   - objective: max sum_c [ sum_{t in c} (1 + edf_s(t)) + gamma * potential_c - lambda * dtank_c / 100 ].
     Defaults gamma 0.3, lambda 0.15 (i.e. 1 target ~ 670 kg; fuel is secondary, as measured), portfolio over both.
   - **edf_s(t)** (earliest deadline first): windows(t) = the target's node events with d < 0.05 AU
     (results/s17/physics/events.json; a copy under results/s18/windows.json is fine). Let w = number of windows
     after T + H. edf = beta if w = 0 (now or never), beta/2 if w = 1, beta/4 if w = 2, else 0. Default beta 1.5.
   - class(t): pin if total windows <= 4, filler if >= 10, medium otherwise (results/s17/physics classes).
3. **Commit.** Each craft adopts its chosen column's twin; regrid to regular 20 d bins (`tools/s15b_regrid.regrid`)
   and `run_ialns.settle` (law 8: honest tank); update `remaining`, counts, T <- T + H. A craft whose column fails
   to re-settle keeps its previous state (the column's targets go back to `remaining`).
4. Stop when T + H > 15 yr - 60 d, or when `remaining` is empty. Final: `fleet/route_<i>.npz` + `fleet.json`
   (covered, misses, per-craft n/tank, twin sum J_i).

### After the loop (existing tools)
closer `tools/s15b_close.py` for the misses -> `tools/s15b_relocate.py` -> regrid all -> exact selection
`tools/s16_select.py` (N <= 8, K = 298, over the run's routes + honest pool) -> export
(`results/s16/export_s16a.sh` pattern: run_ialns export -> combine_submission -> validator2).

## 2. Modules (new files only; never edit existing tools)
- `tools/s18_windows.py`: build `results/s18/windows.json` {ast: [event times (d)], class} from
  results/s17/physics/events.json (or recompute from MEA.txt: torus distance minima < 0.05 AU).
- `tools/s18_segbeam.py`: `propose(route_npz | None, remaining, T, H, L, B, tries, seed) -> columns` (the per-craft
  slice beam; root beam when route is None). Columns are dicts {craft, targets, epochs, twin npz path, tank,
  dtank, potential}. Must be callable as a subprocess job (JSON in / JSON out) for the pool.
- `tools/s18_auction.py`: `choose(columns, state, params) -> {craft: column}` (MILP; quotas; edf; fuel cap).
- `tools/s18_rhfa.py`: driver: `run OUT --N 8 --H 540 --L 360 --B 8 --tries 16 --beta 1.5 --gamma 0.3 --lam 0.15
  --qf 10 --qm 16 --m0max 1120 --nproc 8 --seed 0`; resumable from state.json; logs one line per slice:
  `[slice s T=.. ] launched k/N  covered C  remaining R  fuel F  chosen: i:n(+dtank) ...`.
- `tools/s18_eval.py`: fleet report (coverage, misses by class, per-craft cadence, twin sum J_i, honest check:
  0 irregular nodes, every 20 d window <= cap) and the gate verdicts below.
- `tools/s18_portfolio.py`: runs a list of parameter sets sequentially (each uses --nproc 8), skipping finished
  runs; writes `results/s18/portfolio.jsonl`.

## 3. Gates (fixed now; the builders do not move them)
- **G0 smoke** (integrator): N = 2, H = 540, L = 360, 2 slices: both craft launch, every committed twin re-settles
  after regrid, no target taken twice, log and resume work (kill the process mid-slice 2, restart, same result).
- **G0b day-1461 check** (first N = 8 run, stop after T >= 1461 d): >= 72 distinct flybys across the fleet
  (s16a's 9 craft had 70 by then) and mean fuel <= 160 kg per launched craft. KILL if <= 60.
- **G1 first full N = 8 fleet:** covered >= 268 at total fuel <= 3900 kg -> portfolio. 250-267 -> one parameter
  round (H, beta, quota) then re-judge. < 250 -> KILL (history: 228-249).
- **G2 portfolio + closer:** best 8-craft fleet >= 292 covered with lin-priced closure of the misses <= twin sum J_i
  11.50 -> close, relocate, select, export. Fallback: any 9-craft fleet with twin sum J_i <= 11.64 (beats s16a by
  the 09-25 bar) -> export it.
- Clock: raw J must be < 13.511581 / k(t): 13.658 at 09-25 12:00, 13.609 at 09-26 12:00, 13.560 at 09-27 12:00.
  Export + validation takes 3-4 h. Nothing is ever submitted by an agent.

## 4. Rules
- Python `~/.venvs/astro313/bin/python`, `nice -n 5`, OMP/OPENBLAS threads 1; long jobs under nohup with on-disk state.
- Write only under `tools/s18_*.py`, `results/s18/`, `docs/stage18*`. Never modify existing tools/*.py, ctoc14/*.py,
  results/CTOC14_Result_*.txt, validator logs, results/s16/best. Never use `--swap`.
- Honesty: every reported tank is after regrid + settle. Separate MEASURED from estimated.

## 5. Fixer pass (2026-09-24 05:00 CST) -- reviewer findings applied before G1
All inside tools/s18_*.py; defaults in s18_common.DEF changed accordingly (the G1 run uses these).
1. FATAL pool hang: `s18_rhfa.PoolBox` (apply_async + polling under a deadline: wall_job + job_grace 600 s per beam job,
   regrid_timeout + 120 s per commit); a lost / raised job is logged `LOST` / `RAISED`, counted as no proposal, and the
   pool is re-created.  Verified with the reviewer's crash pattern (scratchpad pool_crash2.py: SIGKILLed worker returns
   None after the deadline, a raising worker is reported, the re-created pool serves the next map).
2. Mission tail: `s18_common.slice_H` -- the last slice absorbs the tail (H_eff = T_END - T when <= H + 120 d; H 540 ->
   slice 9 is 616.8 d, ends at 5476.75 d; H 630 -> last slice 436.8 d).  `S['H_eff_d']` is fixed per slice by
   `start_slice` and used by the prizes, jobs, check_column, auction and commit; recorded in jobs/index.json and the log.
3. Fuel bar and pacing: m0max 1070 (= 600 + 3764 / 8; slack 0 in the auction), `ramp_cap` = 601.5 + 468.5 x
   min(1, (T + H_eff) / T_MISSION)^0.9 + 60 kg (caps 720 / 770 / 818 / 909 kg at 540 / 1080 / 1620 / 2700 d) in
   filter_columns (reason `ramp`) and at the commit (+ m0slack 10 kg regrid tolerance); one fleet fuel row per auction:
   sum x_c max(dtank_c, 0) <= 3764 x H_eff / T_MISSION x 1.3 (482 kg per 540 d slice; the retry auction subtracts what
   the slice already committed); lam starts at 1.0 and is the dual of that budget (`update_lam`: x1.5 when the row is
   >= 95 % used or the ramp removed >= 20 % of the admissible columns, /1.25 when < 50 % used, floor prm lam, cap 6).
4. Quotas: qf 12 / qm 18, lifted in the last 3 slices (`quota_active`); the beam is quota-aware (classes at the quota
   are withheld from the craft job's `remaining`, the head-room travels as job['quota'] and `beam.rem_for` excludes a
   class once a parent's new targets fill it).
5. Cached failures: `s18_segbeam.run_job` re-runs a cached ok-False result once (corrupt ckpt.pkl deleted, the old
   result kept as result_failed.json, `retried` flag); the driver logs a cached failure explicitly.
6. Beam stopping: end after `look_stop` 2 consecutive levels with in_H = 0 (`end.kind = lookahead`); wall_job back to
   2400 as a safety; settle tries capped at 8 after a level with < 30 % settle success.
Minor: check_column asserts launch < first flyby and root launch in [T, T + H]; s18_portfolio.finished requires
state.done with stop kind horizon / covered; portfolio BASE follows the new defaults; results/s18/review/g0_cc deleted.
Not done: the root-grid tail for late launches (small while all 8 launch in slices 0-1) and the G0b day_check 1620
(the 1461 measure is comparable with s16a's, see the integrator report).
Tests: s18_auction selftest OK; fake-proposer 10-slice loop (results/s18/_fix_loop: tail slice 616.8 d, quotas off in
slices 7-9, dual moving 1.0-3.9, quota withholding logged); rewind from slice 5 reproduces the state; G0 re-run
results/s18/g0_fix.

## 6. Launch pacing options (2026-09-24 08:00 CST; additive, both OFF by default)
Motivation: g2_plan (planbeam, N 8) launched all 8 craft in slice 0 (days 0-260); craft 7, launched last into lanes the
other 7 had already crowded, collapsed to 11 flybys.  Two driver options, frozen in state['params'] like every other
run parameter (an older state.json without the keys resumes with both off; a CLI value on such a resume is logged and
ignored):
- `--max-launch K` (s18_rhfa.parse_args / params_of; s18_auction.launch_cap): at most K root columns are chosen per
  slice: the auction's root row becomes sum root x_c <= min(K - launches_used, n_unlaunched) (launches_used = launches
  the slice's earlier commit-retry rounds already accepted, passed by s18_rhfa.commit); root columns offered when the
  cap is exhausted are dropped with reason 'max_launch'.  K = 0 (default) = no cap = the previous behaviour.
- `--launch-stagger D` (days; s18_auction.launch_stagger / build_milp rows ('stagger', k)): the k-th launch of a slice
  by epoch (k = 0, 1, ..; the root assignment hands it to the k-th unlaunched craft) must satisfy t_launch >= T + k D,
  written as cumulative rows sum_{root c: t_launch_c < T + (k + used) D} x_c <= k for k = 1 .. cap - 1 (rows that can
  never bind are omitted).  It is an auction constraint over the existing root proposals (the root jobs already cover
  the whole launch window (T, T+H] in contiguous sub-grids), not a restriction of the root grid.  D = 0 (default) = off.
MEASURED (nice 5, <= 2 cores, while p1_ration held 9): s18_auction selftest OK (the original instance + a cap/stagger
block: no cap -> 2 launches, cap 1 -> the best root only, cap 2 + stagger 100 d -> 2nd launch at 200 d, cap exhausted
-> every root dropped 'max_launch', parameters without the keys -> the uncapped choice bit-identical); fixture_fake
driver runs `--N 3 --nproc 2 --proposer fixture_fake --commit none --stop-T 1620` (13 s each): results/s18/_opt_fake
(no cap) launches 3/3 in slice 0 as before; results/s18/_opt_fake_ml1 (`--max-launch 1`) launches one craft per slice
(craft 0 slice 0 @ 310 d, craft 1 slice 1 @ 680 d, craft 2 slice 2 @ 1530 d), same coverage 44 after 3 slices.
CLI: `nice -n 5 ~/.venvs/astro313/bin/python tools/s18_rhfa.py run results/s18/<run> --N 8 --proposer planbeam --nproc 8
--max-launch 2 --launch-stagger 60`.

## 7. Prize feedback and plan claims (2026-09-24 09:40 CST; additive, both OFF by default)
Motivation: g2_plan (planbeam, N 8) covered 224/298 at 3189 kg (p1_ration 216, p2_count 201); the 74 misses are 24 pin /
35 med / 15 fil.  The fleet's per-craft pace collapses late because the slice-wise re-planning against a changing pool
breaks long lanes, and stranded targets are found only at the end.  Two driver options, frozen in state['params']
like every other run parameter (an older state.json resumes with both off; a CLI value on such a resume is logged and
ignored); with both off every file of a run is reproduced bit-identically (MEASURED below).
- **A. `--prize-json PATH`** (s18_rhfa.parse_args / params_of / load_prize_json; s18_auction.prize_bonus_of): a JSON
  {target_id: bonus}.  The bonus is added to the auction prize of those targets in EVERY slice (s18_auction.value_of:
  1 + edf + bonus per new target) and to the proposer jobs' prize dict (s18_rhfa.prizes_for: job['prize'] = 1 + edf +
  bonus, built from 1.0 even when --beam-edf 0), so the beam / planner steers toward them too.  The table is frozen as
  params['prize_bonus'] (the path as params['prize_json']) when the run is created.  Helper
  `s18_rhfa.py misses-prize RUN_DIR OUT.json --bonus 1.0 [--prev PREV.json --decay 0.5]` (s18_rhfa.misses_prize)
  writes {miss: bonus} for the misses of a finished run (RUN_DIR/fleet/fleet.json) merged with a previous prize file
  (prev values x decay, then + bonus for the current misses): the miss-feedback iteration "run k+1 prizes = misses of
  run k".
- **B. `--claims`** (`--claim-bonus` 0.5, `--claim-ttl` 2; s18_rhfa.claims_on / claim_owner / make_jobs / update_claims):
  after each slice's commit every accepted column's `plan_next` (s18_planbeam.candidates / settle_columns / columns_of:
  the plan's flybys beyond the committed part, in plan order, also in stats.plan; a settle-truncated column's plan_next
  starts with the cut legs; segbeam / fixtures emit none) becomes the claims of that craft, state['claims'] =
  {"<craft>": {"<target>": slice claimed}}.  In the next slice's jobs the other craft AND the root jobs get `remaining`
  minus the targets claimed by others (they cannot plan them), the claimant keeps its claims and gets them prized
  +claim_bonus (job['claims'] lists them; jobs/index.json records the claims in force).  Claims are rewritten every
  slice from the chosen columns: a claim not in the claimant's new plan_next (or of a craft without an accepted
  column) is released, a claim with age >= claim_ttl slices is released even when still planned, a target planned by
  two craft goes to the lower index (logged as a conflict).  Logged per slice: `claims in force: c0:3 c1:2 ..` before
  the jobs and `claims: c0:2 .. (kept / new / released / expired / conflicts)` after the commit; the per-slice log
  record carries `claims`.  Persisted in state.json and in every slices/s<kk>/in/state_in.json snapshot (resume and
  rewind see the same claims).  Test switch: `S18_FIXTURE_PLAN_NEXT=K` makes fixture_fake emit a deterministic fake
  plan_next (the K smallest remaining targets outside the column, shifted by one on odd slices).
MEASURED (nice 5, <= 2 cores, while p3_wide held 9; every check against the INSTALLED files):
- s18_auction selftest OK (the previous blocks + prize bonus {6: 3.0}: c0a + c1a + r1 = 7.91 instead of c0b + r1 = 5.51,
  other columns' values unchanged, an empty / absent table bit-identical).
- fixture_fake `--N 3 --nproc 2 --proposer fixture_fake --commit none --stop-T 1620` (12 s): results/s18/_opt2_fake
  without the flags == results/s18/_opt_fake (remaining, taken, craft targets / tanks, log records, the values of all
  three auctions, commit.json accepted entries); the only difference is the new params keys (prize_json None,
  prize_bonus {}, claims 0, claim_bonus 0.5, claim_ttl 2).  results/s18/_opt2_fake_prize with `--prize-json
  results/s18/_opt2_prize_test.json` ({119: 4, 288: 4}): slice-0 column values +8.0 (s00_r1_k0: 5.289 -> 13.289) and
  +4.0 (s00_r1_k3), all others unchanged, objective 14.759 -> 22.700, craft 2's column flips s00_r1_k1 -> s00_r1_k0;
  job s00_r0 prize[119] = 5.0 (1 + edf 0 + 4).  results/s18/_opt2_fake_claims (`--claims`, S18_FIXTURE_PLAN_NEXT=3):
  slice 0 c0 claims {1, 2, 3} (6 conflicts: every root column planned the same three), slice 1 c0 3 own prized 1.5 =
  1 + edf + 0.5 (0 mismatches over every job's prize dict) and 3 withheld from c1 / c2 and from the root pool, slice 2
  target 2 (claimed in slice 0) expired at age 2, 1 released, 2 conflicts; the claims rebuilt from the accepted
  columns' plan_next equal the state at every slice.  Resume (stop at 1080 d, continue to 1620 d) and rewind (slice 2
  auction.json + commit.json deleted, re-run from the cached columns) reproduce the single run's claims / remaining /
  taken / log exactly; a resume of a run created without the options logs `CLI --claims 1 ignored (frozen 0)` /
  `--prize-json .. ignored (frozen None)`.
- planbeam test-slice, craft 4 of g2_plan slice 3 (results/s18/_opt2_pb, the job resumed from g2_plan's own plan.pkl,
  59 s, 32/32 settles, 16/16 regrids): 16 columns bit-identical to g2_plan's s03_c4 on every pre-existing key, all with
  `plan_next` (12 non-empty, e.g. k14 [189] -> [151, 286, 92, 94, 129], k0 [189, 151, 286, 92, 94] -> [129]).
- `misses-prize results/s18/g2_plan results/s18/prizes_g2.json --bonus 1.0`: 74 entries (= the 74 misses, all 1.0);
  with `--prev {119: 4, 288: 4} --decay 0.5`: 119 -> 3.0 (4 x 0.5 + 1).
CLI:  `nice -n 5 ~/.venvs/astro313/bin/python tools/s18_rhfa.py run results/s18/<run> --N 8 --proposer planbeam --nproc 8
--prize-json results/s18/prizes_g2.json`;  `... --claims --claim-bonus 0.5 --claim-ttl 2`;  iteration:
`tools/s18_rhfa.py misses-prize results/s18/<run_k> results/s18/prizes_k.json --bonus 1.0 --prev results/s18/prizes_k-1.json
--decay 0.5` then `run results/s18/<run_k+1> ... --prize-json results/s18/prizes_k.json`.  Nothing has been run at N 8
with either option yet; the calibration of the bonus (1.0 = one extra target per prized flyby) and of claim_ttl is open.
