# Stage 18 PLANBEAM: plan long, commit short (2026-09-24 06:45 CST; review fixes 07:05, section 5)

`tools/s18_planbeam.py` is a drop-in PROPOSER for the RHFA driver (`s18_rhfa.py run ... --proposer planbeam`; the
driver's `resolve_proposer` got one additive branch). Same JOB / RESULT contract as `s18_segbeam.run_job` (result.json
cache, cached failure re-run once, ok False on error, single-threaded, `col_<k>.npz` under the job dir, every column
passes `s18_common.check_column`). Test dir: `results/s18/planbeam_test/` (everything below is MEASURED there, nice 5,
2 cores, while the G1 run held 8 cores; load average 8-11).

## 1. Diagnosis it answers
Segbeam (twin beam, width 8, every child settled, horizon T + H + L) on `g1_base` slice 3 (T = 1620 d, H 540, L 360):
per craft the cheapest column costs 2-21 kg for 1 target, 16-60 kg for 2, 52-96 kg for 3, 90-220 kg for 4. s16a's deep
routes fly 4.2 flybys per 540 d at ~8 kg each. The short, narrow, settle-every-child beam cannot set up cheap sequences.

## 2. What planbeam does (per job)
1. **Plan** with the production planner (`ctoc14.search`: Lambert 15-400 d + tof refinement, LinLeg 20-600 d, m0 1600,
   `s14_twinbeam.params(dv 1.2, drmax 0.15, mc 150)`, n_per_target 6) as a WIDE depth-synchronous beam (width 150) over
   the horizon `[T + 5 d, min(T_MISSION, T + H + L + 720 d)]` (`Params.win_lo/win_hi` on every target with
   `mandatory = []`: a horizon, not a deadline), prizes `1 + edf` from the job (everything not remaining excluded; classes
   whose quota head-room is 0 excluded), `w_t 3.0` (cadence matters inside a 540 d window; see 4), wall cap 240 s
   (SIGALRM; the beams collected so far are used). Every beam of every depth is collected (`P.collect`).
   - craft job: start = `s14_twinbeam.search_state` of the prefix twin (`s18_segbeam.prefix_state`). **The first level
     is the bottleneck**: from a twin end state the production caps give 0-29 children (c0 stranded at 1363 d: 0; c5 at
     1130 d: 1) because the twin's end velocity is not tuned for a continuation. The first level is therefore the UNION
     of the relaxed planner legs (dv 2.5 km/s, 0.30 AU: 26-109 children) and the segbeam's own level-1 candidates
     (`s18_segbeam.price_parent_full`: whole-prefix linearised twin price; the child state is one integration of the
     linear step: 59-97 children); later levels are production.
   - root job: `search.roots` over the job's launch grid (launches in (t_lo, t_hi], >= T), first flyby <= T + H, deduped
     on (first target, 20 d launch bin), min(`roots_max` 450, 3 x width) kept (`search._beam_loop` keeps at most
     3 x width of its initial front, so a larger roots_max would be silently cut).
2. **Truncate** every planned state at T + H: for every prefix length k of its new legs with epoch <= T + H (cut where a
   class quota head-room would be exceeded) one candidate keyed by the committed target SET; `potential` = the state's
   planned flybys in (T + H, T + H + L] (segbeam's definition, emitted in the column: `pot_mode 'lookahead'`);
   `pot_cont` = the planned flybys in (t_last_committed, T + H + L], i.e. also the in-window flybys the truncation
   dropped: the proposer RANKS by it (dedupe: best (pot_cont, -fuel); frontier per depth), so a shallow truncation of a
   good plan ranks above a dead end, while the auction's gamma keeps segbeam's calibration (`pot_mode 'continuation'`
   emits pot_cont instead; `stats.plan[col].plan_cont / plan_lookahead` show both). Committed planner fuel = the fuel
   of the beam state with exactly those legs AND that launch (an estimate from the leg dv when that ancestor was not
   kept). Every candidate carries the launch (t_launch, v_inf) of the plan that won its dedupe.
3. **Frontier**: per depth 1..8 the best-continuation candidate + the 2 cheapest, then the rest by (continuation, fuel)
   round-robin over the first new target, 16 per job (`per_depth 3`, `max_cols 16`). When the per-depth pass alone
   exceeds 16 (more than 5 in-window depths, i.e. root / early slices) the cut keeps every depth's rank-0 column, then
   rank 1, ..., DEEPEST first within a rank (6 depths x 5: kept 2/2/3/3/3/3; the first version cut by (n_new, fuel)
   and lost the deepest columns: reviewer finding).
4. **Settle** each committed prefix by chaining `s14_twinbeam.lin_price` (first seed) + `settle_child` from the prefix
   twin (root: `root_state` of the candidate's OWN planner launch; settled once if |v_inf| needed thrust), with a
   cache over leg prefixes (columns share their prefixes: 16-29 settles for 16 columns; settle and regrid results are
   both cached, so a failed prefix is never re-tried). A failed leg truncates the column; a settled flyby that drifted
   out of (T, T + H] rejects it.
5. **Regrid** (`regrid_cols 1`): every column twin goes through `C.regrid_settle` and the HONEST twin is emitted,
   re-packed from the regridded problem (`C.pack_state`: tank, miss, t_end, r_end, v_end, tfs all from the same
   twin; the commit's regrid is then idempotent, -0.04..0 kg). A column whose regrid fails is dropped: on c0's
   stranded prefix 6 of 7 fresh first legs (junction 'lam' / 'lin' seeds, 1-33 irregular nodes on a 270 d coast)
   failed the regrid law with misses of 10^4-10^5 km; c1 / c4 / c7 / root: 0 failures.
6. **Budget** (`job_wall 300 s`): one budget for the whole job. The plan gets min(`plan_wall` 240, job_wall - elapsed
   - `settle_reserve` 90) = 210 s on a fresh job; every settle / regrid SIGALRM is capped at what is left; candidates
   beyond the budget are dropped (`stats.budget = {job_wall, plan_cap, settle_left_s, hit, n_unsettled}`, `hit` in
   {None, 'settle_wall', 'job_wall', 'regrid_streak'}). After `regrid_fail_max` 4 consecutive regrid failures only
   candidates that EXTEND a regrid-ok prefix are still settled (c0: its 2-target column extends the 1-target one that
   passed; the other 1-target roots are skipped). MEASURED c0: 288 -> 220 s with the same 2 columns.

Checkpoint `plan.pkl` (a killed job does not re-plan). Parameters: `PB_DEF` in the module, overridable per run by the
env var `S18_PLANBEAM='{"width": 300}'` (inherited by the driver's workers) or per job by `job['planbeam']`.

## 3. MEASURED frontier, slice 3 of g1_base (T 1620 d), crafts 0, 1, 4, 7 (`results/s18/planbeam_test/s3/frontier.txt`)
cheapest dtank [kg] per number of committed targets; potential of that column / max potential at that depth; the
planbeam tanks are HONEST (regridded), segbeam's are planner-twin tanks (the commit regrids them, -0.3..+0.1 kg).

| craft | n_new | planbeam kg (pot) | segbeam kg (pot) | planbeam wall | segbeam wall |
|:--|:--|:--|:--|:--|:--|
| c0 (stranded, t_end 1363 d) | 1 | 10.3 (2/2) | 9.8 (0/1) | 288 s (2 columns: the regrid law dropped the rest) | 715 s |
| | 2 | **15.5** (0/0) | 60.4 (0/1) | | |
| | 3 | - | 82.8 (0/0) | | |
| c1 (13 fb, t_end 1619 d) | 1 | **7.1** (0/3) | 20.9 (0/3) | 138 s (15 columns) | 368 s |
| | 2 | **21.4** (0/3) | 37.0 (0/3) | | |
| | 3 | **48.5** (1/2) | 72.6 (2/2) | | |
| | 4 | - | 90.5 (2/2) | | |
| c4 (8 fb, t_end 1438 d) | 1 | 1.1 (2/3) | 2.6 (3/3) | 130-156 s (16 columns; three runs) | 482 s |
| | 2 | **6.0** (0/3) | 16.1 (3/3) | | |
| | 3 | **11.9** (0/3) | 57.6 (3/3) | | |
| | 4 | **30.0** (1/1) | 184.0 (0/1) | | |
| c7 (8 fb, t_end 1596 d) | 1 | 1.0 (0/3) | 1.8 (2/3) | 111 s (16 columns) | 377 s |
| | 2 | **17.8** (1/3) | 23.2 (0/3) | | |
| | 3 | **25.8** (0/1) | 96.1 (2/2) | | |
| | 4 | - | 219.1 (0/0) | | |

Wall per job (single-threaded, 2 jobs at a time): 111-288 s, mean ~170 s; plan 68-96 s (9-12 levels, 1060-1470
states), settles 28-170 s (c0's failing 'lam' settles take 17 s each), regrids ~1 s per column. Segbeam: 368-715 s.
c4's 61.9 s in `s3/frontier.json` is not a measurement (a stray check ran that job concurrently); the three other runs
of the same job gave 130, 135, 156 s.

Reading: at 2-4 committed targets planbeam's columns cost 1.4-6x less than segbeam's (c4: 3 fb 11.9 vs 57.6 kg, 4 fb
30 vs 184; c7: 3 fb 25.8 vs 96.1; c1: 3 fb 48.5 vs 72.6), at 1 target they are equal (both ~1-10 kg), and the wall
is 2.5-4x shorter. The kill-line pace (~28 targets per slice at ~8 kg) is NOT reached: the deepest honest columns
are 3-4 targets per 540 d (c1 / c7 offer no 4-target column at all; c4's 4-target column is 30 kg, i.e. 7.5 kg per
flyby -- s16a's pace -- but its potential is 1). The plans themselves do reach 9-12 flybys by 3240 d, i.e. the
planner's own cadence from these prefixes is 135-180 d per flyby; the committed window simply holds 3-4 of them.

## 4. Variants (crafts 4, 7; `results/s18/planbeam_test/var/*.out`)
| variant | c4: 1 / 2 / 3 / 4 / 5 fb [kg] | c7: 1 / 2 / 3 fb [kg] | plan wall |
|:--|:--|:--|:--|
| base (w_t 1.0 production, width 150) | 1.1 / 6.0 / 11.9 / 67.0 / - | 1.0 / 17.8 / 25.8 | 70-71 s |
| **w_t 3.0, width 150 (DEFAULT)** | 1.1 / 6.0 / 11.9 / **30.0** / - | 1.0 / 17.8 / 25.8 | 79-92 s |
| w_t 1.0, width 300 | 2.3 / 6.0 / 11.9 / 30.0 / 58.8 | 1.0 / 17.8 / 30.3 | 157-170 s |
| w_t 3.0, width 300 | 2.3 / 6.0 / 11.9 / 24.4 / - | 1.0 / 17.8 / 25.8 | 136-149 s |

The frontier at 1-3 targets is the same in every variant (it is set by the first level, i.e. by the prefix); the
4-5 target columns depend on the beam's cadence pressure (w_t) and width. Width 300 stays under the 240 s plan cap
but doubles the wall; `S18_PLANBEAM='{"width": 300}'` switches it on for a run.

Root job (slice 0, `s00_r0`: launches 0-70 d, `results/s18/planbeam_test/root/`): 16 honest columns, 1-4 targets
(0 / 0 / 3.6 / 23.5 kg), 82 s; segbeam's r0 offered 31 columns up to 6 targets (47 kg) in 226 s -- at launch the
segbeam's settle-every-child beam is deeper (all its states are exact zero-thrust Lambert legs early on); planbeam's
committed depth is capped by the planner's own in-window cadence (4 of 150 depth-4 states end before 540 d).

## 5. Review fixes (07:05; `results/s18/planbeam_test/fix`, `fix_root`, `fix/check_contract.txt`)
Applied: (major) frontier cut deepest-first (section 2.3, unit test `fix/unit_frontier.py`); (major) one `job_wall`
budget shared by plan / settles / regrids + regrid failure streak (2.6); (minor) root candidates settled from their own
plan's launch (2.2/2.4; the first version re-derived the launch from the first-leg key: 40 of 153 keys mapped to 2-3
launches, 9 of 23 root settles failed; now 0 of 23 fail, 16 columns, launches 0 / 40 / 50 / 60 d, 81 s); (minor)
regridded twin re-packed (2.5); (minor) potential of a truncation documented, `pot_cont` ranking, `pot_mode` (2.2);
(minor) `roots_max` 450 = the `_beam_loop` cap (2.1). Not done: a fleet run with `--proposer planbeam` (8 cores).

Re-measured slice 3 crafts 0 and 4 (defaults, 1-2 cores at nice 5 while G1 held 8; `fix/frontier.txt`):

| craft | n_new | planbeam kg (pot) | segbeam kg (pot) | planbeam wall | budget |
|:--|:--|:--|:--|:--|:--|
| c0 (stranded) | 1 | 10.3 (2/2) | 9.8 (0/1) | 220 s (plan 58, settle 163; 2 columns) | hit regrid_streak (1 candidate skipped) |
| | 2 | **15.5** (0/0) | 60.4 (0/1) | | |
| c4 | 1 | 2.2 (2/3) | 2.6 (3/3) | 142 s (plan 80, settle 62; 16 columns) | not hit |
| | 2 | **6.0** (0/3) | 16.1 (3/3) | | |
| | 3 | **11.9** (0/3) | 57.6 (3/3) | | |
| | 4 | **30.0** (1/1) | 184.0 (0/1) | | |

c4's 1-target column moved 1.1 -> 2.2 kg (the dedupe now keeps, per target set, the plan with the longer
continuation; its committed epoch differs); everything else is bit-identical to section 3. Root job s00_r0: 1 / 2 /
3 / 4 targets at 0 / 0 / 3.6 / 23.5 kg, 16 columns, 81 s (unchanged frontier, no settle failures). Contract check
(reviewer's `review/check_contract.py`) on all 34 columns of the three jobs: 0 violations; re-regrid of 3 columns:
-0.02 / -0.04 / 0.00 kg. `s18_auction.py selftest` OK; `resolve_proposer('planbeam')` is `PB.run_job`.
Cheapest-n3 set over the four crafts (c1 / c7 from section 3, unchanged code path): {11.9, 25.8, 48.5, none} vs
segbeam {57.6, 96.1, 72.6, 82.8}: median 25.8 kg over the crafts that offer one (bar: <= 45 kg).

## 6. What was not done / next levers
- A G1-style fleet run with `--proposer planbeam` (needs 8 cores; G1 holds them until ~07:15). Expected per-slice
  wall: max over 8 jobs ~5 min (vs ~12 min for segbeam), i.e. a full 10-slice run in ~1 h.
- The first level from a stranded prefix (c0 / c5) is still poor: the relaxed and lintwin children exist but most of
  their twins fail the regrid law. A `restart` option (a long coast re-seeded by the LinLeg profile at the regular
  nodes rather than one junction impulse) is the obvious fix.
- Column potential is the planner's count in (T + H, T + H + L]; the auction's gamma weighs it. With planbeam the
  potential is informative (it comes from a 9-12 flyby plan), so gamma could be raised.
- `per_depth` / `max_cols` trade settle wall for auction choice: 16 columns cost 28-64 s of settling on a normal prefix.

## 7. pot_mode 'full' + pot_cap (2026-09-24 08:00 CST; additive, default unchanged)
g2_plan (N 8, planbeam, gamma 0.3) covered 224 at 3189 kg with routes of 27-33 flybys at 158-203 d cadence: the
lookahead potential (planned flybys in (T + H, T + H + L]) is 0-3 on every column and cannot tell the auction which
column leads into a long continuation.  `PB_DEF.pot_mode 'full'` (settable per run by `S18_PLANBEAM='{"pot_mode":
"full"}'`, per job by job['planbeam'], on the CLI by `--pot-mode full`) emits as the column `potential` the planned
flybys in the WHOLE continuation (t_last_committed, plan end] (plan end = min(T_MISSION, T + H + L + plan_extra_d), i.e.
3240 d at T 1620); in that mode the candidate dedupe and the frontier ranking use the same full count (`candidates`
sets `rank` = pot_full, `_rk`; the default modes rank by pot_cont exactly as before).  `pot_cap` (default None; e.g.
`{"pot_mode": "full", "pot_cap": 8}`) caps the emitted potential (`column_potential`).  stats.plan carries plan_lookahead
/ plan_cont / plan_full / potential / pot_mode per column.  The auction multiplies `potential` by gamma, so gamma must be
re-calibrated for this mode (a full potential of 10 at gamma 0.3 is worth 3 targets).
MEASURED (`results/s18/_potfull_test`, craft 4 of g2_plan slice 3, T 1620 d, nice 5, 1 core while p1_ration held 9):
same plan as g2_plan's own job (1546 states, depth 12, 185 candidates), 16 columns, 31/31 settles, 16/16 regrids ok,
135 s (plan 80 s); the dtank frontier is bit-identical to g2_plan's (1 / 2 / 3 / 4 / 5 fb: 3.6 / 21.6 / 30.3 / 38.7 /
67.1 kg) while the emitted potentials are 0-11 (mean 4.5; per depth max 11 / 10 / 9 / 8 / 7) against 0-2 (mean 0.94)
in g2_plan's lookahead columns.  Default-mode regression on the same cached plan (scratch pb_equiv.py): candidates and
frontier of the new code are identical to the installed original (185 / 16, same legs, pot, pot_cont, fuel).

## 8. Files
`tools/s18_planbeam.py` (new), `tools/s18_rhfa.py` (resolve_proposer: `planbeam`; docstring line), this file,
`results/s18/planbeam_test/{s3, s3_v1 (before the regrid step), var/*, root, w150, smoke, _wp2_fake, review (the
reviewer's checks), fix (section 5: crafts 0, 4; `c0_streak_v1` = the hard-stop streak variant, 168 s, 1 column),
fix_root}`.
Existing tests re-run after the driver change: `s18_auction.py selftest` OK, `s18_rhfa.py run ... --proposer
fixture_fake --commit none --stop-T 1620` OK (3 slices, done at T = 1620 d).
