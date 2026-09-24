# Stage 18 BUILD: three parallel work packages for the RHFA solver (2026-09-24 01:00 CST)

Design: docs/stage18_rhfa.md (architecture, gates). Rules and laws: docs/stage17_brief.md. This file is the build
contract: who writes which file, the interfaces between the files, and the acceptance test each package must pass
alone. Nothing here moves a gate.

## 0. The contract everybody builds against (read before writing a line)

**Files.** Each package writes ONLY its own files. Every package imports `tools/s18_common.py` (done, tested) and may
import any existing tool (`run_ialns`, `s14_twinbeam`, `s15b_regrid`, `s15b_close`, `ctoc14.*`) but NEVER edits one:
if a function must change, copy it into your s18 file (the two horizon copies in s18_segbeam are already named).
Never write outside `tools/s18_*.py`, `results/s18/**`, `docs/stage18*`. Never touch `results/CTOC14_Result_*.txt`,
validator logs, `results/s16/best`. Never use `--swap`. Never submit.

| package | files (owner writes, others only import) | depends on |
|:--|:--|:--|
| WP1 beam | `tools/s18_windows.py`, `tools/s18_segbeam.py` | s18_common |
| WP2 fleet | `tools/s18_auction.py`, `tools/s18_rhfa.py` | s18_common; at test time s18_fixture (WP3) or its own stub |
| WP3 eval | `tools/s18_eval.py`, `tools/s18_portfolio.py`, `tools/s18_fixture.py` | s18_common |

The skeletons on disk carry every signature, docstring and CLI. Keep the signatures (the other packages call them);
fill the bodies. `raise NotImplementedError` marks what is yours.

**Units (binding, s18_common docstring).** Seconds of mission time in every state, column, job and npz; days ONLY on
the CLI (`--H 540 --L 360 --stop-T 1461`), in log lines, in JSON keys ending `_d`, and in `results/s18/windows.json`
(`t_d`). `C.d2s / C.s2d` convert. **Target ids** are 1-based ints (131, 144 never appear; `C.TARGETS` = the 298);
the `id - 1` conversion lives only next to `eph.ast_states_at`. **Craft** are 0-based ints; a root column has
`craft = None`.

**The regrid law (law 8).** A tank is honest only after `C.regrid_settle(ip)` (= `s15b_regrid.regrid` to regular
20 d bins + `run_ialns.settle`) with `C.honest_check` reporting 0 irregular nodes and `max_window_cap <= 1.02`.
Column tanks are planner-twin tanks (irregular nodes: junction impulse 1 d after a flyby, LinLeg profile every
10 d) and are used ONLY to rank columns. The commit step (WP2) regrids every adopted twin and stores the honest tank;
the state, `fleet.json`, every log line and every eval number use honest tanks; `--commit none` exists for the fake
fixture only and flags every tank `honest=False`.

**Schemas.** `C.COLUMN_SCHEMA`, `C.STATE_SCHEMA` (docstrings in s18_common), the JOB / RESULT schemas in the
s18_segbeam docstring, the CHOICE schema in s18_auction. `C.new_column`, `C.check_column(c, T, H)`, `C.new_state`,
`C.check_state`, `C.save_state` (atomic, validates) are the only way to make or write them.

**Process model.** The driver (WP2) owns ONE `multiprocessing.get_context('fork').Pool(nproc, maxtasksperchild=1)`,
created after `run_ialns.eph()` has been loaded in the parent. It maps `<proposer>.run_job(job)` over the slice's
JOB dicts (`imap_unordered`) and later `s18_rhfa.w_commit` over the chosen columns. A job runs in ONE process,
single-threaded (BLAS threads = 1, `TB._pool_map(fn, jobs, 1)` inside the beam, no nested pools). The root beam is
split into `n_root = max(1, nproc - n_launched)` jobs over contiguous launch sub-windows so the first slices also use
every core. A job is idempotent (its `result.json` is returned unchanged), so a killed slice re-runs only the
unfinished jobs and a restart is bit-for-bit the same except where a settle hit its SIGALRM timeout.

**Timings measured 09-24 00:50 on this machine** (`results/s18/measure_primitives.json`, `measure_deep.json`, one
core, nice 5, load average 4): see section 4. Use them for budgets; mark your own numbers MEASURED with the file.

**Run conventions.** `nice -n 5 ~/.venvs/astro313/bin/python`, `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=VECLIB_MAXIMUM_THREADS=1`
(every s18 file sets them at import). Anything over 5 min: `nohup ... > results/s18/<x>/nohup.out 2>&1 &` with
on-disk state, and poll the log. Round-1 budget per package: <= 2 cores.

## 1. WP1: windows + segment beam (`s18_windows.py`, `s18_segbeam.py`) -- the hard one

### 1.1 s18_windows.py (30 min)
`C.build_windows` already produces `results/s18/windows.json` (95 pin / 123 med / 80 fil, 19 fallback targets whose
windows are the events < 0.08 AU because they have none < 0.05 AU: ids 1, 31, 39, 50, 66, 70, 74, 75, 96, 107, 118,
140, 149, 179, 186, 269, 270, 280, 295). Implement `build` (wrapper), `report`, `check` (class sizes == 95/123/80,
every target >= 1 window, fallback ids exactly those 19, epochs inside the mission) and `recompute_events` (the
0.5 d torus-distance scan of `results/s17/physics/events.py`, so the table can be rebuilt without results/s17).

### 1.2 s18_segbeam.py (the work)
Reuse from `s14_twinbeam` as is: `params`, `root_state`, `search_state`, `seeds`, `child_problem`, `settle_child`
(takes `(par, ast, t, iters, timeout, lin)`), `lin_price`, `twin_score`, `pack`, `save_ckpt`, `best_of`. Copy with
the horizon: `epoch_candidates` -> `epoch_candidates_h(par, remaining, t_lo, t_hi, K, sep)` (mask the 15-400 d and
20-600 d tof grids so `t0 + tof` lies in `(t_lo, t_hi]` BEFORE the Lambert / LinLeg batches) and `_price_parent` ->
`price_parent(job)` with a prize vector (score = `w_t t/T + w_fuel F/1400 - (prize(parent) + prize[ast])`, the
`s15_twintail` convention; `pscore` for settled states). Beam mechanics = `s15_twintail.twin_beam` (ranking,
(visited, 5 d) dedupe, settle in rank order until `B` succeed or `tries` spent, checkpoint every level) with:

1. **Start states.** `prefix_state(route_npz)`: `C.load_twin -> C.twin_of -> C.pack_state` (nreg from
   `C.nreg_of`; the prefix is already regridded+settled by the driver; assert miss <= 150 km, never settle here).
   Root: `root_states(grid, remaining, t_first_max=T+H+L, n_keep, ...)` = `search.roots` over the job's launch
   sub-grid (`excluded` = everything not in `remaining`), keep first flyby <= T+H+L, sort by score, dedupe
   `(first target, 20 d launch bin)`, keep `n_keep`, `TB.root_state` each (exact zero-thrust Lambert arc).
   MEASURED: `search.roots` over 54 launch dates (10 d step over (0, 540]) -> 13189 states in 14.2 s, so a root job
   over 1/7 of the window costs ~2 s + `n_keep` root_states at ~0.05 s each.
2. **Horizon.** Children only at epochs in `(T, T+H+L]`; a beam state with `t_end > T+H+L - 15 d` is not expanded
   (end kind `horizon`). Track `id`, `parent`, `level` on every settled state (extend the pack dict; it is pickled in
   the checkpoint anyway).
3. **Columns** (`columns_of`): settled states with last flyby <= T+H and >= 1 new target; `potential` = max over
   descendants of extra flybys in `(T+H, T+H+L]`; dedupe on the new-target set (lightest); twins to
   `run_dir/<out>/col_<k>.npz` (`C.save_twin`); `C.check_column(c, T, H)` must pass for every column; `counts` from
   `C.count_classes(W, new)`. The prefix itself is not a column.
4. **Diverse root front.** `diverse_select(states, B, key = launch epoch // root_bin, score)`: round-robin over
   launch bins, best first inside a bin. Craft beams keep the plain best `B`.
5. **Budget.** `wall_budget` per job (default 2400 s): stop at the level that exceeds it, still emit columns.
   `max_levels` 14. Every settled state counts (not only the kept `B`): the auction wants many small columns.
6. **Prizes.** `job['prize']` maps id -> weight (default 1 on `remaining`); the driver passes `1 + edf` so the beam
   and the auction want the same targets. Excluded (not remaining) targets have prize 0 AND are in the exclusion
   mask (never expanded), exactly as `s15_twintail` does with `excl_mask`.
7. **run_job** is the contract: dict in, RESULT dict out, `result.json` on disk, idempotent, exceptions caught into
   `ok False / error`. `make_job` builds the JOB for the driver (WP2 calls it; keep the signature).

Budget per craft-slice, ESTIMATE from the measured pieces (section 4): pricing 8 parents x <= 250 candidates
(`price_cap`) x 15-74 ms (depth 4-33) = 30-150 s per level; 16 settles at 0.3-2 s = 5-30 s per level; 6-8 levels per
slice (900 d at ~130 d per flyby) -> 5-25 min per job single-threaded. With 8 jobs in parallel a slice is bounded by
the deepest craft, i.e. ~25 min late in the mission; 10 slices -> a full N = 8 run is ESTIMATED 2-4 h (the design's
14 h figure predates these measurements). If a level's pricing still exceeds 180 s, cut `k_epochs` to 2 for parents
deeper than 25 (log it).

### 1.3 WP1 acceptance test (< 10 min, 1 core; both commands print MEASURED wall)
```
nice -n 5 ~/.venvs/astro313/bin/python tools/s18_segbeam.py test-craft results/s18/_wp1_craft --T 540 --H 540 --L 360 --B 3 --tries 4 --levels 3
nice -n 5 ~/.venvs/astro313/bin/python tools/s18_segbeam.py test-root  results/s18/_wp1_root  --T 0 --H 540 --L 360 --step 30 --roots 10 --B 3 --tries 4 --levels 3
```
Pass: test-craft (r1 truncated at 540 d = 4 flybys, tank 620.8 kg, MEASURED 0.1 s) yields >= 1 column with
`n_new >= 1`, every column npz reloads with `C.miss_of <= 150`, epochs in (540, 1080] d, at least one column with
`potential > 0` at 3 levels (the third level reaches ~1000-1300 d), columns deduped, second `run_job` call returns
the cached result in < 1 s. test-root yields >= 1 root column with `craft None`, `t_launch` in (0, 540] d, first
epoch <= 540 d, and >= 2 distinct 60 d launch bins among the columns. Expected wall: craft 1-3 min (settles at depth
4-7 take 0.2-3 s), root 2-4 min (roots 14 s + level settles). Report the measured numbers and the column tables.

## 2. WP2: auction + driver (`s18_auction.py`, `s18_rhfa.py`)

### 2.1 s18_auction.py (2 h)
Pre-filter (quota / fuel / stale / horizon / empty, each dropped column returned with its reason), objective
`sum_{t in new} (1 + edf(t)) + gamma potential - lam dtank/100` (`C.edf(W, t, T+H, beta)`), constraints: one
column per launched craft, `sum root x <= n_unlaunched`, one column per target (rows only for targets in >= 2 kept
columns), `scipy.optimize.milp` (scipy 1.17.1, HiGHS; < 2000 x 300, seconds), `time_limit` = tlim; accept a
time-limited incumbent and flag `status`. Root columns chosen are assigned to unlaunched craft in index order. Empty
input -> `status 'empty'`, never an exception. `selftest` is the hand-made instance in the docstring with a known
optimum (set `beta = 0` to kill edf ties): `chosen == {'0': 'c0b', '2': 'r1'}`, `dropped == {c1q: quota, c0f: fuel}`.

### 2.2 s18_rhfa.py (4 h)
The three phases with their on-disk files (`slices/s<kk>/jobs/*.json`, `<job>/result.json`, `columns.json`,
`auction.json`, `commit.json`), the resume rule (skip every phase whose file exists; CLI values other than
`--nproc --stop-T --tlim` are frozen in `state['params']` on creation), the commit rule (`C.regrid_settle` in the
pool worker `w_commit`; accept iff `ok` and honest tank <= `m0max + m0slack`; write `craft_<i>.npz` BEFORE
`state.json`; rejected columns' targets stay in `remaining`), the stop rule, `finalize` (`fleet/route_<i>.npz` +
`fleet/fleet.json` with `covered, misses, sumJi, fuel, honest, routes{tank, J_i, flybys, t_launch_d, t_end_d,
counts, cadence_d}`), `status`, and the one-line-per-slice log. The proposer is chosen by name
(`--proposer segbeam | fixture_fake | fixture_replay`) and called as `module.run_job(job)`; the driver never inspects
columns beyond `C.check_column`. `n_root = max(1, nproc - n_launched)` root jobs over contiguous sub-windows of
`(T, T+H]` (step `grid_step` days), `roots` per job = `max(10, ceil(roots / n_root))`.

Fuel and quota facts to keep in mind: 8 x 10 fillers = 80 = the whole filler class, 8 x 16 = 128 >= 123 med, so the
default quotas are exactly tight; `m0max 1120` per craft x 8 = 4160 kg tank > the 3764 kg fuel bar only in sum, the
portfolio (WP3) is what enforces the bar.

### 2.3 WP2 acceptance test (< 10 min, <= 2 cores)
```
nice -n 5 ~/.venvs/astro313/bin/python tools/s18_auction.py selftest                                   # < 5 s
nice -n 5 ~/.venvs/astro313/bin/python tools/s18_rhfa.py run results/s18/_wp2_fake --N 3 --nproc 2 --proposer fixture_fake --commit none --stop-T 1620
nice -n 5 ~/.venvs/astro313/bin/python tools/s18_rhfa.py run results/s18/_wp2_replay --N 2 --nproc 2 --proposer fixture_replay --commit regrid --stop-T 1080
```
Pass (fake, 3 slices, seconds): every slice picks <= 1 column per craft, no target twice (`C.check_state` holds
after every commit), quota / fuel violators appear in `auction.json.dropped` with reasons, `state.json` says
`done` with `T = 1620 d`, and the resume test: delete `slices/s02/auction.json` and `commit.json`, re-run the same
command -> identical `state.json` (compare `remaining`, `taken`, `craft[*].targets`). Kill test: start the replay
run, `kill -TERM` it during slice 1's propose, restart -> the run continues from the cached job results and ends
with the same `state.json` as an uninterrupted run. Pass (replay, 2 slices, both craft launched from r1 / r2
truncations): every committed twin passes `C.honest_check` (0 irregular, window cap <= 1.02) and `fleet/fleet.json`
lists honest tanks (r1 at 1080 d ~ 9-10 flybys; expected wall 2-5 min: truncations 0.1-2 s, regrid_settle 0.6-5 s).
If WP3's fixture is not there yet, a 40-line stub `s18_fixture.run_job` returning fake columns is acceptable for the
fake test only; the replay test waits for WP3.

## 3. WP3: eval + portfolio + fixture (`s18_eval.py`, `s18_portfolio.py`, `s18_fixture.py`)

### 3.1 s18_fixture.py first (2 h; WP2 is blocked on it)
`run_job(job)` with the JOB / RESULT contract of s18_segbeam. `fake`: no twins, deterministic in `seed + slice +
craft`, deliberate conflicts between craft (shared targets), one quota violator (11 fillers) and one fuel violator
(tank 1300) per craft job, `npz = None`. `replay`: real twins = `C.truncate_twin` of `results/s16/best/fleet/
route_r<i+1>.npz` at `T+H` and at the flyby before; `potential` = the route's flybys in `(T+H, T+H+L]`; new targets
not in `remaining` -> truncate earlier; root jobs replay the routes of the unlaunched indices (r9 launches at day
944, so it produces nothing before slice 1 at H = 540: that is a feature, the driver must survive an empty root
result). MEASURED: `truncate_twin(r1, 540 d)` 0.1 s, `regrid_settle(r1 full)` 0.7 s.

### 3.2 s18_eval.py (2 h)
`load_fleet` (run dir or bare fleet dir), `route_report` (one integration per route: `C.twin_summary`,
`C.honest_check`, `C.count_classes`, cadence), `flybys_by(routes, 1461)` (distinct targets + fuel spent by then),
`closure_estimate` (`s15b_close._price_host` + greedy; ESTIMATE, labelled), `gates` with the fixed numbers in
`GATES`, `table`. Reference: `s18_eval.py results/s16/best/fleet` must print covered 298, sumJi 11.7227, fuel 2879
(MEASURED by `C.fuel_of` of the honest tanks: 944.4+908.9+986.3+913.1+958.1+901.5+1017.1+927.2+902.8 - 9x600),
every route honest with the max window caps of `results/s16/best/honest_check.txt` (0.978, 0.933, 0.988, 0.961,
0.995, 0.961, 0.980, 0.997, 0.998), flybys_by_1461d = 70 (the number quoted in the design for s16a).

### 3.3 s18_portfolio.py (1 h)
`command`, `run_one` (subprocess, waits, evaluates, appends the row), `finished`, `run` (skip finished, resume
unfinished), `table`. `--dry` prints the nine commands and writes nothing.

### 3.4 WP3 acceptance test (< 10 min, <= 2 cores)
```
nice -n 5 ~/.venvs/astro313/bin/python tools/s18_fixture.py selftest                     # < 3 min
nice -n 5 ~/.venvs/astro313/bin/python tools/s18_eval.py results/s16/best/fleet --no-g2   # < 1 min; numbers of 3.2
nice -n 5 ~/.venvs/astro313/bin/python tools/s18_eval.py results/s16/best/fleet --nproc 2  # G2 estimate; < 5 min
nice -n 5 ~/.venvs/astro313/bin/python tools/s18_portfolio.py write-default results/s18/portfolio.json
nice -n 5 ~/.venvs/astro313/bin/python tools/s18_portfolio.py run --dry                   # prints 9 commands
```
Pass: fixture selftest asserts the schemas (fake: violators present; replay: twins reload with miss <= 150 and
epochs in (T, T+H]); eval reproduces the reference numbers above (G0b/G1 print `N/A` without a state.json, G2
prints the estimate with its label); portfolio dry run prints nine `nice -n 5 ... s18_rhfa.py run results/s18/p0x_...`
commands and creates no run directory.

## 4. Measured primitive costs (09-24 00:50, one core, nice 5; files under results/s18/)
| primitive | MEASURED | file |
|:--|:--|:--|
| `regrid_settle` of an honest 21-fb route (r9) | 0.6 s, tank 902.83 -> 902.83, miss 31.6 km, ok | measure_primitives.json |
| `regrid_settle` of an honest 42-fb route (r1) | 0.7 s, tank 944.38 -> 944.38, miss 35.5 km, ok | measure_primitives.json |
| `truncate_twin(r1, 540 d)` | 4 flybys kept, tank 620.8, miss 3.6 km, 0.1 s | measure_primitives.json |
| `epoch_candidates` from that prefix over 294 targets | 872 candidates, 0.4 s | measure_primitives.json |
| `lin_price` (4-fb prefix) | 14.8 ms per call | measure_primitives.json |
| `settle_child` (4-fb prefix, lintwin seed) | 0.2-0.4 s, +1.9 / +3.8 kg | measure_primitives.json |
| `search.roots` 54 launch dates (0-540 d, 10 d) | 13189 root states, 14.2 s | measure_primitives.json |
| `truncate_twin(r1, 2000 d)` / `(r1, 4000 d)` | 16 fb, 702.3 kg, 0.4 s / 33 fb, 867.1 kg, miss 55 km, 2.2 s | measure_deep.json |
| `epoch_candidates` from the 16-fb / 33-fb prefix | 841 / 772 candidates, 0.4 s each | measure_deep.json |
| `lin_price` at depth 16 / 33 | 40.9 / 73.9 ms per call | measure_deep.json |
| `settle_child` at depth 16 / 33 (lintwin seed, 3 each) | 0.3-0.4 s / 1.1-1.9 s, all settled, +4.5 to +21 kg | measure_deep.json |
| `regrid_settle` of those children (honest - planner) | 0.2-0.5 s; -0.01 to -0.31 kg; 0 irregular nodes (the lintwin seed adds no irregular node; `lam` / `lin` seeds do) | measure_deep.json |

Consequences for WP1 (from the table): (a) at depth 30 the cost is PRICING (8 parents x 750 candidates x 74 ms = 7 min
per level), not settling (16 x 2 s); cap `lin_price` calls per parent at the `price_cap` = 250 cheapest single-leg
candidates inside the horizon (log the cut) -> ~2.5 min per level, ~20 min per job at 7 levels; (b) the horizon did
not cut anything in these measurements because `epoch_candidates` only looks 600 d past `t_end` and `t_end` was
before T; it binds once the beam approaches T+H+L, which is exactly where it must stop expanding; (c) the regrid law
is nearly free for lintwin children and mandatory for the `lam` / `lin` seeded ones (irregular nodes).

## 5. Integration order and the gates (the integrator = the architect, after the three packages report)
1. WP3 fixture + WP2 driver: fake run (seconds), replay run (minutes) -> the driver's resume and regrid law proven.
2. WP1 columns into the driver: **G0 smoke** `run results/s18/g0 --N 2 --nproc 2 --stop-T 1080` (2 slices; kill
   mid-slice 2 and restart; same result), expected 20-40 min.
3. **G0b** `run results/s18/g0b --N 8 --nproc 8 --stop-T 1461` (3 slices, ESTIMATE 30-60 min): >= 72 distinct
   flybys by day 1461 and mean fuel <= 160 kg per launched craft; KILL at <= 60.
4. **G1** the same run continued to the end: `--stop-T` is a runtime knob, so `run results/s18/g0b --nproc 8`
   (no `--stop-T`) resumes the g0b run past day 1461 to the end of the mission (ESTIMATE 2-4 h for the full run):
   covered >= 268 at <= 3900 kg -> portfolio; 250-267 -> one parameter round; < 250 KILL.
5. **G2** portfolio best >= 292 with closure estimate <= 11.50 sum J_i -> `s15b_close`, `s15b_relocate`, regrid,
   `s16_select` (N <= 8, K = 298), export by the `results/s16/export_s16a.sh` pattern (3-4 h). Nothing is submitted
   by an agent.
