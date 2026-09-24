# CTOC14 Problem A — multi-spacecraft low-thrust NEA flyby tour solver

China Trajectory Optimization Competition 14 (第十四届全国空间轨道设计竞赛), Problem A: fly by as many of 300
catalogued near-Earth asteroids as possible with a fleet of low-thrust spacecraft inside a 15-year window, minimising
J = Σ(1 + x + x²) + misses. Team result on the board: 9 craft, 298/298 targets, raw J 13.747 (shown 13.512).

**Live 3D viewer (GitHub Pages):** <https://exoplanet5.github.io/ctoc14/results/viz/> — Sun, Earth–Moon, all 300
targets and nine fleet results animated over the mission. Local: `cd results/viz && python -m http.server 8765`.
Results catalogue: [`docs/results_catalog.md`](docs/results_catalog.md); solver history: `docs/stage*.md`,
[`docs/stage18_rhfa.md`](docs/stage18_rhfa.md) (the final concurrent fleet solver).

Python: `~/.venvs/astro313/bin/python` (numpy, scipy). Run everything from this directory.

## Layout
- `MEA.txt`, `CTOC14_problem.pdf` — official inputs.
- `docs/` — problem digest (`problem_summary.md`), machine-readable rules (`rules_spec.json`), rules review & trap list
  (`rules_review.md`, `thrust_interpolation_rule.md`), target identification (`nea_identification.md`), population analysis
  (`population_analysis.md`), solver design (`solver_design.md`), final report (`report.md`).
- `ctoc14/` — solver package: `constants`, `kepler` (ephemeris + two-body propagation), `lambert`, `thrust` (validator thrust
  model), `validator` (replica of the official checker), `submission` (writer), `search` (beam search planner),
  `lowthrust` (plan → continuous-thrust conversion with re-planning).
- `tools/` — CLIs; `tests/` — checks against the PDF example; `analysis/` — population / ballistic studies; `data/` — catalogues;
  `results/` — tours, fragments, submissions, logs.

## Pipeline
```bash
P=~/.venvs/astro313/bin/python
$P tests/test_pdf_example.py && $P tests/test_thrust_example.py && $P tests/test_submission_roundtrip.py   # sanity
$P tools/identify_nea.py && $P tools/identify_nea_mpc.py          # (optional) JPL/MPC identification of MEA.txt
$P tools/run_fleet_search.py results/fleet --beam 300 --wfuel 0.7 --max-sc 12 --mmargin 100   # sequential full-tank tours
$P tools/run_tail_search.py results/fleet/tail results/fleet/tail_targets.json 700 850 1000     # cheap mop-up chains
$P tools/convert_fleet.py results/fleet results/CTOC14_Result_TEAM.txt --nproc 8                # low-thrust conversion + validation
$P -m ctoc14.validator results/CTOC14_Result_TEAM.txt                                           # independent re-check
```
Single-spacecraft experiments: `tools/run_search.py`, `tools/convert_one.py`, `tools/combine_submission.py`.

## Resuming the optimisation (state on 2026-09-07 17:00)
- Best validated file so far: `results/CTOC14_Result_TEAM.txt` (= v4: 12 spacecraft, 298 covered, J = 32.259).
- A re-conversion of the same 12 tours with the improved converter is/was running:
  `python tools/convert_fleet.py results/fleet_b300 results/CTOC14_Result_v2.txt --nproc 6 --tag q`
  (log: `results/fleet_b300/convert_fleet_v2.log`; on completion it writes `results/CTOC14_Result_v2.txt` + `_report.json`, validated).
- If some targets are still missed afterwards: compute the leftovers, run `tools/run_tail_search.py` on them, convert the best chain with
  `tools/convert_one.py`, and merge with `tools/combine_submission.py` (see `tools/run_campaign.py` for the automated version).
- Converter defaults now: adaptive_iter, widen_on_fail, expensive_drop on; smart_replan, fuel_guard off (docs/conversion_losses.md).
- 2026-09-09 feasibility study for J < 20: `docs/j20_research.md` (verdict: not reachable; floor ~21–22; realistic 24–27).
  Experiments `tools/exp_j20.py` (tank sweep, fuel-free ceiling, beam/eta, small-tank fleets → `results/exp_j20/`),
  `tools/exp_reach.py` (continuous-thrust reach vs Lambert rule), `tools/select_tours_milp.py` (set-packing tour selection, HiGHS).
  Linearised continuous-thrust leg model (ctoc14/linleg.py, Params.lin_tofs, default off) tested: no gain in flyby count or fuel.
  2026-09-10: patient-swarm campaign (ctoc14/jointsearch.py, tools/run_joint.py, improve_joint.py, select_frags_milp.py, docs/j20_research.md §6):
  planned J 24.7 but converted 34.98; MILP over all converted fragments -> results/CTOC14_Result_milp1.txt J 32.111.
  2026-09-15: phasing campaign (docs/phasing_campaign.md): patient route pools (tools/run_pool_campaign.py, select_routes_milp.py with
  column generation, --tof-refine) + fragment MILP over 134 fragments -> results/CTOC14_Result_milp3.txt = TEAM, J 30.789 (best so far).
  2026-09-16: column generation + neighbourhood search campaign (docs/colgen_campaign.md): ctoc14/colgen.py, tools/run_colgen.py,
  tools/run_lns.py (tail/full re-plans with covered targets as waypoints, spawn, pair, dissolve), tools/convert_pool.py (convert only
  new routes into results/colgen/frags), tools/tighten_frags.py (re-convert at 600 + used + margin), select_frags_milp.py --dry.
  TEAM = results/CTOC14_Result_cgF.txt, J 18.4917 (13 craft, 298 covered, both validators). Previous: cgX 18.5989, cg5 19.5511, milp3 30.789.
  2026-09-16/17: whole-trajectory minimum-propellant SCP + fast impulsive fleet search (docs/globalopt_campaign.md):
  ctoc14/globalopt.py (exact discrete sensitivities, B-plane IRLS, mass-free mode, exact tank scaling, aim-point homotopy
  insertion), ctoc14/impulsive.py (impulses every 20 d + Kepler arcs, same interface, ~500x faster), tools/run_globalopt.py,
  tools/run_gfleet.py (init/dedupe/dissolve/write), tools/run_ialns.py (import/search: relocate + dissolve/export).
  2026-09-17 pm: block moves and the reduced cost model (docs/globalopt_campaign.md sections 8.1-8.6):
  ctoc14/globalopt.insert_block (adaptive backtracking aim-point homotopy, one or many targets at once),
  run_ialns.py `search --lns` (ruin-and-recreate dissolve with block insertion + ejection) and `retime` (every other
  pass of each asteroid from the event catalogue), ctoc14/phasemodel.py (dv = 0.82 [0.028 TV(drift deg/yr) +
  29.8 TV(i_vec rad)], accurate to 2 % on all ten routes). Fleet-size economics: removing a craft pays only if the
  fleet's total Delta-v grows < 23.5 %, while a marginal target move costs 1.6-9 km/s vs a 0.49 km/s average, so
  10 -> 9 by redistribution is ~4x over budget; ctoc14/phasedp.py (phase-plane DP) is NOT a valid generator.
  TEAM = results/CTOC14_Result_t10d.txt, J 14.366055 (10 craft, 298 covered, both validators). Previous: t10c 14.395002, t10b 14.495927, t10a 14.845435, r14 14.950134, r15 15.037440, d11 15.809736 (12 craft),
  dd2 16.720044, go1 16.823063, cgF 18.4917.
  2026-09-17 evening: docs/strategy_stage2.md. The leaderboard's N column is craft + misses: the top five are 8-craft
  fleets at 0.48-0.54 km/s per flyby, i.e. OUR efficiency (0.489); the gap is the partition (10 vs 8 craft). Stage 2 =
  column generation at N = 8 with columns costed by the impulsive twin; 9 craft is not worth the time coefficient.
  Evening: colgen, dive and skeleton growth all stall at ~10 craft (campaign doc 9.1-9.5); plain beams on ~50-target
  neighbourhoods give 33-36-flyby hard-rich routes -> step-by-step plan in docs/stage2_instructions.md.

**2026-09-18.** Two free levers found and measured: the beam's planner mass `m0` caps the craft
(m0 1000 -> 805 kg -> ~34 flybys, which is why every method equilibrated at ten craft; m0 1600 -> one route flies 46
targets at a real 1004-1067 kg), and `vinf_cap` was 2.0 km/s although the contest allows 4.0 free (-63 kg of tank per
craft). t10d was submitted on 09-17 and is banked at 13.81 shown. The open problem is the PARTITION: a beam confined
to a 48-target group captures only ~36, so seven balanced routes reach ~250 of 298. New tools: `tools/greedy_cover.py`,
`dual_greedy.py`, `kroutes.py`, `partition_ls.py`, `partition_mip.py`. Read `docs/stage3_summary.md` first.

**2026-09-18 (evening):** next-stage plan for an eight-craft fleet (raw J < 13) as a gated search tree:
`docs/stage4_n8_search_tree.md`. Instructions only, nothing run yet.
