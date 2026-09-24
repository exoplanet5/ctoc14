# Stage 17 brief: GTOC-style brainstorm for a lower J (2026-09-23 21:30 CST)

Shared briefing for the stage-17 brainstorm agents. Read this first, then the memory files in section 7.

## 1. Problem in five lines (full text: CTOC14_problem.txt, docs/problem_summary.md)
- 300 NEAs (MEA.txt, Keplerian). 131 and 144 are unreachable, so 298 is the ceiling and J always carries +2.
- Low-thrust craft: Tmax 0.5 N, Isp 4000 s, dry mass 600 kg, m0 <= 2000 kg, launch v_inf <= 4 km/s.
  Mission window 2030-01-01 + 15 yr. A flyby means passing within 1000 km. No gravity assists.
- J = sum_i (1 + x_i + x_i^2) + N_miss, where x_i = (m0_i - 600)/1400. Each craft costs 1 J before fuel.
- Shown score = raw J x k(t_submit), with k = 1 - 0.1 (2026-09-28 12:00 - t)/28 d, taken from the LAST valid submission.
- The user submits by hand; agents never submit.

## 2. Where we stand
- Board, snapshot 09-23 ~20:44 (shown J, craft = board N - 2):
  | rank | team | shown J | craft |
  |:--|:--|--:|--:|
  | 1 | HIT | 11.936 | 8 |
  | 2 | THU | 11.973 | 7 |
  | 3 | NUAA | 12.296 | 8 |
  | 4 | eStar | 12.996 | 8 |
  | 5 | 银河远征 | 13.026 | 8 |
  | 6 | 小熊熊 | 13.286 | 9 |
  | 7 | 砺剑尖兵 | 13.406 | 10 |
  | 8 | 未来之星 | 13.493 | 9 |
  | 9 | 地月探索者 | 13.503 | 8 |
  | **10** | **us** | **13.5116** | **9** (s16E, submitted 17:05) |
- Our best file: results/CTOC14_Result_s16a.txt, raw 13.739947, 9 craft, 298 covered, VALID. Twin fleet:
  results/s16/best/fleet (route_*.npz, twin sum J_i 11.7227). The exact conversion costs about +17 kg.
- s16a anatomy (m0 kg / flybys / kg per flyby):
  | craft | m0 | flybys | kg/fb | launch day |
  |:--|--:|--:|--:|--:|
  | sc1 | 948 | 42 | 8.3 | 41 |
  | sc2 | 911 | 41 | 7.6 | 0 |
  | sc3 | 989 | 41 | 9.5 | 88 |
  | sc4 | 916 | 38 | 8.3 | 0 |
  | sc5 | 960 | 31 | 11.6 | 80 |
  | sc6 | 902 | 30 | 10.1 | 0 |
  | sc7 | 1019 | 30 | 14.0 | 15 |
  | sc8 | 928 | 24 | 13.7 | 19 |
  | sc9 | 904 | 21 | 14.5 | 944 |

  Total fuel 3076 kg (10.3 kg/fb). The four deep routes run at 8 kg/fb. The five tails run at 10-14.5 kg/fb.

### Economics against the clock
A new file must beat 13.511581 / k(t) or it LOWERS our shown score:

| submit at | raw J needed |
|:--|--:|
| 09-24 12:00 | 13.708 |
| 09-25 12:00 | 13.658 |
| 09-26 12:00 | 13.609 |
| 09-27 12:00 | 13.560 |

Each day costs about 0.049 J. Export plus validation of a finished twin fleet takes about 3-4 h, so the last useful fleet must exist by roughly 09-27 06:00. That leaves about 3.3 days of 10-core compute.

What a result is worth:
- **8 craft** covering 298 beats the 09-27 bar at up to 3858 kg total fuel (mean m0 1082 kg, 12.9 kg/fb). Our deep routes do 8 kg/fb, so there is roughly 780 kg of slack for "forced" targets.
  - J < 13.0 needs <= 3391 kg (m0 1024).
  - Reaching rank 6 needs raw < 13.33 at 09-27, i.e. <= 3683 kg.
- **9 craft** needs <= 2912 kg (-164 kg) just to break even on 09-27. It needs 2679 kg (t10d efficiency, 9.0 kg/fb) to reach raw 13.33.
- **Conclusion:** 9-craft polishing (0.01 J) is worthless now. Only a structural change pays: 8 craft, or a 9-craft redesign worth -200 kg or more.

## 3. Measured laws (each re-measured by several independent machineries; do not re-derive)
1. **Forced-target premium.** A target forced onto a finished route costs 0.115-0.148 J (100-140 kg). One the route is designed around costs 0.032-0.042 J, about 3x cheaper. Measured by prizes, assigned pools, waypoints, dissolution, ACO realisation and linchpin.
2. **Disjointness, not route quality, binds.**
   - 76 of 759 pool routes individually beat 0.0369 J/fb, yet the largest mutually disjoint set among them is 4 routes covering 157.
   - Free optima want the same targets: 8 free fills gave 206 easy slots but only 132 distinct, and 99 targets were wanted by no route.
3. **Sequential construction starves.** Late routes get scattered leftovers (our sc5-sc9). Sequential 8-craft attempts end at 228-249 covered. The coverage ceilings are:
   - N <= 8: 282 over the pool of about 2650 honest routes (MILP-proven).
   - N <= 9: 290 without the closer.
4. **Insertion into deep finished routes:**
   - Cost is kg ~ d^1.13 depth^0.91 margin^-0.74.
   - P(settle) falls 0.74, 0.60, 0.31, 0.21 by depth band.
   - Near-miss (< 0.06 AU) insertion has a median of 14 kg when it converges.
5. **The planner's leg model is a strict subset of what the twin flies.** A single-leg Lambert or LinLeg with a fixed departure state cannot re-fly 16% of its own routes' legs. Linearised whole-prefix twin pricing (tools/s14_twinbeam.lin_price) admits 286/288 legs but is slower and pricier (12.5 kg/fb).
6. **Impulsive Lambert cost is singular in flyby epoch.** A 37-flyby route that settles at 16.6 km/s prices at 29.8 km/s as an exact-epoch Lambert chain, and 108 km/s at 2 d jitter. Transfer databases on an epoch grid are therefore dead. The state must be the craft ORBIT.
7. **Candidate density sets depth.** Depth is set by the number of targets inside the ~0.035 AU affordable radius of the craft's track, not by the base orbit. Drifting carriers rotate the track but do not lengthen it.
8. **Honest thrust.** Planner-born twins must be regridded to regular 20 d bins before trusting tanks (tools/s15b_regrid.py). Otherwise the per-impulse cap over-credits thrust by 1.1-1.6x.

## 4. Dead ends (do not propose again without a NEW mechanism that addresses the cause)
- Selection over the existing pool (exact set cover, LP, CG; N <= 8 max 282).
- Pool scaling.
- Scarcity or hardness prizes: they trade depth for hard density.
- Stepping stones; wider waypoint windows.
- Skeleton-then-fill: 241.
- Carrier tiling: 249.
- Lloyd co-design: bases collapse.
- ACO membership-before-trajectory: 51% air. This is forcing.
- MILP facility location with closest-approach costs: R^2 0.01 for assigned targets.
- Transfer and encounter databases (law 6).
- Drifting carrier.
- Re-planning inside assigned pools: depth collapses 24-29 of 37.
- Joint K=2 over a 65 pool: worse than sequential.
- Twin fuel polish of unchanged sets (+-5 kg); it does not survive conversion anyway.
- 10 -> 8 dissolution by insertion (over budget).
- `run_ialns.py relocate_pass --swap`: buggy, never use it.

**Listed as "ONLY UNBUILT" in earlier stages:**
- A joint k-route DP / concurrent construction (all craft advanced together in time over one shared target set, with assignment inside the recursion).
- Leftover-specialist craft designed FIRST.
- An orbit-state (twin-velocity) planner.

## 5. Tool inventory (reuse; do not edit shared files — copy into a new s17_*.py if a change is needed)
- **Twin:** ctoc14/impulsive.py
  - ImpulsiveProblem has launch v_inf + 20 d impulse nodes + Kepler arcs, with flyby epochs free.
  - It provides tank(), integrate(), misses().
  - run_ialns.py provides settle, ipr/ist (npz <-> twin), w_remove, w_cands (approach minima), relocate_pass and export.
- **Planners:**
  - ctoc14/search.py: beam; Lambert 15-400 d + LinLeg 20-600 d legs; win_lo/win_hi waypoints; prizes.
  - tools/s14_twinbeam.py: twin-state beam; lin_price, the whole-prefix linearised insertion price, 5-35 ms.
  - tools/s15_isogen.py: premium-steered deep columns.
  - tools/s15b_prefix.py: keep a settled prefix of k flybys and re-plan the tail with keep-own/trading prizes.
  - tools/s15b_twintail.py, tools/s15_twintail.py.
- **Fleet layer:**
  - tools/s15_select.py and s16_select.py: exact MILP selection with coverage frontier.
  - tools/s15b_close.py: cap-ladder closer.
  - tools/s15b_relocate.py.
  - tools/s16_iter.py: iterated MILP over single-target move columns.
  - tools/s16_swap.py, tools/s16_polish.py.
  - tools/fleet_dedupe.py, tools/deep_guided.py (surrogate-guided deepener), tools/insert_model.py.
- **Carriers:** tools/carrier_dp.py, carrier_grow.py, exp_carrier.py; ctoc14/phasemodel.py.
- **Honest pool:** results/s16/honest/route_<key>.npz + index.jsonl, 2490 regridded or regular routes. Selection: `s16_select.py`.
- **Export:** results/s16/export_s16a.sh pattern (run_ialns.py export -> combine_submission.py -> tools/validator2.py). About 3 h on 9 procs.

## 6. Rules for agents
- Python: `~/.venvs/astro313/bin/python`. Prefix runs with `nice -n 5`. Set OMP/OPENBLAS threads = 1.
- Round-1 budget: at most 2 cores and 45 min of compute per agent. The machine has 10 cores shared by 4 agents.
- Run anything longer than about 5 min as `nohup ... > log 2>&1 &` with on-disk state. Closing the Mac lid freezes processes and stalls agents; resume from disk.
- Write only under `results/s17/<your-lens>/` and `docs/stage17/`. Never modify:
  - results/CTOC14_Result_*.txt or results/validator2_*.log;
  - results/s16/best;
  - any existing tools/*.py or ctoc14/*.py.
- Never submit. Never use `--swap`.
- Honesty: separate MEASURED numbers (with the file that holds them) from estimates. Mark literature recall you are not sure of as [recall, unverified].

## 7. Read these (compressed history)
- Memory: /Users/mickey/.claude/projects/-Users-mickey-solarsystem-ctoc14/memory/
  - ctoc14-craft-count-reframe.md
  - ctoc14-stage15-isogen.md
  - ctoc14-stage15b-closed.md
  - ctoc14-stage7-carrier-model.md
  - ctoc14-stage9-gtoc-architecture.md
  - ctoc14-stage6-plan.md
  - ctoc14-n8-budget.md
- Docs:
  - docs/stage10_architecture.md
  - docs/stage13_solver_plan.md (sections 0, 1, 8)
  - docs/stage9_gtoc_architecture.md
  - docs/methods_survey.md: the GTOC4/5/11/12 literature notes
  - docs/stage6_results.md
  - results/s16/PLAN.txt
