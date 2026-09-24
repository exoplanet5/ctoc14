# Stage 13: route-planner capability and cost map (2026-09-21)

This file answers one question: to plan ONE route freely inside a prescribed target pool, or 2 to 3 routes jointly over a shared pool, what should the fleet builder call, how long does it take, how deep do the routes go, what do they cost, and do they settle?

**How it was measured.** All numbers are new measurements unless marked [cited]. Scripts and raw outputs are in `results/s13/planner_probe/`. Every run used the production leg model of `tools/skel_fill.py` / `tools/greedy_cover.py`:
- search mass m0 1600, v_inf cap 4 km/s, per-leg dv cap 1.2 km/s
- Lambert legs 15-400 d in 5 d steps, plus linear legs 20-600 d in 10 d steps
- lin_drmax 0.15 AU, tof_refine on, m_margin 40 kg
- launch grid 0-800 d in 20 d steps
- w_fuel = CM.w_fuel(480, 1600) = 0.520
- impulsive.LIN_FALLBACK = 2.5

Every run used 1 core, except one timing run on 2 cores. The machine load rose from 5 to 250 (10 cores) during the session because other agents were running. So **CPU seconds are the reliable cost figure**, and wall times after about 20:50 are inflated.

"Pool" means the targets the route is allowed to score. The **coherent pools** are unions of t10d route target sets (`results/s12/cover_clean/fleet`): 08+09 = 50, 07+08+09 = 77, 06-09 = 104, 04-09 = 162. The **random pools** are uniform samples of the 298.

## 0. Headline numbers

- **A free route is deep and cheap.** Over the whole pool, beam 100 gives 46 flybys, planner tank 956 kg, twin 1004 kg (8.8 kg/flyby). The Pareto point at 44 flybys settles at 941 kg (7.75 kg/flyby, beam 300). A second route over the 252 targets left over gives 43 flybys at 975 kg (8.7 kg/flyby).
- **Depth falls steeply with pool size.**
  - Random pools: depth ≈ 4.9 · N^0.39 (100→30, 150→35, 200→39, 250→41, 298→46).
  - Coherent pools do worse (50→22, 77→26, 104→31, 162→33).
  - Pool = t10d route 00's own 37 targets: the planner flies 25 of them.
  - **37 flybys per craft needs about 170 random (about 200 coherent) targets available to that craft.** Across 8 craft that is about 4.6 pools per target, so pools cannot be disjoint.
- **Per-flyby efficiency does NOT degrade in restricted pools.** Twin cost is 6.4-10.0 kg/flyby (median 8.35) in every pool, against 9.28 for t10d as a whole. Restriction costs depth, not efficiency. This agrees with the stage-6 P06 note, but P06 honoured its assigned pool only 19/32 (section 3).
- **Three levers do nothing for in-pool depth:**
  - stepping stones at prize 0 (same in-pool depth, +50% CPU);
  - wider linear legs (lin_drmax 0.25: identical route);
  - beam 100→300 (identical depth; only cheaper at fixed depth: 44 flybys 968→941 kg).
  - dvmax 1.6 adds 2 flybys in the 50-pool, but at 110 kg for those 2.
- **The joint planner (jointsearch) is no better than sequential single beams.**
  - Full pool, K=2: 82 covered (twin sum J_i 2.743) vs 46+43 = 89 (2.712) sequential.
  - Pool 77: 41 vs 44.
  - Pool 50: 33-34 vs 22+8 = 30.
  - At K=3, **63% of the stock beam is duplicate/permuted states** (the dedup key depends on craft order). An order-free key recovers part of this (+3 targets).
- **Settling works.**
  - Single-route beam tours: 37/38 settled (97%). Median 13 s, p90 57 s, max 121 s on 1 core.
  - Twin/planner tank ratio: median 1.018 (0.954-1.134).
  - Stored deep beam-1000 columns (38-40 flybys): 4/4 settled in 7-9 s.
  - Joint tours that contain waits: the stock seed settles 2/7. A wait-aware seed (new, `planner_probe/twin_wait.py`) settles 4/5; the other 2 waited tours were not saved, so the wait-aware seed was not tried on them. At least one of the two seeds settles all 5 saved waited tours.

## 1. The table

The wall time per route column gives CPU seconds on 1 core (≈ wall on an idle machine). 2 cores = about CPU × 0.56, since measured parallel overhead is +12%.

| planner | call | wall time/route | typical depth | kg/flyby (twin) | settle rate | caveats |
|:--|:--|:--|:--|:--|:--|:--|
| `search.beam_search`, free, all 298 | `beam_search(eph, excluded=(), m0=1600, t_launch_grid=np.arange(0,801,20)*DAY, P=BP(...), n_proc=1)` | beam 100: 253 s CPU (4709 expansions × 54 ms). Beam 300: about 760 s CPU (14 095 expansions). **Beam 1000: about 2540 CPU s, i.e. about 42 min on 1 core, about 24 min on 2 cores** (extrapolated; expansions are linear in beam) | 46 (Pareto 38-46) | 7.75-8.8 | 6/6 | Default `n_proc=8`, `verbose=True`: pass both. Beam 300 = same deepest route as beam 100 |
| `beam_search`, **strict pool** N=37-252 | same, with `excluded = ALL - pool` and prize 1 on pool | beam 100: 56-231 s CPU (35-41 ms per expansion; 1573-3600 expansions) | 22-43 (curve in section 2.4) | 6.4-13.3, median 8.3 | 24/25 | Depth is set by pool size and geometry, not by beam width (pool 50: beam 100 = beam 300 = 22) |
| `beam_search`, **soft pool** (non-pool targets at prize 0) | `excluded=[]`, prize 0 outside the pool | beam 100: 155 s CPU at N=77 vs 103 strict (55 ms per expansion) | same in-pool depth as strict (22 / 26); deepest-total state adds 1-4 stepping stones | same as strict | 4/4 | No gain; +50% CPU |
| `tools/skel_fill.py`, skeleton-free (`--drop-wp all`, windowed, `--pools`) | section 3.1 | beam 100: 219 s wall / 158 s CPU (pool 50, pad-prize 0) | 26 total = 22 in pool + 4 pad | 9.4 (843 kg / 26) | 2/2 | Reproduces `beam_search` soft. No strict mode. The Pareto key and the "new" count include pad targets |
| `skel_fill.py` with skeleton, segmented, beam 1000, seg-k 50 [cited, stage-6 G] | `results/n8s6/runG.sh` | 175-234 s wall on 8 procs (about 1400-1900 CPU s) | 35-41 with 8-9 waypoints, full pool | 7.6-12.4 (01: 41 @ 918 = 7.8) | 26/29 (90%) | Needs `route_NN.npz` skeletons. Waypoints are ±45 d windows |
| `skel_fill.py --pools` with skeleton, pad 0.35 [cited, P06] | pools32.json | 645 s on 8 procs (contended) | 32 = **19 of the 32 assigned** + 8 waypoints + 5 pad | 9.9 (917 kg) | 3/3 | Pool honoured only 59% |
| `jointsearch.joint_beam_search`, K=2 (order-free key, wait 60 d) | section 4.1 | pool 50: 281 s (beam 100), 1143 s (beam 300). Pool 77: 444 s. All 298: 715 s (beam 100, 1 core) | per craft 15-22 in pools of 50-77; 39+43 in the full pool | 7.2-12.0 | stock seed 7/12, wait-aware 7/8, at least one seed 10/10 (of the saved tours) | ≤ sequential single beams. Waits break the stock seed. Run_joint.py CLI defaults are the 09-10 swarm settings |
| `joint_beam_search`, K=3, pool 77 | same | 530 s (order-free key) / 629 s (stock key), beam 100 | 13+19+19 = 51 (stock key 48) vs t10d's 77 | planner 690-794 kg tanks | not attempted | Stock key: 63% of the beam are duplicate/permuted states |

## 2. `ctoc14/search.beam_search` (single route)

### 2.1 Call

```python
import numpy as np, run_ialns as RI, ctoc14.impulsive as IM
from greedy_cover import BP, CM, DeepCollect, ALL, _timeout   # BP = Params subclass that drops P.collect when pickled
from ctoc14.search import beam_search, tour_json
from ctoc14.constants import DAY, AU
IM.LIN_FALLBACK = 2.5
prize = np.zeros(300); prize[[t - 1 for t in pool]] = 1.0      # J-unit prize per target
P = BP(beam=300, w_fuel=CM.w_fuel(480, 1600), m_margin=40.0, dv_max=1.2,
       tofs=np.arange(15, 401, 5) * DAY, lin_tofs=np.arange(20, 601, 10) * DAY, lin_drmax=0.15 * AU,
       vinf_cap=4.0, tof_refine=True, max_depth=60)
P.prize = prize; P.collect = DeepCollect(10)                    # harvest every level's beam states with >= 10 flybys
best, beam = beam_search(RI.eph(), excluded=sorted(set(ALL) - set(pool)), m0=1600.0,
                         t_launch_grid=np.arange(0, 801, 20) * DAY, P=P, n_proc=2, verbose=False)
# Pareto: per depth (or per in-pool count), keep the state with minimum s.fuel over list(P.collect) + list(beam)
ip, miss, lag = IM.settle_tour(RI.eph(), tour_json(s), lambda ip: RI.settle(ip, 100))   # wrap in SIGALRM 240 s
ok = ip is not None and miss <= 150.0            # km (the skel_fill / greedy_cover criterion)
tank = ip.tank(); J_i = RI.cost(tank); np.savez(out / 'route_NN.npz', **RI.ist(ip))     # loads with RI.IFleet(out)
```

- **Return value.** `beam_search(eph, excluded=(), m0=2000.0, t_launch_grid=None, P=None, n_proc=8, verbose=True)` returns `(best, final_beam)`. `best` is the deepest state (then the lowest score). `final_beam` is only the last level; a Pareto front needs `P.collect`.
- **Use BP, not plain Params.** If `n_proc > 1` and `P.collect` is set, use `greedy_cover.BP`. Plain `Params` would pickle the whole growing collect list into every worker at every level.
- **Continuing from a state.** `beam_search_from_state(eph, t, r, v, m, visited_ids, m0, t_launch, vinf, P, n_proc, seq)` continues a tour from any state (the tail re-plan of `tools/run_lns.py`). Not timed here.

### 2.2 `excluded` vs prize 0

- **`excluded`** is OR-ed into the root visited mask (together with 131 and 144). Those targets are never a leg end, never a child, never a relay. Expansion is cheaper: 35-41 ms vs 54 ms per expansion.
- **prize 0** leaves the target a legal leg end.
  - The score is `w_t·t/T + w_fuel·fuel/1400 − Σprize`, so a detour through a 0-prize target costs time and fuel with no reward at the same depth.
  - Child pruning (`max_children` 60, `n_per_target` 2) also ranks by `cost − prize`.
  - **Measured:** in-pool depth is identical to strict (N=50: 22 vs 22; N=77: 26 vs 26, the same route and tank). The deepest-total state carries 4 and 1 stepping stones. CPU is +50%.
  - Stage 6 Gate A found the same at beam 300. Stepping stones are not a lever.

### 2.3 Launch grid and m0

- **Launch grid.** Production uses 0-800 d in 20 d steps (41 dates). If `t_launch_grid=None`, the default is 0-3 yr in 20 d steps.
  - Roots are all (date, first target, tof) combinations with |v∞| ≤ `vinf_cap`; any excess over 4 km/s is paid by thrust.
  - Roots are deduplicated to the best score per (first target, 20-d launch bin).
  - Launch dates chosen in the probes: 0-260 d.
  - `tools/run_lns.py` full re-plans use [t_L − 200, t_L + 400] d around the old launch.
- **m0 is the SEARCH mass, not the tank.**
  - It sets the thrust acceleration a = 0.5 N / m, which drives feasibility (dv ≤ 0.6·a·tof).
  - It sets the propellant budget (m − dm ≥ 600 + m_margin).
  - The planner's tank estimate is `CM.tank(fuel, m0) = 602 / (1 − 0.70·fuel/m0)`. The real tank is `ip.tank()` after settling (twin/planner ratio: median 1.018, range 0.954-1.134).
  - m0 1600 is the production value. m0 1000 caps the tank at about 805 kg / 34 flybys (memory note, stage 3).

### 2.4 Depth, cost and time vs pool (beam 100, strict, 1 core)

| pool | N | CPU s | expansions | deepest | settled twin (deepest) | kg/fb | 2nd Pareto point |
|:--|--:|--:|--:|--:|:--|--:|:--|
| t10d route 00 (own set) | 37 | 92 | — | 25 | 760 kg | 6.42 | — |
| t10d 08+09 | 50 | 56 | 1573 | 22 | 789 | 8.58 | 21 @ 764 (7.83) |
| 08+09 minus route A (22) | 28 | 20 | 234 | 8 | 654 | 6.74 | — |
| t10d 07+08+09 | 77 | 103 | 2647 | 26 | 821 | 8.50 | 25 @ 803 (8.13) |
| 07-09 minus route A (26) | 51 | 58 | — | 18 | **settle failed** (miss 3.3e6 km) | — | — |
| t10d 06-09 | 104 | 115 | 3052 | 31 | 858 | 8.31 | 30 @ 837 (7.89) |
| t10d 04-09 | 162 | 141 | 3420 | 33 | 930 | 10.01 | 32 @ 910 (9.68) |
| random 100 | 100 | 122 | 3092 | 30 | 999 (lag 0.25 d) | 13.31 | 29 @ 932 (11.46) |
| random 150 | 150 | 147 | 3600 | 35 | 886 | 8.16 | 34 @ 859 (7.61) |
| random 200 | 200 | 204 | — | 39 | 964 | 9.33 | — |
| random 250 | 250 | 230 | — | 41 | 942 | 8.33 | — |
| all minus route A (46) | 252 | 231 | — | 43 | 975 (lag 0.05 d, 121 s) | 8.73 | 42 @ 945 (8.21) |
| **all** | 298 | 253 | 4709 | **46** | 1004 | 8.79 | 45 @ 972 (8.27), 44 @ 968 |
| all, **beam 300** | 298 | about 760 | 14 095 | 46 (same route) | 1004 | 8.79 | **44 @ 941 (7.75)**, 40 @ 859 planner |
| 08+09, beam 300 | 50 | 157 | — | 22 (same route) | 789 | 8.58 | — |
| 08+09, lin_drmax 0.25 | 50 | 77 | — | 22 (same route) | 789 | 8.58 | — |
| 08+09, dvmax 1.6 | 50 | 93 | — | 24 | 899 | 12.46 | 23 @ 800 (8.69) |
| all, **2 cores** | 298 | 283 (wall 253 at load 250) | — | 46 (identical) | 1004 | 8.79 | — |

- **Fit (random pools):** depth ≈ 4.9·N^0.39. 37 flybys needs N ≈ 170.
- **Coherent pools** lose 2-5 flybys against random pools of the same size. A t10d route's own targets are phased for that route's trajectory, which the planner's leg model cannot reproduce (stage 6 R1: 0.15-0.25 AU coasts and >300 d legs).
- **Reference:** t10d flies the 08+09 pool with 2 craft (26+24 at 871/877 kg). The planner gets 22, plus 8 from a second sequential route (30/50), or 33-34 with the joint planner.
- **Time is not the binding limit.** Every restricted-pool route still runs to 14.0-14.9 yr. Depth is lost to longer phasing legs, not to fuel.

## 3. `tools/skel_fill.py`

### 3.1 Skeleton-free single route (verified)

```
~/.venvs/astro313/bin/python tools/skel_fill.py results/n8/skel8b OUT --only 09 --drop-wp all --mode windowed \
   --portfolio 300:0.15:45 --retry '' --pools POOLS.json --pad-prize 0.0 \
   --m0 1600 --vinf 4.0 --dvmax 1.2 --try 3 --lam 0.4 --max-tank 1150 --nproc 2
```

**It needs a skeleton directory, but not waypoints.**
- The positional `src` must be an `IFleet` directory: route names come from its `route_NN.npz` files, and `--pools` / `--only` are keyed by those names.
- `--drop-wp all` empties every skeleton, which leaves a free beam with a prize array.
- Measured at beam 100 on 1 core over the 50-pool: 219 s wall / 158 s CPU. Result: 26 flybys, the same state as `beam_search` soft (22 in the pool + 4 pad), twin 843 kg, 2/2 settled.
- Output: `OUT/route_09.npz`, `fleet.json`, `result.json`, `log.txt`.

**Caveats**
- **Do not use `--mode segmented` without waypoints.** The free tail segment runs with `max_depth = --seg-depth` (default 14; stage-6 runs used 20), so the route is capped at that many flybys. Use `windowed`.
- **The default `--portfolio` is 3 × beam 1000**: `1000:0.15:45, 1000:0.15:90, 1000:0.20:60`. Without waypoints the half-width does nothing, so the first two entries are identical runs. Pass a single entry.
- **`--pools` has no strict mode.**
  - The route's pool targets get prize 1.0 (or `--easy-prizes`).
  - Every other target not yet used stays legal at `--pad-prize` (default 0.5).
  - "Not yet used" excludes targets used by earlier routes of the same run or by the `--start` fleet.
  - `--pad-prize 0` makes the other targets stepping stones: the same in-pool depth as strict, +50% CPU.
  - For a truly strict pool, call `beam_search` directly (section 2.1).
- **The Pareto key and the selection value in `--pools` mode count pad targets as ordinary flybys.** The Pareto key is `len(s.seq)`, and `nnew` counts everything not covered by other routes. The selector can therefore prefer a route padded with out-of-pool targets. P06: 32 flybys = 19 of 32 assigned + 8 waypoints + 5 pad.
- **`--pad` only acts with `--assign`.** It sets the number of extra nearest non-group targets per route: nearest by closest approach to the SKELETON trajectory, from `RI.w_cands` within `--assign-dmax` 2 AU, default 24. They join the route's pool at `--pad-prize`. So `--assign` needs real skeleton trajectories, while `--pools` ignores `--pad`.
- **Other defaults:** `--m0 1600`, `--vinf 4`, `--dvmax 1.2`, `--max-depth 60`, `--try 4` (Pareto depths settled, deepest first), `--lam 0.05`, `--max-tank 1100`, `--settle-timeout 240`, `--lin-fallback 2.5`, and `--prizes results/n8/prizes_H.json`. The hard set derived from `--prizes` is never used afterwards.

### 3.2 With a skeleton [cited, stage 6]

`results/n8s6/runG.sh`:

```
--mode segmented --portfolio 1000:0.15:45 --seg-k 50 --seg-depth 20 --m0 1600 --vinf 4.0 --dvmax 1.2 --wp-dvmax 2.5 --max-tank 1150 --try 4 --lam 0.4 --nproc 8
```

- About 190 s per route on 8 procs.
- Full-pool depth per skeleton: 01 41, 05 40, 12 39, 09 38, 13 38, 03 37, 10 36, 06 35.
- Settled 26/29.
- Failures: route 12 at 39 flybys (miss 4.5e8 km) and 38 flybys (miss 158 km, just over the 150 km threshold); route 13 at 36 flybys (miss 2.1e6 km).

## 4. `ctoc14/jointsearch.py` / `tools/run_joint.py`

### 4.1 Call and state

```python
import ctoc14.jointsearch as JS          # or results/s13/planner_probe/jointsym.py (same code, order-free dedup key)
P = BP(... same as 2.1 ..., prize=prize); P.w_rare = 1.0   # prize enters the joint score only through w_rare
P.wait = 60 * DAY; P.collect_joint = []                    # coast step when the expanding craft has no leg; harvest
best = JS.joint_beam_search(RI.eph(), [1600.0] * K, np.arange(0, 801, 20) * DAY, P=P,
                            excluded=sorted(set(ALL) - set(pool)), n_proc=1, verbose=False)
tours = JS.joint_tours(best)             # list of tour_json; also evaluate every Joint in P.collect_joint
```

**State.** `Joint(craft = K × search.State or None (unlaunched), active flags, ONE shared visited bitmask, launch_t, parent)`.

**Expansion.** Only the active craft with the earliest current time is expanded (via `search.expand`, or through a precomputed root pool if it is unlaunched), so each child has exactly one more flyby.
- If that craft has no feasible leg, it coasts `P.wait` and tries again; with no fuel or time left, it retires.
- The target-to-craft assignment is therefore decided inside the beam, in time order.

**Beam.**
- Children are sorted by `w_fuel·Σfuel/1400 + Σ(t_i − t_launch,i)/T − w_rare·Σprize`. Elapsed time, not absolute time.
- They are deduplicated on `(visited, per-craft 10-d time bins IN CRAFT ORDER)`.
- The width is `P.beam`.

**Do not use the `tools/run_joint.py` CLI for pool re-planning.** It has no prize input, and its defaults are the 09-10 swarm settings: dvmax 2.5, linear legs off, tof_refine off, wait 120 d, mmargin 20, launch-max 400 d. Call the function as above; this is what `tools/run_lns.py` `replan_pair` / `replan_group` do, with `w_rare = 1` and `wait = 60 d`.

**Waiting.** `P.wait` must be greater than 0: with `wait = 0` the coast loop never ends. `--wait 100000` (days) effectively disables waiting.

### 4.2 Measurements (beam 100 unless noted, 1 core)

| run | pool | K | covered (per craft) | planner Σ J_i | CPU s | order-free distinct share of beam (mean / min) | median children per state |
|:--|--:|--:|:--|--:|--:|:--|--:|
| stock key | 50 | 2 | 30 (14+16) | 2.197 | 259 | — | about 2 |
| order-free key | 50 | 2 | **33** (18+15) | 2.223 | 281 | 0.90 / 0.55 | 2.2 |
| order-free, beam 300 | 50 | 2 | 34 (12+22) | 2.243 | 1143 | 0.89 / 0.72 | 2.1 |
| order-free, **no wait** | 50 | 2 | 21 (6+15) | 2.123 | 64 | 0.93 / 0.74 | 1.5 |
| sequential single beams (reference) | 50 | 2 | 30 (22+8) | — | 76 | — | — |
| order-free key | 77 | 2 | 41 (21+20) | 2.283 | 444 | 0.85 / 0.71 | 2.4 |
| sequential single beams (reference) | 77 | 2 | **44** (26+18) | — | 161 | — | — |
| stock key | 77 | 3 | 48 (15+12+21) | 3.350 | 629 | **0.37 / 0.11** | 2.3 |
| order-free key | 77 | 3 | **51** (13+19+19) | 3.340 | 530 | 0.74 / 0.50 | 2.2 |
| order-free key | 298 | 2 | 82 (39+43), twin 982+1025 kg, Σ J_i **2.743** | 2.682 | 715 | 0.95 / 0.60 | 10.7 |
| sequential single beams (reference) | 298 | 2 | **89** (46+43), twin 1004+975 kg, Σ J_i **2.712** | — | 484 | — | — |

t10d references: pool 50 = 2 craft for 50 targets (Σ J_i 2.467); pool 77 = 3 craft for 77 (3.645).

### 4.3 Why 3-4 craft are "weak (beam too thin)" (09-16 note), now measured

1. **Permutation duplicates.**
   - The dedup key is ordered by craft index, and the craft are identical. So the same fleet with craft labels swapped survives as separate beam entries with identical score.
   - At K=3 only 37% of the kept beam is order-free distinct (11% at the worst level). About 63% of the width is wasted. At K=2 the loss is only 10%.
   - The fix is one line: sort the per-craft time-bin tuple in the key, as in `planner_probe/jointsym.py`. It lifts the order-free share to 74% and coverage 48 → 51 (K=3), and 30 → 33 (K=2).
   - It is necessary for K ≥ 3 but not sufficient.
2. **Only the earliest craft branches.**
   - Each joint state expands one craft, so there are about 2 children per state in pools of 50-77 and about 11 in the full pool, against up to 60 per state in a single beam.
   - The width goes into assignment and timing variants, not route alternatives.
   - A 3× beam (K=2, pool 50) bought +1 target for 4× the time.
3. **Wait or retire.**
   - With waiting disabled, a craft retires at its first dead end: K=2, pool 50 covered 21, fewer than ONE free route (22).
   - With waits, tours contain 0-7 coasts of 60-360 d. `impulsive.from_tour` spans each coast plus the following leg with ONE Lambert arc (a different, often degenerate transfer). The stock seed settles only 2/7 waited tours; it fails by miss (3.9e7 km), by the 240 s timeout, or by the 2000 kg screen (seed tank 2955 / 8550 kg).
   - `planner_probe/twin_wait.py::settle_tour_wait` propagates the coast first and seeds the leg from the coast end: 4/5. Using either seed: 5/5.
4. **Even when the beam is not thin, the joint score loses to sequential planning.** At K=2 over the full pool the beam is 95% distinct, and it still finishes with 82 targets at Σ J_i 2.743, against 89 at 2.712 for two sequential free beams.

**Verdict on "make it work with a bigger beam".**
- Settling is no longer the blocker: use both seeds.
- The joint planner does not beat sequential single-route beams at K=2, at any beam width tried (100, 300).
- At K=3 the order-free key is mandatory, and the result is still 51 of the 77 that t10d's three craft fly.
- Its only measured advantage is in very sparse pools (50: 33-34 vs 30 sequential).
- A bigger beam is not the missing piece. Beam 1000 at K=2 would cost about 1 h CPU per call (pool 50, extrapolated from 281 s → 1143 s for 100 → 300).

## 5. Settle success (item 4)

### 5.1 Stored deep planner tours

Source: `planner_probe/settle_deep.json`. Each tour was settled with `LIN_FALLBACK = 2.5` and with `LIN_FALLBACK = None`.

| tour | flybys | planner fuel | seed tank (2.5 / None) | twin tank | kg/fb | time |
|:--|--:|--:|:--|--:|--:|--:|
| n8s6/B/cols_01 #10 | 40 | 786 | 1241 / 1241 | 904 | 7.60 | 7.5 s |
| n8/columns #97 (stage 5, beam 1000) | 40 | 777 | 1231 / 1231 | 910 | 7.76 | 8.6 s |
| n8s6/B/cols_13 #6 | 38 | 887 | **1301 / 1618** | 980 | 9.99 | 7.3 s |
| n8s6/B/cols_12 #0 | 38 | 958 | 1368 / 1368 | 1072 | 12.43 | 8.0 s |

- All 4 settle with both settings.
- LIN_FALLBACK re-seeded a leg on only 1 of the 4 tours (cols_13), and both settings converge there too (980.0 vs 979.7 kg).
- So "deep tours settle" is confirmed, but on these tours LIN_FALLBACK is not what decides success.
- The 40-flyby stage-5 column that was noted as "does not settle" now settles in 8.6 s at 910 kg.

### 5.2 Rates and times

- **Single-route beam tours in this probe:** 37/38 settled, counting skel_fill. The one failure was the 18-flyby route in a 51-target leftover pool (miss 3.3e6 km, 108 s).
  - Settle time on 1 core: median 13 s, p90 57 s, max 121 s. The ladder falls back to lag 0.25 / 0.05 d on some tours at 2-9× the time.
- **Stage-6 G** (beam 1000, with skeleton) [cited]: 26/29 = 90%.
- **Joint tours:** stock seed 7/12 (0/2 on the first stock-key run, whose tours were not saved), wait-aware 7/8, at least one seed 10/10 of the saved craft-tours.
- **Twin/planner tank:** median 1.018, mean 1.022, range 0.954-1.134 (p90 1.056). The outlier, random pool 100 at 30 flybys, came in 13% heavy.

## 6. What this means for the builder (plain consequences of the numbers above)

1. **Re-planning inside a disjoint pool of about 37 targets gives about 25 flybys, not 37.** Measured with t10d route 00's own 37 targets. Per-flyby efficiency is kept (6.4 kg/fb), but depth is not. The stage-6 P06 point (32 flybys, efficiency unchanged) had only 19 of its 32 assigned targets.
2. **A route reaches 37+ flybys only when about 170-200 targets are available to it.** The pools of an 8-craft fleet must then overlap about 4.6×, and choosing which route keeps which target is the packing / set-cover step. Stage 6 measured the cost of forcing a target at about 2 free flybys; `route_cover.py` is exact at this step.
3. **The strongest single-route generator is a free beam over all targets not yet covered.** Measured: 46 @ 1004 then 43 @ 975. `beam_search` over `ALL − covered` costs about 4 CPU-min at beam 100 and about 12 at beam 300 (depth unchanged, cheaper Pareto). Use beam 300 or more for the Pareto front.
4. **For joint planning, `jointsearch` needs the order-free key, and its tours must be settled with both seeds.** Even then it does not beat sequential beams at K=2.

## 7. Files (all new; nothing shared was modified)

- `results/s13/planner_probe/probe_beam.py`: single-route probe (strict / soft pool, beam, nproc, Pareto, settle).
- `results/s13/planner_probe/probe_joint.py`: joint probe (K, pool, order-free key, wait, diagnostics, tours saved to `<stem>_tours.json`).
- `results/s13/planner_probe/jointsym.py`: copy of `joint_beam_search` with the order-free dedup key and per-level diversity stats.
- `results/s13/planner_probe/twin_wait.py`: `from_tour_wait` / `settle_tour_wait`, the wait-aware impulsive-twin seed.
- `results/s13/planner_probe/settle_probe.py`, `settle_joint_tours.py`: settle tests.
- Raw results:
  - `b*.json` / `b*.out`: beam probes. Settled twins are in `<stem>_k<depth>.npz`, in `route_NN.npz` format, loadable with `RI.ipr`.
  - `j*.json`, `j*_tours.json`: joint probes.
  - `settle_deep.json`, `settle_joint_*.json`: settle tests.
  - `sf_free09/`: skeleton-free skel_fill run.
