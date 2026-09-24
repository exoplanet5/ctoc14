# Design: priced column generation with route-fixing dives ("CG-dive") for a 10-craft cover

Author: designer "colgen" (2026-09-15). Status: design + one proof-of-concept measurement (section 15), no production code yet.
Companion docs: `docs/phasing_campaign.md` (section 4: the pool campaigns this design replaces), `docs/insertion.md`,
`docs/conversion_losses.md`, `CTOC14FullSolver.md` sections 17-20 (route pool / set cover), 41-43 (marginal cost, local search, ALNS).

## 0. Summary

The 335k-route pool cannot give a 10-craft cover because it is a dump of beam states from about 15 single-craft searches:
the good routes form a handful of families of near-duplicates, and *complementary* good routes (disjoint target sets) do not
exist in it. Set cover over it therefore either overlaps heavily (12 craft, 284-287 covered) or, at N = 10, is hopeless
(LP value 39.5, 91 targets priced at the miss value). The fix is not a bigger pool but the standard operations-research one:
generate columns *against the dual prices of the master problem* (each new route is asked to cover what the current fleet
covers badly), stabilise the prices, and turn the fractional LP into an integer fleet by *diving* (fix the most-used route,
re-price the residual problem with one craft fewer, repeat), followed by a local repair of the few uncovered targets with
`ctoc14/insertion.py`. Everything needed already exists in pieces: the LP with duals (`select_routes_milp.py --colgen`,
`price_targets.py`), a beam search that accepts per-target rewards (`Params.rarity/w_rare`), warm starts from arbitrary
states (`beam_search_from_state`), exact insertion/re-timing (`insertion.py`), exact conversion and fragment MILP.
What is missing is the loop that connects them with the *right scaling* (rewards in J units, child pruning that sees the
reward, re-pricing after every batch of columns) and the dive.

Expected result (arithmetic in section 12): planned J 17-18.5 with 10 craft (or 17.3 with 11), validated J 18-19.5 after
exact conversion. Reaching ~15 additionally needs 8 kg per flyby (encounter-time refinement everywhere), which this design
enables but does not by itself deliver.

## 1. Diagnosis of the pool (measured on all 13 jsonl files, 334 998 routes)

Script: scratchpad `colgen/pool_diag.py` (reads the pools, packs masks, computes the numbers below; not part of the repo).

| fact | number |
|---|---|
| routes with n >= 25 / >= 28 / >= 30 | 40 176 / 21 626 / 11 776 (camp5 alone has 15 036 / 10 453 / 7 385: `tof_refine` family) |
| union coverage of n >= 25 routes | 298/298 (every target is in >= 49 good routes; 64 and 218 only 49 each, 236: 153, 169: 225, 219: 275) |
| greedy **disjoint** packing of n >= 28 routes (best cost per target first) | **2 routes, 66 targets** |
| greedy near-disjoint (<= 2 shared) packing of n >= 22 routes | 4 routes, 122 targets |
| top-300 routes with n >= 30: pairwise Jaccard | mean 0.55, median max-Jaccard to another route **0.95**; only 150 of 44 850 pairs below 0.05 |
| distinct first targets among n >= 28 routes / launch epochs | 17 / 12 |
| t_end of n >= 28 routes | 10 % end before 10.8 yr, median 12.7 yr (prefix states, tails unused) |
| pool LP, N = 12, margin 1.3 (`colgen/lp10.py`) | **19.55**, fleet constraint slack (sum y = 11.65), 250 routes fractional, max y = 0.24, pi median 0.06, max 0.17 |
| pool LP, N = 11 | 25.18, mu = 13.0, 46 targets at pi = 1 |
| pool LP, N = 10 | **39.50**, mu = 15.6, **91 targets at pi = 1** |
| integer cover N = 12 over union pool (`milp_u123cg_n12`) | 33.86, 284 covered |

Reading:

1. **Correlated families, not missing targets.** The LP at N = 12 reaches 19.55 only by taking 0.1-0.24 of 250 sibling
   routes each; no integer selection can imitate a fractional mixture of near-duplicates, hence the 14-point integrality gap.
   The dominance/dedup pass in `select_routes_milp.py` removes identical sets only; siblings differing in one target survive
   and are useless.
2. **Hard targets live in one family each.** Choosing the family that contains 64/218 fixes 30 other targets that then get
   duplicated by the other families. That is why the greedy exclusion campaign decays 29,31,31,30,28,26,25,22,18,8,3
   (camp1): each round's search is a *dive without re-pricing*, and by round 8 the leftovers are the targets that were cheap
   only in combination with targets already spent.
3. **The rarity campaigns used fixed, mis-scaled prices.** camp2/camp3 gave w_rare 0.3-0.5 x (0/1 or LP dual) per target,
   while a leg changes the beam score by only ~0.046 (1/30 mission time + 1.5 x 12 kg/1400). A 0.35-0.5 reward is 10 legs'
   worth, so the beam sacrificed count for priced targets: 22-24 flybys on 450 kg, J_i 1.62. And the prices were never
   re-solved after new routes appeared, so all 12 rounds chased the same targets.
4. **Child pruning ignores the reward.** `search.expand` keeps `max_children = 60` children by `tof/T + w_fuel dm/FUEL_MAX`
   *before* the rarity term enters the score, so a priced target reachable only by a mid-cost leg is discarded at the parent.
5. **Per-route cost is flat in n.** With tanks 600 + 1.3 fuel + 20, c_r spans only 1.36-1.47 for 18-31 targets; the LP prefers
   many mid-size routes (sum y = 11.65 when N = 12 is not binding). With N = 10 fixed, the only lever is *distinct* targets
   per route, so pricing must maximise collected price, not flyby count.
6. **Conversion factor depends on leg length more than on dv** (71 `*_info.json` files with `dm_plan`, 602 legs, scratchpad
   `colgen/calib.py`): actual/planned propellant summed per bin is 1.05 for legs < 120 d, 1.3 for 120-300 d, **1.6 for
   300-450 d**; per-dv bins are noisy (the 2-leg window shifts cost between neighbours: dv < 0.5 legs show 1.9-2.2 x, dv
   0.9-2.0 legs 0.8-0.96 x); launch legs cost 12 kg each on average though planned at 0. The flat 1.3 margin is right on
   average and wrong per route.

## 2. What a 10-craft cover must look like (target arithmetic)

J = sum_i (1 + x_i + x_i^2) + N_miss, x = (m0 - 600)/1400, 131 and 144 always missed.

| fleet | distinct targets per craft | planned fuel | tank rule | m0 | J_i | sum J_i | misses | planned J |
|---|---|---|---|---|---|---|---|---|
| 10 | 29.6 (296) | 350 | 600 + 1.3 f + 20 | 1075 | 1.455 | 14.55 | 2 + 2 | 18.6 |
| 10 | 29.6 | 350 | 600 + 1.15 f_cal + 20 (calibrated, section 4.2) | 1023 | 1.393 | 13.9 | 2 + 2 | 17.9 |
| 10 | 29.3 (293) | 350 | calibrated | 1023 | 1.393 | 13.9 | 2 + 5 | 20.9 |
| 11 | 27 (297) | 310 | 600 + 1.3 f + 20 | 1023 | 1.393 | 15.3 | 2 + 1 | 18.3 |
| 10 | 29.6 | 250 (8.4 kg/flyby) | 600 + 1.3 f + 20 | 945 | 1.307 | 13.1 | 2 + 2 | 17.1 |

So "planned J <= 18" means: 10-11 craft, **<= 3 non-physical misses**, planned fuel <= 350 kg per craft with a tank rule
no worse than 1.15 x calibrated. Every extra miss costs as much as 60 kg of tank on every craft; the whole design is therefore
organised around coverage (complementary columns), with fuel a secondary term priced in J units.

## 3. Algorithm overview

```
init      columns <- pool routes (n >= 12, dedup) + calibrated costs
          N_path  <- [12, 11, 10]                       (homotopy in the fleet size, section 6.3)
root CG   for N in N_path:
            repeat <= R rounds:
              (pi, mu, y) <- solve RMP(columns, N)      (LP, HiGHS, 2-20 s)
              pi_s <- smooth(pi)                        (Wentges, section 6.1)
              cols <- PRICE(pi_s, mu, excluded = {}, K parallel diversified beams)   (section 5)
              add cols with rc < -eps, diversity filter; stop when no rc < -eps twice
dive      fixed <- []; excluded <- {131,144}; N_left <- N
          while N_left > 0:
            (pi, mu, y) <- solve RMP(columns, N_left, rows = uncovered targets)
            cand <- top-3 routes by y (then by LP value after fixing, section 7)
            fix r*: fixed += r*; excluded |= S(r*); N_left -= 1
            re-price the residual: 2-4 CG rounds with excluded (beams never visit excluded targets)
            (limited discrepancy: at levels 1-3 also try the runner-up once)
repair    uncovered U: tail extension (beam_search_from_state), insertion/re-timing (insertion.py),
          replace moves, extra small craft by the marginal-cost rule, encounter-time polish  (section 8)
convert   convert_fleet.py -> fragments -> select_frags_milp.py with all older fragments (safety net)
ALNS      optional outer loop: destroy/repair on the fleet, every route seen becomes a column,
          RMP recombination every K iterations, pricing warm-started from ALNS routes (section 9)
```

## 4. Master problem (RMP)

### 4.1 Formulation (unchanged from `select_routes_milp.py`, set cover with misses)

min sum_r c_r y_r + sum_a w_a  s.t.  sum_{r ∋ a} y_r + w_a >= 1 (all a != 131,144),  sum_r y_r <= N,  0 <= y, w <= 1.
Duals: pi_a in [0, 1] (the w_a bound makes pi_a <= 1), mu >= 0. Reduced cost of a route: **rc_r = c_r + mu - sum_{a in S_r} pi_a**.
Overlap stays allowed in the master (section 20 of the brief); the dive removes it in practice (section 7).

Implementation: `ctoc14/colgen.py::RMP` holding packed 38-byte masks (as `select_routes_milp.py` does), costs, source
locations, `add(columns)`, `solve_lp(N, rows)` (returns value, y, w, pi, mu), `solve_mip(N, time_limit)` for incumbents,
`fix(route)` (moves its targets out of the row set), `prune(rc_threshold, age)`. LP over 200k columns takes 16-20 s (measured);
over the 10-20k columns the RMP will hold after pruning, 1-3 s.

### 4.2 Column cost: calibrated tank instead of the flat 1.3 margin

c_r = cost_sc(600 + 20 + s * F_cal(r)), F_cal = 12 kg (launch leg) + sum_legs f(tof_leg) * dm_plan_leg,
f = 1.05 (tof < 120 d), 1.30 (120-300 d), 1.60 (300-450 d), 2.0 (> 450 d, i.e. discouraged), s = 1.10 safety.
The table is the bin ratio from the 71 info files (section 1, item 6) and must be re-fitted after the first conversions of
CG routes (experiment E5). The flat rule stays available (`--cost-model margin`) for comparison with older numbers.
Why it matters: the pricing beam feels the *expected actual* fuel per leg, so 300-450 d linearised legs are charged 1.6 x
and the search stops preferring them; the tank right-sizing after conversion then changes J_i by < 0.03 instead of 0.1-0.2.

## 5. Pricing subproblem

### 5.1 Objective and beam score

Find routes with rc_r < 0: maximise sum_{a in route} pi_a - c(F) - mu subject to the planner's leg model (Lambert junctions,
eta 0.6, dv_max, linear legs <= 450 d, 15-yr window, m0 = 1000 kg for the acceleration model). This is a prize-collecting
orienteering problem with a nonlinear (convex) fuel cost and a time budget; we solve it heuristically with the existing
depth-synchronous beam search, scored in **J units**:

score(s) = w_t * t_s / T_MISSION + w_fuel * F_s / FUEL_MAX - sum_{a in s} pi_a

* w_fuel = 1400 * dc/dF = 1.3 * (1 + 2x) ≈ **2.2** at x ≈ 0.33 (calibrated cost: s * f(tof) * (1 + 2x), i.e. per-leg
  weights 1.15-1.75 -> implement as a per-leg fuel weight `w_fuel * f(tof)` in `expand`).
* w_t is the only heuristic: the value of remaining time, needed because a depth-synchronous beam compares states with the
  same flyby count but different elapsed time. Estimate: net prize rate = (mean pi - fuel J per flyby) / mean tof =
  (0.05 - 0.02) / 190 d -> 0.86 per mission; use **w_t = 0.9**, sweep 0.6-1.2 in E2.
* pi enters exactly as `Params.rarity` with `w_rare = 1` today, but (i) **child pruning must include the prize** (rank
  children by `dtof/T * w_t + w_fuel dm/FUEL_MAX - pi_k` before `max_children`; 3-line change in `expand`), (ii) roots need a
  `first_targets` mask for diversification, (iii) `Params.prize` replaces the `rarity/w_rare` pair (keep the old names as
  aliases; prize = None reproduces today's search bit for bit, regression test).
* Route completeness: a craft may stop anywhere, so **every collected state is a column**; its rc is exact from its own fuel
  and visited set. The pricing returns the best few hundred states by rc, not just the deepest.

### 5.2 Columns per pricing run

From `P.collect` (all beam states of all depths, 5-14k per search): keep rc < -eps (eps 0.02); dedupe by target set (keep
cheapest); per (launch epoch, first target) keep the 5 best; greedy Jaccard filter (< 0.85 against columns already in the
RMP and among themselves); cap 300 per beam. Store rc at generation time in the column record (`tour_json` + `rc`, `pi_id`).

### 5.3 Diversification: K beams per round, each on 2 processes

| beam | constraint | why |
|---|---|---|
| A1-A3 | launch window [0,200], [200,450], [450,750] d (`t_launch_grid`) | the pool's good routes use 12 launch epochs, all < 220 d except two |
| B1-B2 | `first_targets` = 15 highest-pi targets (B1), targets 16-40 (B2) | forces routes that start at hard targets while the craft is heavy and Earth-phased |
| C | pi perturbed: pi * (1 + 0.15 N(0,1)), clipped to [0,1] | breaks ties between sibling routes, cheap diversity |
| D | "must include": pi_a += 0.5 for the 8 pi = 1 targets with the fewest columns | Lagrangian-style boost for targets the LP cannot cover |
| E | warm start: `beam_search_from_state` from the prefixes (k = n/3, n/2) of the 10 routes with the largest y, with the current pi and `excluded` = the prefix's targets | columns that share a good head and re-optimise the tail for the prices; also the ALNS coupling (section 9) |

Per round run 4-6 of these (rotate), 2 processes each -> 8-12 processes for one round; with a 5-process cap, 2-3 beams at
a time. All beams use `dv_max = 1.0`, `lin_tofs` <= 450 d, `lin_drmax = 0.15 AU`, `m_margin = 40`, `vinf_cap = 2`,
beam 300, no `tof_refine` (4 x slower) except in the last two dive levels and the repair.

### 5.4 Forbidding duplicates

In the dive, targets of fixed routes are `excluded` in the beam (`roots`/`expand` masks): zero overlap by construction, as
in `run_pool_campaign.py`. Root-level pricing excludes nothing (the LP wants overlapping fractional columns). Option
`--allow-passthrough`: covered targets get pi = 0 but stay visitable (a duplicate flyby can be the cheapest connector);
off by default because conversion time and the validator are indifferent but the beam wastes time on them.

## 6. Stabilisation and bounds

### 6.1 Dual smoothing (Wentges)

pi_used = alpha * pi_best + (1 - alpha) * pi_LP with alpha = 0.5; pi_best = the duals of the best Lagrangian value seen so
far (section 6.2). If pricing with pi_used finds no rc < -eps column, re-price with pi_LP (alpha = 0) before declaring
convergence ("mis-pricing" rule). Duals are already boxed in [0, 1] by the model, mu >= 0.

### 6.2 Bound

For the LP with a cardinality constraint, LP* >= z_RMP + N * min(0, rc*) where rc* is the minimum reduced cost over *all*
routes (Lasdon). The beam is heuristic, so rc* from it is an upper bound on the true minimum and the number is a
**pseudo-bound**: used as a stopping/monitoring quantity (stop the root when z_RMP - LB_pseudo < 0.3 for two rounds), not
as a proof. Two honest bounds remain: the LP over all generated columns (an upper bound on LP*) and the arithmetic of
section 2 (10 x 1.39 + 2 = 15.9 with zero extra misses is the floor for this tank rule). An exact pricing (label-setting
over (asteroid, time bin) with velocity-independent leg costs) is not affordable: the junction cost depends on the incoming
Lambert velocity, so a precomputed edge cost would need 300 x 300 x 550^2 / 2 Lambert solves. Do not attempt.

### 6.3 Homotopy in N and initial duals

Start the RMP at N = 12 (pool duals are informative: median 0.06, max 0.17, mu = 0), run 4 CG rounds, then N = 11 (4 rounds),
then N = 10 (until convergence). At N = 10 the pool duals are degenerate (91 at 1.0, mu 15.6): pricing against them chases
91 "misses" at once and the beam collapses to 20-flyby routes, which is exactly camp3's failure. The N = 12 columns make
the N = 11 and N = 10 duals graded, and mu falls as good columns appear (mu is the marginal value of an 11th craft; the
target is mu < 1.5, i.e. an 11th craft would not pay for itself).

### 6.4 Column management

Keep at most 20k columns: every 5 rounds drop columns with rc > 1.0 for 5 consecutive solves unless they are the best
column of some target. Never drop fixed routes or incumbent routes. Write every added column to `outdir/columns.jsonl`
(same schema as the pools plus rc / round / beam tag) so later runs start from it.

## 7. From fractional to integer: dive with re-pricing (heuristic branch-and-price)

Branching on single route variables (y_r = 1) is the natural rule here because a route is a craft: fixing y_r = 1 is a bulk
"target-to-craft assignment" branch (all targets of r go to one craft, the other craft see them as excluded).

* **Candidate choice at each level**: top-3 routes by y in the current LP; for each, solve the LP with it fixed (2 s each)
  and take the one with the lowest LP value *after* fixing (one-step look-ahead). Tie-break by number of targets with
  pi >= 0.5 (hard-target content: fix those first, while the residual has the most freedom).
* **Residual re-pricing**: after fixing, the row set shrinks to the uncovered targets, N_left -= 1, and 2-4 CG rounds run with
  `excluded` = all fixed targets (beam CPU drops roughly with the number of remaining targets: camp1's rounds went 207 -> 75 s
  from 298 to 76 targets). The pricing beams here are the greedy campaign's searches, but *steered by the prices of the
  residual LP* instead of by "as many as possible", which is what flattens the 29,31,...,18,8,3 decay.
* **Limited discrepancy search**: at levels 1-3 keep the runner-up candidate as an alternative branch and explore it once
  (3 extra dives at most); at deeper levels no backtracking. Each dive is scored by planned J after repair (section 8).
* **Last two levels**: pricing with `tof_refine` (the tail craft must be as efficient as possible; +10 % flybys measured).
* **Ryan-Foster (target pairs together/apart)** is the exact alternative compatible with the beam (apart: mask; together:
  drop columns containing exactly one of the pair, and in the beam give the partner a bonus once one is visited). It gives a
  balanced tree but each node needs a full re-pricing; not for v1. **Route-pair branching** (y_r + y_s <= 1 for two sibling
  routes) is a cheap cut against family correlation and is applied implicitly by the Jaccard filter in 5.2.
* **Forbidding the fixed route's siblings**: after fixing r*, delete columns with Jaccard > 0.5 to r* (they can only duplicate).

Output of a dive: `outdir/dive_k/tour_sc1..N.json` (planner schema with `m0` from the cost model) + `summary.json`
(covered, missed, planned J, LP values per level, pseudo-bounds).

## 8. Repair of the uncovered targets

Expected after 10 fixed routes: 5-15 uncovered (U). Ordered by cost per target gained:

1. **Tail extension** (new `colgen.extend_tail`): for each fixed route with t_end < 14.6 yr, `beam_search_from_state` from
   its final state with prize 1 for U, 0 otherwise, w_fuel 2.2, `excluded` = everything covered, beam 200, 1-2 processes,
   20-60 s per route. Half the pool's good routes end before 12.7 yr, so this is the cheapest source of extra targets.
2. **Insertion with re-timing**: `insertion.insert_targets` over the 10 routes with `InsParams(margin=10)`; the routes are
   planned with `m_margin 40`, so 30 kg is available before the tank grows; then `grow_m0=True` (the tank grows by the
   inserted propellant x 1.15: +12 kg -> +0.016 J, always worth a target).
3. **Replace** (new `insertion.evaluate_replacements`): for u with no feasible insertion, for every route r and every
   b in r that is either duplicated or insertable elsewhere for < 0.3 J: remove b (merge legs k-1 -> k+1 at the existing
   times, one Lambert call, re-run `mass_chain`), then `evaluate_insertions(u)` on the shortened route; accept the cheapest
   (delta J < 1). Removal needs `insertion.remove_target(TM, k)` (build_tour on the shortened sequence).
4. **Extra craft by the marginal-cost rule** (brief section 41): if |U| >= 3 after 1-3, run a tail beam (launch grid over
   the whole 15 yr, m0 700-1000, `excluded` = covered) and keep the craft iff it covers K > c = cost_sc(m0) targets
   (c = 1.08 at m0 700: 2 targets already pay).
5. **Encounter-time polish**: coordinate descent on all flyby times, +-6 d at 1 d then 0.25 d (the measured -29 % Lambert-chain
   dv on a 32-leg tour), on every route; then right-size tanks with the calibrated cost. This is `tof_refine` applied a
   posteriori and costs 30-60 s per route (one `build_tour` per trial time).
6. **Convert** (`tools/convert_fleet.py --nproc 5`), validate, and run `tools/select_frags_milp.py` over the new fragments plus
   the 134 existing ones: the fragment MILP is the safety net that can never be worse than 30.789.

## 9. Coupling with the fleet-level ALNS

The fleet solution is a list of `TourModel`s ((asteroid, time) sequences with launch epoch and tank). ALNS on it is the
VRP-style local search of brief sections 42-43 with the planner's leg model as the surrogate objective (calibrated planned
J). All operators are cheap because junctions are local (`insertion.py`: 0.4 s per (tour, target) evaluation, one
`build_tour` per removal).

Destroy: random removal (5-10 targets), worst removal (largest propellant saved by removal), route removal (whole craft ->
U), **price-guided removal** (targets with the highest current pi: the ones the master finds hard, so that repair may place
them elsewhere), time-window removal (all flybys in a 1-yr window across all craft). Repair: cheapest insertion, regret-2/3
(from the (tour, target) table `insert_targets` already keeps), tail extension, new small craft (marginal-cost rule),
replace. Acceptance: simulated annealing on planned J; every 200 iterations the incumbent is converted for a true J check.

Two-way coupling with CG (the "ILS-SP" hybrid of the VRP literature):

1. **Every route ALNS produces is a column** (n >= 15, after the Jaccard filter) -> `RMP.add`. Every K = 100 iterations the
   RMP is solved as a MILP (time limit 120 s) over all columns; if it beats the ALNS incumbent, ALNS restarts from it. This
   is where families of correlated routes become *useful*: recombination picks the best member of each family.
2. **Pricing warm-started from ALNS routes**: beam E in section 5.3 seeds `beam_search_from_state` with ALNS routes' prefixes
   and the current pi, and ALNS's price-guided destroy uses the same pi; the two searches share the dual vector through
   `outdir/duals_latest.json`.
3. Schedule: root CG + dive first (sections 5-7), then ALNS from the repaired dive fleet, with the RMP recombination every
   K iterations and one CG round (2 beams) every 500 iterations if the LP duals changed by more than 0.1 on 20 targets.

## 10. Data structures

* `Column`: packed target mask (38 bytes), cost (float32), n_distinct, fuel_cal, launch epoch, first target, source
  (file, byte offset) or inline tour JSON for generated ones, rc at creation, round id. Held in numpy arrays (masks as
  (n, 38) uint8, unpacked lazily per LP build as today).
* `RMP`: arrays above + `rows` (active target indices), `fixed` (route ids), `N`, HiGHS LP via `scipy.optimize.linprog`
  (as now) and `milp` for incumbents; duals returned as a (300,) array with zeros for inactive rows.
* `PricingSpec`: dict(tag, launch_window, first_targets, pi_override, warm_states, beam, tof_refine, nproc).
* `PricingResult`: list of columns + best rc + wall/CPU time + the pi used (for the Lagrangian record).
* Fleet (dive/ALNS): dict name -> `insertion.TourModel`; serialised with `TourModel.to_json` (drop-in for `convert_fleet.py`).
* Logs: `outdir/cg.log` (round, N, z_RMP, sum y, mu, #pi=1, best rc per beam, columns added, pseudo-bound), `columns.jsonl`,
  `duals_round_k.json`, `dive_k/summary.json`.

## 11. Runtime (from measured numbers)

* One full pricing beam (beam 300-400, 298 targets, no `tof_refine`) = 200 s on 10 processes = **~2000 CPU-s**; camp5 with
  `tof_refine` = 975 s x 10 = 9750 CPU-s. CPU falls with the residual target count (207 s -> 75 s for 298 -> 76 targets).
* Root CG: 12 rounds x 4 beams x 2000 CPU-s = 96k CPU-s = **2.7 h on 10 CPUs, 5.3 h on 5**. LP solves: 12 x 20 s.
* Dive: 10 levels x 3 rounds x 3 beams x ~1000 CPU-s (shrinking residual) = 90k CPU-s = **2.5 h on 10 CPUs**; last two
  levels with `tof_refine` add ~1 h. Limited discrepancy (3 extra partial dives from levels 1-3): +2-4 h if used.
* Repair: tail extensions 10 x 60 s, insertion table 10 routes x 15 targets x 0.4 s = 60 s per pass, polish 10 x 60 s:
  **< 30 min**. Conversion: 10 tours x 5-15 min on 5 processes: 30 min. Fragment MILP: 1-2 min.
* Total: **6-9 h per full run at 10 processes**; a "fast mode" (beam 200, 3 beams per round, 6 root rounds, no LDS) fits in
  2-3 h and is the mode for E2/E3. All of it is embarrassingly parallel at the beam level, so the process cap is the only
  constraint; never run two conversions concurrently (memory).

## 12. Expected J

Root LP: expect the N = 10 LP to fall from 39.5 (pool) to 15-17 once complementary columns exist (10 x 1.39-1.45 with
fractional overlap savings + 2 physical misses). Dive: integrality typically costs 1-2 J here (2-4 extra misses before
repair, 0-2 after). Planned J after repair, 10 craft: 13.9-14.6 (sum J_i at tanks 1023-1075) + 2 + (1-2) = **17-18.5**; with
11 craft (27 targets, 310 kg each): 15.3 + 2 + 1 = 18.3. Validated J: patient legs convert at 1.0-1.1 x (dv < 0.9; the
calibrated cost already contains the leg-length factors), tanks re-sized after conversion move J_i by +0.03-0.05 each, and
one lost flyby per two craft: **18-19.5 validated**, i.e. 11-12 points better than the current 30.789. Getting to ~15 needs
8-9 kg per flyby (250-270 kg per craft: `tof_refine` in every pricing beam, +4 x CPU, plus the polish), which turns the
same 10-craft cover into sum J_i 13.1-13.4 -> **planned 15.1-15.4, validated ~16-17**.

## 13. Implementation plan (ordered; files and functions)

1. `ctoc14/search.py` (extend, backward compatible): `Params.prize` (array 300, J units) with `w_t`; `State.score` uses
   `w_t * t/T - prize`; `expand` ranks children with the prize before `max_children` and applies a per-leg fuel weight
   `w_fuel * f(tof)` when `Params.fuel_factor` is given; `roots(..., first_targets=None)`; keep `rarity/w_rare` as aliases.
   Regression: `prize=None` reproduces `results/phasing/camp1/round1_best.json` bit for bit (test).
2. `ctoc14/colgen.py` (new): `route_cost(tour, model)` (margin / calibrated, table from section 4.2), `RMP` (section 4.1,
   10), `price(eph, rmp, spec)` (wraps `beam_search`/`beam_search_from_state`, returns `PricingResult`),
   `columns_from_states(states, pi, mu, filters)` (section 5.2), `smooth_duals`, `pseudo_bound`, `dive(...)` (section 7),
   `extend_tail(eph, TM, targets, P)`.
3. `tools/run_colgen.py` (new CLI): `outdir pools... --N 10 --homotopy 12,11,10 --rounds 4,4,8 --beams A1,A2,A3,B1,C,D,E
   --beam 300 --nproc 5 --cost-model calibrated --dive --lds 1 --tof-refine-last 2 --resume`; writes the logs of section 10
   and the dive tour directories; `--resume` restarts from `columns.jsonl` + `duals_latest.json`.
4. `ctoc14/insertion.py` (extend): `remove_target(TM, k)`, `evaluate_replacements(eph, TM, u, P, removable)`,
   `retime_polish(eph, TM, half=6 d, steps=(1, 0.25))`.
5. `tools/repair_fleet.py` (new CLI): steps 1-5 of section 8 on a tour directory -> a new directory for `convert_fleet.py`,
   with `--targets` defaulting to the dive's `missed`.
6. `tools/select_routes_milp.py`: keep as is (the fragment/route MILP for comparisons); its LP code moves into
   `colgen.RMP` and the script imports it.
7. `tests/test_colgen.py`: (a) rc of a column recomputed from its tour JSON equals the RMP value; (b) pi in [0,1], mu >= 0,
   LP value equals `price_targets.py` on `prices3` inputs; (c) prize = 0 regression; (d) `columns_from_states` filters;
   (e) a 3-route toy dive covers the toy target set.
8. Phase 2: `tools/alns_fleet.py` (section 9) after E3 passes.
9. Calibration refresh: `tools/calib_legs.py` (the scratchpad `calib.py` made permanent) re-fits the f(tof) table from all
   `*_info.json` after every fleet conversion; the table lives in `ctoc14/colgen.py::CALIB` with the fit date.

## 14. Experiments (pass/fail)

| id | what | pass criterion |
|---|---|---|
| E1 | one pricing beam with N = 11 pool duals as prize (w 1.0, w_fuel 2.2, beam 200, 2 proc), rc of all collected states | >= 20 columns with rc < 0 and pairwise Jaccard < 0.8; N = 11 LP drops by >= 1.0 after adding them (result: section 15) |
| E2 | root CG, fast mode (homotopy 12/11/10, 3 beams per round, 6 rounds per N), w_t in {0.6, 0.9, 1.2} | N = 10 LP <= 20 (pool: 39.5); MILP incumbent at N = 10 covers >= 285; w_t chosen by the lowest LP |
| E3 | one dive (no LDS) + repair from the E2 columns | 10 routes cover >= 290 before repair, >= 293 after; planned J <= 18.5; zero overlap |
| E4 | exact conversion of the E3 fleet, fragment MILP with the 134 old fragments | >= 8/10 craft lossless; actual fuel <= 1.15 x calibrated for >= 8 craft; validated J <= 21 (any improvement over 30.789 is kept as TEAM) |
| E5 | calibrated cost vs flat 1.3 margin on the 71 info files and the E4 conversions | median |actual - predicted| / actual < 6 % for calibrated vs > 10 % for flat; if calibrated fails, keep 1.3 and widen `f` bins |
| E6 | LDS (3 extra partial dives) and `tof_refine` in the last two levels | planned J improves by >= 0.5 over E3 for <= 4 h extra |
| E7 | ALNS + RMP recombination from the E3 fleet, 2 h | planned J improves by >= 0.5; every recombination MILP solves in < 120 s |

Fast-mode E2+E3 together fit in ~3 h at 5 processes; run them before anything else.

## 15. Proof of concept (E1) — measured 2026-09-15

See the appendix at the end of this file (filled in from scratchpad `colgen/price_N11_w1.log`).

## 16. Risks

1. **Heuristic pricing**: the beam may miss negative-rc routes, so convergence is "no column found", not optimality; the
   pseudo-bound can lie. Mitigation: diversification (5.3), mis-pricing re-check with raw duals, warm starts from incumbents.
2. **w_t mis-set**: too small -> the beam hoards prize and runs out of time early (short routes); too large -> today's
   count-maximising behaviour returns. Sweep in E2; it is one scalar.
3. **Dual degeneracy at N = 10** (0/1 prices): handled by the N homotopy and smoothing; if mu stays > 3 after the root, the
   10-craft structure does not exist in the model and the answer is 11 craft (still planned 18.3).
4. **Correlation creeps back**: beams from the same duals produce siblings; the per-(epoch, first target) cap, the Jaccard
   filter and sibling deletion after fixing are essential, not cosmetic.
5. **Conversion overhead on long legs** (1.6 x at 300-450 d): the calibrated cost prices it in; additionally cap `lin_tofs`
   at 450 d and keep `dv_max = 1.0`. Tanks are sized on calibrated fuel with 10 % safety; a dropped flyby costs 1, a 30 kg
   larger tank 0.04, so err on the larger tank for the last conversion pass.
6. **Runtime**: 6-9 h per full run on 10 processes; with the shared machine, plan on fast mode first and `--resume`.
7. **Model gap to the leaders**: allocation alone gives 17-19; the leaders' 5-10 kg per flyby implies better encounter
   timing and possibly multi-revolution legs the single-rev Lambert model cannot see. This design is the vehicle for those
   improvements (they enter through `expand`), not a substitute.
8. **Dive myopia**: the first fixed route determines much; the look-ahead + LDS at levels 1-3 is the guard, and the fragment
   MILP over all dives' converted fragments is the last resort.

## 17. Reuse (existing code that fits)

* `tools/select_routes_milp.py` LP/MILP construction and packed-mask storage -> `colgen.RMP` (same model, adds rows/fixing).
* `tools/price_targets.py` -> replaced by `RMP.solve_lp`; its dual export becomes `duals_round_k.json`.
* `ctoc14/search.py`: `beam_search`, `beam_search_from_state`, `Params.rarity/w_rare` (prize), `Params.collect` (columns),
  `tof_refine`, `lin_tofs`; `roots` gets the `first_targets` mask. The pricing subproblem *is* this beam search with a
  different score; no new search engine.
* `tools/run_pool_campaign.py`: the round loop and pool dump become the dive's residual pricing; `--exclude` = fixed targets.
* `ctoc14/insertion.py`: `build_tour`, `mass_chain`, `evaluate_insertions`, `insert_targets`, `min_m0_for`, `TourModel.to_json`
  for the repair and every ALNS operator (exact planner model, vectorised, verified to 1e-10 km/s).
* `ctoc14/jointsearch.py` wait action and `tools/improve_joint.py` leave-one-out: the latter is a special case of ALNS route
  removal + tail re-plan and its acceptance rule (planned fleet J with right-sized tanks) is reused verbatim.
* `tools/convert_fleet.py`, `tools/select_frags_milp.py`, validators: unchanged; the fragment MILP over all converted
  fragments (old + new) guarantees monotone progress of the TEAM file.
* Pools `results/phasing/camp*/pool.jsonl`, `results/phasing/pool/*.jsonl`: initial columns (n >= 12) and the N = 12 duals.
* Calibration data: the 71 `*_info.json` files with `dm_plan` (scratchpad `colgen/calib.py` shows the extraction).
