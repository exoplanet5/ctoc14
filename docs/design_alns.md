# Fleet-level Adaptive Large Neighbourhood Search (ALNS) for CTOC14-A — design

Author: designer "alns", 2026-09-15. Status: design only (nothing implemented; the numbers quoted were measured on this
machine with throw-away scripts in the session scratchpad, single core, load average 7.6).

## 0. Thesis

The gap between our best validated file (J = 30.789: 16 craft, 294 covered, 88 duplicated visits, 25 visits per craft)
and the leaders (raw J 13-14: 9-10 craft, 30-33 flybys each, no overlap, <= 3 misses) is **allocation**, not per-craft
capability: single-craft beam searches already produce 31-37 planned flybys on ~360 kg. The route-pool + set-cover
architecture cannot close it because the pool only contains tours the myopic single-craft beam happened to generate
(LP bound of the pool 21.1, best integer 33.9 planned with 12 craft, 16 misses and 19 duplicates). A fleet-level ALNS
attacks allocation directly: it moves targets between craft, opens time windows on one craft to absorb the targets of
another, dissolves whole craft, and re-times encounters at 0.25 d — every move priced with the planner's own leg model,
so that accepted fleets are exactly the kind of plan the exact converter has been calibrated on. The objective is the
true planned J with right-sized tanks (tank = fixed point of 600 + reserve + margin x planned propellant) plus misses.
With 10-11 craft x 28-32 flybys, <= 5 misses and <= 12 kg planned propellant per flyby the arithmetic (section 15)
gives planned J 17.5-19.5 and validated J 19-22; J < 20 validated requires lossless conversion (ratio <= 1.1), which
the patient-leg calibration (median 1.01, p75 1.13) says is the normal case.

## 1. Facts that shape the design (all measured this session)

1. **The planned fleets are not Lambert chains.** In the 12-craft fleets of `results/phasing/{milp_u123cg_n12,
   camp1/milp_n12, camp2/milp_n12}` 334 of 896 legs (37 %) are *linearised* legs (`ctoc14/linleg.py`, "coast + small
   corrections"): their single-revolution Lambert reconstruction is meaningless (junction dv 20-58 km/s, Lambert-chain
   propellant 30 418 kg vs 13 836 kg planned for the 36 tours). `insertion.py`'s pure-Lambert `build_tour` therefore
   cannot even evaluate the starting fleets; the ALNS needs the hybrid leg model of `search.expand`.
2. **A 2-state DP over the hybrid model reproduces the planner.** Per leg two variants exist — `L` (Lambert junction,
   `dv <= eta a tof`, propellant `m(1-exp(-kappa dv/VE))`) and `N` (LinLeg from the incoming state, coast miss
   <= 0.15 AU, `|u_k| <= 0.8 a`, propellant `m(1-exp(-cost/VE))`) — and the `N` variant's arrival velocity depends on
   the incoming velocity. Choosing the variant per leg by a 2-state dynamic programme (state = variant of the previous
   leg) gives, at the plan's m0, propellant <= `fuel_est` for 23 of 24 tours (equal to 1e-9 kg where the planner's own
   choice was already optimal, e.g. camp2 sc4 / sc12; 1-10 % lower elsewhere because the beam's choice was greedy),
   382 ms per 30-leg tour (236 ms for the greedy per-leg minimum). One tour (camp2 sc11) is infeasible at leg 10 in the
   model and needs a re-time repair before it can be a restart.
3. **Tank size and the model are coupled.** `fuel_est` in the pools was computed at m0 = 1000 kg; the MILP wrote
   m0 = 1086 (600 + 1.3 fuel + 20). Re-evaluating the same sequence at 1086 kg gives 1.086 x the propellant (mass chain
   is ~linear in m0) and a tighter eta rule (a = T/m): camp2 sc7 is feasible at 1100 kg and infeasible at 1211 kg. The
   ALNS must size tanks by a fixed point and check feasibility at the sized tank.
4. **Conversion calibration (82 `*_info.json`, 1 407 legs).** Non-launch legs, all regimes: sum actual / sum planned
   = 1.05. Patient tours (m0 <= 1300, >= 10 legs, 27 tours): per-tour ratio median 1.01, p75 1.13, p90 1.23 (the
   K12_m1200w_imp "s" family with dv > 1.2 legs is the 1.1-1.25 tail). Small planned legs convert worst (plan < 5 kg:
   3.3x median; 5-20 kg: 1.2-1.4x; > 20 kg: 0.85-0.95x), i.e. there is a per-flyby floor of roughly 5-8 kg. The launch
   leg costs 22.8 kg actual vs 2.3 planned (systematic +20 kg). Surrogate for tank sizing: `F_act = 20 + rho x F_plan(m0)`
   with rho = 1.15 initially (p75), later per-tour from conversion feedback.
5. **The TEAM file as a restart.** The 16 converted chains re-evaluated in the hybrid model with relaxed limits
   (dv_max 3 km/s, coast miss 0.3 AU) are feasible for 13/16 and priced within +-10 % of the actual propellant for 9 of
   those (SC3 350 vs 378 actual, SC9 472 vs 475, SC12 514 vs 512, SC16 461 vs 461); with the planner's limits
   (1.5 / 0.15 AU) only 2/16 are feasible. The pure Lambert chain was infeasible for 9/10 (docs/insertion.md 3.4).
   The TEAM fleet has 398 visits for 294 distinct targets: **88 duplicates = 3.5 craft worth of flybys wasted**.
6. **Timings** (single core): `evaluate_insertions` full grid (1 d, 9 shifts) on a 31-leg tour: 43 353 candidates in
   0.30-0.39 s; coarse (3 d, no shift): 1 511 candidates in 0.049 s; Lambert 14 us per solve (40k in 0.57 s);
   `build_tour` 4.4 ms; `mass_chain` 400 x 31 in 0.37 ms; LinLeg build 3.8 ms (one TOF), solve 0.15 ms per target,
   build over 77 TOFs 19 ms, solve 77 TOF x 300 targets 249 ms; all-asteroid daily position table (5 479 d x 300)
   0.7 s to build, 40 MB float64 for r (another 40 MB for v).
7. **Target prices.** LP duals of the pool set cover (`results/phasing/prices3.json`): 138, 47, 132, 177 at 1.0 (never
   covered cheaply), 41, 174, 11, 137 at 0.9-0.98; 56 targets above 0.5, 40 at 0. These are the natural weights of the
   target-related destroy operator and the shortlist of the 4-5 misses we should budget for.

## 2. The problem as a vehicle-routing problem

- Vehicles: K craft (K free, 9-16), each with a launch epoch t_L (free in [0, ~600 d]; later launches are allowed and
  used by tail craft), free launch v_inf <= 4 km/s, tank m0 <= 2000 sized by the plan, dry 600 kg.
- Customers: 298 targets (131/144 excluded) with *moving* positions; a visit is a (target, time) pair; a target may be
  visited by several craft (only the first counts) — overlap is allowed but wasteful.
- Route cost: propellant of the leg chain through the planned (target, time) sequence in the hybrid model; the cost
  enters J through the tank: `J_i = 1 + x + x^2`, `x = (m0_i - 600)/1400`, `m0_i = min(2000, 600 + R + rho F_i(m0_i))`
  (R = reserve 30 kg = 10 margin + 20 launch-leg offset).
- Objective: `J = sum_i J_i + (298 - |covered|) + 2` (the +2 are 131/144).
- Hard constraints (model validity, checked at the sized tank): per leg `dv <= min(dv_max, u_max eta a tof)` for
  Lambert legs (eta 0.6, dv_max 1.2 km/s, u_max 0.85 for legs created or modified by the ALNS — the planner's tours
  have legs at u = 1.0 and those are the ones the converter loses; existing legs keep their slack), LinLeg legs
  `|u_k| <= 0.8 a`, coast miss <= 0.15 AU, cost <= dv_max; leg duration >= 15 d (>= 20 d for N legs); every flyby
  <= T_MISSION - 5 d; m0 <= 2000 (else the move is infeasible); no target twice on the same craft.
- The J landscape is dominated by misses: dJ/dm0 = (1 + 2x)/1400 = 0.0011 per kg at x = 0.3, so **one miss is worth
  ~900 kg of tank**. In the model any *feasible* insertion beats a miss; feasibility (time windows, eta rule, dv cap,
  2000 kg cap) is the real constraint, and the fuel term decides *which* craft absorbs a target and how big its tank
  is. Consequently the operators that matter are the ones that *create feasibility*: opening windows (cluster
  destroy), re-timing neighbours, swapping a duplicated visit for an uncovered target, exchanging tails between craft.

## 3. Solution representation and evaluation

### 3.1 Data structures (`ctoc14/fleet.py`, new)

```
Leg   : ast (1..300), t (s), typ ('L'|'N'), tof, dv, dm, m_start, slack, u, r (3), v_in (3), v_out (3), coast_id
Craft : id, t_launch, vinf (3, derived), m0 (derived, fixed point), legs [Leg], fuel, J_i, feasible, locked (bool),
        rho (tour-specific conversion ratio, default global), n_unique (targets no other craft visits)
Fleet : crafts [Craft], cover (int16[300] visit counts), uncovered (set), J, n_miss, version counter,
        caches: insertion table {(craft_id, target): CandidateRecord | None, valid_version}, removal savings per visit
```
`Craft.to_json()` writes the `search.tour_json` format plus `leg_types` so the result directory is a drop-in input for
`tools/convert_fleet.py`; `Craft.from_json()` accepts pool routes / MILP selections (types by DP) and converted chains
(`validator.parse` rows: launch state row -> t_launch, v_inf, m0; Event=3 rows -> (ast, t)).

### 3.2 Leg evaluation and locality

`leg_options(r_prev, v_prev, r, tof, m, first)` returns the feasible variants `{L: (dv, dm, v_out), N: (cost, dm, v_out)}`
exactly as `search.roots/expand` price them (Lambert: `ctoc14/lambert.py`; N: one `LinLeg(r_prev, v_prev, [tof])` +
`solve`, i.e. ~4 ms; launch leg: Lambert only with the 4 km/s free v_inf).

Locality: an `L` leg's `v_out` does not depend on `v_in`; an `N` leg's does. So a change at flyby k (target, time or
removal) invalidates leg k (arrival), leg k+1 (departure) and then every consecutive `N` leg after k+1 up to and
including the junction into the next `L` leg (the "N-run"). With 37 % N legs the expected N-run is 0.6 legs, i.e.
3-4 legs re-evaluated per point change (~15 ms). `Craft.reeval_from(k)` re-runs the 2-state DP from leg k with the DP
state (fuel, v_out, m) of leg k-1 as the boundary until the DP path re-joins the stored path at an `L` leg with the same
variant *and* the mass difference is only a downstream rescaling; the downstream legs are then re-accumulated with
`insertion.mass_chain` (0.4 ms) instead of re-solved. `Craft.reeval_full()` (the DP of section 1.2) is used after
loading, after re-typing and as the assertion in tests (`reeval_from` == `reeval_full` to 1e-9 kg).

### 3.3 Tank fixed point and J_i

`size_tank(craft)`: `m0 <- 600 + R + rho * F(m0)` with `F(m0) ~ F(m_ref) * m0/m_ref` for the first guess, then two exact
mass-chain passes at the new m0 (0.4 ms each; converges to < 0.1 kg because dF/dm0 ~ F/m0 ~ 0.35 < 1 / rho). Leg
feasibility (eta rule, u_max) is re-checked at the sized m0; if a leg fails, the craft is infeasible (the move is
rejected). Fixed point > 2000 kg -> infeasible. `J_i = cost_sc(m0)`.

### 3.4 Move evaluation contract

Every move returns `(dJ, patch)`; `patch` lists (craft_id, new legs, new t_launch) for 1-2 craft. dJ is computed from
the re-evaluated craft(s) with sized tanks plus the change in misses. Applying a patch bumps `Fleet.version` and the
version of the touched craft; caches keyed by (craft version, gap) invalidate themselves.

## 4. Caches

1. **Daily ephemeris table** (`ctoc14/ephcache.py`, new): `R[d, k, :]`, `V[d, k, :]` for d = 0..5478 days, k = 0..299
   (2 x 40 MB float64, built in 1.4 s, shared by fork). Used only by prefilters and relatedness; candidate evaluation
   always calls `Body.state` at the exact times (vectorised, ~10 us per state).
2. **Coast tables per flyby state**: for each craft leg k the coast from `(r_k, v_out_k)` sampled daily over the next
   600 d (`propagate_twobody` over 600 dt: ~1 ms; 30 per craft; invalidated with the leg). The coast-miss matrix of one
   craft to *all* targets, `D[k, d, j] = |R[t_k + d, j] - C_k[d]|` (30 x 600 x 300 = 5.4 M distances, ~20 ms), is the
   single most useful prefilter: it tells, for every gap and every target at once, when the target passes within
   x AU of the craft's coast.
3. **Insertion table**: for each (craft, uncovered target) the best verified candidate (gap k, t_j, delta, dJ, variants)
   and the second-best craft's value (for regret). Entry valid while the craft's legs k-1..k+2 (and the N-run) are
   unchanged — implemented with a per-leg version stamp; a move touching 3 gaps of one craft invalidates ~10 % of that
   craft's column and nothing on other craft.
4. **Removal savings** per visit: dJ of deleting the visit (legs k-1 -> k+1 merged, N-run re-evaluated, tank re-sized);
   same invalidation rule.

## 5. The insertion screen: one vectorised pass per (craft, target)

Goal: find the cheapest feasible insertion of target j into craft c over all gaps k, intermediate times t_j and
neighbour shifts delta, at 10-20x the throughput of `insertion.evaluate_insertions` (0.35 s), then verify exactly.

**Stage A — gap prefilter (table lookups, per craft per move, shared by all targets).** From the coast-miss matrix: gap k
is a candidate for target j at day d if `D[k, d, j] <= min(0.3 AU, 2 dv_max (d - 0) )` (a transfer moving the craft by
distance D over duration tof needs at least ~D/tof of dv; the factor 2 is slack) and `t_k + d` lies inside the gap
(`t_k + 15 d <= t_j <= t_{k+1} + delta_max - 15 d`). Typically 3-6 of 30 gaps survive per target, with a few tens of
candidate days each.

**Stage B — Lambert screen ('L' variants for legs A = k->j, B = j->k+1, C = k+1->k+2 shortened, junction D).** This is
`evaluate_insertions` restricted to the surviving (gap, t_j) pairs on a coarse grid (3 d) and delta = 0 first: the
three Lambert calls and one `mass_chain` over ~300 candidates take ~10 ms. Reuse the function body (the assembly of
the modified DV/TOF sequences, the "unchanged legs may not get worse" rule and `u_mod`) with two changes: (i) the
candidate list comes from stage A instead of the full grid, (ii) legs downstream of the modified block that are `N`
are flagged "needs verify" instead of assumed unchanged.

**Stage C — coast screen ('N' variants).** From stage A the days where `D[k, d, j] <= 0.15 AU` are exactly the
candidates where leg A can be an `N` leg (the target comes to the craft): build one `LinLeg(r_k, v_out_k, tofs)` over
those few TOFs (5-20 ms) and solve for j. For each such A candidate, leg B as `N` requires the *continued* coast from
j (arrival velocity ~ coast velocity) to pass near k+1: `D_B = |R[t_{k+1}+delta, k+1] - coast_k(t_{k+1} + delta)|`
(again from the same coast table because an N-A leg barely changes the orbit) -> if <= 0.15 AU solve one more LinLeg
(4 ms). For A = 'L', B = 'N' (Lambert to j, then coast to k+1): needs a coast from `(r_j, v2A)`; use a batched
universal-variable propagator `propagate_batch(r0[N,3], v0[N,3], dt[N])` (new, ~1 ms per 1 000 states) on the
top-30 Lambert candidates only. Combined the N-screen costs ~30-60 ms per (craft, target) when candidates exist,
~0 otherwise.

**Stage D — refine + verify (top-3 per (craft, target)).** Around each of the best three (k, t_j) from B/C: 1 d grid
+-4 d, then 0.25 d +-0.75 d (as `search.expand`'s `tof_refine`), delta in {-10..+10 d step 2.5} on the shifted
neighbour — 250 candidates per refinement, one vectorised Lambert pass (~5 ms) — then the exact sequential hybrid
re-evaluation `reeval_from(k)` of the best candidate with sized tank (~20 ms). The verified dJ goes into the insertion
table. Total per fresh (craft, target) pair: **~50-80 ms** vs 350 ms today, and only invalidated pairs are recomputed.

Removal (`remove_visit(c, k)`): legs k-1 -> k+1 merged into one leg of tof = sum; both variants tried; N-run
re-evaluated; if the merged leg is infeasible the removal is still allowed *inside a destroy* (the repair must fix it)
but not as a stand-alone improvement.

## 6. Destroy operators (q visits removed, q ~ U[4, min(40, 0.12 x visits)])

| name | selection | purpose |
|---|---|---|
| random | q visits uniformly over all craft | diversification |
| worst-cost | rank visits by removal saving dJ (cache), pick with the Ropke-Pisinger rank^p rule (p = 3) | frees the most propellant |
| time-cluster | pick a craft (prob ∝ 1/n_unique) and a window [t0, t0 + W], W ~ U[0.7, 2] yr; remove every visit of that craft inside it | opens a free window so the repair can re-sequence 2-5 targets in a different order; the only way to insert targets that need a different order |
| related (Shaw) | seed = random uncovered or visited target; relatedness r(i, j) = min over i's visit times of |Δt|/120 d + |Δr|/0.15 AU (daily table); remove the q most related visits across craft | regroups targets that are close in space-time onto one craft |
| whole-craft | remove *all* visits of one craft (prob ∝ 1/n_unique, never a locked craft); the craft is deleted, K -= 1 | consolidation: the repair reinserts what it can, the rest become misses; SA on J decides |
| make-room | pick an uncovered target j; from the coast-miss matrix find the (craft, gap) pairs with the smallest miss to j at any day; remove the 1-2 visits bounding those gaps | the operator that hunts the 16 current misses directly |
| overlap | remove every duplicated visit from the craft where it costs the most | deterministic improvement; also run after every accepted move |

## 7. Repair operators

1. **Greedy cheapest insertion**: repeat: for every uncovered target take the best verified insertion over all craft
   (insertion table; fresh pairs computed with section 5); apply the globally cheapest by dJ; with noise (cost x U[0.9,
   1.1]) as a second variant. Targets with no feasible insertion remain uncovered (miss).
2. **Regret-2 / regret-3**: regret_j = sum over the 2nd (3rd) best craft of (dJ_best_other - dJ_best); insert the target
   with the largest regret first; a target feasible on exactly one craft gets regret = +inf and is inserted first —
   this is what protects the rare-window targets (duals ~1.0) from being crowded out.
3. **Window re-sequencing** (after a time-cluster destroy): the free window [t_a, t_b] between fixed flybys a and b on
   craft c is refilled by a small beam search from the state at a over the uncovered targets *plus* the removed ones
   (`search.beam_search_from_state`, beam 60, depth <= 6, tofs <= 400 d, `visited_ids` = everything the craft
   already visits and, with a rarity bonus, targets other craft cover), truncated at t_b - 15 d, then the junction into
   b is checked with `leg_options`; the best complete refill by dJ is applied. 1-3 s per call — used with lower base
   weight than 1-2.
4. **Create craft from the pool**: pack the 335 k pool routes as 300-bit masks once (`select_routes_milp.py` does
   this in 3 s); score = |mask ∧ uncovered| - cost_sc(tank(fuel_est x 1.15 x 1.086)); take the best 5, load them
   (`Craft.from_json` + DP), delete their already-covered visits (removal moves), re-time, size the tank and add the
   craft if dJ < 0 (the marginal-cost criterion of brief §41: K_new > c). Also used at restart from a 10-craft MILP.
5. **Tail exchange (2-opt*)** between craft A and B: for all (k, l) with |t_k - t_l| <= 300 d evaluate the junction
   A[:k] -> B[l+1:] and B[:l] -> A[k+1:] (two Lambert calls over the (k, l) grid, ~900 pairs = 30 ms) then verify the
   best 3 with the N-run rule; apply if dJ < 0 or by SA. Merges craft (route merge of brief §42) and moves whole tails
   between them.
6. **Post-repair local improvement** (always): overlap removal -> re-time pass on the modified segments (coordinate
   descent of the flyby times k-2..k+2, +-4 d at 1 d then +-0.75 d at 0.25 d, two Lambert arcs per shift vectorised,
   N legs re-solved; ~0.1 s per flyby) -> launch-epoch shift (t_L +- 30 d grid 5 d: only the first junction changes,
   one vectorised Lambert call) -> re-type DP of the touched craft -> tank fixed point. Full re-time sweeps of every
   craft (2-3 s each) run every 50 accepted moves and on every new best.

## 8. Craft deletion (the consolidation move)

`dissolve(c)`: remove craft c (whole-craft destroy), reinsert its targets with regret-2 into the other craft (their
insertion tables already contain most of the answers), then run the local improvement. dJ = -J_c + sum of the
insertion dJ's + (targets left uncovered) x 1. Accept by the SA rule. The move is proposed with probability
proportional to `max(0, c_J - n_unique x p_reinsert)` where p_reinsert is the running fraction of reinsertion successes,
so it concentrates on craft whose unique-target count is close to their cost (currently sc12: 18 flybys, J_i 1.36, and
the TEAM craft with 3, 16, 16 flybys). Conversely `spawn` (repair 4) adds a craft; K floats between 9 and 16.

## 9. Acceptance, adaptive weights, schedule

- Simulated annealing on J: accept if `dJ < 0` or with `exp(-dJ/T)`. Temperature relative to one miss: T0 chosen so a
  +0.25 move (a quarter miss, or +220 kg of tank) is accepted with p = 0.5 (T0 = 0.36); exponential cooling to
  T_end = 0.01 over the chain budget; reheat to 0.5 T0 on restarts. Because misses dominate, also record the incumbent
  under the lexicographic order (misses, then sum J_i) — the chain reports both.
- Operator weights (Ropke & Pisinger): scores 33 (new global best), 9 (better than current), 13 (accepted, worse),
  0 (rejected); weights updated every 100 iterations with reaction factor 0.1; roulette wheel over destroy and repair
  separately; minimum weight 0.05 so no operator dies.
- Feasibility is hard except inside a destroy-repair pair (an infeasible merged leg after a removal is allowed until
  the repair ends; the pair is rejected as a whole if the repair cannot restore feasibility).

## 10. Restarts (initial fleets)

1. `results/phasing/milp_u123cg_n12` (12 craft, 284 covered, 19 duplicates, planned J 33.9 at the MILP tanks;
   re-priced by the DP with fixed-point tanks) — the main start.
2. `camp1/milp_n12`, `camp2/milp_n12` (12 craft each; camp2 uses dv_max 1.5, so its legs above 1.2 km/s are marked
   "fragile": allowed to stay, never created).
3. TEAM-16 converted chains (`results/CTOC14_Result_TEAM.txt` via `validator.parse`): loaded with the relaxed limits
   (dv_max 3, 0.3 AU) as *locked-converted* craft (sequence fixed, only their coverage counts, J_i from the file's
   m0), the 3 model-infeasible ones included as locked. The ALNS then works only on unlocked craft and on unlocking:
   an `unlock` move re-times a locked craft into the planner limits (re-time sweep with dv_max 1.5 / 0.15 AU); if it
   succeeds the craft becomes an ordinary one. This start has 88 duplicates to harvest.
4. `select_routes_milp.py --max-craft 10` over the union pool (10 craft, ~25-30 misses): the "fill" start — the ALNS
   must insert misses, not delete craft.
5. Random: greedy pick of 10-11 pool routes with >= 28 flybys minimising overlap (12 558 routes with >= 25 flybys,
   618 with >= 30).
Every restart first runs `repair_model()`: DP typing, full re-time sweep, tank fixed point; legs still infeasible are
removed (the target goes to the uncovered set) so that the chain starts from a feasible fleet.

## 11. Parallelisation

Independent chains (`tools/run_alns.py --seed s --init X --hours h`), one process each, 4-6 concurrently on this
machine (10 CPUs, load 7.6 shared with ten agents; never more than 6; never detached). Chains differ in seed, start
fleet, T0 (x0.5, x1, x2) and destroy-size range. Exchange: every 500 accepted moves a chain writes its best fleet to
`results/alns/exchange/<chain>.json`; every 2 000 moves it reads the global best and, with p = 0.5, restarts from it
(reheat); with p = 0.25 it imports only the best craft from the global best for the target subset it covers worst
(craft-level path relinking: the imported craft's targets are removed from the chain's own craft, which are then
re-timed). A `results/alns/leaderboard.jsonl` records every new best (planned and, when known, validated J).

## 12. Interface to exact conversion (`ctoc14/lowthrust.py::convert_tour`, `tools/convert_fleet.py`)

- **When**: (a) a new planned global best that improves the last converted best by >= 0.5 and is >= 60 min old
  (avoid converting every transient); (b) on demand at the end of a chain; (c) per craft: any craft whose sequence
  has been stable for 2 000 moves and whose rho is still the global default. Conversion runs in a separate pool of
  <= 3 processes (5-15 min per 30-leg tour), never inside the chain process.
- **Feedback** (`tools/alns_feedback.py`, reads `<frag>_info.json`): per craft rho_c = actual / F_plan(m0) (launch leg
  excluded), stored with the craft's sequence hash; per-leg ratios update a global table (bins by planned dm, dv, tof,
  leg type) that replaces the constant rho by `rho(dv, tof, typ)` for tank sizing — this is the "re-fit the fuel
  surrogate" loop, and it should be rerun after every fleet conversion. Dropped flybys mark their (craft, target)
  pair as *fragile*: the ALNS adds a penalty of 0.5 to keeping that visit unchanged (so it prefers to move or re-time
  it) and the next conversion of that craft uses `widen_on_fail` and `adaptive_iter` (already default) plus
  `allowed_extra` = the chain's current uncovered targets (the converter's re-plan may pick them up).
- **Locking**: a craft whose conversion kept every flyby with rho_c <= 1.05 is locked (only insertions with verify are
  allowed on it, no destroy); locked craft keep their converted fragment file. The final submission is assembled with
  `tools/select_frags_milp.py` over *all* fragments (new + the 110 old ones + TEAM), so every conversion adds a
  column and the validated J can only improve. Unlocking happens when the chain finds a planned improvement worth
  >= 0.3 on that craft.

## 13. Runtime budget (single core, measured components)

Per move (destroy q = 8, repair by regret-2 into K = 11 craft):
- destroy: 8 removals x 15 ms (re-evaluation) + cache updates = 0.15 s;
- repair: pairs to evaluate = fresh (craft, target) pairs (~30 after cache reuse) + re-evaluation of the modified craft
  for the remaining targets (7 + 6 + ... = 28) = ~60 pairs x 50-80 ms = 3-5 s; with the pure
  `evaluate_insertions` at 0.35 s per pair and no cache it would be 116 x 0.35 = 40 s;
- verify + refine of the 8 applied insertions: 8 x 30 ms = 0.25 s;
- local improvement on modified segments: 8 flybys x 0.1 s = 0.8 s; tank fixed points negligible.
Total ~4-6 s per move -> **1 000 moves in 70-100 min per chain**; with the coarser setting (5 d grid, 0.2 AU
prefilter, verify top-1) ~2 s per move -> 35 min. Six chains in parallel: ~8-10 k moves per hour fleet-wide, i.e.
25-50 k moves per chain overnight, which is a normal ALNS budget. Reduction levers in order of value: insertion cache
(10x), gap prefilter (5x), coarse-then-refine grid (7x), delta = 0 in the screen (9x), verify top-3 only.
Memory: 80 MB tables + ~1 MB per fleet; fork-shared.

## 14. Expected J (arithmetic)

Tank fixed point with R = 30 kg, rho = 1.15, F(m0) = F_1000 m0/1000: m0 = 630 / (1 - 1.15 F_1000/1000).
- F_1000 = 360 kg (31 flybys, today's best single-craft plans): m0 = 1 075, J_i = 1.454.
- F_1000 = 330 kg (30 flybys x 11 kg after 0.25 d re-timing): m0 = 1 015, J_i = 1.384.
- F_1000 = 300 kg: m0 = 962, J_i = 1.325.
Fleet scenarios (298 targets, 131/144 always missed, i.e. +2):
- **Conservative** (fill only, 12 craft kept, 16 -> 6 misses at +30 kg per absorbed target): sum J_i 17.9 + 12 x 0.035
  = 18.3, J_plan = 18.3 + 6 + 2 = 26.3; validated (rho 1.1 on tanks already sized with 1.15, 1-2 lost flybys) ~27-28.
  Already beats 30.789.
- **Target** (consolidation to 10 craft x ~30 visits = 300 visits, <= 6 duplicates, 4 misses): 10 x 1.45 + 4 + 2 =
  20.5 at 360 kg; 10 x 1.38 + 6 = 19.8 at 330 kg. Validated 20-22 if conversion is at the median (1.01) and no flyby is
  lost; each lost flyby that is not a duplicate costs +1.
- **Stretch** (11 craft x 27 visits at 300 kg, 3 misses: the leaders' 9-10 craft regime is 33 per craft, beyond our
  planner): 11 x 1.325 + 3 + 2 = 19.6 planned; or 10 craft at 330 kg with 3 misses: 18.8 planned, ~19.5-20.5 validated.
- The pool LP bound (21.1 at 12 fractional craft, with fuel_est at 1000 kg and margin 1.3) is *not* a bound for the
  ALNS because the ALNS creates tours outside the pool; the true floor is the 4 hard targets (138, 47, 132, 177 with
  dual 1.0) + 2 = 6 misses unless a craft is dedicated to them (K > c criterion: 4 targets for 1.3 is worth it only if
  a 4-target tour exists — the TEAM sc1 with 3 flybys for J_i 1.35 shows that such tail craft are currently bad deals).
Honest expectation: planned J 19-22 after a full campaign (E4), validated 20-23; < 20 validated needs the stretch
allocation *and* lossless conversion. Every step of the way produces validated files that beat 30.789.

## 15. Implementation plan (new files only; existing modules are imported, not modified)

1. `ctoc14/ephcache.py`: `DailyTable` (R, V), `coast_table(r, v, t0, ndays)`, `propagate_batch(r0, v0, dt)` (vectorised
   universal-variable Kepler over N states; validated against `kepler.propagate_twobody` to 1e-6 km).
2. `ctoc14/fleet.py`: `Leg`, `Craft`, `Fleet`, `leg_options`, `Craft.reeval_full` (2-state DP — the prototype in the
   scratchpad `alns/hybrid_dp.py` is the reference), `Craft.reeval_from(k)`, `size_tank`, `Fleet.J`, JSON/submission
   I/O (`Craft.from_json`, `Craft.from_rows`, `to_json` in `search.tour_json` format), `repair_model`.
3. `ctoc14/moves.py`: `gap_prefilter` (coast-miss matrix), `screen_lambert` (refactor of
   `insertion.evaluate_insertions` taking an explicit candidate list; keep `mass_chain`, the DV/TOF assembly and the
   `u_mod` rule), `screen_coast` (N variants), `refine_verify`, `remove_visit`, `retime_segment`, `retime_sweep`,
   `shift_launch`, `tail_exchange`, `insertion_table` (cache with version stamps), `removal_savings`.
4. `ctoc14/alns.py`: operators of sections 6-8, `Chain` (SA, weights, schedule, checkpoints every 100 moves to
   `results/alns/<chain>/state.json`), `exchange` (section 11), pool loader for repair 4 (reuse the bitmask packing
   of `tools/select_routes_milp.py`), window re-sequencing via `search.beam_search_from_state`.
5. `tools/run_alns.py`: CLI (`--init {u123cg|camp1|camp2|team|milp10|random} --seed --hours --moves --T0 --exchange-dir
   --convert-every --nproc-convert --lock-ratio 1.05 --coarse`), logging of per-move timings and operator stats.
6. `tools/alns_feedback.py`: convert the incumbent with `convert_fleet.py`, read `*_info.json`, update
   `results/alns/calibration.json` (per-craft rho, per-bin rho), mark fragile visits, lock craft; final
   `select_frags_milp.py` over `results/alns/**/frag_*` + old fragment dirs.
7. `tests/test_fleet.py`: (a) DP fuel <= `fuel_est` + 1e-6 at the plan m0 on the 24 camp tours and == to 1e-9 on camp2
   sc4/sc12; (b) `reeval_from` == `reeval_full` after 1 000 random insert/remove/re-time moves; (c) screen candidates
   agree with `insertion.evaluate_insertions` on Lambert-only tours (`results/fleet_b300`, the existing
   `tests/test_insertion.py` cases, e.g. 137 into sc4 for 30.26 kg); (d) `propagate_batch` vs `propagate_twobody`;
   (e) tank fixed point monotone and < 0.1 kg.
Order: 1-2 (+ test a, b, e) -> 3 (+ test c, d) -> E1 -> 4-5 -> E2 -> 6 -> E3/E4.

## 16. First experiments (pass/fail)

- **E0 evaluator**: tests (a)-(e) pass; full DP <= 0.5 s per 30-leg craft; `reeval_from` <= 20 ms median.
- **E1 fill** (no craft deletion, start = u123cg 12 craft, 2 000 moves, 1 chain, ~2 h): pass if planned misses
  16 -> <= 8 and sum J_i <= 19.5 (planned J <= 29.5 from 33.9); also record how many of the 16 have *any* feasible
  insertion in the hybrid model (expect >= 12; the pure Lambert model found detours for 9 of 22 tail targets).
- **E2 consolidate** (craft deletion + tail exchange on, 10 k moves, 3 chains x 4 h): pass if a fleet with <= 11 craft,
  <= 6 misses and planned J <= 24 appears; report duplicates (target <= 8) and per-craft visits (target >= 27).
- **E3 conversion** (convert the E2 fleet, nproc 3, ~1.5 h): pass if >= 95 % of planned flybys survive, actual
  propellant <= 1.15 x F_plan(m0) for >= 9/10 craft, and `select_frags_milp` over all fragments validates at
  J <= 26 (< 30.789 is the minimum bar; the TEAM fragments remain in the MILP as a hedge).
- **E4 campaign** (6 chains x 8 h with exchange, conversion feedback every new best >= 0.5, locking): target validated
  J <= 22; stretch < 20. Report the operator-weight history (which neighbourhoods produced the bests).
- **E5 runtime**: per-move timing histogram; pass if the median move <= 6 s in the default setting and <= 2.5 s in
  `--coarse`; insertion-cache hit rate >= 70 %.

## 17. Risks and mitigations

1. Hybrid model optimism on N legs near the 0.15 AU limit (0.35 AU coast miss cost 27 kg actual vs 7 planned):
   keep lin_drmax 0.15 AU for new legs and rho(typ) from feedback; fragile marking after conversion.
2. Tank/feasibility coupling: a move that grows the tank can break a tight leg elsewhere on the craft — always
   re-check at the sized tank (section 3.3); prefer u_max 0.85 for new legs.
3. The Lambert screen cannot see insertions whose feasibility relies on re-ordering neighbours; the time-cluster
   destroy + window re-sequencing repair covers that, at 1-3 s per call — keep its weight adaptive.
4. Miss-dominated objective makes the SA either frozen (T small: only insertions accepted) or noisy (T ~ 1: misses
   accepted freely): the T0 = 0.36 / lexicographic incumbent choice above is a guess; E2 must tune it (try 0.15, 0.36, 0.7).
5. Conversion throughput: one fleet conversion is 1-2 h on 3 processes; converting every planned best is impossible,
   hence the >= 0.5 / 60 min rule and the per-craft stable-sequence rule.
6. Hard targets (138, 47, 132, 177, 41, 174, 11, 137): a floor of 4-6 misses is likely; a dedicated 4-6 target tail craft
   is only worth it if K > 1.3 — leave to the create-craft repair and the SA.
7. Data inconsistencies in the pool (camp2 sc11 infeasible at leg 10 in the model; legs above dv 1.2 in camp2): the
   restart `repair_model` step drops or re-times them; report what was dropped.
8. CPU contention with ten other agents: chains are single-process; use `--nproc-convert <= 3`; checkpoint every 100
   moves so a killed chain resumes.
9. Duplicate accounting in the final MILP: `select_frags_milp` counts a target once; the ALNS's overlap operator keeps
   planned duplicates near zero, so planned J and the MILP's J agree.

## 18. Existing code reused

- `ctoc14/search.py`: leg formulas (`roots`, `expand`), `beam_search_from_state` (window re-sequencing repair),
  `tour_json` format, `Params` (eta, dv_max, lin_* parameters must be identical).
- `ctoc14/insertion.py`: `mass_chain` (vectorised mass chain), the DV/TOF assembly and rules of `evaluate_insertions`
  (refactored into `screen_lambert` with an explicit candidate list), `min_m0_for` (tank bisection fallback),
  `apply_insertion`/`build_tour` for the Lambert-only unit tests.
- `ctoc14/linleg.py` (`LinLeg`), `ctoc14/lambert.py`, `ctoc14/kepler.py` (`Body.state`, `all_ast_states`,
  `propagate_twobody`), `ctoc14/constants.py` (`cost_sc`).
- `tools/select_routes_milp.py`: bitmask packing of the pool (create-craft repair), `--max-craft 10` restart.
- `tools/select_frags_milp.py`, `tools/convert_fleet.py`, `tools/convert_one.py`, `ctoc14/validator.py` (`parse` for the
  TEAM restart; `validate`), `ctoc14/lowthrust.py::convert_tour` (`allowed_extra`, `global_excluded`).
- `tools/price_targets.py` / `results/phasing/prices3.json`: target weights for the related and make-room destroys.
- Scratchpad prototypes (reference implementations, not project files): `alns/hybrid_eval.py` (greedy hybrid chain),
  `alns/hybrid_dp.py` (2-state DP; produced the numbers of section 1.2 and 1.5).
