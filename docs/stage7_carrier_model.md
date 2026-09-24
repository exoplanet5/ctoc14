# Stage 7 — a new model: carrier orbits + exact phase-plane DP (2026-09-20 evening)

t10d (`results/CTOC14_Result_t10d.txt`, raw J 14.366, VALID) and all stage 1-6 tools are untouched.
New code: `tools/exp_carrier.py`, `tools/carrier_dp.py`; results under `results/s7/`.

## 1. What the leaderboard really says (CORRECTED 2026-09-20 late)

The board's N column counts the two permanent misses (131 at >= 10.7 AU and 144 at >= 33.8 AU all mission long) as
craft: our own row shows N = 12 for t10d, which has 10 craft. With real N = board N - 2 and the time coefficient
removed (k = 1 - 0.1 (t_end - t_submit)/28 d):

| team | shown | board N | real N | raw J | tank | flybys per craft | km/s per flyby |
|:--|--:|--:|--:|--:|--:|--:|--:|
| HIT | 11.936 (09-15) | 10 | **8** | 12.52 | 952 kg | 37.2 | 0.48 |
| THU-LAD | 12.016 (09-20) | 9 | **7** | 12.37 | 1097 kg | 42.6 | 0.55 |
| NUAA | 12.296 (09-13) | 10 | **8** | 12.97 | 1004 kg | 37.2 | 0.54 |
| us, t10d | 13.802 (09-17) | 12 | **10** | 14.37 | 877 kg | 29.8 | 0.50 |

(The check: our row reproduces t10d's 877 kg / 0.50 km/s exactly.)  An earlier version of this section read the board
N literally and concluded that the leaders fly 10 craft on ~70 kg each. That was WRONG. The leaders fly 8 and 7 craft
at OUR per-flyby price; what they have and we lack is 37-43 disjoint flybys per craft. Stage 6's target (N = 8,
<= 1098 kg) was the right one, and its wall -- routes that want the same targets, 2:1 forcing cost -- is the real
problem. The carrier model below therefore has to earn its place as a way through THAT wall, not as a new regime:

- it supplies ~20 flybys per craft for ~5 km/s, leaving ~13 of an 18 km/s budget for the other ~17;
- its routes are built from a different principle (orbit geometry, not time-greedy Lambert legs), so they prefer
  different targets than the beam's routes -- the "99 targets nobody wants" are exactly what a second generator is for;
- its pricing is exact and instant, so dual prices can steer it without the beam-pruning distortion of stages 4-6.

## 2. The model

**Carrier orbit.** Launch v_inf (4 km/s, free) buys each craft an orbit: a ~ 1 AU, e <= 0.13, i <= 7.7 deg, any
orientation. Call its fixed e-vector and i-vector the carrier. A craft that never changes them only steers its PHASE,
and phase is nearly free: a period offset costs v/3 per unit of lag slope = 0.028 km/s per deg/yr (the K_DRIFT
measured in `ctoc14/phasemodel.py`).

**Which asteroids are cheap for a carrier is pure geometry.** Two orbits meet only at their mutual node line; the
asteroid is cheap iff the radial gap there, r_ast(u) - r_car(u), is ~0. That is one number per (carrier, asteroid,
node), computed analytically for 100 000 carriers x 298 asteroids in seconds (`tools/exp_carrier.py`):

| eps (node gap) | Earth's own orbit | best of 100k launch-reachable carriers | 10 carriers (greedy) cover |
|--:|--:|--:|--:|
| 0.003 AU | 13 | 24 | 161 / 298 |
| 0.010 AU | 30 | 53 | 254 / 298 |
| 0.020 AU | 62 | 87 | - |

**Events.** For a candidate asteroid every passage through that node is an event (t, lag): the craft must be at the
node `lag` days later than its ballistic schedule. A route is a piecewise-linear lag curve from (t_launch, 0) through
events; cost = (v/3) x total variation of the slope.

**Gap-slope coupling (the piece `phasedp.py` missed).** A lag slope s IS a semi-major-axis offset, so the whole track
sits (2/3) s a further out. The gap an event pays is gap0 - (2/3) a s, a function of the arriving slope. This is
Markov on the line graph, so it costs nothing to model, and it turns a nuisance into a control.

**Exact DP, no grid, no beam.** Nodes = events (1000-2500 per carrier), arcs = event pairs with |slope| <= 0.08,
state = incoming arc (so the slope is exact), value = flybys - lam x dv. 0.02-0.05 s per carrier; 60 000 carriers
ranked by true route value in 50 s on 8 cores. The old phase DP failed because it granted eccentricity slack
independently at every event and let the phase slip on a 10 d grid; here the carrier's e/i are fixed by construction
and there is no grid.

**Settling.** `from_tour`'s Lambert chain is degenerate on these routes (legs are near-integer revolutions). New
seed (`carrier_dp.seed_ip`): ballistic carrier + the DP's tangential phasing impulses, no flybys; then ALL flybys
enter through one joint aim-point homotopy (`globalopt.insert_block`) and the usual `run_ialns.settle`. 40 of 40
calibration routes settled, 4-15 s each.

## 3. Evidence (impulsive twin, the same judge as every earlier stage)

| route | flybys | twin tank | dv | km/s per flyby | J_i |
|:--|--:|--:|--:|--:|--:|
| cheap end (lam 1.0), typical | 11-15 | 619-642 kg | 1.1-2.5 | **0.10-0.17** | 1.014-1.03 |
| best DP-ranked carrier | **21** | **689 kg** | 5.35 | **0.255** | 1.068 |
| other DP-ranked | 15-19 | 659-711 kg | 3.6-6.6 | 0.24-0.39 | 1.04-1.09 |
| deep end (lam 0.3) | 22-24 | 855-1033 kg | 14-21 | 0.58-0.88 | 1.22-1.41 |
| t10d average, for reference | 29.8 | 875 kg | 14.6 | 0.49 | 1.237 |

Calibration over 40 settled routes: twin dv = 1.48 x phase term + 6 km/s/AU x gap term, corr 0.82. The model is
honest for gaps <= 0.015 AU; at eps 0.03 it is 2x optimistic (large gaps are not linear) -- keep eps small.

What this shows: very cheap flybys (0.10-0.25 km/s) exist in this catalogue and we can now generate them at will:
nothing in six stages produced 21 flybys for 89 kg of propellant. What it does not show yet: depth. One fixed carrier saturates near
1.4 cheap flybys per year; past ~21 the DP starts buying expensive phase swings. In J per flyby a 21 @ 689 route
(0.051) is still behind t10d's 30 @ 875 (0.0415), and the real target is 37 @ <= 1000 kg (0.035). A carrier route is
half a route; section 4 is about the other half.

## 4. Plan (each step has a number that kills or keeps it)

1. **Carrier optimisation, continuous.** So far carriers are random samples, ranked. The DP is a 0.03 s objective
   in 6 variables (launch date, v_inf vector) plus launch year: run CMA-ES / pattern search from the top 200 samples.
   Gate: a single-carrier route with >= 25 flybys at <= 720 kg (leaves >= 11 km/s for the remaining 12).
2. **Full encounter geometry instead of the node only.** Low-inclination asteroids run alongside the carrier track for
   weeks; the node-gap test sees one event where there are several. Replace the node by the local minima of the
   track-to-orbit distance (MOID-style, still analytic/vectorised). Expect +20-40 % events on exactly the asteroids
   (i < 5 deg) that dominate the cheap lists.
3. **Carrier hop.** Route = 2-3 carrier segments; a hop is one e/i-vector move (phasemodel metric: (v/2)|de| + v|di|,
   ~1-2 km/s for a neighbouring carrier). DP state gains a carrier index from a library of the ~500 best carriers;
   hop arcs carry the metric cost. Two 7-year segments at 0.15-0.25 km/s per flyby plus one hop = ~30 flybys for
   ~7-8 km/s = ~730 kg; the N = 8 target needs a third segment or step 4 on top (37 flybys within ~18 km/s).
4. **Deepen by insertion, not by search.** Carrier routes are sparse in time (one flyby per 8-10 months), which is
   the opposite of the time-full routes on which insertion cost 1.6-9 km/s. Re-measure `insert_block` / ialns relocate
   on a 21 @ 689 route. Gate: 21 -> 37 flybys at <= 1000 kg on one route (the N = 8 line); +8 targets for <= +100 kg is the first checkpoint.
5. **Fleet layer = real column generation.** Pricing is now EXACT and takes 0.03 s, so duals can be used as node
   prizes without the beam-pruning distortion that killed every prize scheme in stages 4-6 (there is no beam).
   Master = set partitioning over (carrier, route) columns with N fixed (existing `colgen.py` RMP). This also answers
   stage 6's wall directly: "8 routes want the same 132 targets" is a statement about greedy independent fills; an
   LP with exact pricing either finds the disjoint cover or proves (dual bound) that it does not exist in this model.
6. **Hybrid bank.** Even before 1-5 land: replace t10d's most expensive routes by carrier routes wherever a set-cover
   MILP over {t10d routes, carrier routes, fragments} lowers J (`tools/select_routes_milp.py` exists). Export through
   the existing twin -> exact SCP -> validator2 path. Anything that validates below raw 13.86 is a new bank.

Clock: the coefficient costs 0.05 J per day; freeze 09-27 12:00. Steps 1, 2, 4 are hours each and decide whether
3 and 5 are worth building.

## 5. Reproduce

    tools/exp_carrier.py --n 100000 --eps 0.01                       # geometry screen, greedy cover
    tools/carrier_dp.py --n 60000 --eps 0.012 --pool 60000 --top 8 \
        --lam 1.0 --kphase 1.5 --kgap 10 --settle --out results/s7/x.json   # DP-rank all carriers, settle the best

`--kphase 1.5 --kgap 10` are the twin-calibrated cost weights; `--static-gap` switches the slope coupling off;
`--lambert-seed` tries the old seed first.

## 6. Results of step 4 and the first fleets (2026-09-20, 21:15-22:55)

**Step 4 gate: PASSED.** `tools/carrier_grow.py` (closest approaches -> parallel `insert_homotopy` trials -> accept the
cheapest, 250 d apart -> joint re-settle) took the 21 @ 689 kg carrier route to **39 flybys @ 975 kg in 215 s** and to
**40 @ 992 kg** (J_i 1.358, 0.034 J per flyby; t10d 0.0415, HIT 0.0353). Insertion into a SPARSE route costs 4-20 kg
per target for the first dozen and 10-50 kg later, against 100-280 kg on the time-full routes of stage 5. It even took
hard target 199 cheaply. A single N = 8-grade route now costs five minutes, from nothing, with no beam.

**The fleet wall is still there.** Three fleet builders, all 8 craft, all at the N = 8 tank level:

| builder | covered | sum J_i | per craft |
|:--|--:|--:|:--|
| sequential (`carrier_fleet.py`): DP base on uncovered targets, grow, next | 228 | 10.13 | 40 37 32 32 26 20 21 20 |
| + leftover pass (`carrier_pass.py`, step limit 90 kg) | 249 | 11.56 | leftovers cost 55-70 kg each |
| parallel regret growth (`carrier_regret.py`) from the 8 disjoint bases | 237 | 11.21 | 35 33 28 32 25 27 25 32 |

The regret run stalls at 210 under +45 kg per insertion (usable trials fall from 37/48 to 1/48) and reaches 237 at
+90 kg. That is the stage 6 number again (fill2 241, PART 234), reached in 22 minutes instead of days, by a generator
that shares no code with the beam. So the wall is a property of the PROBLEM AS POSED TO INDEPENDENT ROUTES, not of any
search: about 7 % of the catalogue is cheap for a given route, 8 routes chosen for their own value give each target a
1 - 0.93^8 = 44 % chance of a cheap home, and 100 + 0.44 x 200 = 188-210 is what we measure.

**What follows.** Coverage must be designed into the choice of the 8 carriers, before any growth:
1. library of ~2000 DP carrier routes (no exclusions); for each, its NEIGHBOURHOOD = targets whose closest approach to
   the seed trajectory is < ~0.06 AU (cheap to insert; to be calibrated from logged trial costs vs distance);
2. choose 8 carriers by a transportation MILP: y_c binary, x_tc <= y_c continuous, sum_t x_tc <= 40, sum_c y_c = 8,
   maximise covered targets (then minimise distance) -- 8 routes whose neighbourhoods TILE the catalogue;
3. build the bases with the MILP's assignment as each route's pool, then `carrier_regret.py`.
The cost of this experiment is ~2 h; its kill number is the MILP optimum itself: if no 8 neighbourhoods of <= 40 tile
>= 290 targets at 0.06 AU, this model cannot reach N = 8 either and the honest conclusion is that the leaders hold
something structurally different (or accept tanks well above ours on a few craft).

## 7. The tiling experiment (2026-09-20, 23:00-23:45): NEGATIVE

Tools: `tools/carrier_tile.py` (library + MILP), `tools/carrier_build.py`, `carrier_regret.py --pools`.

**Calibration (insertion cost vs closest approach to the route, 58 trials on two bases).** < 0.03 AU: 3-6 kg;
0.03-0.06: median 14 kg; 0.06-0.10: median 41 kg; beyond 0.10: median 35-50 kg with a 100-290 kg tail. So "cheap"
= within 0.06 AU, as assumed.

**Library.** 60 000 carriers DP-valued in 61 s; the top 3944 kept with their base (median 13, max 22) and their
neighbourhood from the unsettled seed trajectory (median 27 targets within 0.06 AU, max 50).

**Tiling MILP (8 carriers, cap 40 each).**

| delta | greedy union, no cap | MILP tiling |
|--:|--:|--:|
| 0.04 AU | 203 | 205 |
| 0.06 AU | 236 | **244** |
| 0.08 AU | 263 | 269 |

The pre-stated kill line was 290 at 0.06 AU. It is missed by 46.

**Real test anyway** (tiling at 0.08, bases settled, regret growth restricted to each route's assigned pool, then open):

| phase | covered | sum J_i |
|:--|--:|--:|
| tiled bases | 95 | 8.3 |
| pooled growth until pools exhausted | 196 | 9.86 |
| open growth, +45 then +90 kg steps | **249** | 11.61 |
| reference: untiled regret growth (R1) | 237 | 11.21 |
| reference: sequential + leftover pass | 249 | 11.56 |

Jointly chosen carriers buy +12 targets over independent ones at +0.4 J: nothing. The static neighbourhood is only
valid for the first few insertions; each insertion bends the route and the usable-trial rate collapses exactly as in
the untiled run (40/48 -> 1/39 in eight rounds).

## 8. The leftovers are NOT a family -- the wall is packing (2026-09-21 00:15)

Section 7 ended by asserting that the ~50 targets left over "are the same family in every run". That was a hypothesis,
not a measurement, and measuring it (`results/s7/leftovers.py`, six fleets: F1, F1p, R1, T1, stage 6 PART/best, fill2)
falsifies it:

- each fleet leaves 49-70 targets, but only **3** targets (138, 275, 152) are left by ALL six, and 10 by five of six;
- the union of the 32 carrier routes of four fleets covers **290 of 298**.

So there is no hard family to design specialists for. Every construction drops a DIFFERENT ~50: the wall is packing,
i.e. which routes are combined, not which targets exist. And routes are INDEPENDENT spacecraft -- J = sum_i J_i +
N_miss couples nothing -- so any subset of routes ever built is a legal fleet, and the fleet problem is a plain set
cover with N as an OUTPUT:

    min sum_r J_i(tank_r) y_r + sum_t z_t    s.t.   sum_{r: t in r} y_r + z_t >= 1

`tools/route_cover.py` solves it over every `route_*.npz` on disk: 899 files -> 387 distinct routes, union 298/298,
HiGHS optimal in under a second. Result: **10 craft, 298 covered, J 14.355** -- it re-derives t10d (14.366) from
scratch. The pool is the limit, not the packing algorithm.

### What the pool says about the target
J_i = 1 + x + x^2 has a FIXED cost of 1 per craft, so J per flyby is dominated by depth:

| route | flybys | tank | J per flyby |
|:--|--:|--:|--:|
| best carrier route (`grow_c00b`) | 40 | 992 kg | **0.0339** |
| best t10d route | 37 | 918 kg | 0.0346 |
| t10d fleet average | 29.8 | 875 kg | 0.0415 |
| farm routes at 24-29 flybys | 24-29 | 765-857 kg | 0.042-0.049 |

A disjoint fleet of routes at 0.034 J per flyby covering 298 costs J = 298 x 0.034 + 2 = **12.1**; at t10d's 0.0415 it
costs 14.4. The carrier generator already produces ONE route at 0.0339. Stage 7's remaining question is therefore
narrow and concrete: can it produce ENOUGH deep (35-40 flyby) routes, diverse enough that eight or nine of them tile
the catalogue? `tools/carrier_farm.py` mass-produces them (one every 3-9 min per core) and `route_cover.py` decides.
First lesson from the farm: a 40 % random target exclusion starves growth (24-29 flybys, 0.042-0.049 J/fb) -- diversity
must come from the carriers themselves (different geometry, different natural targets), not from hiding targets.

## 8. Where stage 7 stands

Established:
- a second, independent route generator (carrier geometry + exact phase DP + sparse-route insertion) that builds ONE
  N = 8-grade route (40 flybys @ 992 kg) in five minutes;
- the fleet wall at ~240-250 of 298 is reproduced by it with no shared search code, in 20 minutes, under four
  different fleet constructions. With stage 6 that makes ten constructions converging on the same number.

240-250 is what eight routes GROWN GREEDILY IN ONE PASS cover, whichever generator grows them: the greedy fleet keeps
whatever each route happened to take. It is not a property of the catalogue (section 8: the leftovers differ from run
to run, and 32 carrier routes already cover 290 of 298 between them). The fleet layer must therefore be selection over
a large pool of routes, never a single pass -- `route_cover.py` over `carrier_farm.py` output is that experiment, and
its outcome decides whether stage 7 beats 13.86.

## 10. The fleet layer, measured end to end (2026-09-21, 01:00-10:00)

### 10.1 Fleet assembly is a set cover, and it is solved
Routes are independent spacecraft, so `tools/route_cover.py` (min sum J_i y_r + sum z_t, N an OUTPUT) is exact.
Over 489 distinct routes: **10 craft, 298 covered, J 14.355**, LP bound = MIP optimum. The algorithm is not the
limit: the pool is. (14.355 vs t10d's 14.366 is 0.011 J, far below the cost of re-exporting and re-validating.)

### 10.2 Column generation is unusable here: the duals are degenerate
Priced against the LP duals, new routes looked strongly improving (reduced cost +1.0 to +2.3). Adding them moved the
bound by nothing. One route's top targets had duals 0.336 / 0.293 / 0.248; after the column entered, all three were
0.000 and its sum pi fell 3.82 -> 1.47. Set-cover duals here have many optimal vertices, so pricing against one is
meaningless. (Stage 6 saw the same thing as "N=11 degenerate, mu 13"; this is the mechanism.)

### 10.3 Deepening works, and is not enough
`tools/fleet_deepen.py` grows every chosen route (insertion is cheap into insertion-built routes). All 10 improved:

| route | flybys | tank | J per flyby |
|:--|:--|:--|:--|
| 01 | 36 -> **48** | 904 -> 1143 kg | 0.0351 -> **0.0321** |
| 02 | 33 -> 40 | 904 -> 994 kg | 0.0383 -> 0.0340 |
| 00 | 37 -> 42 | 918 -> 999 kg | 0.0346 -> 0.0325 |
| 06 | 27 -> 36 | 792 -> 948 kg | 0.0428 -> 0.0364 |

Route 01 at 48 flybys / 1143 kg beats THU-LAD's per-craft average (42.6 @ 1097). Fleet surplus went 0 -> 67 flybys.
But the cover still returns the ORIGINAL ten: deepening raises each J_i (sum 12.355 -> 13.430) and only pays if a
craft can be dropped -- and drop-one shows every route uniquely owns **19-29** targets. The surplus is spread, not
concentrated.

### 10.4 The invariant, now measured a third time
`tools/fleet_dissolve.py`, 6.5 h, five victims,each abandoned over budget:

| victim | saved | placed | spent | J per target | budget |
|:--|--:|--:|--:|--:|--:|
| 09 | 1.322 | 10 of 19 | 1.483 | **0.148** | 0.068 |
| 05 | 1.406 | 14 of 24 | 1.697 | **0.121** | 0.054 |
| 08 | 1.368 | 12 of 25 | 1.385 | **0.115** | 0.052 |

A target FORCED onto a finished route costs 0.115-0.148 J. A target a route is DESIGNED AROUND costs 0.032-0.042 J.
**Freedom is worth about 3x**, and that single ratio is now measured by four independent machineries: stage 6's
prizes / assigned pools / mandatory waypoints, the 2026-09-16 redistribution campaign (0.10-0.25 J), stage 7's
restricted-growth farm (0.045 vs 0.034 J/fb), and this dissolve.

### 10.5 Co-design does not escape it either
`tools/carrier_lloyd.py` is Lloyd's algorithm on carrier orbits: assign each target to its geometrically cheapest
carrier, re-optimise each carrier (4-parameter pattern search on the DP) on its own targets only, repeat. It
converges in seconds and the partition IS disjoint by construction (260 assigned, 260 distinct). But the bases
collapse to **8-14 flybys against 17-22 when free**: respecting an assignment costs exactly the depth that the
assignment was meant to buy. Same 3x law, approached from the opposite side.

## 11. Conclusion of stage 7

`results/CTOC14_Result_t10d.txt` stands: raw J 14.366055, VALID, banked at 13.81 shown.

Stage 7 built a genuinely new generator and a correct fleet layer, and both work:
- a carrier route reaches 40-48 flybys at 992-1143 kg (0.032-0.034 J per flyby) in 5-10 minutes from nothing, which is
  better per craft than any t10d route and better than THU-LAD's 7-craft average;
- the fleet layer is an exact set cover whose LP bound it attains.

They do not combine, for one measured reason: **every route is cheap only on the targets it chose, and 3x more
expensive on any other.** Ten routes chosen freely overlap (union of 82 carrier routes needs 20 craft); eight routes
forced to tile the catalogue go shallow (Lloyd 8-14, tiling 244/298, restricted farm 0.045 J/fb). The leaders' 7-8
craft at 0.035 J per flyby require routes that are simultaneously deep AND disjoint, which no mechanism in this
codebase produces.

The one algorithm that could, and has never been built: a **joint k-route DP** -- all k lag curves advanced together
over one shared target set, a target consumed by whichever route takes it, with the k-way assignment inside the
recursion instead of outside it. The single-route DP is exact and costs 0.03 s; the k-route version is not a small
change (state = k slopes + the consumed set), but it is the only construction that never pays the 3x premium.
Anything short of that will keep landing on 10 craft near J 14.35.
