# Stage 9-10 — the solver architecture, and what the measurements forced us to change (2026-09-21)

## 1. The decisive measurement

Take a REAL 37-flyby route (one of the 10-craft fleet) that settles in the exact model at **16.64 km/s**, keep its
target sequence, and price it with an impulsive Lambert chain — the cost model every transfer-database architecture
is built on:

| what is priced | cost |
|:--|--:|
| the route, settled (truth) | **16.64 km/s** |
| same chain, exact epochs, impulsive Lambert junctions | 29.82 |
| same chain, epochs jittered by 0.25 d | 33.58 |
| same chain, epochs jittered by 1.0 d | 70.44 |
| same chain, epochs jittered by 2.0 d | 108.20 |
| same chain, snapped to a 6 d encounter grid (median error 1.1 d) | 69.90 |

Two independent failures: the impulsive model **overprices by 1.8x even at exact epochs**, and its cost is
**near-singular in the flyby epoch** — a quarter-day of timing error is already visible, one day doubles the route.

**Therefore a fixed-epoch transfer database cannot work for this problem.** GTOC9 and GTOC11 databases work because
those are LEO/rendezvous problems whose transfer cost is smooth in nodal drift and arrival time. Ours is heliocentric
FLYBYS: the impulsive cost is singular in epoch, while the continuous-thrust cost is smooth, because a low-thrust
craft never matches a conic through two points — it stays near one slowly-changing orbit and lets the targets come
to it. This is the mechanical reason our carrier routes (0.032-0.034 J per flyby) beat everything the beam produces.

## 2. The architecture

```
     CATALOGUE (298 reachable NEAs)
              |
   [A] CARRIER  -- the craft's ORBIT is the state; flyby epochs are OUTPUTS, never grid points
       tools/carrier_opt.py     4-parameter pattern search (t_launch, v_inf) on the DP value
       tools/exp_carrier.py     which targets are cheap = analytic node-gap test, 100k carriers x 298 in seconds
              |
   [B] ROUTE   -- exact DP on the (time, lag) event plane of ONE carrier, state = incoming arc, no grid, no beam
       tools/carrier_dp.py      0.03 s per route;  gap coupled to lag slope (gap0 - (2/3) a s)
              |                 ~20 cheap targets per carrier, 0.10-0.25 km/s each
   [C] DEPTH   -- greedy insertion into the SPARSE carrier route (4-50 kg per target)
       tools/carrier_farm.py / fleet_deepen.py       -> 40-48 flybys at 992-1143 kg = 0.032-0.034 J/fb
              |
   [D] FLEET   -- exact set cover, N is an OUTPUT:  min sum J_i y_r + sum z_t
       tools/route_cover.py     LP bound = MIP optimum -> the pool, not the algorithm, is the limit
              |
   [E] EXACT   -- impulsive twin -> continuous-thrust SCP -> validator
       ctoc14/impulsive.py, ctoc14/globalopt.py, tools/run_ialns.py export
```

Measured, and now closed off:
- **[D] cannot be driven by column generation**: set-cover duals are degenerate (a column priced at reduced cost
  +2.3 moved the LP bound by zero; its targets' duals collapsed 0.336/0.293/0.248 -> 0).
- **[D] is pool-limited**: over 489 routes it returns 10 craft at J 14.355, re-deriving t10d from scratch.
- **The wall between [C] and [D]**: a target a route CHOOSES costs 0.032-0.042 J; a target FORCED on a finished route
  costs 0.115-0.148 J. Freedom is worth 3x, measured by four independent machineries.
- **Stage 9's transfer graph is rejected** by section 1, after being built and measured
  (`tools/encounter_db.py`, `tools/transfer_db.py`, `tools/graph_beam.py`, `tools/graph_fleet.py`). It is genuinely
  1000x faster (0.6 s per route, 0.5 s per 12-craft fleet) and the tooling is sound; the COST MODEL under it is not.

## 3. Stage 10: the drifting carrier — right economics, missing term

One fixed carrier can only ever serve ~20 targets: its cheap set is geometry. Moving to a neighbouring carrier costs,
in the calibrated element metric (`ctoc14/phasemodel`, 2 % on the real fleet),

    dv = K_TOTAL [ K_DRIFT |d drift| + K_PLANE |d i_vec| ] = 0.26-1.5 km/s

while ONE extra flyby is worth 1.29 km/s at a 900 kg tank. So a drift that unlocks another ~17 cheap targets should
pay for itself ten times over. `tools/carrier_drift.py` searches exactly that: a route is a SEQUENCE OF CARRIERS over
the mission joined by continuous element drifts — no impulsive hop, no epoch grid — by a DP over (segment, carrier)
with the element metric as the transition cost, duplicates blacklisted and the DP re-run until the target set is clean.

It runs at **1.5 s per route** and does raise the base: **28 flybys over 8 carriers for 3.30 km/s of drift**, against
17-22 for a fixed carrier. But settling it in the twin (`tools/drift_settle.py`) gives 17-22 flybys at 842-865 kg,
i.e. 0.055-0.072 J per flyby — WORSE than a fixed carrier, and 3 of 6 routes do not settle at all.

**The missing term is phase continuity.** Each carrier's event lags are defined from that carrier's own launch. The DP
charges the element change at a segment boundary but lets the craft arrive on the new carrier at an arbitrary phase;
the real trajectory has to pay for that, and the settle reveals it (real dv 13-14 km/s against the DP's ~3).

The fix is structural, not a tuning knob: carry the LAG across the boundary in the DP state — state
(segment, carrier, lag bin) with the drift changing the lag slope, so the entry phase of a segment is the exit phase
of the previous one. That is ~8 x 240 x 40 states, tractable, and it is the next thing to build.

---

## Stage 10 campaign (2026-09-21)

Four tracks ran in parallel, each built new machinery, each was then re-verified adversarially by an independent
agent (every fleet re-settled from its stored state, every headline number recomputed). **The numbers below are the
VERIFIED ones, not the claimed ones.**

### Verified results

| track | machinery built | verified sum J_i | covered | verified J | verdict |
|:--|:--|--:|--:|--:|:--|
| drift-dp — phase-continuous drifting carrier | `tools/drift_dp.py`, `deep_guided.py`, `cand_density.py`, `fleet_partition.py`, `fleet_dedupe.py` | 12.3546 | 298 | **14.3546** | dead end as a route-quality lever; 2 spin-offs live |
| fleet-milp — convex-cost facility location | `tools/fleet_milp.py` (bases / milp / greedy / build) | 18.1516 | 244 | 74.1516 | dead end |
| fleet-aco — ant colony over the whole partition | `tools/fleet_aco.py`, `carrier_reach.py` | 8.2358 | 151 | 157.2358 | dead end |
| pool-scale — scaled route pool + exact cover | `cover_curve.py`, `pack_limit.py`, `carrier_capacity.py`, `cap_greedy.py`, `cover_delta.py`, `carrier_optw.py`, `overlap_prune.py`, `fleet_curve.py` | 12.3546 | 298 | 14.3546 | dead end (premise falsified) |
| **FINAL EXACT COVER** (this run, `results/s12/cover_clean`) | `tools/route_cover.py` over **1731 route files → 756 distinct**, union = 298 | **12.354601** | **298** | **14.354601** | LP bound 14.355 = MIP optimum |

References: banked **t10d raw J 14.366055** (validator-CONFIRMED continuous thrust, `results/validator2_t10d.log`);
exact-cover baseline **14.355**; target **J < 13.0**.

**The final cover reproduces the banked t10d partition exactly** — all 10 target sets identical (`A==B`, 10/10),
298 slots for 298 distinct targets, zero overlap, missing = {131, 144}. It differs only in that the pool held
marginally lighter settled copies of 3 of the same target sets (919.318→918.154, 904.783→904.426,
877.791→876.521 kg): **Δ = −0.002841 J, and it is an IMPULSIVE-TWIN number.** t10d's own twin reads 12.357442 and
its validated export reads 12.366055, so the twin→continuous conversion costs **+0.008613 J**. Converted, this
fleet would land at J ≈ 14.363 — inside conversion noise of the banked 14.366, not progress.

Adding the entire stage-10/11 output (180 further files, 4 new tracks, 42 brand-new farm routes, 2 new carrier
families, 680 optimised carriers) moved the cover result by **exactly 0.000000**.

### What is now ESTABLISHED

1. **The fleet layer is not pool-limited in the count sense, and not algorithm-limited.** J is flat at 14.355 from
   46 routes through 756; LP bound = MIP optimum at every pool size. "LP = MIP" was being read as "pool-limited";
   it only means the MIP is easy.
2. **A zero-overlap, full-coverage 10-craft fleet costs 0.0415–0.0418 J per flyby, however it is built.** Exact set
   cover over 756 routes → 12.3546 / 298 fb = **0.041459**. Redundancy-stripping the deepest fleet ever built
   (`results/s7/deep1`, 365 slots at 0.0368 J/fb) → 12.445 / 298 fb = **0.041761**. Two constructions from opposite
   directions agree to 0.7 %. The contest needs **0.0369**. That 12 % is the partition premium, now measured by six
   independent machineries.
3. **Individual route quality is NOT the constraint.** 76 of the 759 distinct routes in the corpus beat the 0.0369
   breakeven on their own (best 0.0312 at 41 flybys / 918 kg; deepest 48 flybys / 1143 kg at 0.0321).
4. **Disjointness is the constraint, and it is hard.** Re-measured here on the enlarged corpus: among routes under
   0.0369 J/fb the largest **mutually disjoint** set is **4 routes covering 157/298**; under 0.042 it is 5/186;
   under 0.045 it is 7/224. A J<13 fleet needs 8 disjoint routes covering 298.
5. **Free growth destroys disjointness.** P(two efficient routes disjoint) = 0.0023 for free-choice farm routes vs
   0.1118 for routes taken from previous partition constructions (random null 0.009–0.012) — a 49× ratio. 30 new
   farm routes gave 0 disjoint pairs out of 153. Random-graph clique law: an 8-clique needs ~3.3e10 farm routes.
6. **Carrier geometry ceiling.** At 0.037 AU (the radius where an extra target is worth its price on a deep host),
   8 carriers reach only 185/298 targets and 10 reach 210, growing **+6.4 targets per DOUBLING** of the carrier
   library (measured 60 → 680 carriers). Full coverage therefore forces 29–65 targets bought at 0.10–0.15 AU for
   102–160 kg (0.11–0.18 J) each.
7. **The insertion law, on 1225 real trials** (`tools/insert_model.py`, `results/s10/insert/samples.jsonl`):
   `kg ~ d^1.13 · depth^0.91 · margin^-0.74`, and P(settle) falls 0.74 → 0.60 → 0.31 → 0.21 across depth bands
   10-19 / 20-29 / 30-37 / 38-48. The fraction of attempts that both settle and beat breakeven is
   **53 % / 29 % / 6 % / 1 %**. All 653 failures were aim-point HOMOTOPY divergences, not settle failures — so
   feasibility is a reachability question and is not fixable with more iterations. **The old 5-bucket insertion
   table was a shallow-host table** (its 58 samples came from depth-12 and depth-15 hosts) and underprices
   insertions into the deep routes a J<13 fleet is made of by 3–12×.
8. **The refund asymmetry.** Stripping the 67 duplicate flybys out of `deep1` refunded 0.985 J = **0.0147 J per
   removed flyby**, against the 0.0415 J the fleet pays per flyby kept. A route's cost is dominated by its base
   structure, not by its marginal flybys, so overlap is cheap to add and nearly worthless to remove.

### What is newly FALSIFIED (with mechanism)

- **Drifting carriers are not a route-quality lever.** The phase-continuity fix of §3 WORKS: carrying the lag
  through the DP state with a common longitude reference raised the settled base from 19 to 24 flybys and 12 of 16
  seeds now settle (against 3 of 6 for the broken version). It buys nothing after deepening. Paired over 12 seed
  carriers, drift 0.0406 vs fixed 0.0394 J/fb, drift better on **4 of 12**, Wilcoxon **p = 0.42** (the track
  reported 6/12 and p = 0.85; the verifier's recomputation makes drifting look *worse*, and one drift "route" was a
  byte-identical duplicate, so the arm is 11 distinct routes, not 12). **Mechanism:** final depth is set by
  candidate density inside the affordable insertion radius (~0.035 AU at 950 kg), not by the base. Targets within
  0.02/0.035/0.06/0.10 AU of the settled trajectory at matched depth: drift 3/11/21/51, fixed 0–5/8–14/21–30/44–56
  — identical. A craft's track is a 1-D curve through a 3-D region; changing which carrier it rides **rotates** the
  curve, it does not lengthen it, and the drift spends ~1.3× its modelled dv on plane and eccentricity changes to
  arrive at the same ~950–1000 kg ceiling. Keep drifting only as a **supply/diversity generator** for the pool
  (65 drift routes cover 261 of 298; 13 fixed routes cover 117).
- **Facility-location MILP over carriers is correctly posed and dies on its cost matrix.** Closest approach to the
  host trajectory prices a target a route CHOOSES (R² 0.651) and is uninformative for a target it is GIVEN
  (R² 0.010, barely better than a constant 37 kg), while underpricing by **×1.94** aggregate — worst (×2.4) exactly
  where the model is most confident. Predicted 15.53 covering 287 → settled **18.15 covering 244**. The margin term
  explains only ×1.24 of the ×1.94; insertion index and regression-to-the-mean explain ~16 % of the log-residual.
  **Do not re-fit this surrogate on assigned insertions — the informative feature is not in the
  (d, depth, tank, margin, v_rel) family.**
  ⚠ **CITATION CORRECTION:** the MILP's B&B dual bounds (J ≥ 12.73 at N=8, 13.39 at N=9, 13.98 at N=10) bound
  **that formulation**, not the contest problem. Its free-N bound of 14.630 is *beaten by the existing VALID t10d at
  14.366*, which proves the formulation cannot express the incumbent. These bounds must never be quoted as evidence
  that J<13 is unreachable.
- **ACO over the whole partition is the same wall in new clothes.** The premise ("an ant assigns each target exactly
  once by construction, so nothing is forced") is false: **fixing membership before the trajectory exists IS
  forcing**. Surrogate plan 8 craft × 38 flybys covering 298 at J 12.259 realised as **7 craft, 151 covered,
  J 157.24** — 51 % air. The colony searches a space whose feasibility it structurally cannot see, because the real
  gate is P(settle) on the route's own settled trajectory, which does not exist until after the assignment.
  Verification sharpened two points: (a) the stops are **majority PRICE, not infeasibility** — of 18 probe
  insertions on the shallowest stopped routes, 12 settled cleanly (miss 11–82 km) and were rejected only for
  costing +64 … +296 kg, i.e. **0.059–0.272 J, median ≈ 0.13 J**, an independent fifth measurement of the
  forced-target premium; (b) a **matched-carrier control** the track did not run: on carrier 216 the free route flew
  **29 flybys for 789 kg** against the assigned route's **22 flybys for 828 kg**, with the free build handicapped on
  every acceptance parameter. The assignment, not the builder, is the cause.
  ⚠ The ACO loop itself is reusable; its fitted selection-curve surrogate is **not** — its marginal insertion price
  is 0.5–8.3 kg against 52–303× more when measured (median 143×), which is why every unconstrained campaign drives
  to 7 craft × 45 flybys, past the deepest route ever settled.
- **Scaling the route pool is finished.** +42 new routes from two new carrier families, 99 fleet directories, 680
  optimised carriers, 574 → 756 distinct columns: J = 14.355 at every step, LP = MIP at every step. Carrier
  price-steering works on geometry but not economics (rare/popular base ratio 0.45 → 1.08, base depth 17 → 11, but
  the resulting routes came out 19 % WORSE at 0.0486 J/fb).
- **Overlap-then-prune is negative twice** (13.430 → 12.541; 13.459 → 12.400), for the reason in *established* §8.

### The single highest-value next step

**Reframe the target.** J < 13 needs sum J_i < 11.0 = 0.0369 J/flyby with zero overlap — an **11 %** reduction of an
invariant that two independent constructions reproduce to 0.7 %. Nothing measured in this campaign comes near it,
and the four architectures that attacked it all failed at the same wall. The *scoring* bar is much closer: beating
raw **14.153** needs sum J_i **12.153**, i.e. **−0.20 J = −1.6 %** on a partition we already hold.

So the one thing worth the remaining days is **the orphan-targeted partition generator**, because it is the only
lever aimed at the measured cause rather than a symptom:

> Over the 36 stage-10 base trajectories, **124 of 298 targets lie on no route**, and at a depth-30 / 950 kg host
> only 50 of them are affordable. The other **74 cost sum ΔJ 7.33 against 2.73 if they were free choices — a
> +4.60 J surcharge**, which is the entire gap between 12.35 and the 11.0 that J<13 needs.

Concretely: run `tools/exp_carrier.py` (100k carriers × 298 targets in seconds, analytic node-gap test) and
`carrier_opt.py` with the objective set to **"how many of those 74 orphans have a track gap under 0.035 AU"** —
it has only ever been run with DP value as the objective. Build routes on the winners with `carrier_dp` +
`deep_guided.py` (now 50–190 s per route at a 90 % settle rate, 15–25× faster than `fleet_deepen`), generate the
complements with `fleet_partition.py`'s avoid-lists so the products are **partition-derived** (P(disjoint) 0.112,
49× the farm's 0.0023), strip each candidate fleet with `fleet_dedupe.py --bulk`, and feed everything to
`route_cover.py`. **Gate it on one number by 09-25: does any NEW fleet reach sum J_i < 12.30 in the twin?** If not,
the invariant in *established* §2 is physical and t10d stands.

Two cheap, low-risk items that do not need a research bet: the `results/s12/cover_clean` fleet is −0.0028 J in the
twin and has never been exported, and the twin→continuous conversion costs +0.0086 J that a longer SCP pass may
partly recover. Together ~0.011 J — worth doing as a *parallel* submission file, never by touching t10d.
