# Stage 13: the J<13 solver plan (PCC-8G), 2026-09-22

Goal (user): review the stage-10 solver and result; build a new solver that reaches J < 13.
Hard rule: `results/CTOC14_Result_t10d.txt` and `results/validator2_t10d.log` are never modified. Every new artefact
goes into a new file under `tools/s13_*`, `results/s13/...` or `results/CTOC14_Result_s13*.txt`.

This plan takes the judges' winner (FREE, "PCC-8": covering passes + set-cover selection + twin closing), grafts the
best runner-up ideas, and settles the judges' disagreements explicitly (section 3.2). Before the plan was fixed, the
two quantities every judge called decisive were measured on the best settled 8-route pass that exists: the closing
price, and the tail capacity (section 2, 02:40-03:55, all files under `results/s13/plan13/`). The measurements
change the verdict, so this document is **capacity-first**, with early kill gates.

---

## 0. Verdict

**J < 13 is not reachable with any route generator or closer we own. The plan is built so that this is confirmed or
refuted by 09-23 08:00, at a cost of about one day of 8 cores.** The expected outcome is that t10d (raw 14.366,
shown 13.80) stands.

What was measured tonight, on FREE's N6a pass (263 covered, sum J_i 10.207, all 8 routes settled):

1. **Closing is expensive, and half the leftovers cannot be closed at all.** A sequential twin closer placed 17 of
   the 35 leftovers for +1.726 J:
   - the cheapest 8 cost a mean 0.026 J;
   - the next 9 cost a mean 0.168 J;
   - **17 could not be placed at all**, even with a 2x longer homotopy.

   FREE's "~13.5 with point insertion" assumed ~0.04 J for all 35. The judges' expected ~13.3 inherits that assumption.
2. **Re-planning cannot close.** Pinned segment re-plans of all 8 hosts inserted **0 of 35**.
3. **A 9th route cannot sweep.** It flies 5 of the 35 leftovers (dv 1.2). At dv 2.0 it plans 21, but 0 of 8 depths
   settle. Over the 18 that closing left, it flies 2.
4. **The tail portfolio works, but only at the frontier price.** Keeping N6a's first 5 routes and re-planning the tail
   with six leg-model variants gives the **best settled 8-route pass on record: 273 covered at sum J_i 10.476**. That is
   +10 targets for +0.269 J, or 0.027 J per target.
5. **`run_ialns.relocate_pass --swap` is broken.** Two of the three applied swaps each dropped a third target
   (coverage 280 -> 277). The shared tool must not be used with `--swap`, least of all on the banked fleet (section 2, M-C).

**Consequence.** J is decided by settled pass capacity, not by closing. The closer can finish about 8 leftovers at
~0.026 J; each further leftover costs 0.17 J to 1 J.
- J < 13 needs a settled 8-route pass of about **290 covered at sum J_i <= 10.75** with every leftover cheap.
- Beating the banked 13.80 shown (raw < 13.90 on 09-26) needs about **282 at <= 10.6** with no unplaceable leftover.
- The best settled pass is 273 at 10.48.

**Architecture: PCC-8G, capacity-first.**
- [P] pass generator with the tail portfolio, bisect repair and a tank guard;
- [G] generation loop with closability-aware premiums and pilot choice;
- [M] exact master;
- [C] finisher closer;
- [X] polish and export, only if the capacity gate passes.

**Expected J (raw):**

| outcome | raw J | probability |
|---|---|--:|
| t10d stands | 14.366 | ~0.90 |
| improvement | ~13.7-13.9 | ~0.08 |
| J < 13 | ~12.9 | ~0.02 |

Probability-weighted raw J is about 14.3.

---

## 1. Review findings

### 1.1 The reframe is right; three numbers need correcting

J = sum_i J_i + N_miss, J_i = 1 + x + x^2, x = (tank - 600)/1400. For 298 covered at N craft, J < 13 is equivalent to
sum(x + x^2) < 11 - N. With equal tanks (the optimum, since x + x^2 is convex):

| N | sum(x+x^2) budget | equal-split tank | kg per flyby at 298/N per craft |
|--:|--:|--:|--:|
| 7 | 4.0 | 1169 kg | 13.36 |
| **8** | **3.0** | **1006.8 kg** | **10.92** |
| 9 | 2.0 | 862 kg | 7.92 |
| 10 | 1.0 | 728 kg | 4.30 |

- The note "J<13 with 8 craft allows ~11.5 kg/flyby (tank 1030 kg)" is a *shown*-J figure: 8 x 1030 kg gives raw J
  13.212. The raw limit is 1006.8 kg, i.e. 10.92 kg/fb.
- t10d's real propellant is **2750.2 kg = 9.23 kg/fb** (validator2_t10d.log, m0 per craft), not 2766. Spread over
  8 craft that gives raw J 12.447.
- Leader propellant is an equal-tank upper bound and assumes exactly 2 misses.
- N=9 is *not* an easier J<13. It needs 7.92 kg/fb fleet-wide, 14% better than t10d; N=8 allows 18% more than t10d.
  N=7 is out of reach: sequential passes cannot give 7 routes of 42.6.

The reframe stands. The gap to the leaders is the craft-count term; J<13 means 8 disjoint routes averaging 37.25
flybys at 10.92 kg/fb or less. **What the reframe missed is that the binding quantity is 8-route disjoint *capacity*,
not per-route efficiency.** Section 2 shows that capacity cannot be bought after the fact.

### 1.2 What stage 10 got wrong (`docs/stage10_architecture.md`)

1. **The "0.0415 J/fb invariant" (lines 115-123, 194).** It is 0.0336 of fixed craft cost (10/298) plus 0.0079 of
   fuel. It describes a 10-craft fleet only. "An 11% reduction of an invariant" treats a number that is 81% craft
   count as an efficiency.
2. **"Two independent constructions agree to 0.7%" (lines 119-123).** False. `s7/deep1` is the t10d partition
   (coverO == cover_clean == newgen/ifleet10, 10/10 identical), deepened and then stripped. The 0.7% agreement comes
   from the shared fixed term; the fuel terms differ by 3.8%.
3. **The insertion depth law (lines 136-142) is confounded with host provenance.** Every host of depth 36 or more in
   `s10/insert/samples.jsonl` was grown by insertion to saturation. At matched depth 28-36 and d < 0.05 AU:
   - beam-born hosts settle 80% of insertions at a median 17.5 kg;
   - insertion-grown hosts settle 36-44% at 65-72 kg.
4. **"Forcing costs 3x whether first or last" (lines 53-54, 176-177).** Every forced measurement was an insertion into
   a finished route. It also compares an average J/fb that includes the fixed cost with a marginal fuel dJ.
5. **The carrier ceiling and the "+4.60 J partition-wall surcharge" (lines 132-135, 202-204)** are larger than HIT's
   entire fuel term. They are priced by fixed carriers plus insertion.
   - *But* section 2 shows that the order of magnitude of an orphan's price (~0.1-0.2 J) is right for every realiser
     we own.
   - The error was calling it a property of the problem. The leaders prove otherwise.
6. **The ACO and MILP verdicts judge the realiser and the surrogate, not the algorithm class.**
   - `fleet_milp` prices assignment with `predict_kg`, whose own docstring gives R^2 0.010.
   - ACO stopped for price, not infeasibility.
7. **"[D] cannot be driven by column generation" (lines 50-51)** is over-general. Stage-5 CG raised LP coverage from
   240.8 to 254.3 in 4 rounds. Every CG/forcing run used prizes of 0.35-2.5, i.e. 2-14 years of mission time in beam
   units. Premiums of at most 0.02 are free.
8. **"Not pool-limited" (lines 52, 116-118).** `picked.json` shows the cover picks t10d's own route files. A flat
   14.355 only says the incumbent is optimal in pools that contain it.
9. **Tool defects.**
   - `route_cover.py --exact` keeps a failed pick.
   - `w_load` never checks the stored miss.
   - New tonight: **`run_ialns.relocate_pass` with `--swap` keys its swap results by (host, inserted target) and not
     by the removed target.** So `min(sw[(B, X)])` can pick "B minus Y' plus X" for the pair (X, Y): Y' is silently
     dropped and Y is duplicated. Two corrupt swaps lost 278, 75 and 56. t10d is unaffected: it is validated with
     298 covered, and the 09-17 swap passes applied no swap.
10. **gc16 (source of "the 67 orphans chain into 6 flybys") ran at the default `--vinf 2.0`.** At v_inf 4 the plain
    greedy covers 251.

### 1.3 What stage 10 got right

- **The last ~25-35 targets of any 8-route construction are expensive for every realiser we own.** Section 2 measures
  0.03-0.25 J sequentially, with half of them impossible.
- **Re-planning cannot add capacity:**
  - linchpin: re-plans fly 0.65-0.85 of a pool;
  - PAR M4: 0 of 15 absorbed;
  - SR probe: 0 of 35.
- **The 10-craft fleet is not the thing to optimise**: t10d is saturated (swap, retime and relocate find nothing).

### 1.4 The linchpin experiment (10 -> 8 by re-planning inside assigned pools)

Setup: t10d partition; two routes dissolved; capacitated min-distance pools (cap 38); every survivor re-planned from
scratch; production legs; beam 1000 / seg-k 50; every candidate settled.

| fleet | covered | per craft | kg/fb | sum J_i | J |
|---|--:|--:|--:|--:|--:|
| A strict (dissolve 05+09) | 193 | 24.1 | 7.08 | 9.100 | 116.10 |
| B strict (dissolve 06+07) | 194 | 24.3 | 7.80 | 9.235 | 115.24 |
| A best (dv 2.0, drmax 0.25) | 234 | 29.25 | 9.39 | 9.881 | 75.88 |
| the 8 survivors before re-planning | 246 | 30.75 | 8.48 | 9.777 | n/a |

- **Control.** A route re-planned over its own t10d set re-flies only 84%.
- **Depth tracks pool size.** About 37 needs a pool about 2.5x the route, half of it other routes' targets.
- **No effect:** beam 3000, stepping stones.
- **The stage-6 "P06 efficiency unchanged" data point** was a soft pool that flew 19 of its 32 assigned targets.
- **Verdict.** kg/fb is never the obstacle (6-10); depth per disjoint pool is.

### 1.5 Planner facts the design relies on (`results/s13/planner_capability.md`, PAR/FREE probes)

- **Free beam over all 298: 46 flybys at 1004 kg twin.** The second route flies 43 at 975 kg. Beam 100, 300 and 1000
  give the same depth; beam 300 is about 3% cheaper at a fixed depth.
- **Fixed point.** A beam-born route re-planned over its own set returns itself (A46 46/46, B43 43/43, fill2/05 36/36).
- **Prize regime.** Score = t/T + 0.520 fuel/1400 - sum(prize), T = 5478.75 d, so a premium of 0.01 = 55 d = 27 kg.
  - Premium 0.01 on the hard set: 46 flybys with 13 hard, against 4 hard at no premium.
  - At 0.03 or more the route collapses (37, then 30, then 21).
  - So every prize must stay in [1.00, 1.02].
- **Pinning.** Denial elasticity is 0.66 for hard waypoints pinned at +-45 d, against 0.065 for free routes. Never pin.
- **Settling.**
  - Deep dv-1.2 tours settle 37 of 38 times.
  - In M-D, 12 of 19 dv-2.0 tail candidates failed to settle or came out at absurd tanks (1284-2781 kg against
    806-953 planned). A 1.15x tank guard is mandatory.
  - A failing tour breaks at a single leg: bisection finds it in about 2 min.
- **Joint beam.** With the order-free key, K=2 over a 50-pool covers 33-34 against 30 sequential; over a 77-pool, 41
  against 44.

### 1.6 The pass frontier (FREE passes, now in `results/s13/seed_passes/`, plus this plan's tail probe; m0 1600, v_inf 4, beam 100)

| pass | levers | covered | sum J_i | unsettled |
|---|---|--:|--:|---|
| V0 | plain | 251 | 9.766 | none |
| G4a | premium 0.01 on the hard set, dv 1.8 on those legs | 259 | 9.962 | none |
| N6a | G4a + npt 6 / mc 150 | 263 | 10.207 | none |
| **T-N6a (this plan)** | N6a routes 1-5 + tail portfolio | **273** | **10.476** | none |
| V2 | V1 + relaxed tail from route 6 | 272 | 10.396 (planner x1.02) | routes 7, 8 |
| P20 | all relaxed | 278 | 10.955 (planner x1.02) | route 7 |

- All the loss is in the last three routes: N6a's depths are 48, 39, 42, 36, 34, then 26, 23, 15.
- Coverage beyond 263 costs 0.02-0.05 J per target along this frontier.
- Seven premium variants land at 249-258.
- MILP recombination over 117 settled pass routes gives 259.

---

## 2. Measurements taken for this plan (2026-09-22 02:40-03:55, `results/s13/plan13/`)

All four probes use FREE's N6a:
- 263 covered, sum J_i 10.2068;
- depths 48, 39, 42, 36, 34, 26, 23, 15 at 1048, 1039, 993, 930, 881, 836, 848, 699 kg;
- `results/s13/seed_passes/N6a` (copied from FREE's scratchpad).

It leaves 35 leftovers. Distance to the nearest host (`N6a_leftovers.json`):

| nearest host within | 0.03 AU | 0.05 AU | 0.08 AU | 0.12 AU | 0.20 AU | 0.35 AU |
|---|--:|--:|--:|--:|--:|--:|
| leftovers | 4 | 7 | 13 | 19 | 30 | 34 |

### M-A. Planner segment re-plan as a closer (FREE's "segment re-sequencing" / J3's "strict windowed re-plan"): 0 of 35

**Setup** (`sr_fill.py`, `run_sr_N6a.sh`, `sr_N6a/NN/log.txt`):
- every host is re-planned segment by segment with all its own flybys pinned as waypoints (+-45 d);
- the 35 leftovers are available at prize 1;
- beam 300, npt 6, max_children 150, timing variants kept at each waypoint;
- `skel_fill` is used through a wrapper with one fix: segment 0 collects the root states too.

**Result**
- **No host inserted a single leftover.**
- Four hosts re-fly exactly (01 48/48 @ 1048, 04 36/36 @ 931, 05 34/34 @ 881, 06 26/26 @ 834).
- 02 re-flies 38/39, 03 35/42, 07 21/23 @ 890 kg; 08 does not settle.

With PAR M4 (a free re-plan absorbs 0 of 15), this settles it: a beam-born host is at the planner's own capacity, so
only the twin can add targets to it.

### M-B. A 9th "sweeper" route over the leftovers

`sweeper_probe.py`, `sweeper_probe2.py`; strict pool; launch grid 0-1500 d.

| pool | dv 1.2 / 0.15 AU / npt 6 | dv 2.0 / 0.25 AU / npt 6 |
|---|---|---|
| all 35 leftovers | **5 flybys @ 634.5 kg** (settled) | planned up to 21 @ 924 kg; **depths 14-21: 0 of 8 settle** |
| the 18 left after closing (M-C) | n/a | **2 flybys** @ 602 kg |

The leftovers of a good pass are anti-chainable. They cannot seed a route of their own, so an N=9 track is pointless.

### M-C. Sequential twin closing (the quantity FREE's Gate 2 and all three judges called decisive)

**Setup** (`close_probe.py`, `close_N6a/closing.json`, `close_N6a.out`). Rounds of:
1. `RI.candidates`: 3 per target, radius 0.08, widened to 0.15 and then 0.30 AU on failure;
2. `RI.w_insert` trials;
3. cheapest insertion per host, **one per host per round**, re-priced every round, accepted if dJ < 0.25;
4. a 2x longer homotopy (stages 20, iters 10) on the 4th try.

**Closing result**
- **17 of 35 placed for +1.7256 J.** Covered 263 -> 280; sum J_i 10.2068 -> 11.9325.
- Sorted dJ: 0.002, 0.0036, 0.013, 0.019, 0.032, 0.037, 0.049, 0.056 | 0.084, 0.103, 0.119, 0.145, 0.171, 0.201,
  0.203, 0.243, 0.245.
- **The cheapest 8 average 0.0265 J; the next 9 average 0.168 J.**
- **17 targets have no settling insertion at all** (3, 15, 33, 48, 53, 58, 114, 118, 132, 137, 139, 158, 196, 219,
  236, 243, 247). A further one (17) prices above 0.25.
- Several of them are only 0.06-0.09 AU from a *deep* host (33 -> 02, 53 -> 01, 118 -> 01/05, 114 -> 03). The deep
  hosts (39-48 flybys) accepted almost nothing cheap.
- The shallow tail (08: 15 flybys, 699 kg) accepted far targets cheaply: 110 at 0.25 AU for 0.019 J, 216 at 0.17 AU
  for 0.049 J.

**Relocate afterwards** (3 passes, `RI.relocate_pass` with swap=True, about 14 min per pass on 8 contended cores):
- **Legitimate moves recovered 0.350 J**: 117: 02->08 0.150; 293@04 <-> 206@06 swap 0.064; 191: 04->06 0.061;
  195: 03->05 0.038; 28: 03->05 0.036.
- **The other two applied swaps were corrupt.** Each dropped a third target and duplicated the swapped one (lost 278,
  75, 56; coverage 280 -> 279 -> 277). This is the keying bug of section 1.2, item 9.
- The honest post-relocate state is about **280 covered at sum J_i ~11.58**, with 18 targets unplaced.

**The closer's price curve:**
- about 8 leftovers at ~0.026 J;
- about 9 at ~0.17 J;
- then roughly half the leftovers unplaceable. Block moves with ejection on the t10 fleet placed such targets at
  0.10-0.25 J, or stalled (`docs/globalopt_campaign.md` section 8.2).

On N6a that means:
- **J = 11.58 + 18 + 2 = 31.6** as it stands;
- **about 15.3-16** even if block moves placed all 18 at 0.2 J.

Both are worse than t10d's 14.366.

### M-D. Tail portfolio on N6a's head: the capacity test

**Setup** (`tail_probe.py`, `tail_N6a_k5/log.txt`, fleet in `tail_N6a_k5/fleet/`):
- keep N6a routes 1-5 (199 covered);
- re-plan routes 6-8 over what is left with six members, in parallel:

  | member | settings |
  |---|---|
  | T1 | dv 1.2, npt 6 |
  | T2 | dv 1.6, drmax 0.20 |
  | T3 | dv 2.0, drmax 0.25, npt 2 |
  | T4 | T1 + launch grid 0-1500 d |
  | T5 | T2 + w_t 0.5 |
  | T6 | dv 2.0, npt 6, grid 0-1500 d |

- for each route, take the deepest candidate that settles with twin tank <= min(1.15 x planner, 1050 kg).

**Planner fronts (deepest point per member):**

| route | pool | T1 | T2 | T3 | T4 | T5 | T6 |
|---|--:|---|---|---|---|---|---|
| 6 | 99 | 26 @ 874 | 29 @ 904 | 32 @ 1008 | = T1 | 30 @ 945 | 31 @ 995 |
| 7 | 68 | 22 @ 800 | 23 @ 821 | 25 @ 997 | = T1 | 23 @ 817 | 27 @ 1015 |
| 8 | 43 | **8** @ 660 | 18 @ 790 | 20 @ 879 | = T1 | 19 @ 822 | 22 @ 953 |

**Settled:**

| route | taken | N6a's route |
|---|---|---|
| 6 | 31 @ 938 kg (T3) | 26 @ 836 |
| 7 | 25 @ 949 (T6; 27 @ 1067 also settled but was over the cap) | 23 @ 848 |
| 8 | 18 @ 772 (T2) | 15 @ 699 |

On route 8, 6 of 6 dv-2.0 candidates of 19-22 flybys failed or settled at 1284-2781 kg.

**Fleet: 273 covered, sum J_i 10.4757.** That is +10 targets for +0.269 J (0.027 J per target), the best settled
8-route pass so far.

Findings:
- The per-leg dv cap is the only lever in sparse pools: 8 -> 18 -> 20-22 for dv 1.2 -> 1.6 -> 2.0 in the 43-pool.
- The long launch grid adds nothing (T4 = T1).
- The time weight adds nothing beyond its dv setting (T5 about equal to T2).

### What the four probes mean

J is decided by **settled pass capacity**:

| outcome | settled 8-route pass must reach | at sum J_i | leftovers left for the closer |
|---|--:|--:|---|
| J < 13 | ~290 | <= 10.75 | <= 8, all within a shallow host's reach (~0.026 J each) |
| beat banked 13.80 shown (raw < 13.90 on 09-26) | ~282 | <= 10.6 | <= 16 (8 at 0.026 + 8 at 0.17), none unplaceable |

Today the best settled pass is 273 at 10.48, with 25 leftovers, about half of which will be unplaceable.

---

## 3. Architecture: PCC-8G, capacity-first

### 3.1 Pipeline

```
 [P] PASS      tools/s13_pass.py     8 sequential FREE beams over the uncovered pool (never pinned); premiums in
                                      [1.00,1.02]; dv 1.8 on premium legs; npt 6 / mc 150; TAIL PORTFOLIO (M-D) for
                                      routes >= 6; deepest-settled choice with a kappa tie-break; bisect-ban repair;
                                      1.15x tank guard
        |
 [G] GENERATE  tools/s13_gen.py      8 slots per generation: 4 full passes (premium jitter), 2 tail re-runs from the
                                      best head, 2 pilot passes (rollout choice at routes 3-6); after each generation a
                                      1-round PRICE PROBE of the best fleet's leftovers drives CLOSABILITY-AWARE
                                      premiums (+0.006 unplaceable, +0.002 dear, 0 cheap, -0.002 multiply covered)
        |
 [M] MASTER    tools/s13_master.py   exact set cover over the archive, sum y <= 8 (and <= 9, reported); LP bound;
                                      miss-checked columns; the --exact bug fixed
        |
   G1: settled N=8 >= 282 at <= 10.6 (improvement) or >= 290 at <= 10.75 (J<13)?  NO -> STOP, t10d stands
        |
 [C] CLOSE     tools/s13_close.py    FINISHER only (<= ~16 leftovers): sequential twin point insertion, cap ladder
                                      0.03 -> 0.20, long homotopy, block-with-ejection last; interleaved relocate
                                      from a FIXED copy (swap keyed by removed target; coverage asserted)
        |
 [X] POLISH +  fixed relocate/retime to saturation, then export -> combine_submission -> validator2 into
     EXPORT    results/CTOC14_Result_s13<x>.txt; submit only if VALID and raw J < the dated bar (section 5.3)
```

The closer, the polish and the export are not built until G0.5 passes. Until then the only cores spent are on
capacity.

### 3.2 The judges' disagreements, resolved

| # | Question | Positions | Resolution | Basis |
|--:|---|---|---|---|
| 1 | Which N is the deliverable | J1, J2: N=9 is the realistic deliverable; FREE: N=9 via a 9th leftover-specialist route; J3: N=8 main | **N=8 only; no N=9 track.** The master reports N <= 9 for free | M-B: a 9th route flies 5 of 35 leftovers (dv 1.2), 0 settled at dv 2.0, and 2 of the 18 hard ones. N=9 J<13 needs 7.92 kg/fb |
| 2 | Closing operator | FREE: point insertion + planner re-sequencing; J1: sequentially priced point insertion; J2: distance-ordered; J3: replace point insertion by strict windowed re-plans | **Twin-level only**, and only as a *finisher*: sequential point insertion with a cap ladder, long homotopy, block-with-ejection last, fixed relocate between | M-A: 0/35 by planner re-plans. M-C: the closer's curve (8 at 0.026 J, 9 at 0.17 J, the rest impossible) |
| 3 | Tank policy | FREE: fixed planner cap 990; J3: budget-aware | **Deepest settled** Pareto point, tie-break value = n - (tank - 600)/kappa with **kappa = 90 kg/flyby**; twin tank cap 1080 on routes >= 6, 1050 on routes 1-5 | A pass flyby costs 0.027 J at the margin (M-D), a closer flyby >= 0.17 J (M-C), so depth is worth far more than fuel until the J<13 budget binds (8 x 1000 kg is still sum J_i 10.94) |
| 4 | Two-leg look-back leg model | J3: in parallel now; J2: after Gate 2 | **Not built in this stage.** Optional only if a spare agent exists after G0.5 passes, with gate R | Its cheap proxy (dv relaxation) is measured at 0.027 J/target (M-D), i.e. no better than the frontier. The 09-09 leg-model swap gave identical sparse-pool counts. dv-2.0 legs already fail to settle in the tail (M-D), and a look-back model makes such legs *more* available |
| 5 | PAR's PULL/PUSH archive loop | J1: graft for tails; J3: use as the closer | **Rejected as a coverage engine.** Only the archive master is kept | M4, the Lloyd step (259 -> 237), milp_mix (259), M-A (0/35) |
| 6 | Prize band | all agree | [1.00, 1.02]; hard premium 0.01; any dual/urgency price rescaled into the band | PAR M4/M6, FREE's delta table |
| 7 | Hard waypoints | JOINT: never pin | Never pin | e = 0.66 vs 0.065 |
| 8 | Tail generation | FREE: dv 1.6/2.0 tails; capability review / JOINT: joint K=2 | **Tail portfolio** T1-T3, T5, T6 (drop T4); joint K=2 as member T7 for the last pair only | M-D: +10 settled; K=2 +3-4 on a 50-pool |
| 9 | Settle failures | PAR: bisect and ban; FREE: settle-until-success + 1.15x guard | **Both**: walk the candidate list (up to 18 settles); bisect the deepest failure, ban that leg's target for this route, re-run once | M-D: 6/6 dv-2.0 failures on route 8, absurd tanks 1284-2781 kg |
| 10 | Schedule | J2: dry-run export 09-24 12:00, freeze 09-26 12:00 | Adopted, but **conditional on G1**. Submit at the first validated improvement | Section 5.3 |
| 11 | *(new)* Premium update | FREE: +0.004 on all leftovers | **Closability-aware**: +0.006 on leftovers with no settling insertion, +0.002 on those priced > 0.08 J, 0 on cheap ones, -0.002 on multiply covered; clipped to [-0.01, 0.02] | M-C: cheap leftovers cost ~0.026 J to close and need no premium; only the unplaceable ones must be absorbed by a pass |

---

## 4. Components to build (new files only; shared tools imported, never edited)

All Python runs use `~/.venvs/astro313/bin/python` with `OMP/OPENBLAS/VECLIB/MKL_NUM_THREADS=1`. Every multiprocessing
pool uses the `fork` context. Every settle call runs under a 240 s SIGALRM (`greedy_cover._timeout`).
**Settled** means all three of: miss <= 150 km, twin tank <= 1.15 x planner tank, and twin tank <= the route's cap.

### 4.1 `tools/s13_pass.py`: one covering pass

Built from `results/s13/seed_passes/scripts/probe_pass5.py` (head routes) and `results/s13/plan13/tail_probe.py` (tail portfolio).

```
s13_pass.py OUT [--start FLEET --keep 01,02,03,04,05] [--prem PREM.json] [--delta 0.01] [--dvt 1.8]
            [--beam 100] [--npt 6] [--mc 150] [--tail-from 6] [--members T1,T2,T3,T5,T6,T7] [--kappa 90]
            [--cap-head 1050] [--cap-tail 1080] [--grid 0,800,20] [--pilot 0|M] [--seed S] [--nproc 6]
```

- **`plan_route(avail, prize, member) -> list[(n, planner_tank, tour_json)]`**
  - Calls `search.beam_search(RI.eph(), excluded=ALL-avail, m0=1600, t_launch_grid, P, n_proc=1, verbose=False)`.
  - `P = greedy_cover.BP(beam, w_fuel=CM.w_fuel(480,1600), m_margin=40, dv_max, tofs=15..400/5 d,
    lin_tofs=20..600/10 d, lin_drmax, vinf_cap=4.0, tof_refine=True, max_depth=60, n_per_target, max_children, w_t)`.
  - `P.prize = 1 + lambda_t` on avail; `P.dv_max_t = max(dv, 1.8)` on the premium set; `P.collect = DeepCollect(4)`.
  - Returns the min-fuel state per depth, deepest 8.
- **Members** (M-D): T1 (1.2/0.15/npt 6), T2 (1.6/0.20), T3 (2.0/0.25/npt 2), T5 (1.6/0.20/w_t 0.5),
  T6 (2.0/0.25/npt 6/grid 0-1500). Head routes (< tail_from) use T1 only.
  - T7, the last pair only: `results/s13/planner_probe/jointsym.joint_beam_search(eph, [1600,1600], grid, P, excluded,
    canon=True)` with `P.wait = 60 d` and `P.w_rare = 1`.
  - T7's tours settle through `results/s13/planner_probe/twin_wait.settle_tour_wait`. It is taken only if its two
    routes together beat the sequential best of routes 7+8.
- **`settle(job)`**: `impulsive.settle_tour(RI.eph(), tour, lambda ip: RI.settle(ip, 100))` with
  `IM.LIN_FALLBACK = 2.5`. Candidates are tried deepest first, 6 in parallel, at most 18.
- **`repair(tour)`**: when the deepest candidate fails, `s13_bisect.failing_leg(tour)` finds the leg. Its target is
  banned for this route (excluded), and the route is re-planned once.
- **Choice**: the deepest settled candidate; ties broken by value = n - (tank - 600)/kappa.
- **Library**: other settled candidates go to `OUT/lib/route_<k>_<member>_<n>.npz` as master columns.
- **Pilot** (`--pilot M`, routes 3-6): take the M best diverse candidates (at most 85% target overlap). Roll each out
  with beam 30, T1, npt 2, planner only. Keep the one with the most rollout coverage (tie: the lower planner sum J_i).
- **Outputs**: `OUT/route_NN.npz` (`RI.ist`), `OUT/pass.json` (per route: targets, member, planner and twin tank,
  miss, settle tries, lambda digest), `OUT/cols.jsonl`, `OUT/lib/`.
- **CPU** (measured in M-D, where 14 workers shared 10 cores):
  - head routes: ~20 CPU-min;
  - tail portfolio: ~105 CPU-min of beams + ~60 CPU-min of settles, i.e. 46 min wall on 6 contended workers;
  - **total about 2-3 CPU-h per full pass**; a pilot pass adds ~1.5 CPU-h of beam-30 rollouts.

### 4.2 `tools/s13_bisect.py`

- `failing_leg(tour, max_settles=6) -> (k_ok, k_fail, leg)`: binary search on prefixes `legs[:k]` with `settle_tour`
  (from `results/s13/seed_passes/scripts/bisect_settle.py`).
- `settled_prefix(tour) -> ip | None`: the longest settling prefix, kept as a library column.

### 4.3 `tools/s13_gen.py`: the generation loop

```
s13_gen.py results/s13/gen --gens 6 --slots 8 --seed-fleets results/s13/plan13/tail_N6a_k5/fleet,results/s13/seed_passes/N6a
```

- Each generation fills 8 slots, one `s13_pass.py` process each, with `nice -n 5` and `--nproc 1` except the tail
  portfolio pool:
  - 4 full passes with premiums jittered by +-0.004 (seeded);
  - 2 tail re-runs of the current best fleet (`--keep 01..05`, and once `--keep 01..04`);
  - 2 pilot passes (`--pilot 4`) from generation 2 on.
- **After each generation:**
  - (a) run `s13_master.py`;
  - (b) run a **price probe** of the best N=8 fleet's leftovers: one `RI.candidates` + `RI.w_insert` round at radius
    0.30, 6 per target, *nothing applied*, about 5-10 min on 8 cores. It classifies each leftover as cheap
    (<= 0.08 J), dear, or unplaceable;
  - (c) update `PREM.json` with the closability-aware rule of section 3.2 #11;
  - (d) log the generation, the best N=8 (covered, sum J_i), the projected closed J (ledger), and the LP bound.
- Stops at G1.

### 4.4 `tools/s13_master.py`: exact selection (copy of `route_cover.solve`, bugs fixed)

```
s13_master.py OUT 'results/s13/gen/**/route_*.npz' 'results/s13/plan13/**/route_*.npz' 'results/s13/seed_passes/**/route_*.npz' --N 8,9
```

- `load(f)`: `RI.ipr(st)`, tank, and miss (`Yf,_ = ip.integrate(); dm,_ = ip.misses(Yf)`). The column is dropped if
  miss > 150 km or tank is outside 600-1400. Dedupe by target set, lightest first.
- `solve(cols, N)`: HiGHS `milp`, minimise sum cost(tank) y + sum z, coverage rows >= 1, **sum y <= N**, 300 s, gap
  1e-3, plus the LP value.
- `--exact` re-settles the picks and **removes** any that fail before re-solving.
- `--strip` runs a strip pass (`fleet_dedupe` bulk rule via `RI.w_remove`) when the picked fleet overlaps.
- Outputs: `OUT/N8/fleet`, `OUT/N9/fleet`, `OUT/master.json`.

### 4.5 `tools/s13_close.py`: the finisher (grown from `results/s13/plan13/close_probe.py`)

```
s13_close.py FLEET OUT [--caps 0.03,0.05,0.08,0.12,0.20] [--radii 0.08,0.15,0.30] [--per-target 6]
             [--relocate-every 6] [--kill-projected 11.9] [--nproc 8]
```

- **Round.**
  - `RI.candidates` (per-target radius widening on failure), then `RI.w_insert` jobs.
  - A target on its 3rd failure uses `w_insert_long`: stages 20, iters 10.
  - Best price per target; place cheapest first, **one placement per host per round**, only if dJ <= the current cap.
- **Cap ladder.** Start at 0.03. Raise one step when a round places nothing; step back down after any placement. This
  fixes the probe's fault of letting 0.2 J placements escalate a host early.
- **`relocate_pass_fixed(fleet, pool, a, say)`**, every `relocate-every` placements. It is a copy of
  `RI.relocate_pass` with two fixes:
  - (1) swap results are keyed by **(host, inserted, removed)**;
  - (2) every applied move asserts that the union coverage does not shrink and that the two touched routes differ
    from before only by the moved targets. A violating move is reverted and logged.
- **Ledger.** projected sum J_i = current + sum over remaining targets of the best known price (1.0 if none). Stop if
  projected > `--kill-projected`.
- **Last resort** (a target with no point insertion <= 0.20): `RI.w_block(host, st, [(u, t_u)], 30 d margin, tag)`.
  Ejected targets re-enter the queue under a 3-round tabu.
- Outputs: `OUT/closed/`, `OUT/closing.json`, `OUT/log.txt`.

### 4.6 Polish and export (existing tools, new names)

```
s13_close.py ... (final relocate_pass_fixed passes to saturation)
run_ialns.py retime  <closed> results/s13/final/iretime --events results/newgen/scratch/events2.npz --retime-top 60
run_ialns.py export  results/s13/final/iretime results/s13/final/efleet --nproc 6
combine_submission.py results/CTOC14_Result_s13a.txt $(ls results/s13/final/efleet/frags/frag_*.txt | sort)
validator2.py results/CTOC14_Result_s13a.txt > results/validator2_s13a.log 2>&1
```

- **Never run `run_ialns.py search --swap`** (section 1.2, item 9). `search --relocate` without `--swap` is safe but
  should also be replaced by `relocate_pass_fixed`.
- Wrap the chain in `results/s13/export_s13.sh` with the INVALID-fragment check of `results/newgen/export_t10d.sh`.
  **Never** copy export_t10d.sh's output paths.

### 4.7 Seeds (already preserved, 2026-09-22 04:00)

The scratchpad is session-local, so the seeds were copied to `results/s13/seed_passes/`:
- FREE's passes `N6a`, `V2`, `P20`, `G4a`, `V0`, `V1`, `N10a`, `P16`, each with its `.out` log;
- `closing_V1.json`;
- `scripts/`: `probe_pass5.py`, `closing_trial.py`, `bisect_settle.py`, `milp_mix.py`, `probe_vdiv.py`.

The scripts `chdir` to the project root and read `results/newgen/gc16`, so they run from there unchanged.

---

## 5. Run schedule and CPU budget

### 5.1 Budget

The machine has 10 cores. This plan uses **8 cores**; 2 are left for other agents and for export/validation. Tonight's
probes showed that 14 workers on 10 cores slow everything about 1.6x.

| item | CPU-h |
|---|--:|
| One generation: 4 full passes (~10 CPU-h) + 2 tail re-runs (~4.5) + 2 pilot passes (~8) | ~22 (about 2.8 h wall on 8 cores) |
| Price probe | ~1 |
| Master | minutes |
| Closing run (finisher, <= 16 leftovers) | ~2 CPU-h (M-C: 7 rounds in 11 min wall on 8 cores) |
| One `relocate_pass_fixed` | ~1.5-2 CPU-h (M-C: 14 min wall on 8 contended cores) |
| Export | ~1.5 CPU-h (t10d: 10 routes, 6 procs, ~15-25 min) |

### 5.2 Timeline (CST)

| when | what | cores |
|---|---|--:|
| 09-22 04:00-05:00 | Seeds already copied (4.7). Zero-code hedge: run the existing `tools/route_cover.py` (N free) over `results/**/route_*.npz` + the seeds + `plan13/tail_N6a_k5/fleet`, to see whether any hybrid beats 14.35 in the twin | 8 |
| 09-22 05:00-10:00 | Build `s13_pass` (merge `probe_pass5` + `tail_probe`), `s13_bisect`, `s13_master`, `s13_gen`. **G0** | 2-4 |
| 09-22 10:00-20:00 | Generations 1-3 (tail re-runs of T-N6a first) | 8 |
| 09-22 20:00 | **G0.5** early kill | - |
| 09-22 20:00-09-23 08:00 | Generations 4-7 (only if G0.5 passed); build `s13_close` + `relocate_pass_fixed` in parallel | 8 |
| 09-23 08:00 | **G1** capacity decision | - |
| 09-23 08:00-20:00 | Finisher closing + fixed relocate on the top-3 G1 fleets in parallel | 8 |
| 09-23 20:00 | **G2** | - |
| 09-24 00:00-12:00 | Polish; export dry run (`export_s13.sh` -> `s13a`). **G3** | 6 |
| 09-24 12:00 | Submit s13a if VALID and raw J < 14.00 | - |
| 09-24 12:00-09-26 12:00 | If J<13 is still open: more generations on the G1 winner's head, then close/export `s13b`. Otherwise polish only | 8 |
| 09-26 12:00 | **G4** | - |
| 09-27 12:00 | **Freeze** | - |
| 09-28 06:00 | Latest submission of anything | - |

### 5.3 Break-even bars against the banked t10d

Shown J = k x raw J, with k = 1 - 0.1 (09-28 12:00 - t_submit)/28 d. t10d is banked at 13.801851 shown. A new
submission helps only if its raw J is below:

| submitted | 09-23 12:00 | 09-24 12:00 | 09-25 12:00 | 09-26 12:00 | 09-27 12:00 | 09-28 11:00 |
|---|--:|--:|--:|--:|--:|--:|
| raw J must be < | 14.053 | 14.003 | 13.952 | 13.901 | 13.851 | 13.804 |

Each day of delay costs about 0.049 J of raw margin.

---

## 6. Checkpoints and kill criteria

All quantities are settled twin numbers; "fleet" means the master's best N=8 pick.

| gate | when | test | pass -> | fail -> |
|---|---|---|---|---|
| **G0** build | 09-22 10:00 | `s13_pass --start seed/N6a --keep 01..05` reproduces T-N6a (273 +- 2 at sum J_i 10.48 +- 0.05); `s13_master` over N6a's own routes returns N6a | start generations | fix; do not start |
| **G0.5** early kill | 09-22 20:00 (~3 generations) | fleet >= **278 at sum J_i <= 10.55**, *or* the projected closed J of the best fleet (price-probe ledger) < 13.9 | continue to G1 | **STOP the N=8 line. t10d stands.** Release the cores |
| **G1** capacity | 09-23 08:00 | (a) fleet >= **290 at <= 10.75**, with every leftover priced <= 0.08 J -> J<13 attempt; (b) fleet >= **282 at <= 10.6**, no unplaceable leftover -> improvement attempt | build and run the finisher | **STOP. t10d stands** |
| **G2** closing | 09-23 20:00 | after the finisher + fixed relocate: all 298 covered and twin J < 12.95 (a) or < 13.85 (b) | polish + export | STOP; t10d stands |
| **G3** export | 09-24 12:00 | validator2 VALID; raw J < 14.003 | submit s13a | retry export at 0.05/0.025 d steps; resubmit by 09-25 12:00 against 13.952 |
| **G4** | 09-26 12:00 | a further twin gain larger than the time decay (0.049 J/day) | export s13b | keep s13a |
| Freeze | 09-27 12:00 | no new optimisation; only a validated export may still be submitted | - | - |
| Track R (optional) | 09-24 18:00 | 87-target gc16 remainder (`probe_vdiv.py p87`) reaches >= 27 settled at <= 800 kg twin (today: 25 @ 775) | one generation with `ctoc14/search_lb.py`; re-test G1 by 09-25 12:00 | drop R |

**Kill rules that override everything:**
- Never submit anything that validator2 does not call VALID.
- Never touch the t10d files.
- If at any gate the projected J cannot beat the next dated bar, stop spending cores on this line.

---

## 7. Expected J

| outcome | condition | raw J | shown (submit 09-24/25) | probability |
|---|---|--:|--:|--:|
| t10d stands | G0.5 or G1 fails | 14.366 | 13.80 | ~0.90 |
| improvement | G1(b) and G2(b) pass | ~13.7-13.9 | ~13.55-13.75 | ~0.08 |
| J < 13 | G1(a) and G2(a) pass | ~12.9 | ~12.75 | ~0.02 |

**Probability-weighted raw J is about 14.3.** Why so low, against the judges' ~13.3:

1. **Capacity.**
   - The best settled 8-route pass is 273 at 10.48 (M-D).
   - The remaining levers each add a few targets at the 0.027 J/target frontier price: pilot choice, closability-aware
     premiums, more generations, joint K=2 for the last pair. None is measured above +3-8.
   - G1(b) needs +9 at no extra cost; G1(a) needs +17.
2. **Closing.**
   - The judges assumed ~0.04 J for every leftover.
   - Measured: 8 at 0.026 J, 9 at 0.17 J, and about half impossible (M-C).
   - A 9th route cannot take them (M-B), and re-plans cannot (M-A).
3. **Shown-J arithmetic.** An improvement must also beat the time decay: each day costs 0.049 J of raw margin.

**What would change the verdict.** A generator that flies the orphans *inside* deep routes at the pass stage without
losing depth. The leaders do this: HIT flies 8 x 37.2 at 9.45 kg/fb or less. The data say where to look:
- routes designed with the orphan set as a gentle-premium target (section 3.2 #11), so that the price probe's
  "unplaceable" list shrinks generation by generation;
- a planner that keeps sparse-pool legs settleable at dv 1.6-2.0. In M-D the deepest dv-2.0 tours failed to settle:
  route 6's 32, and route 8's 19-22 (six of six).

If the unplaceable list does not shrink over the first 3 generations, G0.5 will fail, and it should.

---

## 8. Do not retry (measured dead ends)

- **Planner re-plans as closers:**
  - pinned segment re-plan 0/35 (M-A);
  - free re-plan of an at-capacity route 0/15 (PAR M4);
  - assigned-pool re-plans 0.65-0.85 of the pool (linchpin);
  - Lloyd best response: union 259 -> 237.
- **A 9th leftover route:** 5/35 at dv 1.2; 0/8 settle at dv 2.0; 2/18 on the hard ones (M-B).
- **`run_ialns.py search --swap` / `relocate_pass(..., swap=True)`:** it corrupts coverage (M-C).
- **Premiums of 0.03 or more; pinned hard waypoints** (section 1.5).
- **In sparse pools:** stepping stones, longer TOFs, drmax 0.25, launch grid 0-1500 d, w_t 0.5, and beam width beyond
  100-300. Only the per-leg dv cap moves depth there (M-D).
- **MILP recombination as a coverage engine:** 259 = best single pass.
- **Insertion-cost surrogates as assignment prices** (R^2 0.010).
- **Anything touching `results/CTOC14_Result_t10d.txt` or `results/validator2_t10d.log`.**

---

## 9. Files

- This plan: `docs/stage13_solver_plan.md`.
- Tonight's probes, all new, with no shared tool modified, in `results/s13/plan13/`:
  - `leftover_dist.py` -> `N6a_leftovers.json`
  - `sr_fill.py` (wrapper around `tools/skel_fill.py`, fixed `fill_segmented` copy) + `run_sr_N6a.sh` -> `sr_N6a/NN/`
  - `sweeper_probe.py`, `sweeper_probe2.py` -> `sweep_dv12.json`, `sweep_dv20.json`, `sweep_all_dv20b.json`,
    `sweep_dead_dv20.json`, `N6a_after_close_leftovers.json`
  - `close_probe.py` -> `close_N6a/{closing.json, closed/, relocated/ (coverage-corrupted by --swap; do not reuse), log.txt}`, `close_N6a.out`
  - `tail_probe.py` -> `tail_N6a_k5/{log.txt, fleet/}` (the 273 @ 10.476 fleet, settled)
- Inputs:
  - FREE passes, preserved in `results/s13/seed_passes/` (originals in the session scratchpad `.../scratchpad/passes/`);
  - `results/s13/planner_probe/{jointsym.py, twin_wait.py}`;
  - `results/s13/linchpin/linchpin_numbers.json`;
  - `results/s10/insert/samples.jsonl`.
