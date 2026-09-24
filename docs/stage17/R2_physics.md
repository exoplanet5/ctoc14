# Stage 17 R2: PHYSICS lens — cross-examination (route physics, nodal events)

Agent PHYSICS, 2026-09-23 23:35-00:10 CST (machine clock). MEASURED numbers name their file (all under
`results/s17/`); everything else is est. Compute this round: one check, 2 cores, 3 min wall, plus offline reads.

## 0. My round-2 check: do pin lanes realise in the TWIN at the model price? PASS
- Pre-registered in `physics/r2/THRESHOLDS_R2.txt` before the run. Script: `physics/r2/lane_twin.py`.
  Results: `physics/r2/lane_twin.json` and `.out`.
- Pipeline for the 8 T1 carriers: pin-only events (eps 0.03) -> carrier DP (kphase 1.5) -> `carrier_dp.seed_ip`
  -> `insert_block` -> `settle` -> `s15b_regrid` (20 d bins, law 8) -> settle.

| lane | pins | first fb (d) | model dv | twin dv | twin tank | dv/pin |
|:--|--:|--:|--:|--:|--:|--:|
| 0 | 16 | 356 | 11.1 | 10.8 | 792 | 0.68 |
| 1 | 13 | 542 | 7.2 | 5.4 | 689 | 0.41 |
| 2 | 13 | 478 | 7.3 | 3.5 | 658 | 0.27 |
| 3 | 15 | 486 | 9.8 | 8.7 | 751 | 0.58 |
| 4 | 12 | 758 | 5.7 | 2.9 | 647 | 0.24 |
| 5 | 13 | 576 | 7.5 | 6.0 | 702 | 0.46 |
| 6 | 14 | 546 | 9.6 | 4.2 | 670 | 0.30 |
| 7 | 15 | 1189 | 7.9 | 14.6 | 873 | 0.97 |

- **Verdict: PASS.** All 8/8 settled (miss 5-68 km). Median twin/model is 0.78, and median twin dv is **0.44 km/s
  per pin**. A 13-pin skeleton costs a median 95 kg of fuel.
- **Meaning.** In a lane built for them, pins cost what any deep-route flyby costs (0.43 km/s/fb). R1's worry ("40%
  of a craft's budget before any fill") is withdrawn: it is about 24% (5.4 of 22.4 km/s).
- **Caveats.** The lanes use w=1 over all 95 pins, so they overlap (204 and 88 each sit in 4 lanes); the disjoint
  T1b set stored no carriers and was not realised. Lane 7 (cruise to day 1189) is the outlier at 1.8x the model.
  **The fill half of P1 is still untested**, and its history is the stage-5 result of 241.

## 1. F4 resolved: the three readings hold for different subsets. My proposals depend on "pin AND off-track"
Sources:
- `physics/events.json` and `pool_popularity.npz` (honest pool; "deep" means >= 36 fb);
- `campaign/miss_windows.json`;
- decomposition of the closest approaches: `physics/r2/decomp.py` -> `decomp.json` (s16a twins, 1 d step).

| hard core | 17 | 36 | 75 | 108 | 118 | 180 | 212 | 216 | 219 | 220 | 243 | 247 | 258 | 288 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| node events < 0.05 AU | 5 | 3 | 0 | 3 | 0 | 5 | 4 | 3 | 3 | 6 | 4 | 7 | **10** | **9** |
| all pool routes | 329 | 271 | 227 | 220 | 425 | 412 | 194 | 371 | 246 | 325 | 273 | 303 | 204 | 447 |
| DEEP pool routes | 19 | 12 | 12 | 10 | 41 | 36 | 52 | 37 | 10 | 35 | 13 | 35 | 104 | 53 |

- **CAMPAIGN D3 ("popular") counts shallow, steered routes.** Among deep honest routes the hard core has a median of
  35 carriers, against 58 for all targets. Five of the 13 (36, 75, 108, 219, 243) sit at or below the 4th percentile.
- **PHYSICS D2 ("pins") holds for 12 of 14.** The median is 4 events, against 6 for all targets. **258 and 288 are
  not pins** (10 and 9 events, and 104 and 53 deep carriers). They are pure packing misses; I concede them to CAMPAIGN.
- **REDTEAM 1.1 ("isolated") is a 3-D fact, not a phase fact.**
  - Median offsets of the 18 approaches below 0.15 AU: along-track 0.071, radial 0.034, **vertical 0.042 AU**.
  - Only 38% of hard-core windows lie within 15° of Earth-relative lag of a survivor lane. The baseline is 56%.
  - Several match in lag and still miss in 3-D. Example: 17@4959 vs r1 is 4° in lag, but its closest approach is
    0.11 AU (radial -0.079).
- **Price of bending a track through one hard-core window (est.).** Δi ≈ 2.4° is ~1.25 km/s, Δe 0.034 ~0.5 km/s,
  phasing 0.2-0.7 km/s: ≈ 2 km/s = 45-50 kg at 950 kg. That is the whole 46 kg budget (F5), and only if the bend is
  **never undone**. If the track must return to its own targets, the price doubles: why tails trade pin-for-pin.
- **What my proposals assume.** A hard-core target is cheap only on a track DESIGNED through its window. That is
  measured now (§0), not assumed.

## 2. Cross-examination

### CAMPAIGN
- **P1 (stitched children). FATAL as a coverage engine; CAMPAIGN already killed it.** A junction bridges tracks that
  overlap (their D5); the hard core is 0.04 AU off in z and 0.03 AU in r from every track (§1), which no ≤ 1 km/s
  junction crosses. Rescue: none for N=8.
- **P2 (joint prefix re-plan + exact tail MILP). FATAL for its own 292 gate.**
  - Depth arithmetic: the base-286 depths are 47+46+44+38+31+31+29+24 = **290**. That leaves 4 duplicates.
  - F1 and LNS measured that re-planned tails conserve depth: r6 30→30 and r8 24→24 at most. The 1-for-2 trades of
    s15b show the same.
  - So re-planned tails at conserved depth reach at most 290 < 292.
  - Rescue: (a) LNS Part B shows a tail gains ≥ 3 flybys; or (b) the base is chosen with Σ depth ≥ 302, so that
    swaps drop DUPLICATES. Rescue (b) is my E1.
- **P3 (concurrent beam + congestion V). MAJOR.** It needs a 12-20 h build before its first gate; joint beams already
  lost to sequential (K=2: 35 vs 40); the congestion it wants to price is mostly fillers (57% of the targets shared by
  two deep routes), so V must model the fillers' remaining windows. Rescue: none on this clock.
- **D3 as evidence against scarcity prizes: MINOR.** It is right for 258 and 288. For the other 12, the popularity
  comes from shallow steered routes (table §1).

### LNS
- **P1 (suffix-column set partitioning). FATAL as posed.** The physics adds WHY the chain screen fails.
  - The own targets a tail drops are themselves pins: 53 (5 events), 138 (2), 223 (3), 150 (5).
  - Tails hold 9-14 pins each at a fixed ~180 d pace, so a tail can only swap one pin for another.
  - Rescue: make the dropped target a duplicate (E1). Re-read that way, F1's links cost: r6 30_1 +2 coverage at
    +49.5 kg (25 kg each); r8 23_1 +2 at -70 kg; r6 29_2 +1 at +9.6 kg; r8 22_x +1 at -83 kg; r8 24 +2 at +306 kg.
    So 4 of 5 links come in at ≤ 30 kg per +1 coverage, but only IF another chosen route flies the dropped targets.
- **P2 (ramped emptying). FATAL.** With single-move columns it IS forced dissolution, 3.6x over budget. LNS agrees (~0).
- **P3 (prefix re-plan with a free launch). MAJOR.**
  - Physics is on its side in one respect: the launch v_inf (4 km/s ≈ 7.7° of free tilt) buys exactly the vertical
    offset (0.042 AU ≈ 2.4°) that makes the hard core dear.
  - But the kept suffix pins the end-of-prefix plane, so the tilt must be undone at ~1.25 km/s. The 6-8 h build does
    not fit the clock.
  - Rescue: only as the launch stage of a whole new lane (E2 already has free launch).

### REDTEAM
- **P1 (closure-aware MILP with lin arcs). MINOR.** A confirmatory kill, nearly done.
  - lin prices PINNED insertion: a 0.1 AU 3-D offset over a ~120 d gap cannot come in under 0.15 J (bump ~4Δr/T), so
    the arcs cannot reach the hard core by construction. It can only re-prove 285-287.
  - `full.out` crashed on its output path (its one line: N8/K294, time limit, no incumbent); restarted 23:34
    (`RESTART_NOTE.txt`). Let it finish; do not extend it.
- **P2 (13-target contract). FATAL.** F1 is its first kill test, already run. The one hard-core swap on record (247)
  pushes out 53, whose re-home costs 1617 kg. The mean dJ including re-homing is far above 0.08 J.
  - Rescue: the same as LNS P1, drops of duplicates (E1).
- **P3 (early lane separation). MAJOR.** Opening a lag lane Δφ in time T costs (2/3)·v·Δφ/(nT) in total: 90° in
  1 yr ≈ 5.0 km/s, v_inf pays half, leaving ~60 kg. While drifting the craft sits Δa/a = 2Δv/v ≈ 0.17 AU off the node
  torus and flies almost nothing. The "15 lost flybys" are mostly launch geometry every team shares (our first
  flybys: days 86-383). Rescue: only as a parameter of a new-lane generator (my carriers already launch over days
  0-730 with free v_inf).
- **Their #4 against me (leftover specialist) is correct, and I had killed it in R1.** One lane holds ≤ 16 pins, and
  the realised lanes carry 12-16 at 647-873 kg. A craft that catches 22 of the 30 isolated targets does not exist.

## 3. What the others' evidence changes in my proposals
- **P1 (8-craft pins-and-grout from scratch): 4% → ≤ 1%.**
  - The pin half passes (§0).
  - But LNS measured that depth is set by track capacity. Fill history is 241. The full rebuild (~30 h of pipeline)
    does not fit before 09-27 06:00 with polish.
  - Cut to scope as E2: heads kept, 4 tails.
- **P2 (EDF prize): dropped.** It is a prize inside sequential generators; law 3 and depth conservation make it moot.
- **P3 (surrogate): dropped this stage.** It is enabling only.
- **Concessions.** 258 and 288 are packing misses, not pins (to CAMPAIGN). Tails are not an orbit family (to
  REDTEAM); my classes are temporal, and their AUC 0.63 on #node passes agrees with my ρ 0.55.
  - "Pace binds" (my D7) and LNS's depth conservation are the same law. I merge them:
    **a route's depth is set by the density of its track; re-planning changes which targets it flies, not how many.**
- **Registered before LNS Part B finishes (23:43: depth 22, min t_end 3977 d).** Prediction: n* ≤ 32 (WEAK or
  KILL). Even a PASS is zero-sum: the 38 offered targets live on the r5/r7/r8 tracks, 0.1-0.25 AU away in 3-D.

## 4. Ranked experiments (max 3)

### E1. Dedup-swap columns (LNS P1 x CAMPAIGN P2 x PHYSICS D4). Rank 1
**Mechanism.**
- Pick an 8-route BASE by MILP with deliberate overlap (Σ depth ≥ 302), so that duplicates exist, mostly fillers.
- On the base routes whose suffix passes a miss window (approach ≤ 0.12 AU after the cut), run trading suffix
  re-plans (`s15b_prefix`) with prizes by multiplicity: unique own 2.0, duplicated own 0.2, misses 1.5.
- The exact master is `s16_select` over base + columns. Then `s15b_close`, `s16_iter`, relocate and export.

**Laws escaped.**
- Depth conservation: it needs only 1-for-1 swaps, which F1 showed exist at 25 kg.
- Isolation conservation: the dropped targets are covered elsewhere, so their re-home cost is 0.
- Law 1: every miss is chosen by a re-planned suffix.
- Not escaped: law 7. If the base's suffixes do not pass the misses, the gate says so.

**Gate (thresholds fixed now).**
- Step 0: MILP, 15 min, 2 cores. Choose ≤ 8 routes (honest pool + s16 fleet + F1 columns), maximising C0, subject to
  D0 ≥ 302 and fuel + 30·(298 - C0) ≤ 3723 kg.
  - PASS-0: feasible with C0 ≥ 280.
  - KILL: infeasible even at 46 kg/swap with C0 ≥ 270.
- Step 1: 4 trading re-plans, beam 12, 90 min, 4 cores.
  - PASS: ≥ 6 distinct misses gained through drops of duplicated targets, at ≤ 30 kg net per +1 coverage (net = tank
    delta).
  - KILL: < 3 misses at ≤ 46 kg.

**Resources.** Cores 4-8. Build 3 h (job generator + prize map; tools exist). Run ~14 h (3 rounds) + 8 h
close/polish + 4 h export. On the clock it could be ready by about 09-25 06:00, against the lenient 13.658 bar.

**Probabilities.** P(gate passes) **0.12**. P(submittable by 09-27 06:00) **0.02**.

### E2. Heads-fixed K=4 pin-lane tail rebuild (my P1 cut to scope). Rank 2
**Mechanism.**
- Keep r1-r4, which have 162 targets and 1364 kg of fuel.
- Rebuild the 136-target residual with 4 craft. First, 4 jointly chosen pin lanes over the residual's 58 pins,
  realised in the twin with the §0 pipeline (20 s per lane). Then a class-quota fill (fillers ≤ 5 per tail) with
  `skel_fill` segmented mode.
- Budget at the 09-27 bar: each tail ≤ 1177 kg twin, i.e. 34 flybys at ≤ 17 kg/fb, 0.78 km/s/fb.
- **Fuel is slack; pace is not:** 34 flybys means ≤ 150 d per leg.

**Laws.**
- Escapes law 1: pins are designed in at a measured 0.44 km/s.
- Escapes law 3: the 4 lanes are chosen jointly before any fill.
- Uses law 7: the extra fuel per leg widens the affordable radius.
- Does not escape depth conservation, if the residual is too sparse.

**Gate (thresholds fixed now).**
- Step A: 15 min, 2 cores. 4-lane exact coverage of the residual pins, realised in the twin.
  - PASS-A: ≥ 50/58 pins covered, mean lane tank ≤ 760 kg.
- Step B: build 3 h, run 90 min, 2 cores. Fill the richest lane over the residual minus the other lanes' pins.
  - PASS: ≥ 34 flybys at a regridded twin tank ≤ 1177 kg.
  - KILL: ≤ 30 flybys, or a tank > 1250 kg.

**Resources.** Cores 2 (gate), 8 (run). Build 4 h. Run ~24 h + export.

**Probabilities.** P(gate passes) **0.10**. P(submittable) **0.015**.

**Warning.** REDTEAM measured that the best single pool route over the tail set flies 33-35. Here each of 4
DISJOINT routes must reach 34. That is why P is low.

### E3. Stop rule (not a compute job)
- Run E1 step 0 first; it takes minutes and kills cheaply. Run E2 step A in parallel.
- **If neither gate passes by 09-24 20:00, end stage 17.** s16E stays banked at 13.5116. Never resubmit a raw above
  13.5116/k(t).
- If one passes, no new builds start after 09-25 18:00.

## 5. The run-nothing case
**The strongest argument for stopping.**
- Four machineries now agree that existing tracks give 285-286 at 8 craft, and that they conserve both depth and
  isolation.
- Every remaining path needs new track, and our last 16 stages of new-track generators landed at 228-249.
- The combined P(submittable) of E1 and E2 is about 3.5%. With a gain of ~0.15 J shown, the EV is about 0.005 J.
- A 9-craft rebalance needs -76 kg by 09-25 12:00 just to break even. No lever for it is on the table.

**Whether I endorse it: partly.**
- The cores are otherwise idle, and both gates cost ≤ 2 h.
- So: run the gates, and **stop at 09-24 20:00 without a pass**. Beyond the gates, the EV does not justify the user's
  attention.
