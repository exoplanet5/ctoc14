# Stage 17 R2: RED TEAM cross-examination (2026-09-23 23:30-23:58 CST, machine clock)

MEASURED numbers name their file (under `results/s17/`). EST marks estimates. Compute: one check on 2 cores for
~20 min (`redteam/r2/forced_cov.py`); the P1 frontier finished in its own 4-core slot (now released).

## §0 Data first

**P1 frontier, final** (appended as §4 of R1_redteam.md).
- `full.json` crashed at 23:33 on a path bug (the out argument was doubled).
- Its restart (`full2_*`) showed the full 1382-column / 24-27k-arc MILP is intractable. After 900 s the N=8
  incumbents covered 26 and 81, and N=9 had no incumbent at all.
- On the tractable F2 set (`p1/first455.json`, all runs optimal):
  - N=8 max coverage is **285 with arcs ≤ 0.08 J AND with arcs ≤ 0.15 J**, and it is the same solution both times;
  - K ≥ 286 is infeasible;
  - N=9 at K=298 gives **11.7227, which is s16a exactly**, at both thresholds.
- My pre-registered kill fires, so P1 is dead, both for 8 craft and as a 9-craft lever.
- CAMPAIGN's reopen triggers (b) and (c) cannot fire.

**My one check: is the hard core SPECIFIC, or does it ROTATE?** (`redteam/r2/forced_455.json`)
- Setup: the F2 set (455 columns, 3671 arcs, lin ≤ 0.08 J, cap 150 kg per host).
- Method: max coverage at N ≤ 8 with y_X forced to 1, for each of the 13 misses of the 285-optimum.
- All 14 runs reached the HiGHS optimum.

| X forced | covered (Δ vs 285) | X carried by | newly missed |
|--:|--:|:--|:--|
| 247 | 284 (-1) | pool 657888c0 | 35, 212 |
| 17 | 283 (-2) | pool 2fa05bfc | 35, 55, 212 |
| 118 / 243 | 282 (-3) | pool 4fcccd74 / 79f07202 | 45, 65, 200, 285 / 21, 35, 45, 55, 200, 212 |
| 108 / 288 | 281 (-4) | pool f783998b / e37b2b67 | 35, 41, 55, 212, 214, 260 / 21, 45, 55, 198, 200, 212 |
| 219, 258 | 280 (-5) | pool 400d6af8 (219 by a 40 kg arc) | 47, 153, 171, 179, 187, 196, 205, 212 |
| 180, 216 | 279 (-6) | pool 8bb829fc | 45, 92, 103, 123, 197, 198, 200, 212, 250, 286 |
| 36, 75, 220 | 278 (-7) | s16a r9 itself | 20 of r8's 24 (r8 drops out) |

**None of the 13 misses enters an 8-craft cover without a NET loss of 1-7 targets.**
- For 36, 75 and 220, the only way in is r9, and then r8 falls out.
- For the other ten, the carrier pushes out a second tier of isolated targets; 35, 45, 55, 200 and 212 recur.
- So isolation is conserved at FLEET level, as LNS found for single links (`lns/r2/chain_screen.json`: 247 in, 53
  out, and 53's forced price is 1617 kg).
- (CAMPAIGN R2 quotes this as "costs 3-8 others". The net is -1 to -7.)

**F4 resolved (honest pool, 2349 distinct sets).** All three readings are true. The DEEP-route count is what binds.
- **All routes:** the misses are flown by 182-436 sets, against a median target's 149. "Popular" (CAMPAIGN D3) holds.
- **Deep routes (≥ 36 fb):** only 10-51 sets fly them, the 2nd-47th percentile (median target 54). The one
  exception is 258, with 99 sets.
- The popularity comes from s15/s15b premium steering: 1444 of the 2349 sets.
- 8 of the 14 are pins (PHYSICS `events.json`). PHYSICS adds the geometry: the median offset is 0.042 AU VERTICAL.
- **Reading used by all my proposals:** isolated, and fleet-incompatible at depth. Taking a miss costs a route depth
  (law 7), so every 8-cover leaves about 13 out.

## §1 Objections to the other lenses

### CAMPAIGN
**P1 (fragment stitching): FATAL. It declared its own FAIL, and I agree.**
- M1 stays at 286 under exact k ≤ 3 refill.
- Its own D5 explains why: junctions join only tracks that already overlap, and just 2.3% of joined pairs are disjoint.
- Rescue: none for coverage. As a 9-craft column source it is now also dead: children at K=298 still give s16a,
  11.7227 (their r2).

**P2 (joint prefix re-plan + exact tail MILP): MAJOR, near fatal. Now conceded.**
- F1 already ran P2's unit operation on two of its four tails, with 50x its 1.02 premium. Every r9-carrying state
  is a trade, at a saving of -118 kg (best) to -567 kg.
- An exact MILP over 6^4 tail products is a subset of my column MILP. That MILP loses 1-7 targets per forced miss.
- Its 3858 kg PASS budget is really the 09-24 bar; at 09-27 the budget is 3723 kg.
- Rescue: only a SINK tail (keeps all own targets, adds a miss, costs ≤ +46 kg). LNS Part B is the last test of that.

**P3 (concurrent fleet beam): FATAL on the clock.**
- It needs a 12-20 h build.
- Any reachability value V rates the misses "safe", because shallow routes make them popular.

### PHYSICS
**P1 (pins and grout): MAJOR. Its gate is decoupled from submittability.**
- Credit where due: `physics/r2/lane_twin.out` settles 8/8 lanes at a median of 0.44 km/s per pin (twin/model 0.78).
  The pin half is real.
- But those 8 lanes are not a partition: 111 pin slots cover only 68 distinct targets (lanes 2 and 4 share 7). The
  disjoint T1b set (87/95) is an unproven incumbent and was never twin-checked.
- PASS at 262 leaves 36 leftovers. At the forced price that is 3.6-5.0 t, against a 961 kg budget.
  - Stage 5 set the bar itself: "the fill itself must reach ≥ ~280" (docs/stage5_skeleton_fill.md §6).
  - The best constructive 8-craft result is 249.
- The quota redistributes depth; it does not create it. Fill2 went 36/36/33/33/30/25/25/23 = 241, and balancing that
  still gives 241. D4 (fillers buy depth, corr 0.56) cuts both ways.
- Rescue: PASS ≥ 275 at ≤ 3200 kg regridded fuel, KILL ≤ 258. The heads-fixed scope (their E2, my E3) fits the clock
  better.

**P2 (EDF deadline prize): MINOR, and it yields no fleet.**
- It is a prize (the dead family).
- A pin's last window is the same calendar window for every route, so the prize synchronises routes onto it: the
  cross-epoch conflict of D4, on purpose.

**P3 (leg surrogate): MINOR, and irrelevant to the clock.**
- It is trained on CHOSEN legs, so it cannot price forced legs (law 1; stage 10 measured R² 0.01).

### LNS
**P1 (suffix-column set partitioning): FATAL for 8 craft. LNS now says ≤ 1%.**
- There is no sink and no positive link (F1, chain_screen).
- The forced-X check is the fleet-level proof: even the global MILP over 455 columns and 3671 arcs cannot add a miss
  without a net loss.
- Remaining gap: r5 and r7 are untested, and so are cuts before day 2500. But the depth law (30→30, 24→24 with 21
  extras on offer) predicts the same result.

**P2 (ramped emptying): FATAL.**
- Its own census predicts net +0.10 and +0.40, against PASS lines of 0.06 and 0.20.
- Without P1's columns it IS forced dissolution, 3.6x over budget.

**P3 (free-launch prefix): FATAL on the clock.**
- It needs a 6-8 h build.
- 17 lies 0.111 AU from every track, 36 lies 0.246 AU, and 108 has no approach within 0.25 AU.
- The kept suffix pins the plane, so the launch tilt must be undone (PHYSICS: about 1.25 km/s).

## §2 What the others' evidence does to MY proposals
- **P1 (closure-aware MILP): DEAD**, by its own kill (§0). As a 9-craft lever it returns s16a exactly.
- **P2 (13-target contract): DEAD.**
  - F1 plus chain_screen refute its premise: re-planned suffixes take isolated targets only by trading out equally
    isolated own ones.
  - **Concession to CAMPAIGN:** my kill line (mean > 0.08 J) was mis-set. 13 × 0.08 = 1.04 J, against a 0.48 J
    budget. The pass line should have been ≤ 0.037 J.
- **P3 (early lanes): DEAD.**
  - **Concession to PHYSICS:** opening a 90° lane in 1 yr costs about 5 km/s, and a drifting craft sits about 0.17 AU
    off the node torus. The 28-flyby early deficit is shared launch geometry.
- **Merged law (MEASURED four ways): ISOLATION IS CONSERVED.**
  - Each tail owns 10-16 hard-core targets that no other track reaches cheaply (census_r5..r9).
  - The four measurements: single links (LNS), fleet MILP with forced misses (mine), stitched corridors (CAMPAIGN D5),
    and slot census (CAMPAIGN chain_slots: 6 of 13 have any drop slot).
  - Existing machinery only chooses WHICH tail's set is orphaned. r9's 13 is the smallest, hence 285-286 everywhere.
  - PHYSICS's "depth is set by track density" is the same law in flyby counts.

## §3 Ranked next experiments (≤ 3)

The twin ΣJ_i bar is the same at 8 or 9 craft (raw = twin + 2.0172):

| submit | twin bar | vs s16a (11.7227) |
|:--|--:|--:|
| 09-24 12:00 | 11.690 | -0.033 |
| 09-25 12:00 | 11.641 | -0.082 |
| 09-27 12:00 | 11.543 | -0.180 |

Rank 8 (13.493) on 09-24 would need -0.051 J (about -49 kg). The 9-craft frontiers (mine and CAMPAIGN's) show none.

### E1 (rank 1): bank s16E and close the in-flight gates on pre-fixed triggers. No new build.
- **Mechanism:** LNS `L2r6k17O` and CAMPAIGN `kx9` finish on the cores they already hold.
- **Law escaped:** none. This rung confirms the law at about zero cost.
- **Gate:**
  - PASS if (a) n* ≥ 33 AND a settled column carries a hard-core target at net (Δtank + census re-home) ≤ 46 kg, or
    (c) kx9 finds a 298 cover at ≤ 11.64. A PASS re-opens "realise the columns" (CAMPAIGN E1).
  - KILL if neither fires by 09-24 01:00.
- **Resources:** 0 extra cores, 0 h build, ≤ 1 h run.
- **Odds:** P(gate) 0.04, P(submittable) 0.005.
- **Prediction:** n* ≤ 32. At 23:51 Part B stood at depth 26 with min t_end 4501 d, against F1's 4697. It is
  faster, but not by 3 flybys, and it is zero-sum over r5/r7/r8.

### E2 (rank 2): PHYSICS E1 "dedup-swap", step 0 only, with a SLOT criterion
- **Mechanism:**
  - Base MILP: ≤ 8 routes from F2 + F1 columns, arcs ≤ 0.06 J. Maximise C0 subject to Σ depth ≥ 302 and
    fuel + 30·(298 - C0) ≤ 3723 kg.
  - Then a join: for each base miss X, count hosts that pass ≤ 0.12 AU of X at t_X AND fly a DUPLICATED target within
    ±250 d of t_X.
- **Law escaped:** isolation conservation. A dropped duplicate costs 0 to re-home, but only if such a slot exists.
- **Objection priced in (EST):**
  - CAMPAIGN chain_slots gives any cheap slot for only 6 of 13 targets; 36, 108 and 216 have no host within 0.2 AU.
  - About 16 duplicates over 8 routes gives p ≈ 0.17 of one falling in a ±250 d window.
  - So E[slots] ≈ 2, against the ~11 needed.
- **Gate:**
  - PASS: C0 ≥ 280 AND ≥ 9 base misses with a duplicate slot.
  - KILL: C0 < 275, or ≤ 5 slots.
- **Resources:** 2 cores, 1 h build (p1_milp + chain_slots join), 0.5 h run.
- **Odds:** P(gate) 0.05, P(submittable) 0.008.

### E3 (rank 3, a 3-day lottery; only on the user's explicit go): PHYSICS E2, heads-fixed pin-lane tail rebuild, gated on the LAST lane
- **Mechanism:**
  - Keep r1-r4: 162 targets, 1364 kg.
  - Rebuild the 136-target residual with 4 jointly chosen pin lanes (0.44 km/s per pin, measured), then a quota fill
    (fillers ≤ 5) in sequence.
- **Law escaped:** isolation capacity is designed in from launch. It is the only proposal that attacks the conserved
  quantity rather than reshuffling it.
- **Objection to PHYSICS's own gate:** it tests the RICHEST lane. The pool's best single tail-set route already flies
  33-35 (R1 §1.3), so the first fill always passes. Law 3 bites on lanes 3-4, and the pool's best 4-cover of these 136
  targets is 116 (`milp4_*.txt`).
- **Gate:** fill all 4 lanes, regridded twin.
  - PASS: residual ≥ 128/136 AND tail fuel ≤ 2120 kg (leaving 240 kg for 8 closer insertions) AND last lane ≥ 30 fb.
  - KILL: ≤ 120 covered, or > 2400 kg, or last lane ≤ 26 fb.
- **Resources:** 2 cores for the gate; 4 h build, 6 h gate run; then about 24 h of production on 8 cores.
- **Odds:** P(gate) 0.05, P(submittable) 0.01.

**Order and stop rule.** E1 is running; start E2 now; E3 only if the user asks. If nothing passes by **09-24 12:00**,
stage 17 ends: s16E stays banked, and no file with raw > 13.5116/k(t) is ever submitted.

## §4 Run-nothing case: I ENDORSE it once E1 and E2 close

**The five measurements agree.**
- Selection gives 286, stitching 286, and arcs ≤ 0.15 J give 285.
- Forcing any miss costs -1 to -7.
- Suffix re-plans conserve both depth and isolation.
- s16a is optimal at 9 craft in the arc-extended and the stitched column spaces.

**The odds do not pay.**
- Every remaining path is a new-track build. Constructive history for such builds is 228-249.
- Summed P(submittable) over E1-E3 ≈ 0.023, for a gain of 0.1-0.2 J shown. EV ≈ 0.003 J, a third of the gap to
  rank 9.
- The submission rule makes the score downside zero. The real cost is the user's attention over 3 days, and that is
  not worth 0.003 J beyond the two cheap gates.
