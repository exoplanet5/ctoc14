# Stage 17 judge ruling: RUN NOTHING. Stage 17 closes now (2026-09-24 00:30 CST)

Inputs: stage17_brief.md, R1_* and R2_* of CAMPAIGN / PHYSICS / LNS / REDTEAM, facts F1-F6, and the result files
named below (re-read on disk by the judge). The judge ran no new compute beyond reading and arithmetic.

## 0. Decision
- **run_nothing.** No new build, no new gate. s16E stays banked (shown 13.511581). Nothing is submitted.
- **Why, in one line:** every candidate left has a blocker that is already on disk. Their summed EV is about
  **0.002 J shown**, against 0.009 J to reach rank 9. And the dominant failure mode of every one of them is the
  3x forced premium measured again.
- **In-flight work is finished; no process is running** (ps at 00:10 shows only Jupyter kernels):
  - LNS `L2r6k17O` done: n* = 31 at 937.4 kg (`lns/r2/L2r6k17O.out`, `open_r6.json`);
  - CAMPAIGN `kx9` done: gain 0.000 on 129/129 exchanges (`campaign/r2/kx9.out`);
  - REDTEAM `full2` done: no incumbents (`redteam/p1/full2_*.out`).
- **Every pre-registered reopen trigger has resolved negative** (CAMPAIGN E0, REDTEAM E1):
  - (a) needs n* >= 33 AND a hard-core column at <= 46 kg net. Both conjuncts fail: n* = 31, and the best link
    saves -266 kg (`lns/r2/chain_screen_B.json`).
  - (b) and (c) have no incumbent, or kx9 shows a gain of 0.

## 1. Disagreements, resolved on evidence
1. **F4 (popular vs pins vs isolated).** All three are true, each on a different measure. The binding measure is
   the number of fleet-compatible DEEP carriers.
   - The misses are flown by 182-436 of all 2349 sets, but by only 10-52 deep sets (2nd-47th percentile). 1444 of
     the sets come from s15/s15b steering.
   - The offset from every track is 3-D (vertical 0.042 AU, `physics/r2/decomp.json`).
   - 8 of the 13 are pins. 258 and 288 are packing misses.
   - Every surviving proposal depends on the reading "isolated from fleet-compatible deep tracks". So does this
     ruling.
2. **L2r6k17O: WEAK (LNS) or KILL?**
   - KILL on the pre-registered primary metric (n* <= 31).
   - The wall projection (t26 = 4501 d) says WEAK. But the depth-31 front ends on days 5303-5406, and the window
     closes about day 5479, so the ceiling is 32 < 33 under either reading.
   - Depth responds to a 2x pool by +1 flyby at most. Law 9 stands.
3. **PHYSICS P1 gate (262 on the fill2 A/B).** It is decoupled from the objective (CAMPAIGN, LNS and REDTEAM are
   right).
   - At the forced closer price, 10.5K + 120(298 - K) <= 3723 kg needs K >= 293.
   - A PASS at 262 is 24 below the 285-286 already held.
   - Pin lanes used as waypoints are the skeleton-then-fill dead end (241).
4. **Dedup-swap (PHYSICS E1) vs E[slots] ~ 2 (REDTEAM).** Settled by existing data, with no new run needed.
   - `redteam/r2/forced_455.out` shows that the best 8-cover containing 36, 75 or 220 is **278**: only s16a r9
     carries them in the 455-column + 3671-arc space.
   - So every base with C0 >= 280 misses all three, and each must enter by a swap on a host passing <= 0.12 AU.
   - `campaign/r2/chain_slots.json`: 36 has 0 hosts within 0.2 AU; 75 has one (r8, 0.114 AU) with no slot; 108 and
     216 have 0 hosts.
   - Closure is therefore unreachable by dedup swaps, whatever the step-1 trade price. Not run.
5. **LNS E1 (9-craft tail cyclic exchange) vs CAMPAIGN kx9.** LNS is right that kx9 never saw prefix-anchored
   re-plan columns. The columns that do exist show the cycle has no raw material:
   - every column cheaper by more than 15 kg (-51/-54/-58 kg) drops 166 (no other route within 0.25 AU) or 53
     (1617 kg);
   - own-only re-plans re-fly r6 at +1.8 kg and cannot re-fly r8 (best 9/11);
   - removal refunds have a median of 1.7-3.2 kg (`redteam/removals_all.json`).
   - The clock also works against it: >= 77 kg by a 09-25 12:00 submit only breaks even; -100 kg earns 0.024 J.
     Not run.
6. **Heads-fixed 4-lane tail rebuild (PHYSICS E2 / REDTEAM E3).** REDTEAM's objection holds: gating on the richest
   lane is uninformative.
   - The EXACT best pool 4-cover of the tails' 136 is **116** (`redteam/milp4_all.txt`, HiGHS optimal, over the
     steered pool).
   - 4 new lanes must beat that optimum by +20. Each must fly about 34, which is the pool's single-route maximum
     (33-35). The fill is sequential (law 3), and its history is 241/298.
   - The pin half is real: 8/8 lanes settle at 0.44 km/s per pin (`physics/r2/lane_twin.out`). But those lanes
     overlap: 111 slots cover 68 distinct pins.
   - Not run. It is kept as the only user-override lottery (section 5).
7. **Budget.**
   - F5 holds: 3723 kg of fleet fuel at the 09-27 bar, about 46 kg per hard target. The brief's 3858 kg is the
     09-24 bar.
   - REDTEAM P2's 0.08 J kill line was mis-set: a pass must be <= 0.037 J per target (conceded).
8. **REDTEAM P1 full model (intractable) vs LNS E2 (warm start).** Not needed.
   - The 455-column and forced-X runs are HiGHS-optimal and settle the existing-tracks question: 285, and a net
     -1 to -7 for every forced miss.
   - Even a hypothetical K = 290 leaves 8 hard targets at >= 140 kg each (REDTEAM lin calibration 0.73), about
     1100 kg against a 46 kg-per-target budget.
9. **REDTEAM P3 (early lanes).** PHYSICS is right: opening a 90 deg lane in 1 yr costs about 5 km/s, and the
   drifting craft sits 0.17 AU off the torus. Dead (conceded).
10. **CAMPAIGN E2 concurrent build.** Dead on the clock:
    - a 14 h build, and a day-1461 gate that is only necessary, not sufficient;
    - a tuple beam is a weaker master (joint K=2 scored 35 vs 40);
    - 73% of conflicts are cross-epoch.

**Law 9 (adopted, measured 5 ways): re-planning an existing track conserves the isolation of what it flies and
nearly conserves its depth (+0/+1).**
- The five measurements: F1, L2r6k17O, chain_screen, forced_455 and CAMPAIGN's X3 detours.
- Existing machinery only chooses which tail's hard set is orphaned. r9's 13 is the smallest such set, hence
  285-286 everywhere.

## 2. EV against the clock
Bars: raw < 13.511581 / k(t), where k = 1 - 0.1 (09-28 12:00 - t)/28 d. Twin bar = raw bar - 2.0172.

| submit | k | raw bar | twin sum J_i bar | 9-craft saving needed vs s16a |
|:--|--:|--:|--:|--:|
| 09-24 18:00 | 0.98661 | 13.695 | 11.678 | 42 kg |
| 09-25 12:00 | 0.98929 | 13.658 | 11.641 | 77 kg |
| 09-26 12:00 | 0.99286 | 13.609 | 11.592 | 123 kg |
| 09-27 12:00 | 0.99643 | 13.560 | 11.543 | 169 kg |

What a success would be worth (arithmetic):
- **8 craft at 298:**
  - 3600 kg of fleet fuel, submitted 09-26 12:00: shown 13.319, a gain of 0.19;
  - 3723 kg at 09-27: a gain of 0.
- **9 craft:** -120 kg submitted 09-25 12:00 gains 0.045.

| best remaining candidate | blocker already on disk | P(gate) | P(submittable by 09-27 06:00) | gain if yes | EV (J shown) |
|:--|:--|--:|--:|--:|--:|
| heads-fixed 4-lane tail rebuild | pool 4-cover 116/136 is exact; law 3 on lanes 3-4 | 0.03 | 0.006 | ~0.15 | 0.0009 |
| concurrent from-launch build | 14 h build; joint beams lost; gate only necessary | 0.10 | 0.004 | ~0.2 | 0.0008 |
| tail cyclic exchange (9 craft) | cheaper columns drop 166/53; kx9 3-exchange-optimal | 0.03 | 0.012 | ~0.03 | 0.0004 |
| dedup-swap columns | 36/75/220 only via r9 (278); 36 has no host within 0.2 AU | 0.05 | <= 0.002 | ~0.15 | <= 0.0003 |

- **Sum about 0.002 J shown.** P(any rank gain by 09-28) is about 1%.
  - Rank 9 needs -0.009 J shown (13.503), and the board is also moving.
  - With 9 cores for 3 days and zero score downside (the never-submit-worse rule), the only real cost is the
    user's supervision.
- **What running would cost.** Each candidate's modal outcome is a fourth measurement of law 1 / law 9, which is
  exactly the standing risk:
  - dedup and cyclic exchange are re-plans of existing tracks;
  - the lane rebuild and the concurrent build use waypoint-then-fill.

## 3. Timeline to 09-27 06:00 (run_nothing)
- **09-24 00:30.** Stage 17 is CLOSED. No jobs are running or may be started by agents. s16E stays banked at
  13.511581.
- **Submission rule (all dates).**
  - Never submit a file whose raw J >= 13.511581 / k(t_submit) (table above).
  - s16a (raw 13.7399) and any 9-craft re-export of it fail every bar from now on: 13.7074 at 09-24 12:00,
    falling.
  - The user submits by hand; agents never submit.
- **09-24 12:00.** Last point at which the user-override lottery (section 5) can start and still land before
  09-27 06:00 with closer and export. After this, no structural work of any kind.
- **09-27 06:00.** Absolute last twin fleet (export plus validation takes 3-4 h). With no override, nothing happens
  at any checkpoint.

## 4. Global stop rule
- **Stage 17 ends with this ruling.** No agent reopens it.
- It reopens only if the user explicitly says so before 09-24 12:00, and then only for section 5's track, with its
  thresholds exactly as fixed below. No extensions, and no second track.
- Any KILL on that track ends the contest campaign. s16E stays banked.

## 5. User override only (NOT recommended; it runs only on the user's explicit go)
**Track: the heads-fixed 4-lane tail rebuild.** It is the only proposal that designs new capacity for the
isolated set instead of reshuffling it.

**Build (new files only):**
- `tools/s17_pinlanes4.py`:
  - residual = 298 minus the targets of s16a r1-r4 (136; 58 pins by `physics/events.json`);
  - pins-only carrier events via `carrier_dp` (dp_route, eps 0.03, kphase 1.5), column generation, then an exact
    4-lane DISJOINT max-coverage MILP (scipy/HiGHS);
  - each lane realised as in `results/s17/physics/r2/lane_twin.py`: `carrier_dp.seed_ip` ->
    `globalopt.insert_block` -> `run_ialns.settle` -> `s15b_regrid.regrid` (20 d) -> settle.
- `tools/s17_quotafill.py`:
  - a copy of `tools/skel_fill.py` plus a copied `ctoc14/search.py` beam carrying a class counter (fillers <= 5 per
    tail);
  - segmented mode, lanes filled pin-richest first, each over the residual minus the other lanes' pins;
  - every tail regridded (law 8).

**Gate (<= 3 h wall, 8 cores).**
- Step A (15 min): >= 50/58 pins in 4 disjoint lanes, at a mean twin tank <= 760 kg. KILL if < 45.
- Step B (all 4 lanes filled, gated on the LAST lane): PASS if all three hold:
  - residual covered >= 128/136;
  - 4-tail regridded fuel <= 2120 kg (the 09-26 bar allows 2346 including the closer);
  - last lane >= 30 fb.
- KILL if any one holds: <= 120 covered, > 2400 kg, or last lane <= 26 fb.

**On PASS:**
- `s15b_close` (cap ladder) for the <= 8 leftovers;
- `s15b_relocate`, then `s16_iter` single moves;
- regrid all 8 routes;
- `s16_select` exact at N=8, K=298;
- export with the `results/s16/export_s16a.sh` pattern (run_ialns export -> combine_submission -> validator2);
- submit only under the bar at the moment of submission.

**Odds:** P(gate) about 0.03; P(submittable) about 0.006; EV about 0.001 J shown.

## 6. For the record (memory)
- **Stage 17 CLOSED 09-24 00:30.** s16E banked at 13.5116. The 8-craft ceiling on existing tracks is 285-286,
  measured by 3 machineries.
- **Law 9** (isolation and depth are conserved under re-planning).
- **Pool exact 4-cover of the tails is 116/136.**
- **s16a is 3-exchange-optimal** over 42,650 columns, children included.
- **The only unbuilt lever is new capacity for the isolated set, designed from launch.** Its history (228-249) and
  the clock rule it out for this contest.
