# Stage 17 R2: LNS lens (set partitioning, ejection chains), cross-examination

Agent LNS, 2026-09-23 23:30 to 09-24 00:10 CST (machine clock).
MEASURED numbers name their file (under `results/s17/`); the rest is est. Twin bar at 09-27 12:00: sum J_i <= 11.543
(F5); dJ/dkg ~0.0011 at m0 ~1000, so the designed law (0.032-0.042 J) is 28-37 kg per target, the forced law 100-140 kg.

## 0. My test and my round-2 check (they decide most of what follows)

**F1: my §3 test FAILED** (`lns/THRESHOLDS.txt`; now written into R1 §3). Keep-own suffix re-plans of r6 and r8, with
r9's 21 targets prized 1.5, produced no state that keeps every own target and adds an r9 target. Neither host ever
flew more flybys than its original route: r6 stayed at 30, r8 at 24.

**The one check, part A: do the F1 trades chain?** (`lns/r2/chain_screen.py` -> `chain_screen.json`; offline.) Each
saved F1 column is an ejection-chain link; a dropped own target is re-homed at its census price (best lin host, not
the column's host, not r9). Saving = forced price of the added r9 targets - (tank change + re-homing).

| link | tank change | in | out (forced re-home) | net | saving |
|:--|--:|:--|:--|--:|--:|
| r6 `30_1` | +49.5 | 184, **247** | **53 (1617 kg, 0.168 AU)**, 170 (55) | +504 | **-118** |
| r6 `29_2` | +9.6 | 142 | 53, 170 | +464 | -446 |
| r8 `23_1` | -70.2 | 142, **212** | **138 (no approach <= 0.25 AU)**, 150 (269), 223 (65) | +664 | -450 |
| r8 `22_1/2` | -83 | **212** | 150, 173 (412), 223 / 138, 150, 223 | +650 | -455 |
| r8 `24` | +306 | **118, 288** | 138, 172 (293) | +999 | -567 |

**Part B: does a tail get DEEPER when offered the joint residual?**
- Pre-registered in `lns/r2/THRESHOLDS.txt` before the start. Job `job_L2r6k17O.json`; analysis `analyse_open.py`
  -> `open_r6.json`, `chain_screen_B.py` -> `chain_screen_B.json`. 2 cores, 30 min.
- Same r6 prefix as F1 (17 flybys, day 3328). Pool = CAMPAIGN P2's residual, 72 targets (F1: 34): own suffix 13 +
  r5/r7/r8 targets flown after day 3186 (38) + r9's 21. CAMPAIGN's prizes (own 1.0, r9 1.02). Beam 12/24 (F1 16/32).
- **n\* = 31 flybys at 937.4 kg**, against the original 30 at 901.5. The wall budget stopped the beam at depth 31. Every
  depth-31 state ends on day 5303-5406, at most 176 d before the window closes, so 33 was out of reach.
- The earliest end at depth 26 is day 4501, against 4697 in F1. That 196-d lead in mid-route turned into +1 flyby.
- **Pre-registered verdict: WEAK** under the wall rule (t26 falls between 4300 and 4600). On the primary metric, n\*
  = 31 is on the KILL line. PASS (>= 33) was unreachable.
- **Every cheaper or deeper state again pushes out an isolated own target:**
  - depth 30 at 850.5 kg (-51) drops **166, which no other route passes within 0.25 AU**, plus 170;
  - depth 31 takes 219 and drops 53 (1617 kg), 256, 40 and 170;
  - column `29_2` carries 5 r9 targets, 17 and 247 among them, at -5.9 kg, but drops 6 own targets, 53 and 166 among
    them.
  - Link savings are **-266 to -677 kg** (`chain_screen_B.json`).

**What this settles.**
- Depth responds weakly to the pool: +1 flyby over 13 suffix flybys, from 2x the candidates.
- **Isolation is conserved exactly.** No link has a positive saving, and no host is a sink (a column that absorbs a
  target without dropping one). An exact master can only chain links, and chains with no sink and no positive link
  do not close.
- On the orchestrator's F1 question: the r6 trade {184, 247} for {53, 170} cannot be chained cheaply. 53 is the
  problem, not the solution.

## 1. Cross-examination (strongest objection per proposal)

### CAMPAIGN
- **P1 (stitching): fatal as a coverage engine (conceded; exact refill k <= 3 gives 286).** I add D5: s16a's 72 route
  pairs have no junction at all, so children cannot cheapen a tail across corridors; the 9-craft salvage has no
  mechanism. Rescue: none.
- **P2 (congestion LNS): its budget is below the best design-in price ever measured, and I have now measured its
  per-tail step twice. Fatal (major without part B).**
  - (i) Arithmetic from P2's own base m0 list (1256/1047/1053/947/918/953/1015/991):
    - sum J_i = 11.183, so the slack is 0.360 J = 314 kg for 12 misses;
    - that is **26 kg (0.030 J) per miss, below the 0.032 J floor** of the designed law;
    - P2's "478 kg" uses the 09-24 bar (F5).
  - (ii) The per-tail step is `s15b_prefix` with trading prizes, which is what F1 and part B ran.
    - Every one of the 13 has 1-4 node events < 0.1 AU after day 3328 (`physics/events.json`), so all were reachable.
    - Tails took them only by pushing out equally isolated own targets (§0).
    - With P2's exact residual and prizes, one tail gained +1 flyby, not the +3 that 298 needs (or +1.5 for 292).
  - (iii) At 1.01-1.02 the "diversified" candidates will be near-copies of the own-only re-plan. F1's r6 re-flew
    itself at +1.8 kg.
  - Rescue: a sink, which neither test found. Prediction for P2's own gate: <= 288, i.e. FAIL.
- **P3 (concurrent beam): fatal on the clock.**
  - Its own D3/D4 falsify the premise, and the build takes 12-20 h.
  - Assignment inside the recursion is what an exact master over columns already does. A tuple beam is a weaker
    master with a diluted beam (joint K=2: 35 vs 40).

### PHYSICS
- **P1 (pins and grout): the gate is not tied to the objective, and the clock cannot carry it. Fatal for 09-27.**
  - (i) If the closer pays the forced price (about 120 kg) and the fill covers K at 10.5 kg/fb, then
    10.5 K + 120 (298 - K) <= 3723 needs **K >= 293** (>= 297 at 12 kg/fb).
    - PASS at 262 therefore still implies a fleet of about 7000 kg.
    - Every sequential 8-craft fill so far ended at 228-249.
  - (ii) Pin lanes become skeletons, i.e. waypoints, one of the six instruments that measured law 1. Stage-5
    skeleton fill is the same design and made 241.
  - (iii) T1b's 87/95 is a carrier-MODEL incumbent. P3 itself says the model dv is calibrated 1.5-2x. Its 8 lanes
    miss 36, one of our 13.
  - (iv) A 12-16 h build plus 1.5 d on 6-8 procs lands about 09-26 midday. That is before the closer, the regrid
    (law 8: x1.1-1.6 on planner tours) and a 3-4 h export.
  - Rescue: post-deadline research only.
- **P2 (remaining-window prize): F1 already measured what a prize does at a tail. Major.**
  - A 1.5 premium on targets with live late windows (r9 flies 7 of them after r6's cut, days 3489-5382) pulled
    them in only by trades.
  - A deadline premium is a weaker, time-gated version of that. Its gate ("pins +8 on 4 free routes") never touches
    closure.
  - Rescue: none within prize space.
- **P3 (surrogate): minor.** Enabling only, no fleet by 09-27, law 1 unescaped (conceded).

### REDTEAM
- **P1 (closure-aware MILP): its own runs decide it. Fatal as a route to a fleet; the best closure proof of the
  stage.**
  - `first455.out` at thr 0.15 still gives 285 at 11.007, with the same 13 misses.
  - K = 298 therefore costs >= 11.007 + 13 x 0.08 = **12.05 > 11.543**.
  - The census puts 10 of the 13 at 195-1617 kg on s16a's tracks, i.e. 0.22-1.9 J.
  - The full 1382-column model found no incumbent in 900 s. Its dual bound of 297.99 is the model's LP (REDTEAM's
    own trap 1), not evidence of feasibility.
  - Rescue: none from pool plus insertion.
- **P2 (13-target contract): F1 plus part A is its kill test, run with a stronger prize. Fatal.**
  - The kill line is > 0.08 J mean including re-homing. Measured net per link: +459 to +999 kg = 0.5-1.1 J.
  - Heads are not a rescue. R1 1.4 records that keep-own appends to heads failed and that M-A absorbed 0 of 35.
- **P3 (early lanes): minor.** Every craft leaves Earth at v_inf <= 4 km/s, so year 0-2 crowding is the launch
  constraint (the 28 early flybys include the outbound transfer); REDTEAM concedes it cannot remove a craft.
- **REDTEAM's stop rule (nothing passes its first kill test by 09-25 06:00, then stop): endorsed.** I would tighten
  it for 9 craft (§3).

## 2. What the others' evidence changes in my proposals

- **All three of my proposals are retired.** P1: the master is sound but no column is a sink; P(8 craft by 09-27)
  3% -> **<= 0.5%**. P2: without P1's columns it is the measured forced dissolution (3.6x budget; REDTEAM 4.73 J).
  P3: law 7 unescaped (17's nearest track 0.111 AU, 36's 0.246 AU) and a 6-8 h build.
- **Proposed law 9: re-planning conserves the isolation of what a track flies, and nearly conserves its depth**
  (+0 to +1 flyby). Four machineries give it independently:
  - F1 and part B: 30 -> 30 -> 31, and 24 -> 24;
  - parts A and B: isolated in, isolated out (247 for 53; 219 for 53 and 256);
  - CAMPAIGN X3: a detour to another corridor takes >= 530 d and costs >= 3 own flybys; 0 detours are net-positive;
  - REDTEAM 1.1: removal refunds have a median of 1.7-3.2 kg, so a tail's cost is its track.
- **Concessions.**
  - REDTEAM is right that 8 craft is exactly the 13-target contract. My census (3453 kg against a 961 kg budget) and
    their 4.73 J screen are the same fact.
  - CAMPAIGN D2 is right that the misses of the N<=8 optimum are a subset of r9's targets.
  - REDTEAM's lin calibration (true/lin 0.73 above 0.2 J) lowers my hard core to >= 140 kg per target. That is still
    3x the budget.
- **F4: which reading my proposal depended on.** It depended on fleet-relative isolation, REDTEAM's reading. All three
  readings hold at once:
  - PHYSICS: 8 of the 13 are pins by PHYSICS's own definition (`events.json`, < 0.05 AU; counts 5, 3, 0, 3, 0, 5,
    4, 3, 3, 6, 4, 7, 9).
  - CAMPAIGN: the popularity is stage 14-16 steering. Routes steered onto these targets leave the fleet's corridors
    (D5: only 2.3% of junction-joined pairs are disjoint).
  - Pin-ness is the cause. Isolation from the fleet is the price every 8-craft plan must pay.

## 3. Ranked experiments (at most 3)

**The clock fact that re-ranks everything (computed from s16a raw 13.7399 and the tails' dJ/dkg of 0.00106).** A
9-craft saving has to beat only the bar on the day it is SUBMITTED:

| submit | raw bar | 9-craft saving needed | shown if -120 kg |
|:--|--:|--:|--:|
| 09-24 18:00 | 13.695 | 42 kg | 13.43 |
| 09-25 12:00 | 13.658 | 77 kg | 13.47 |
| 09-26 12:00 | 13.609 | 123 kg | 13.51 (no gain) |
| 09-27 12:00 | 13.560 | 169 kg | — |

A -120 kg 9-craft fleet submitted before 09-25 12:00 is rank 8 on the 20:44 board (est.). Every 8-craft path needs
10-13 isolated targets at 26-46 kg each, 4-35x their measured prices.

**E1. Tail cyclic-exchange master (9-craft cost; the 8-craft coverage master comes free on the same columns).**
- **Mechanism.**
  - Cut all five tails r5-r9 at a common day (~3000).
  - Run 4 trading re-plans per tail over the tails' JOINT post-cut residual, about 80 targets. Use `s15b_prefix`
    with own 1.0 and others 1.0; beam 12/16 x w_fuel 1.0/1.5.
  - Every settled front state becomes a column, and the five originals are columns too.
  - Master A: exact set partitioning, one column per tail, the union covered, at most 3 leftovers via lin arcs
    <= 0.06 J into heads (REDTEAM's calibrated band), minimising sum J_i.
  - Master B: the same columns with r9 removed, N=8, maximum coverage.
  - Regrid and settle the chosen columns (law 8).
- **What it escapes.**
  - Law 3: no tail is built last.
  - The s16_iter plateau: an exact master closes CYCLIC exchanges (A gives X to B, B gives Y to C, C gives Z to A)
    that no single-link screen can see. Part B's columns show the raw material: 5 r9 targets at -5.9 kg, depth 30
    at -51 kg.
  - It does NOT escape law 9. It wins only if a cycle of isolation-swaps nets cheaper. The blockers are known and
    must be carried by some column: 166 (only r6 passes it) and 53 (only r7, at 0.168 AU).
- **Gate (fixed now), after round 1.**
  - Master A PASS: regridded twin tank sum of the 5 tails <= 4706.7 - 80 = 4626.7 kg, with 298 covered.
  - Master A KILL: improvement < 40 kg, or no partition cheaper than the originals.
  - Master B PASS: >= 292 at <= 3723 kg. KILL: <= 288.
  - **Hard stop 09-24 20:00.** The bar rises by about 45 kg per day.
- **Cost:** 8 cores; build 2.5 h (job generator + 5-block partition MILP reusing s16_select); run 3 h (20 jobs x 30
  min) + 1 h regrid/settle. **P(A gate) 5%, P(B gate) 1%, P(submittable) 3%**, all via A.

**E2. Formal 8-craft closure: REDTEAM P1 on the full 1382 columns, warm-started from the 285 incumbent, at thr 0.15.**
- Exact facility location with arcs; escapes nothing, it is a stop proof. **Gate (fixed now):** PASS if an integer
  K >= 290 exists at optimistic <= 11.40; KILL at 1 h with no such incumbent.
- **Cost:** 2 cores, build 0.3 h, run 1 h. **P(gate) 2%, P(submittable) 0.3%.**

**E3. Nothing else.** CAMPAIGN P2, PHYSICS P1 and REDTEAM P2 are each predicted to fail their gates by §0 and §1. Run
none of them unless E1's Master B passes.

## 4. The case for running nothing

- **The strongest argument.** Five machineries measured 8-craft ceilings of 285-286 or forced prices at 3.6x the
  budget, and two more (F1 and part B) show that re-planning conserves isolation. For 9 craft, s16_iter, polish and
  relocate have converged. So E1 is a bet that a CYCLE exists where every measured link loses 118-677 kg.
- **EV in shown J (est.):**
  - E1: 0.03 x 0.06 = 0.002;
  - E2: 0.003 x 0.1 = 0.0003;
  - PHYSICS P1 and CAMPAIGN P3: about 0 on the clock.
- **My position.**
  - I endorse running nothing for every 8-craft program.
  - I do not endorse it for E1. The cores are otherwise idle, the never-submit-worse rule makes the risk zero, and
    E1 is the only experiment whose clock improves with speed: 42 kg is enough by 09-24 18:00.
  - If E1 misses its gate by 09-24 20:00, stop. s16E (13.5116) stays banked.
