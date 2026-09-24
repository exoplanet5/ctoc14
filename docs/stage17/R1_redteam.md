# Stage 17 R1: RED TEAM (data first), 2026-09-23 21:40 CST

All numbers are MEASURED unless marked EST. Data and scripts are in results/s17/redteam/. Compute used: about 30 core-min on 2 cores.
Twin fleet: results/s16/best/fleet (twin sum J_i 11.7227). The honest pool is results/s16/honest (2490 routes).
"lin" means s14_twinbeam.lin_price, run exactly as in s15b_close._price_host (3 approaches per host-target pair, t and t±4 d).

**lin calibration.** Source: closer_calib_lines.txt, 75 real twin insertions from earlier closers.
| lin dJ band | n | settle rate | true/lin (median) |
|:--|--:|--:|--:|
| < 0.03 | 47 | 100% | 1.04 |
| 0.03-0.06 | 10 | 100% | 0.94 |
| 0.06-0.10 | 5 | 60% | 0.87 |
| 0.1-0.2 | 6 | 67% | — |
| 0.2-0.5 | 6 | 50% | 0.73 |

Below 0.06 J, lin can be trusted. Above 0.15 J it only ranks candidates.

## §1 Measured facts

### 1.1 The 9 -> 8 lower bound by insertion
Sources: q1_dissolve.json, dissolve_*.jsonl, nbr.json, removals_all.json.

Every target of route V was priced into every other route: approaches within 0.3 AU, cheapest host by lin. "Screen" is the sum of the per-target minima. "Convex" is a greedy assignment that lets host tanks grow.

| V | n | J_i saved | screen ΣdJ | convex | ≤0.03 | 0.03-0.06 | 0.06-0.15 | >0.15 | no other track ≤0.06 AU |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| r9 | 21 | 1.263 | **4.73** | 5.19 | 5 | 3 | 4 | 9 | 12 |
| r8 | 24 | 1.288 | **7.44** | 8.90 | 0 | 2 | 9 | 13 | 18 |
| r3 | 41 | 1.352 | 5.10 | 6.39 | 8 | 11 | 11 | 11 | 17 |
| r7 | 30 | 1.387 | 6.19 | 7.64 | 3 | 3 | 8 | 15 | 17 |
| others | | | 6.6-9.1 | | | | | | |

- **Every route costs 3.7-7x its own J_i to dissolve by insertion.** r9 is the cheapest to dissolve.
- **r9 splits cleanly into two groups:**
  - 8 near-track targets: lin ≤ 0.06, Σ 0.22 J.
  - 13 isolated targets: 17, 36, 75, 108, 118, 180, 212, 216, 219, 220, 243, 247, 288. The nearest other track passes 0.024-0.25 AU away (median 0.12 AU). Their lin Σ is about 4.5 J, and 9 of them are above 0.15 J.
- **Where r9's price comes from, against the 09-27 budget.** Dissolving r9 must cost ≤ 1.263 - 0.18 = **1.08 J** to beat 13.56.

  | pricing model | cost of dissolving r9 | vs 1.08 J budget |
  |:--|--:|:--|
  | insertion screen | 4.73 J | 3.7x over |
  | forced-premium law: 0.22 + 13 × 0.115-0.148 | 1.72-2.14 J | net +0.46..+0.88 J (loses) |
  | designed-around law: 0.22 + 13 × 0.032-0.042 | 0.64-0.77 J | fits; raw EST 13.12-13.25 |

  So the 4.73 J screen breaks down as:
  - about 0.7 J intrinsic;
  - about 1.2 J of 3x forced premium;
  - about 2.8 J of isolation surplus. These are targets more than 0.1 AU from every track, where lin exceeds the forced law.
- **The r9 break-even:** each isolated target at ≤ 0.066 J. That is 1.8x the designed-around price and 0.5x the forced price. The only designed-in isolated target on record (64 into h06, s15b) cost about 0.12 J with its re-homings.
- **r8 is worse.**
  - Only 2 of its 24 targets are cheap; 18 have no other track within 0.06 AU.
  - Forced law: about 2.9 J. Designed-around: 0.88 J, against a 1.11 J budget.
- **Removal refunds are tiny.**
  - Every single-target removal from every route settled.
  - The median refund is 1.7-3.2 kg per target, and the maximum is 34 kg.
  - Summed, they return only 39-54% of a tail's fuel term (62-77% for deep routes).
  - A tail's cost is its TRACK, not a few forced targets. That is why moving targets between the existing 9 never pays (s16_iter converged).

### 1.2 The "unwanted" targets
Sources: targets.csv, anatomy.json, auc.txt.

**Tail and deep targets are NOT separate orbit families.**
- Cross-validated logistic AUC for tail vs deep membership on (a, e, i, q, node distance, Earth MOID, flyby v_rel, node longitude, #node passes) is **0.61 ± 0.01**.
- The best single features are a, P and #node passes, all at AUC 0.63. Tails have median P 2.9 yr against 2.4, i.e. 5.2 against 6.3 near-Earth node passages in the window.
- i, node distance, MOID, flyby r and z are all at AUC ≤ 0.58.
- Only pool preference separates the two (AUC 0.82): deep pool routes (n ≥ 37, 439 of them) contain a tail target 30 times against 74.5 for a deep target.
- Only 3 targets (86, 169, 204) are almost never in deep pool routes.
- **Every tail target has been flown by some 43-48-flyby pool route** (median max depth 47). "Unwanted" is an assignment artefact.

**The real subfamily is r9's 13 isolated targets.** Compared with the rest:
| property | r9's 13 isolated | rest |
|:--|--:|--:|
| P | 3.7 yr | 2.6 yr |
| near-Earth node passes in window | 4.1 | 5.8 |
| i | 19.6° | 13.4° |
| flyby v_rel | 21.9 km/s | 16.4 km/s |
| node distance | 0.040 AU | 0.044 AU |

They are rare, fast events, at the same node distance.

**The tracks are already disjoint and saturated** (nbr.json). This counts targets with no non-owner route passing within a given distance:
| within | targets |
|:--|--:|
| 0.035 AU | 241 of 298 |
| 0.06 AU | 163 |
| 0.10 AU | 88 |

Each route passes only 5-10 foreign targets within 0.035 AU.

### 1.3 The 9-craft ceiling (Q3)
Sources: timeline.json, milp4_*.txt, anatomy.json. The tails carry 136 flybys on 1707 kg (12.55 kg/fb).

| tails at | twin ΣJ_i | raw EST (+2, +0.017 conversion) |
|:--|--:|--:|
| 8.3 kg/fb (deep efficiency) | 11.14 | **13.16** |
| 9.0 kg/fb | 11.23 | 13.25 |
| 10.0 kg/fb | 11.37 | 13.38 |
| pool-average law, fuel = 62 + 8.97 n | | 13.555 (only break-even on 09-27) |

**Not reachable by re-planning the tails over their own union:**
- A MILP over the pool shows the only 5-route cover of the 136 is the current tails. Four pool routes cover at most 116 of 136.
- The best single tail-set routes in the pool fly 33-35 targets at 1037-1139 kg, i.e. 11-15 kg/fb.
- The tails are density-limited:
  - median gap between flybys: 176-194 d (deep routes: 111-132 d);
  - median leg dv: 0.47-0.59 km/s (deep routes: 0.30-0.40);
  - every route spends 15.9-20.6 km/s whatever its depth.
- The linchpin law (depth collapse under pool re-planning) points the same way.
- **Verdict:** a 9-craft -200 kg needs NEW tail tracks, not re-planning. There is no evidence it exists.

### 1.4 Reading the leaders (Q4)
Source: anatomy.json.

**Total Δv to cover 298** (ve × Σ ln(m0/600)):
| fleet | total Δv | per flyby |
|:--|--:|--:|
| THU | 165.7 km/s | 0.556 km/s |
| us | 158.2 km/s | 0.531 km/s |
| HIT | 144.9 km/s | 0.486 km/s |

- **Our own 158 km/s spread evenly over 8 craft gives raw EST 12.90; over 7 craft, 12.11.** The leaders' edge is not propellant efficiency.
- What binds is the flyby RATE per craft. 8 × 37.25 needs a mean gap ≤ 142 d on EVERY craft. Our deep routes meet that; our tails run at 176-194 d.

**The free-optimum envelope** is the lightest pool route at each depth:
| depth | fuel | kg/fb |
|:--|--:|--:|
| 30-38 | 161-225 kg | 5.4-6.2 |
| 39 | 322 kg | |
| 39-48 | | 8.2-8.9 |

The jump at 38/39 is a generator artefact: the old m0 cap.

**Everyone pays about 1.1 t of disjointness tax** above that envelope:
| fleet | tax above envelope | how it is spread |
|:--|--:|:--|
| us | 1094 kg | deep routes 183 kg (-1..88 each); tails 911 kg (140-256 each) |
| HIT | 8 × 139 = 1112 kg | even |
| THU | 7 × 155 = 1085 kg | even |

- The leaders spread the tax evenly across balanced routes. Ours is concentrated in 5 starved tails. EST ±100 kg: the envelope above n = 38 comes from other generators.
- THU's routes are "our deepest free optima + about 1 forced target each". HIT's routes are "our r4 + about 40 kg".
- Both leaders must own 7-8 mutually disjoint deep TRACKS. Our pool's largest disjoint set of good routes is 4 (law 2).

**Phase lanes** (phase.json). Measured as Earth-relative heliocentric longitude:
- Our deep routes are lanes:
  - r1 stays within -40..+33° all mission (span 116°);
  - r2 within +20..+127°;
  - r4 within -28..-145°.
- The tails wander (spans 207-866°). r9 laps the Sun 2.4 times relative to Earth to collect the rare events.
- **Early crowding:** in year 1 all 8 early craft sit within -86..+54°.
  - The fleet flew 28 flybys in years 0-2, against 42-45 in every later 2-year bin.
  - The tails flew 0-3 each in years 0-2.

## §2 My proposals

### P1. Closure-aware N = 8 selection: facility location with lin-priced insertion arcs
- **Precedent:** GTOC9 (JPL/NUDT) campaign-level assignment [recall, unverified]; in OR terms, capacitated facility location over columns.
- **Mechanism:**
  - Columns: honest pool + s16 routes with n ≥ 25 (about 700).
  - Arcs: every (target, column) approach within 0.10 AU, priced by lin, kept when lin ≤ 0.08 J.
  - MILP: pick 8 columns; assign every target to a picked column that contains it, or via an arc of a picked column; cap each column's arcs at ≤ 150 kg added.
  - Objective: Σ cost(tank) + Σ arc dJ.
- **Law it escapes:** "N ≤ 8 pool ceiling 282". That ceiling counts only targets a column already flies. But every route passes 14-31 foreign targets within 0.06 AU. My greedy version (r1-r4 + 17 arcs ≤ 0.06 + 4 pool tails) already reaches 285 at twin 11.06 (close8.txt).
  - This is NOT the dead stage-10 facility location. That one priced closest approach (R² 0.01). lin is calibrated: ratio 1.0 and 100% settle at ≤ 0.06 J.
- **Falsification:** the MILP optimum is OPTIMISTIC (additive, no interaction between arcs).
  - If min twin ΣJ_i at 298 with 8 columns is > **11.54** (the 09-27 bar), the whole pool + insertion route to 8 craft is dead, and so is every LNS variant that does not create new tracks.
  - Go on only if ≤ **11.40**.
- **Cost:** lin pricing is about 4 s per column, so about 50 core-min, plus a MILP of minutes. **About 1.5 h on 4 cores.**
- **Payoff:** P(bound ≤ 11.40) EST 20%; P(realised under the bar with closer insertions) EST 6-8%. Raw 13.3-13.5 on 09-26.
- **Its main value is information.** It gates P2, LNS and PHYSICS for 1.5 h.

### P2. "The 13-target contract": gated design-in of r9's isolated set
- **Precedent:** s15b prefix-seeded trading re-plan (64 into h06); LNS destroy-repair on a single target.
- **Mechanism:** the pool-level 8-craft skeleton is 285 covered at twin 11.06:
  - r1-r4,
  - plus the 17 cheap lin arcs (q_tailabsorb.txt),
  - plus pool tails b33a7d00, 007b9f03, 1f61b93d and r8.

  The **whole** 8-craft question is then whether exactly 13 targets can be designed into these 8 tracks. They are rare-event, long-period, inclined and fast. The budgets:
  - ≤ 0.48 J total (**0.037 J each**) to beat the 09-27 bar;
  - ≤ 0.66 J (0.051 each) just to beat the current 9-craft twin.

  Order them by nearest track: 180 (0.024 AU), 118 (0.066), 75 (0.114), 212, 288, ...
- **Law it escapes:** the forced-premium law, but only if the host's suffix is re-planned (keep-own 2.0 + trading prize) rather than inserted into.
- **Kill test:**
  - First 3 targets (180, 118, 75), each by s15b_prefix trading re-plan from the flyby before its epoch, 3 h on 4 cores.
  - Kill if the mean true dJ including re-homing is > **0.08 J**, or if 2 of 3 fail.
- **Payoff:** P EST 5%. It succeeds only if designed-around pricing holds for isolated targets, and no measurement so far says it does.
- **Effort:** tools exist (s15b_prefix, s15b_close, s15b_regrid). About 12-18 h in total if it passes.

### P3. Early lane separation (weak, flagged for PHYSICS/CAMPAIGN)
- Years 0-2 fly 28 instead of about 43 flybys. That is about 15 flybys lost to 8 craft crowding the same ±60° of phase.
- A launch design that sends the craft to separated phase lanes (v_inf 4 km/s ≈ ±100°/yr drift) could recover about 10-15 flybys. That is 0.3-0.4 of a craft.
- **Kill:** count the distinct targets with a cheap early event (≤ 1 km/s from an Earth departure, v_inf ≤ 4, first 700 d) in lanes other than the current deep ones. Kill if < 20.
- Alone it does not remove a craft, so it is only a component.

## §3 Red team: predicted proposals, ranked by EV against the clock

**Frame.**
- A result pays only if raw < 13.56 at 09-27 12:00, or < 13.609 at 09-26 12:00.
- The last useful twin fleet must exist by about **09-27 06:00**, because export takes 3-4 h.
- The current shown score is 13.5116. s16a at 09-23 22:00 would show about 13.515, so it is already stale.
- Every 8-craft path must pass the P2 arithmetic: the pool skeleton leaves 13 rare-event targets and about 0.48 J.
- EV = P(valid fleet under the bar by 09-27 06:00) × (13.51 - shown).

| # | proposal (lens) | P EST | gain if success | EV (J) | cheapest KILL test, fixed threshold |
|:--|:--|--:|--:|--:|:--|
| 1 | **Closure-aware N = 8 MILP with lin arcs** (mine; CAMPAIGN may pose it as "joint assignment") | 7% | 0.1-0.25 | 0.012 | Optimistic MILP bound > 11.54 twin at 298 means dead; go only if ≤ 11.40. 1.5 h. |
| 2 | **Suffix re-partition / dissolve r9 by designed-around re-plans** (LNS) | 5-8% | 0.1-0.3 | 0.010 | r9 has 5 isolated targets before day 2500 (17@1309, 36@1496, 75@1679, 118@1890, 212@2456), so t* must be ≤ about 1300 d. That is 88% of every host, which is the depth-collapse regime. Test one host (r1) re-planned from t* = 1300 with keep-own 2.0 + r9 prizes. Kill if it drops any own target or absorbs < 2 r9 targets at ≤ +80 kg. Then the P2 3-target gate (mean > 0.08 J means kill). |
| 3 | **Concurrent K = 8 time-synchronous beam / joint k-route DP** (CAMPAIGN) | 3-5% | 0.2-0.5 | 0.010 | Truncate the horizon: run to day 1461 on the full set. Kill if < 72 distinct targets flown (our 9 craft had 70), or if mean fuel spent is > 160 kg per craft, or if any target's last event passes unflown with no craft in range. Priors: jointsearch K = 12 reached 227-280 at 18 kg/fb; joint beam 10 craft reached 200; joint K = 2 lost to sequential (35 vs 40). The new code must beat these before day 1461 is worth simulating further. |
| 4 | **Leftover-specialist craft designed first** (PHYSICS) | 3% | 0.1-0.3 | 0.005 | The tail set is not a family (AUC 0.61). Only the isolated ones are: 13 of r9 plus 18 of r8 = 30. Build the specialist with prizes on those 30. Kill if it catches < 22 of them at ≤ 1100 kg; history is a 26-flyby route with 14 hard. Then fix it in s16_select N ≤ 8: kill if the frontier rises by < 6 above 282. |
| 5 | **9-craft tail redesign for -200 kg** (any lens) | 8% | 0.01-0.05 | 0.002 | Needs tails ≤ 1500 kg; the pool's 5-cover is the current tails. Kill if no new single tail track reaches ≤ 10 kg/fb at depth ≥ 27 over the tail residual within 4 h. Hard stop 09-25 12:00. |
| 6 | **Orbit-state / phase-graph path cover** (PHYSICS, GTOC4-style) | 2% | 0.2-0.6 | 0.006 | Kill if the model's own K = 8 optimum, before settling, covers < 295 at ≤ 20 km/s per path. Also kill if its leg price has Spearman < 0.6 against twin leg dv on our 298 flown legs. Carrier tiling gave 244-269; it must beat that on paper first. |
| 7 | Transfer/fragment DB + GA (GTOC11-style); ACO partition; more pool selection | <1% | | 0 | Already killed by laws 6 and 2 and the ACO forcing result. Reject without a new mechanism. |

**Ranking logic.**
1. Run #1 first. It is the cheapest decisive test (1.5 h), and its bound gates #2 and #4. If the optimistic bound is above 11.54, the 8-craft program cannot come from existing tracks plus insertion. Only #3 or #6 (new tracks) remain, and neither is likely to finish before 09-27 06:00.
2. #2 is the only path with existing tools end to end. But it needs designed-around pricing on isolated rare-event targets, which is unmeasured, and the one data point (0.12 J) is 3x too high.
3. #3 has the largest upside but needs new code. Its priors are 200-280 covered. Allow at most 6 h to reach the day-1461 gate.
4. Stop rule for the whole stage: if no route passes its first kill test by **09-25 06:00**, stop. s16E stays banked: shown 13.5116; do not resubmit.

**Traps I expect in the other R1s:**
1. Citing an LP or MILP bound of a MODEL as a bound on the problem (the stage-10 lesson).
2. Unregridded planner tanks: they are 1.1-1.6x over-credited (law 8).
3. "Tails are a physical family": measured AUC 0.61.
4. "Kg/fb of the tails can reach deep efficiency": removal refunds show the tails' cost is the track itself.
5. Counting coverage without pricing the leftovers: the 282 and 285 skeletons both leave rare-event targets at lin 0.1-1.3 J each.

## §4 P1 frontier (appended 2026-09-23 23:55, round 2)

**The full run failed.**
- `p1/full.json` crashed at 23:33 after 1 of 20 jobs, because of a path bug in the out argument.
- Its restart, `p1/full2_*.json` (4 cores, 900 s per job), shows that the full 1382-column, 23.8k-27.1k-arc MILP is
  INTRACTABLE in 15 min:
  - N=8 max coverage: incumbents of 26 and 81 covered, against an LP bound of 297.99;
  - N=9 at K=298: no incumbent.
- The frontier below uses the tractable F2 set instead (`p1/first455.json`, `p1/partial_*.json`, `r2/forced_455.json`):
  the 9 s16a routes plus the 446 honest pool routes with n ≥ 36, cap 150 kg of arcs per host. Every run reached the
  HiGHS optimum.

| arcs | N | K | result |
|:--|--:|:--|:--|
| lin ≤ 0.08 J (3671) | 8 | max | **285**, ΣJ_i 11.007 additive / 11.008 convex. Misses 17 36 75 108 118 180 216 219 220 243 247 258 288 |
| lin ≤ 0.15 J (4467) | 8 | max | **285**, the same solution. Loosening to untrusted arcs adds nothing. |
| ≤ 0.08 | 8 | 286, 290, 294, 298 | infeasible |
| ≤ 0.08 | 8 | max with one miss forced in | 278-284 for every one of the 13 (R2_redteam §0). Isolation is conserved. |
| ≤ 0.08 and ≤ 0.15 | 9 | 298 | **11.7227 = s16a exactly**. No pool column or arc beats the current 9 routes. |

**Verdict against the pre-registered kill ("> 11.54 at 298 means dead"):** at 8 craft, 298 is infeasible, so P1 is
**DEAD**, and so is every "existing tracks + insertion" route to 8 craft.
- Agreeing results: CAMPAIGN 286 (3585 sets, no arcs) and stitched children 286.
- At 9 craft, s16a is optimal in the arc-extended deep column space, and in CAMPAIGN's 3347-parent space.
- Not closed: the full-pool (n 20-35 columns) + arc version. It is intractable as posed. The prior that it helps is
  low, because s16_iter already explored fleet-internal arcs.
