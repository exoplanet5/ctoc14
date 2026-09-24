# Stage 17 R1: LNS lens (route minimisation, set partitioning, ruin and recreate)

Agent LNS, 2026-09-23 22:10 CST. MEASURED numbers name their file (all under `results/s17/lns/`); everything else is
marked est. Literature I am not sure of is marked [recall, unverified].

## §1 Diagnosis: what binds

### 1.1 The 8-craft budget is tighter than the brief's table (computed from `results/s16/best/fleet/fleet.json`)
With r9 gone, the survivors r1..r8 sit at twin sum J_i 10.4596, and conversion adds +0.0172 J (s16a exact 11.7399
minus twin 11.7227). The 21 targets of r9 must then be re-homed for at most:

| submit | raw bar | kg added to r1..r8 (spread evenly) | per r9 target | balanced 8-craft total fuel |
|:--|--:|--:|--:|--:|
| 09-24 12:00 | 13.707 | 1084 | 51.6 | 3846 |
| 09-25 12:00 | 13.658 | 1043 | 49.7 | 3805 |
| 09-26 12:00 | 13.609 | 1002 | 47.7 | 3764 |
| 09-27 12:00 | 13.560 | **961** | **45.8** | **3723** |

- The brief's "8 craft beats the 09-27 bar at up to 3858 kg" is really the 09-24 bar. At 09-27 the limit is 3723 kg,
  135 kg less.
- The brief's "J < 13 needs <= 3391 kg" is also off: that fuel gives raw 13.16. J < 13 needs 8 x 993 kg = 3145 kg.
- The 9-craft number (-171 kg at 09-27) agrees with the brief.

### 1.2 The forced price of dissolving r9 (MEASURED)
Source: `census_r9.json`, made by `census.py`. For every r9 target I took the approach minima within 0.25 AU of the 8
other routes and priced each one with `s14_twinbeam.lin_price` (the `s15b_close._price_host` screen, 3 epoch
offsets). Best host per target:

| class | n | targets (best host, kg) | sum |
|:--|--:|:--|--:|
| cheap, <= 40 kg | 7 | 265 (r2, 13), 176 (r3, 18), 142 (r2, 18), 208 (r6, 22), 130 (r1, 26), 6 (r6, 29), 184 (r1, 30) | 156 |
| mid, 40-150 kg | 4 | 243 (r6, 53), 258 (r5, 91), 180 (r7, 96), 75 (r8, 102) | 342 |
| **hard core**, > 150 kg or no approach | **10** | 17, 36, 108, 118, 212, 216, 219, 220, 247, 288 (best 195-880 kg; 108 has no approach within 0.25 AU) | >= 2955 |

- Sum, capping each target at 400 kg: **3453 kg, 3.6x the 961 kg budget.** This is law 1 again. It also matches the
  earlier dissolve_lns result (globalopt_campaign 8.2: 4x over budget).
- The other tails are worse victims (`census_summary.json`, same cap): r8 5353 kg, r7 5459, r6 5314, r5 6038. Each has
  13-16 hard-core targets. So r9 is the right victim.
- What the budget demands: the cheap and mid targets use about 500 kg, which leaves **46 kg for each of the 10
  hard-core targets.** Their forced prices are 195-880 kg. So 8 craft means **designing the hard core into other
  tracks at the 1x price while those tracks keep their own targets.**
- In general terms: if a fraction f goes in at the designed price (~32 kg) and the rest at the forced price (~120 kg),
  then f >= 0.84.

### 1.3 Time does not bind; track density does (MEASURED, `anatomy.json`)
- Cadence of the heads r1-r4: 128-144 d per flyby. Tails r5-r7: 174-181 d. r8: 211 d. r9: 211 d (it launches on day 944).
- At head cadence, 298 flybys fill only 7.0-7.9 craft lifetimes, so the fleet has spare craft-days. They sit on the
  tails, whose tracks are 0.11-0.25 AU from the hard core (law 7).
- Each hard-core target is near no track at all. That is why the problem is not VRPTW-like slack-shuffling.

### 1.4 The three failures, read through OR
- **s16_iter ran out of moves at 9 craft.** Its acceptance rule is fleet J descent over single moves. Emptying a
  route means crossing a plateau: by the census, every single move out of r9 costs more than it saves, so no J-descent
  ever starts. VRPTW route minimisation handles this with a separate phase and its own objective (1.5). s16_iter had
  no such phase, and its columns were single-target, i.e. forced by construction.
- **Joint K=2 over a 65-pool did worse than sequential.** A two-craft beam splits one width budget across a product
  space, so each craft's effective beam collapses. The failure indicts joint BEAMS, not joint ASSIGNMENT. Assignment
  can be done exactly by a MILP over single-craft columns (s16_select already does this).
- **The insertion depth law.** Point insertion into a pinned deep route pays the "bump" dv ~ 4 dr / T over short
  gaps (globalopt_campaign 8.1). Only a re-planned SUFFIX removes the pins. Evidence so far:
  - s15b prefix jobs traded +1 designed target for -1 or -2 own at +7 to +34 kg (`results/s15b/log.txt`, 02:50-07:18);
  - keep-own "appends" to heads all failed;
  - stage 13 M-A: pinned segment re-plans of 8 hosts absorbed 0 of 35.
  So heads have no 1x absorption capacity. Tails had not been tested before §3.

### 1.5 What OR says about fleet-size reduction
- Nagata & Bräysy (2009, route minimisation for VRPTW): ejection pool, penalty counters, and a "squeeze" through
  temporary infeasibility.
- Bent & Van Hentenryck (2004, two-stage hybrid): stage 1 maximises the sum of squared route sizes so the smallest
  route drains.
- Both work because slack exists somewhere and ejections are cheap. In our problem slack exists (time on the tails),
  but no ejection is cheap for the hard core. [recall: both papers fairly sure; details unverified]

## §2 Proposals

### P1 (top). Suffix-column set partitioning: heuristic branch-and-price over prefix-anchored suffix columns
(a) **Precedent.**
- Set-partitioning recombination of a route pool built by local search: Subramanian, Uchoa & Ochi 2013 ILS-SP;
  Rochat & Taillard 1995 adaptive memory [recall].
- Branch-and-price with heuristic pricers tried before exact ones (Desaulniers, Desrosiers & Solomon 2005) [recall].
- GTOC9 (JPL): a database of independently optimised missions, assembled combinatorially, then improved by
  replacing a few missions at a time [recall, unverified].

(b) **Mechanism.**
- A column of craft i is its settled twin prefix up to cut k plus a suffix re-planned by `s15b_prefix`.
  - Pool: own suffix, plus r9's targets, plus the suffix targets of craft whose tracks pass within 0.1 AU.
  - Prizes: own 1.0, r9 1.5, hard core 2.0. This is TRADING, not keep-own; LP duals were degenerate in stage 7, so
    the prizes are fixed.
  - Cuts: one flyby before each hard-core target's approach to the host (from the census), plus every 4th flyby in
    the second half.
- Every settled front state becomes a column; the tool already writes them.
- Each column is a whole route, so the master IS `s16_select` (N <= 8, K = 298) over the honest pool plus these columns.
- The MILP resolves ejection chains exactly: A takes X and drops Y, and B's column takes Y. Neither single moves nor
  sequential dissolution can do this. Run about 3 rounds; each round cuts around the new incumbent.

(c) **Laws.**
- Escapes law 1: every placement lives in a re-planned suffix, so it is designed, not forced.
- Escapes law 3: the assignment is joint (MILP), while each beam is single-craft, the regime where s15b trades worked.
- Avoids the linchpin depth collapse: the prefix is a settled twin, whereas stage 13 re-planned from launch.
- Does NOT escape laws 2 and 7. If no track can bend through the hard core without dropping own targets, the columns
  trade 1-for-1 and the master cannot close.

(d) **Falsification.**
- Stage 1 is §3 (keep-own absorption into two tails). Thresholds are in `THRESHOLDS.txt`, written before the results.
- Stage 2, if stage 1 is not FAIL, takes 60 min on 2 cores:
  - 6 trading jobs, beam 12, cut one flyby before the approach;
  - hosts nearest the hard core: r2 for 118, 212 and 219; r5 for 216 and 220; r8 for 288; r6 for 247.
- PASS if >= 4 distinct hard-core targets appear in columns whose net cost is <= 46 kg per hard-core target. Net cost
  = tank delta + the census lin price of every own target the column drops.

(e) **Payoff.**
- 8 craft at 298, raw about 13.45-13.56 at 09-27 (est.).
- By-product at 9 craft: multi-target suffix columns for the existing master, about -0.02 to -0.05 J (est.).
- Probabilities: see §3.4.

(f) **Build: about 3 h.** A job generator plus fixed prizes; no new master. Reuses `s15b_prefix`, `s15b_twintail`,
`s15_twintail`, `s16_select`, `s15b_regrid` and the `s16_iter` cache. Compute is about 25 core-h per round
(100 jobs x 15 min).

### P2. Ramped emptying objective inside the s16_iter MILP loop (dissolve r9 gradually)
(a) **Precedent.** The Bent & Van Hentenryck stage-1 objective, and the Nagata & Bräysy ejection pool with penalty
counters (1.5) [recall].

(b) **Mechanism.**
- Objective: sum J_i + lambda x n(r9), with lambda ramped 0 -> 0.02 -> 0.04 -> 0.06 -> 0.08 J per target over s16_iter
  iterations.
- Columns: s16_iter's single moves, P1's suffix columns, and r9 truncations (r9 minus its first or last j targets,
  settled with `w_remove`).
- Once n(r9) <= 4, the rest go in by point insertion. The final fleet is accepted by true J only.

(c) **Laws.**
- It crosses the plateau of 1.4, one target at a time.
- With only single-move columns it IS the measured forced dissolution (3.6x over budget). It escapes nothing without
  P1's columns.

(d) **Test (<= 45 min, 2 cores).**
- Settle 3 r9 truncations (drop the 7 cheap targets; drop the 7 cheap + 4 mid; drop the late 6).
- Point-insert the dropped targets into their census best hosts (`w_insert`, 11 trials).
- PASS if the first 7 leave for <= +0.06 J net and the first 11 for <= +0.20 J net. That is the pace that leaves
  >= 46 kg per hard-core target.
- FAIL otherwise. Prediction from the census: net +0.10 and +0.40.

(e) **Payoff.** It is a driver for P1. Alone: P(8 craft) < 1%.

(f) **Build: 2 h.** Objective term plus truncation columns; reuses `s16_iter`, `s16_select`, `run_ialns.w_remove`
and `w_insert`.

### P3. Prefix re-partition with the launch free (end-anchored twin beam) for r9's early hard core
(a) **Precedent.** Bidirectional labelling in ESPPRC pricing (Righini & Salani 2006) [recall]. Backward search from
late fixed events is common GTOC practice [recall, unverified].

(b) **Mechanism.**
- Keep craft i's settled suffix from flyby k on.
- Re-plan its prefix from a free launch (date free, v_inf <= 4 km/s) over own prefix targets plus r9's early hard core
  (17@1309, 36@1496, 118@1890, 212@2456).
- Children are inserted before the kept suffix and priced by `lin_price` on the whole twin; the suffix pins the end
  state.
- This needs a new `s17_prefixfree.py`: the current beam only appends after t_end.

(c) **Laws.**
- For early targets it escapes law 1: the launch is a route's freest variable, and early insertions into sparse routes
  cost 4-20 kg (stage 7 carrier_grow).
- It does NOT escape law 7. In the census, 17's best approach is 0.111 AU (r1, day 4955) and 36's is 0.246 AU (r4).
  No early track comes near them either.

(d) **Test (60 min, 2 cores, after the build).** Re-plan r2's prefix up to 43@2269 with 118 and 212 prized 1.5. PASS
if a settled prefix keeps all own targets, adds >= 1 of the two, and costs <= 46 kg per added target.

(e) **Payoff.** Only as P1's complement for the first half of the mission. P(8 craft) about 1%.

(f) **Build: 6-8 h** (a new beam direction). Too risky for a 3.3-day window unless P1 passes.

## §3 Measurement: can a tail absorb r9's targets at the designed (1x) price?

Pre-registered in `THRESHOLDS.txt` (21:31, jobs started 21:26, no result read before). Two keep-own suffix re-plans
with `tools/s15b_prefix.py`: own prize 2.0, r9's 21 targets prize 1.5, beam 16 / tries 32, 1 core, 30-33 min each.
Outputs: `pre/L1r6k17S/`, `pre/L1r8k13S/` (logs `pre/*.out`). Final columns scored by the orchestrator (fact F1) and
re-scored in round 2 by `r2/chain_screen.py` -> `r2/chain_screen.json`.

| host, cut | original | own-only best | r9-carrying states (all TRADE own targets) |
|:--|:--|:--|:--|
| r6 after 17 fb (day 3328) | 30 fb, 901.5 kg | all 13 own re-flown, 30 fb, 903.4 kg (+1.8) | +{184, 247} -{53, 170} +49.5 kg; +{142} -{53, 170} +9.6 kg |
| r8 after 13 fb (day 3112) | 24 fb, 927.2 kg | never all own (best 9/11) | +{118, 288} -{138, 172} +306 kg; +{142, 212} -{138, 150, 223} -70 kg; +{212} -3 own -83/-85 kg |

**Verdict: FAIL** (no all-own state carries any r9 target). Two further facts:
- **Depth is conserved.** Offered 21 extra targets, neither host ever exceeded its original flyby count (r6 30, r8
  24; r8's re-plan hit a settle wall at 24). The re-plan changes WHICH targets a tail flies, not HOW MANY.
- **Isolation is conserved** (round-2 addendum, `r2/chain_screen.json`). Every own target a trade drops is itself
  isolated. Its forced re-home price (census_r6 / census_r8, best host other than itself and r9):
  - 53: 1617 kg (r7 at 0.168 AU);
  - 138: no approach within 0.25 AU;
  - 150: 269 kg; 172: 293 kg; 173: 412 kg; 170: 55 kg; 223: 65 kg.

  Net per link = tank delta + re-home of the dropped. Against the forced alternative for the added targets, the
  saving is -118 kg for the best link (r6: 247 in, 53 out) and -441 to -740 kg for the others. No link has a
  positive saving, and no host offered a sink (absorbing without dropping). An exact MILP over such columns can only
  chain links. A chain with no sink and no positive link cannot close.

### §3.4 Probabilities after the test
- P1: P(8-craft submittable by 09-27 06:00) falls from ~3% to **<= 1%**. It survives only if the round-2 open-pool
  test (R2_lns.md) shows that tail depth is pool-bound.
- P2: ~0.
- P3: ~0 on the clock.
