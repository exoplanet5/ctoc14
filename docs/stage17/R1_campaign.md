# Stage 17 R1 — CAMPAIGN lens: build the fleet jointly (2026-09-23 22:10 CST)

Question: how do we build all craft together, so that disjointness and coverage are designed in rather than
selected afterwards? This lens draws on GTOC9 (JPL, NUDT), GTOC11 (Tsinghua), GTOC12 (Tsinghua, ACT) and GTOC7.

**Summary.**
- My top proposal going in was **state-proximal recombination**: stitch settled twin segments at orbit-state
  junctions (GTOC11-style fragments), then pick 8 chains jointly. I built it and ran it (§3).
- **The children are honest.** 5 of 5 settle in 3-8 s, within -15..+1 kg of the estimate.
- **They do not lift the 8-craft coverage ceiling.** It stays at 286, which is the pre-registered FAIL.
- The runs also measured why, and that changes what "joint design" has to mean (§1).
- What I recommend building next is P2, joint prefix re-planning with exact tail selection. Its probability of
  giving a submittable 8-craft fleet by 09-27 06:00 is about 5%.

All numbers marked MEASURED come from files under `results/s17/campaign/` (scripts are in the same directory).

## §1 Diagnosis: what binds

**D0. Fuel does not bind at N=8; coverage does.** (From the brief.) 8 craft at 298 beat the 09-27 bar at up to
3858 kg of fuel (12.9 kg per flyby), against 8 kg per flyby on our deep routes. So the problem is packing 298 targets
into 8 tracks, not paying for them.

**D1. The pool is narrow (MEASURED, `xover_scan.py`, `log.txt`).**
- The honest column cache `results/s16/colcache.pkl` holds 6659 ok columns.
- They reduce to 3585 distinct target sets.
- Clustering at Jaccard < 0.8 leaves only **1220 families**.

**D2. Control (MEASURED, `milp.json`).**
- N <= 8 max coverage over all 3585 sets is **286**, not the 282 of stage 16. The stage-16 move columns
  (iter/swap/polish/colgen) add 4.
- LP bound 290.35 (`lns_N8.json`).
- The 8 picks have depths 47/46/44/38/31/31/29/24 at 1256/1047/1053/947/918/953/1015/991 kg.
- Misses: 17, 36, 108, 118, 180, 216, 219, 220, 243, 247, 258, 288.
- The deep heads are ~27 flybys above 4 x 37, while the four tails are ~33 below. This is law 3 again.
- **All 12 misses are targets of s16a's r9**, the day-944 specialist with 21 flybys. It flies them one every
  ~320 d, at days 1309/1416/1495/1889/2693/3036/3245/3488/4101/4302/4602/5189 (checked on
  results/s16/best/fleet).
- So the N=8 problem is exactly: absorb r9's 12 targets into 8 tracks. Stage 16's move columns already absorbed
  r9's other 9.

**D3. The chronic misses are the MOST contested targets, not rare ones (MEASURED, `scarcity_check.json`,
`miss_windows.json`).**
- Each of the 12 misses is flown by 89-232 of the 1220 families.
- In popularity they sit at the 64th-99th percentile of all targets (median 83rd).
- 190 targets are flown by fewer families than the least-flown miss.
- Every miss has 3-11 distinct encounter windows in the pool, and each has at least one window after day 4000.
- Caveat: stages 14-16 deliberately steered routes onto hard targets, which may inflate these counts. Even so, the
  misses are not isolated, so this is a packing loss (law 2), not a reach loss.
- Consequence: any value function that prices scarcity or "coverage potential" (how many craft could still take
  a target) would rate these 12 as safe. That is the mechanism behind the dead scarcity and hardness prizes.

**D4. Conflicts are mostly across epochs (MEASURED, `overlap_epochs.json`).**
- Pool families overlap on a target 1.60M times. Only **26.8%** of these are the same encounter (flyby epochs < 60 d
  apart); the median gap is **1243 d**.
- The pool uses a median of 6 encounter windows per target (10th percentile 3).
- So a craft that takes a target in year 2 usually removes a different craft's year-6 encounter. Advancing all
  craft together in time does not see this conflict when the decision is made.

**D5. A disjoint fleet lives in separated 6-D corridors, and state proximity predicts overlap (MEASURED,
`fleet_phase.json`, `junction_stats.json`).**
- A junction here is a Lambert transfer from craft A's twin state at t to craft B's state at t+T, with T = 30-150 d
  and total dv <= 1 km/s.
- **s16a:** none of its 72 ordered route pairs has a single junction anywhere in 15 yr. The median pairwise
  heliocentric-longitude separation is 89.7 deg.
- **Pool:** there are 8.77M junctions, joining 58,202 ordered family pairs (3.9% of all pairs).
- Joined pairs share a median of 6 targets, against 2 for random pairs. Only 2.3% of them are disjoint, against
  24.8% of random pairs; at dv <= 0.3 km/s the figure is 0.11%.
- Junctions are early: median day ~950, with 90% before day ~2980. Tracks diverge after that.

**What binds, in one line.** Coverage is limited by the set of corridors (6-D tubes in time) that our generators
ever fly. Recombining or reselecting existing tubes cannot repack them (§3). A craft can only take a target it was
built around, and a target on another tube's window costs an excursion of at least 530 d (§3, X3). New coverage
therefore needs NEW track pieces, designed around the congested targets' late windows for several craft at once and
chosen jointly.

## §2 Proposals

### P1 — State-proximal recombination (fragment stitching). BUILT, RAN: FAIL as a coverage engine

(a) **Precedents.**
- GTOC11, Tsinghua (+SISE): beam search over a database of 3-8-asteroid flyby fragments chained into mothership
  trajectories, then a GA picks 10 motherships with disjoint asteroid sets (docs/methods_survey.md; junction model
  [recall, unverified]).
- VRP 2-opt* tail exchange (Potvin & Rousseau 1995) [recall, unverified].

(b) **Mechanism.** A child is A's flybys before t, then a Lambert junction t -> t+T, then B's flybys after t+T. B's
impulses are unchanged, so the unsettled child is already feasible in the impulsive model. Children then enter the
max-coverage MILP.

(c) **Laws.**
- It escapes law 1: no target is forced, since every segment was flown around its targets.
- It escapes law 6: junctions go orbit-state to orbit-state, not along an epoch chain.
- It escapes law 5: there are no planner legs.
- **It does not escape law 2.** Junctions exist only between tracks that already overlap (D5), so children
  reshuffle targets inside a corridor and never across corridors.

(d) **Test and thresholds.** Pre-registered before the run; see §3.

(e) **Payoff and probability.** Coverage stays at 286. P(8-craft fleet) about 1%.
- **Salvage 1:** children are an honest column source at about 5 s each, usable for the 9-craft cost frontier
  (not run).
- **Salvage 2:** `junctions()` is a cheap 6-D same-encounter overlap predictor for generators.

(f) **Effort.** 2 h, already spent. Reuses run_ialns (ipr, settle), s15b_regrid and the colcache.

### P2 — Joint prefix re-planning with exact tail selection ("congestion LNS"). RECOMMENDED next build

(a) **Precedents.**
- GTOC12, ACT & Friends: many candidate ships per round, then a 0-1 LP maximum-weight independent set
  (methods_survey [V]).
- Multi-ship BIP + SCP campaigns of 37-39 ships (arXiv 2411.11281 [V]).
- GTOC12, Tsinghua winner, ~35-40 ships on one asteroid set, generate-select-regenerate
  [recall, unverified].
- GTOC9, JPL winner: campaign-level resynthesis over a mission database; NUDT: ACO where each ant builds all
  missions (docs/stage9).
- VRP set-partitioning matheuristics (ILS-SP, Subramanian et al. 2013) [recall, unverified].

(b) **Mechanism.** Take the 286 selection (or any 8-route base), and do the following.
- Pick K = 3-4 routes to re-plan (the tails at 24-31), plus optionally one head. Cut each at a common epoch t_c
  near day 2500. Every one of the 12 misses has a pool-flown window after day 2500 (D3); the earliest are
  17@2522, 118@2549, 247@2624, 108/243@~2680.
- Freeze the other routes and every prefix.
- Generate **many candidate tails per cut route in parallel**, from the same residual. Use
  tools/s15b_prefix.py-style twin tails with trading prizes (own 1.0; fleet misses 1.01-1.02, in the safe premium
  range of stage 13 §1.5).
- Each candidate steers to a DIFFERENT subset of the 12 misses, and some candidates also free-take the other cut
  routes' current tail targets. This deliberately diversifies the candidates so that law 2's "everyone wants the
  same easy targets" produces alternatives rather than duplicates.
- Regrid and settle every tail (law 8). Candidates are prefix + tail columns.
- An exact MILP chooses one column per cut route, together with the frozen routes: maximise coverage, then minimise
  sum J_i. Reuse xover_scan.maxcov / s16_select.solve.
- Iterate on a different K-subset and cut epoch. Both stay fixed within one iteration.

(c) **Laws, and why this does not re-measure the 3x premium.**
- **Law 1.** Nothing is inserted into a finished route, and nothing is assigned before its trajectory exists. Each
  miss sits in a tail that its own beam built around it, which is the 0.032-0.042 J regime. The only coupling is
  the MILP choosing among flown columns.
- **The failed joint K=2.** That was a joint BEAM over tuple states, so its width was diluted. Here each tail beam
  runs at full width, alone, and the joint part is an exact selection over the product of candidates.
- **Law 3.** All cut routes are regenerated from the same residual at once. No route is "last".
- **Laws 5 and 8.** Twin-prefix beam, regridded.
- **Admitted risks.**
  - The re-plans could collapse depth, as in linchpin (24-29 of 37). The linchpin re-planned from scratch inside
    assigned pools, while P2 keeps a flown prefix and uses free tails; s15b's prefix re-plan of h06 kept depth.
  - P2 does not escape the corridor limit (D5). It only lets each late tail choose its corridor knowing what the
    others will take.

(d) **Falsification test (<= 60 min, 2 procs).**
- Base: `milp.json:parents_N8` (286).
- Cut the four tails (31/31/29/24) at the flyby nearest day 2500.
- 6 candidate tails each (3 miss-subset prize variants x 2 beam widths), then the exact MILP over 6^4 combinations.
- **PASS:** >= 292 covered at total fuel <= 3858 kg (the 09-27 N=8 bar).
- **FAIL:** <= 288, or any regridded tail above its column's planned tank by > 60 kg.

(e) **Payoff and probability.**
- **Budget arithmetic.** The base's fuel is 3380 kg (tanks sum to 8180 kg). The 09-27 bar leaves 478 kg, i.e.
  about 40 kg for each of the 12. That is exactly the designed-around price (30-40 kg) and far below the forced
  price (100-140 kg). So P2 works only if every one of the 12 is designed in, never inserted.
- If it reaches 292-294, the stage-13 closer (cheapest ~8 at ~0.026 J) could close 298 at sum J_i ~11.3-11.5, i.e.
  raw 13.3-13.5 at 8 craft (rank 6-8 territory).
- P(submittable 8-craft fleet by 09-27 06:00) about **5%**.
- P(the 60-min test passes) about 15%.

(f) **Effort.** 5-7 h to build: a driver around s15b_prefix + s15b_regrid + the MILP. Runs need 6-10 h on 8-10 cores.

### P3 — Event-driven concurrent fleet beam with a coverage-potential value function. NOT recommended as posed

(a) **Precedents.**
- GTOC12 multi-ship beams whose state is the tuple of ships (methods_survey §2 "balanced construction").
- GTOC9 JPL beam search for campaign synthesis [recall, unverified].
- GTOC7 motherships + probes: not transferable. We have no mothership, and s16a's craft launched on the same day
  have zero junctions afterwards (D5) [GTOC7 team details: recall, unverified].

(b) **Mechanism.** All 8 craft advance in time; the craft with the earliest clock extends next; ownership is shared;
the score is flybys - fuel + V(remaining).

(c) **Laws. My data falsify both premises before any build.**
- 73% of conflicts are cross-epoch, a median 3.4 yr apart (D4). Time-ordered expansion resolves only the 27%
  same-encounter conflicts at decision time.
- The natural V (reachability or scarcity) is blind to exactly the targets that get lost (D3: misses at the 83rd
  popularity percentile).
- A tuple beam also dilutes width: joint K=2 was already worse than sequential on 65-77 pools.
- If it were ever built, V must price **congestion**: demand for a corridor's windows against the ~40-flyby time
  capacity of one craft. That makes it a transport or flow problem, not a count of reachable targets.

(d) **Falsification test (30 min, 1 core).**
- Compute, per target, a congestion index from the pool: the number of families flying it within 60 d of each
  other, divided by the capacity of those families' corridors.
- **PASS:** its AUC for "missed by the N<=8 optimum" is >= 0.75.
- **FAIL:** < 0.65. For comparison, the scarcity AUC is below 0.5 by D3.
- Only if it passes does a concurrent beam get built (12-20 h).

(e) **Payoff and probability.** About 1-2% for a submittable fleet by 09-27.

(f) **Effort.** 12-20 h. Reuses tools/s14_twinbeam.lin_price and the ownership logic of ctoc14/jointsearch.py.

## §3 The measurement I ran: X1/X2/X3 (all 1-2 cores, about 35 min wall total, results/s17/campaign/)

Pre-registration (written into this file before any run):
- M1 = N<=8 max coverage over parents plus 2-segment children.
- PASS >= 292, FAIL <= 286.
- X2 PASS if >= 4 of 6 children settle with tank <= estimate + 40 kg.

| step | what | result (file) |
|:--|:--|:--|
| states | 1220 families x 548 grid epochs (10 d), twin states + cumulative dv | 651k states, 1 min (`states.npz`, `reps.pkl`) |
| junctions | ballistic KD-tree prefilter, then vectorised Lambert, leads 30/50/80/110/150 d | 29.4M candidates, **8.77M junctions <= 1 km/s**, 3 min (`pairs.out`) |
| children | A-prefix + junction + B-suffix, tank <= 1250 kg, B impulses <= 1.05 tcap | 81,977 distinct sets (45k rejected on tank, 192k on tcap); 39,790 non-dominated; 42,650 columns in all |
| control | parents only, N<=8 | **286**, LP 290.35 (`milp.json`, `lns_N8.json`) |
| M1 LP | parents + children | **LP bound 293.74**, i.e. children add +3.4 to the LP |
| M1 MILP | HiGHS, 600 s | incumbent 271, dual bound 293.66: no proof (`milp.json`) |
| M1 LNS | from the 286 pick: drop k, EXACT MILP refill over all 42,650 columns | k=1,2 (greedy): 286; **k=2 exact: 286; k=3 exact (56 drops, 404 s): 286**; k=4 stopped unfinished (`lns.out`) |
| X2 | settle the 5 children of the MILP incumbent (regrid + RI.settle) | **5/5 settle**, miss 39-77 km, tank settled - estimate = -14.8/-2.0/-0.6/+0.9/+0.8 kg, 3-8 s each (`settle_with_children_N8.json`) |
| X3 | 3-segment detours A -> B -> A to the 12 misses | 893 families fly a miss, 64,666 / 63,808 junctions picks<->families, but only 25 miss flybys are bracketed by an out-junction and a return-junction; shortest bracket **530 d, costing >= 3 own flybys (median 11)**, so **0 net-positive detours** (`detour.json`) |

**Verdict against the pre-registration.**
- M1 = 286 (exact refill for k <= 3; the MILP bound of 293.7 is not closed). This is **FAIL** as written.
- X2 **PASS**: the children are honest, and settling them costs nothing.

**Interpretation.**
- The LP gain (+3.4) is fractional mixing of overlapping children. Every integer neighbourhood up to 3-route
  exchanges returns the same 286.
- X3 is law 1 seen from the state side. A target on another track is reachable only by leaving your own track for at
  least 530 d, so a "chosen" target on someone else's corridor costs 3-11 of your own flybys.

**This kills P1 as a coverage engine.** It supports P2's premise: the misses have late windows that many tracks do
fly (D3), so they must be designed into new tails, not stitched or detoured in.

**Not claimed.**
- M1 optimality with k >= 4 exchanges is unproven. The 293.7 bound leaves room; my prior on an integer 292 is
  < 10%.
- The children were not regridded before the MILP; only the 5 settled ones were.
- The 9-craft min-cost question (can children cut s16a by 200 kg?) was not run.

**Housekeeping.**
- `junctions.npy` (491 MB) was deleted after `junction_stats.json` was written. It regenerates in 3 min with
  `xover_scan.py pairs`.
- No shared tool or result file was modified. Nothing was submitted.
