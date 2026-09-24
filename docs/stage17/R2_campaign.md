# Stage 17 R2 — CAMPAIGN lens: cross-examination (2026-09-23 23:55 CST)

MEASURED numbers name a file under `results/s17/`. Budgets use F5: 3723 kg of fleet fuel at the 09-27 bar, and about
46 kg for each hard-core target.

**Headline.**
- Four lenses have now measured the same law from different sides. **With the existing tracks, coverage is conserved.**
  - Selection, stitching and insertion arcs give 286, 286 and 285 (F3).
  - Forcing any single hard target into the selection drops coverage from 285 to 278-283, i.e. it costs 3-8 others
    (`redteam/r2/forced_455.out`).
  - Suffix re-plans conserve DEPTH: r6 re-plans to 30 and r8 to 24, their original counts (F1).
  - Trades conserve ISOLATION: every dropped own target is itself isolated (`lns/r2/chain_screen.json`: best link
    saves -118 kg).
- **New zero-compute join (`campaign/r2/chain_slots.json`).**
  - Only 6 of the 13 hard targets have even ONE cheaply droppable own flyby (re-home <= 30 kg) within ±250 d of any
    surviving host's approach (< 0.2 AU). Within ±150 d the count is 3.
  - 36, 108 and 216 have no host within 0.2 AU at all.
  - Closure needs about 11 of 13 at the designed price (f >= 0.84, LNS 1.2). So **every re-plan-existing-tracks
    variant is capped below closure even at zero trade cost.**
- **My P2 is conceded.** Its generator is F1's generator, and LNS is running its gate (`lns/r2/L2r6k17O`, using my
  prizes).
- **My quick check** was the 9-craft min-cost frontier with the stitched children. Result: **FAIL, no 9-craft side bet**
  (§3).
- **Recommendation:** no new structural build. s16E stays banked. The last in-flight gate (LNS open pool) finishes under
  a fixed reopen trigger (§4).

## §1 Cross-examination

### PHYSICS P1 "pins and grout" — objection MAJOR (fatal on the clock)
- **(a) Its gate cannot decide anything.**
  - PASS is 262 covered on the stage-5 fill2 A/B. We already hold 285-286 at 8 craft (F3), so a PASS at 262 is 24
    targets WORSE than the existing ceiling, and it would not justify the 12-16 h build.
  - A decision-relevant gate would need >= 287, i.e. +46 on a 241 baseline. Nothing in stages 5-8 moved a fill by more
    than about 10.
- **(b) The pin class misses 5 of the 13 targets that actually bind** (`physics/events.json`, d < 0.05, my count).
  - 8 of r9's hard 13 are pins: 36, 75, 108, 118, 212, 216, 219, 243.
  - 17, 180, 220, 247 and 288 have 5-9 events, so they are "medium".
  - Pins are everywhere in the fleet: r5 and r7 hold 14 each, r4 12, r1 10. The binding set is REDTEAM's
    fast/inclined subfamily (v_rel 21.9 vs 16.4 km/s), not the window count.
- **(c) Skeleton-then-fill waypoints are the forced regime.**
  - Stage 6 measured it: capacity 304, but forcing coverage costs 2:1, so fills saturate at about 240.
  - The lanes cost about 0.7 km/s per pin in the model, which is about 37% of the 22.5 km/s a craft can spend at the
    09-27 bar, before any fill.
- **What I concede.** PHYSICS' D4 is the ONLY depth-positive lever anyone has put forward, and my data support it.
  - The s16a depth-vs-filler curve saturates: 0 fillers → 24 flybys, 6 → 30-31, 11 → 41, 15 → 38-41, 21 → 42.
  - 33 of the heads' 48 cheaply re-homable targets (<= 46 kg) are fillers (69%, against a catalogue base rate of 27%;
    `chain_slots.py` join).
  - Fillers are the movable currency, and r1's 10 excess fillers buy it about 1 flyby.
- **Rescue.** A from-launch concurrent build that uses the filler quota as its congestion price. That is my R1 P3's
  missing "V"; see E2. It will not fit in the window as a production run.

### PHYSICS P2 (earliest-deadline prize) — MAJOR
- It is still a prize. Prizes are the dead family, because they trade depth for hard density.
- Its test uses four SEQUENTIAL free routes, which is law 3's starving regime.
- For the 5 medium hard targets (5-9 windows), a deadline prize switches on only in the last ~2 windows, i.e. after
  day ~4000, when every track is already full (F1: no r6 or r8 state adds one without a drop).
- **Rescue:** as one term inside E2's value function, never on its own.

### PHYSICS P3 (leg surrogate) — MINOR
- It is enabling only, and it is trained on CHOSEN legs (law-1 selection bias; stage 10's R² was 0.01).
- There is no path from it to a fleet by 09-27.

### LNS P1 (suffix-column set partitioning) — FATAL
- The MILP gains coverage only through columns with net (added - dropped) >= +1. Every measured column is <= 0:
  - F1 on r6 and r8;
  - the s15b history: +1 for -1 or -2 (LNS 1.4).
- LNS's own chain screen finds no link with a positive saving.
- My slot census caps even the IDEAL version, with re-home-aware prizes that drop only cheaply re-homable own targets,
  at 6/13 hard targets. That leaves at least 7 at the forced price (160-1446 kg each; `chain_slots.json`) against a
  46 kg budget.
- **Rescue:** none by re-planning suffixes. It needs a column that adds without dropping (a "sink"), and no host has
  shown one.

### LNS P2 (ramped emptying) — FATAL alone
LNS agrees: with single-move columns it is forced dissolution, 3.6x over budget (`census_r9.json`).

### LNS P3 (free-launch prefix) — FATAL on the clock
- It needs a 6-8 h build.
- 36, 108 and 216 have no surviving host within 0.2 AU at any epoch (`chain_slots.json`), and 17's only approach is on
  day 4955.
- Prefixes are depth-conserved like suffixes: the same re-plan machinery, the same F1 law.

### REDTEAM P1 (closure-aware selection with lin arcs) — MINOR as a gate, FATAL as a route to a fleet
- It has done its job. F2 shows K >= 286 infeasible, and `forced_455` shows 278-283 whenever a hard target is forced.
  Its own kill rule ("min ΣJ_i at 298 > 11.54") fires trivially if K = 298 is infeasible.
- **Process defect:** the first full run crashed after one result (the `full.out` traceback: the output path was
  doubled). It was restarted at 23:34 as `full2_*`.
- **full2 (done 23:52) is uninformative:** N=8 maxcov hit 900 s at incumbents 26/81 (bound 298); N=9 K=298 found
  no incumbent at thr 0.08 or 0.15. The 455-column partial (F2) stays the only usable REDTEAM P1 result.

### REDTEAM P2 (13-target contract) — FATAL
- **The gate is mis-set.** It kills at a mean > 0.08 J including re-homing. But 13 × 0.08 = 1.04 J, against a 0.48 J
  budget, so a PASS would still lose. The pass line must be <= 0.037 J.
- **Its first three targets already fail on data.**
  - 75: no drop slot on its only host, r8, and no F1 state on r8 ever carried it.
  - 118: no slot on r2 or r5; its F1 state costs +306 kg.
  - 180: one slot (r7, own 241 at day 3892), but the forced price there is 96 kg.
- **Rescue:** none. It is LNS P1 with a narrower pool.

### REDTEAM P3 (early lane separation) — MINOR / component
- It rests on a real effect: 28 flybys in years 0-2, against 42-45 in every later 2-year bin.
- My D5 says the tracks are fully separated after year 2 (0 junctions in 72 ordered pairs), so the loss is launch
  crowding only.
- It needs every craft re-launched, i.e. E2's from-launch build. It cannot stand alone.

## §2 What the evidence changes in my own proposals
- **P2 (joint prefix re-planning with exact tail selection): CONCEDED.**
  - It is the same generator as LNS P1, with a common cut. F1, the chain screen and the slot census kill its premise
    that the misses' late windows can be designed into re-planned tails.
  - Re-priced: P(gate) 15% → **3%**; P(submittable) 5% → **< 0.5%**.
  - Its falsification test is now LNS's `L2r6k17O` (open pool of 72, my trading prizes). I do not duplicate it.
  - Even a PASS there (n* >= 33) comes from taking r5, r7 and r8 targets. That is zero-sum across tails, and it still
    needs a sink for the hard core.
- **P1 (stitching):** FAIL stands (286 at N=8). Salvage 1, the 9-craft cost frontier, is now measured dead too (§3).
- **P3 (concurrent beam):**
  - Merged with PHYSICS P1. In R1 I said P3 needed a congestion-pricing V, and the filler quota IS a congestion price
    on the only movable class.
  - D4 still stands: 73% of conflicts are cross-epoch, so time-synchronous expansion alone does not see them.
  - Stays at about 1%.
- **F4, resolved.**
  - All three readings are true:
    - 8/13 are pins;
    - all 13 are isolated from the s16a tracks (0.024-0.25 AU);
    - they are "popular" (89-232 of 1220 families) because stage 15's isogen STEERED deep columns onto them (depth
      46-48).
  - My R1 P2 depended on "popular ⇒ designable into re-planned tails", and that is refuted.
  - What survives of D3 is weaker but decisive: the hard core IS flyable at depth by FROM-LAUNCH routes. So the
    obstacle is the partition, not reach. Every live proposal now depends on that reading.
- **Budget:** my R1 figures (3858 kg of fleet fuel, about 40 kg per miss) are replaced by F5's 3723 kg and about
  46 kg per hard target, after roughly 500 kg for the cheap and mid targets.

## §3 My quick check: does the stitched-children pool cut the 9-craft cost? (2 cores, about 25 min)
- **Pre-registration:** `campaign/r2/PREREG.txt`, written at 23:40 before any result.
  - PASS: an improvement of >= 0.08 J (about -75 kg, the 09-25 12:00 bar).
  - FAIL: an improvement < 0.03 J.
- **Parents** (3347 columns, dominance-pruned; all 9 s16a sets are in the cache at their twin tanks), from
  `r2/frontier_parents.json`:

  | K | min ΣJ_i | misses | status |
  |:--|:--|:--|:--|
  | 298 | **11.7227** | none | optimal, and exactly s16a |
  | 296 | 11.6552 | 148, 212 | optimal |
  | 294 | 11.6052 | 106, 148, 212, 234 | incumbent, bound 11.457 |
  | N=8, 286 | 11.1835 | | reproduces R1 |

  - Each dropped target saves only about 0.03 J, and 212 is hard-core (forced price >= 195 kg).
  - So no "select fewer, close the rest" 9-craft path exists.
- **Parents + 39,790 children, K = 298** (`r2/frontier_children.json`): after 420 s the incumbent is still s16a's
  11.7227, with the bound at 11.426. The gap is the same LP fractional mixing seen at N=8 (LP 293.7 against integer
  286).
- **Exact k-exchange from s16a** (`r2/kx9.py` → `r2/kx9.json`): drop k of the 9 routes and refill exactly from
  parents + children. Result:
  - all 129 exchanges (k = 1, 2, 3) solved to optimality in 642 s;
  - **every optimum is s16a's own routes (gain 0.000), and no optimum uses a child.**
  - So s16a is 3-exchange-optimal at 9 craft over the whole honest cache plus the stitched children.
- **Interim, LNS L2r6k17O (23:51):**
  - depth 26 at t_end 4501-4959 d (F1: 4697 d), max S 2-4 r9 targets: WEAK by LNS' projection rule (4300-4600);
    the final n* is still pending.
- **Verdict: FAIL.** No 9-craft side bet. REDTEAM's §1.3 conclusion stands: a 9-craft -200 kg needs new tail tracks.

## §4 Ranked experiments (at most 3)
The cutoff is 09-27 06:00. Every probability below is a P(submittable fleet) under the bar at that cutoff.

**E0 (rank 1): bank s16E; no new structural build; let the last in-flight gate finish at zero extra cores.**
- **In flight:**
  - LNS `L2r6k17O`, open-pool depth, due about 00:05;
  - (REDTEAM full2 finished without incumbents, and my `kx9` is FAIL, so only trigger (a) is still live.)
- **Reopen triggers, fixed now.** Any one of these re-opens E1:
  - (a) L2r6k17O has n* >= 33 AND some settled column carries >= 1 hard-core target at net (dtank + census re-home of
    dropped own) <= 46 kg;
  - (b) full2 finds N=8, K=298 feasible at optimistic ΣJ_i <= 11.40;
  - (c) full2_N9 finds a 298 cover at optimistic ΣJ_i <= 11.64 (kx9 has already said no).
- **Kill:** no trigger by 09-24 01:00 → stage 17 closes and s16E stays banked. P(a trigger fires) about 3%.

**E1 (rank 2, only on a trigger): realise the triggering columns.**
- Regrid, settle, then s15b_close, then s16_iter/relocate, then export.
- **Gate:** settled, regridded twin ΣJ_i <= 11.54 at 8 craft, or <= 11.64 at 9 craft, by 09-25 12:00.
- **Kill:** > +40 kg drift from the plan on any column, or no gate pass by 09-25 12:00.
- 4-6 cores, 1 h build, 8-12 h run. P(gate | trigger) 30%, so overall P(submit) about 1.5%.

**E2 (rank 3, a lottery, only if the user explicitly wants to spend 3 days): a from-launch balanced concurrent build.**
- **Mechanism:** CAMPAIGN P3 × PHYSICS P1.
  - All 8 craft extend round-robin over one shared set.
  - A per-route filler cap of 11: the s16a curve is flat beyond it.
  - Hard targets as lane waypoints (T1b), with launch dates or phases spread per REDTEAM P3.
- **Law escaped:** law 3, sequential starvation, because no route is last. The quota prices the only movable class.
- **Gate:** REDTEAM's truncated horizon, fixed in its R1, re-stated here.
  - PASS: >= 78 distinct targets by day 1461, mean fuel <= 140 kg per craft, and no hard-core window before day 1461
    passes unflown.
  - KILL: < 72 targets, or > 160 kg per craft.
- 2 cores for the gate; build 14 h, gate run 3 h; the production run would need 8 cores for about 36 h.
- P(gate) 15%. P(submit) about **0.8%**: the build ends about 09-24 16:00, the production run about 09-26 06:00, and
  closer plus export leave no slack.

**Run-nothing case (I ENDORSE it, apart from E0's zero-cost triggers).**
- The best build (E2) has an EV of about 0.008 × 0.15 J ≈ **0.001 J** shown.
- Even one rank needs -0.009 J shown (13.503), and rank 8 needs -0.019 J.
- A 9-craft route to those ranks needs -85 kg by 09-25 12:00, and my §3 check shows the pool has none of it.
- Every lens now reports the same conservation law. Four more days of supervision buy about 0.001 J of expectation.
