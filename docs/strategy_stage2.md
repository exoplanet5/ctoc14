# Stage 2: eight craft, or nothing (2026-09-17 evening, after idea 1 held J = 14.366)

Idea 1 of `docs/strategy_next.md` (block insertion + LNS dissolve + global retiming) is implemented, tested and
refuted: TEAM stays `results/CTOC14_Result_t10d.txt`, J = 14.366055, 10 craft, 298 covered. This document sets the
next stage from two facts that were not on the table before: what the leaderboard's N column means, and the real
time-coefficient rule.

## 1. New facts

**The leaderboard's N is craft + misses, not craft.** Decoded from our own rows (snapshots in `leaderboard/`):

| our submission | craft | misses | shown N |
|---|---|---|---|
| v4 (09-07) | 12 | 2 | 14 |
| milp3 (09-15) | 16 | 6 | 22 |
| cg4 (09-16 09:10) | 14 | 2 | 16 |
| go1 (09-16 22:13) | 13 | 2 | 15 |
| d11 (09-17 01:18) | 12 | 2 | 14 |

So the top five (N = 10) are **eight-craft fleets with the two unreachable misses**, and the N = 9 entry is seven craft.
Undoing k = 1 - 0.1 (t_end - t_submit)/28 d (t_start 08-31 12:00, t_end 09-28 12:00) and the even-split tank model:

| team | shown | raw J | craft | tank | Delta-v per craft | fleet Delta-v | per flyby (298) |
|---|---|---|---|---|---|---|---|
| HHIT-SAT | 11.936 | 12.52 | 8 | 952 kg | 18.0 km/s | 144 km/s | **0.484** |
| TTHU-LAD | 12.070 | 12.65 | 8 | 967 | 18.6 | 149 | 0.500 |
| NNUAA | 12.296 | 12.97 | 8 | 1004 | 20.1 | 161 | 0.539 |
| eeStar | 12.996 | 13.69 | 8 | 1081 | 23.0 | 184 | 0.617 |
| (7-craft entry) | 12.827 | 13.55 | 7 | 1228 | 28.0 | 196 | 0.658 |
| ours t10d | (13.81 if submitted now) | 14.366 | 10 | 873 | 14.6 | 145.8 | 0.489 |

**Our per-flyby cost already equals the top team's.** The whole gap is ten craft against eight (2 x ~0.9 J). The
earlier claim "the leaders' ~12.5 is 8 craft at OUR efficiency" (`strategy_next.md` section 1) is now confirmed by
the leaderboard, not inferred. Route 01 (37 flybys, 16.6 km/s, 0.45 per flyby) is what every one of the leaders'
eight routes looks like; routes 10, 12, 03 are as efficient (0.40, 0.36, 0.44) but short (27, 30, 27 flybys), which
is a partition problem, not a route-quality problem. The 7-craft entry shows 42-43 flybys per craft is feasible.

**t10d is not on the leaderboard yet** (our row is d11, 15.16). Submitted now (k = 0.9615) it shows 13.81, rank 8.
The leaderboard keeps the best shown score, so submitting it now has no downside.

**Resubmission break-even** (a later result must satisfy J_new k_new < 14.366 x 0.9615):

| days from now | k | raw J needed | gain needed |
|---|---|---|---|
| 3 | 0.972 | 14.21 | 0.16 |
| 5 | 0.979 | 14.10 | 0.26 |
| 7 | 0.987 | 14.00 | 0.36 |
| 10 | 0.997 | 13.85 | 0.51 |

**What each fleet size is worth** (even split; shown score at the date it could realistically be submitted):

| fleet | fleet Delta-v | raw J | shown | rank |
|---|---|---|---|---|
| 9 craft | 146 (no growth) | 13.43 | 13.14 @09-22 | 7 |
| 9 craft | 160 (+10 %) | 13.78 | 13.48 @09-22 | 7 |
| 8 craft | 160 (+10 %) | 12.95 | 12.81 @09-25 | 4 |
| 8 craft | 170 (+17 %) | 13.25 | 13.11 @09-25 | 5-6 |
| 8 craft | 180 (+23 %) | 13.56 | 13.42 @09-25 | 7 |

Conclusion: **nine craft is not a worthwhile final target** (it barely beats t10d after the time cost); eight craft
at <= 170 km/s is worth 1.1-1.4 J and rank 4-6. Any 10-craft polish (launch-date sweep, tank balancing, worth
0.02-0.04 J) is below the daily time cost and only makes sense bundled with the 8-craft result.

## 2. Why 8 is within reach of our tooling
- 298/8 = 37.25 targets per craft. Route 01 carries 37 at 0.45 km/s per flyby with our SCP, so the per-route
  regime exists in our own fleet. The budget for 8 craft is 160-170 km/s, i.e. **+10-17 % on today's 145.8**;
  consolidation 11 -> 10 cost +12.4 %, and 10 -> 8 by a fresh partition (not by redistribution) does not have to
  pay the marginal price twice.
- What failed was REDISTRIBUTION (a marginal target costs 11x the average because the hosts are pinned by ~30
  flybys). What is untested is a PARTITION built for 8: which 37 targets fly together. That is a column-generation
  question, and the master (`ctoc14/colgen.py` RMP/MIP, duals, craft cap) already exists; the weak part was the
  pricer's cost (planner estimate 3.2 J pessimistic at N = 10, section 8.8 of `globalopt_campaign.md`).

## 3. Ideas, ranked
F. **Gate experiment: twin-recost the dense pool (~1 h).** Cost the 1 665 pool columns with >= 35 targets (and
   the 12 001 with >= 30) with the impulsive twin (`impulsive.ImpulsiveProblem` + `globalopt.restore/optimise`,
   ~1.5 s settled, ~10 s with the restore of a planner tour), instead of the Lambert-junction estimate. Output: the
   true per-flyby cost distribution of planner-born 35-38-target routes. If the best are <= 0.5 km/s per flyby the
   planner produces 8-craft-grade material and idea G is the pricer; if they are >= 0.7 the planner's sequences are
   structurally worse than SCP-consolidated ones and idea H is the pricer.
G. **Column generation at N = 8 with twin-costed columns (2-3 days).** Pricer = the existing beam search with
   duals as prizes (`colgen.price`, waypoints for covered targets); every column is recosted by the twin before it
   enters the RMP, so the master sees true tanks (J_i = 1 + x + x^2). Iterate until the easy targets' duals are
   ~0 and the columns carry the 101 hard targets; final set partition at N = 8 (MIP over twin costs, or a dive);
   import the 8 routes into `run_ialns` (ifleet), polish (relocate/retime/swap), export with the validator.
   Per-round metric: LP value at N = 8 with twin costs; number of hard targets in columns with negative reduced cost.
H. **Fallback pricer: dual-driven skeleton growth (2 days).** Seed 8 skeletons from the flown routes' hard-target
   subsequences (routes 01, 05, 02, 12, 13, 04, 10, 03 keep theirs; 06 and 09's hard targets go to the best host by
   `host_rank`), then grow all eight simultaneously by cheapest twin insertion (`insert_block`, global event
   candidates from `events2.npz`) with the remaining targets priced by duals. Mechanism: a bump costs ~4 dr / T over
   the free window T; a 13-flyby skeleton has 2-3x longer windows than a 30-flyby host, so early insertions cost near
   the average, not the marginal price. `grow_fleet.py` failed from scratch because the easy targets went first and
   101 hard ones were left; skeletons that already carry the hard targets remove that failure mode.
I. **Polish bundle, applied only to an 8-craft fleet (+0.05-0.1 J).** Launch-date sweep (tL is FIXED in the SCP;
   routes 03/04/06/09/13 launch 0.6-1.5 yr late, so their first phasing legs have less lead than they could),
   relocate/retime/swap passes, tank tightening at export.
J. Note only: the 7-craft entry (0.66 km/s per flyby, 196 km/s) is within our 7-craft break-even of 213 km/s, but
   it is not a target for this stage.

## 4. Workflow with gates and dates (deadline 2026-09-28 12:00; last safe export day 09-27)
0. **Today (09-17):** submit t10d (shown 13.81). Start gate F in the background (~40 min on 8 cores). Meanwhile
   build the twin-recost hook: planner tour -> `ImpulsiveProblem` (launch state from the tour, 20-d impulses
   initialised from the Lambert junction Delta-v, `restore` then `optimise`, `enforce_cap`), returning tank and J_i.
1. **09-18:** pricing rounds at N = 8 (G, or H if F fails). Log per round: LP(N = 8) with twin costs, coverage of
   the hard set, best column per-flyby cost.
2. **Gate 2 (09-19 evening):** LP(N = 8) <= 13.5 and an integer cover of 298 exists -> continue; otherwise stop,
   keep t10d, and spend no more compute (every further day costs 0.05 J-equivalent).
3. **09-20/21:** dive or MIP -> 8 routes -> `run_ialns import` -> polish passes (I) -> export at 0.1/0.05/0.025 d
   with the validator -> submit if shown < 13.81.
4. **09-22 onward:** if 8 craft are in hand at <= 170 km/s, a second pass of the same pipeline at N = 7 is the only
   remaining lever; otherwise stop.
