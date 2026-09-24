# Next strategy: fewer craft, not cheaper flybys (2026-09-17, after J = 14.366 with 10 craft)

> **SUPERSEDED by `docs/strategy_stage2.md` (2026-09-17 evening):** the leaderboard's N column is craft + misses, so the
> top five are EIGHT-craft fleets at our per-flyby cost (0.48 vs our 0.49 km/s). Stage 2 targets 8 craft by column
> generation with twin-costed columns; 9 craft is not worth the time cost.

> **REVISED the same day (see docs/globalopt_campaign.md sections 8.1-8.6).** Idea A below was implemented and tested
> and does NOT work; ideas B and C keep their ranking but for a different reason. The measured economics:
> removing a craft is only worth it if the fleet's total Delta-v grows by less than 23.5 % (34.3 km/s of the current
> 145.8), while a marginal target moved between routes costs 1.6-9 km/s against a fleet average of 0.49 km/s per
> flyby. The 10-craft fleet is converged (400 extra SCP iterations move the tanks by 0.2-0.6 kg) and every marginal
> operator is exhausted. Section 1's claim that "the leaders' ~12.5 is 8 craft at OUR efficiency" holds only if the
> per-flyby cost stays at 0.49 while a craft carries 37 targets; ours degrades with every consolidation
> (0.394 km/s per flyby at 13 craft -> 0.489 at 10), and that degradation, not the craft count, is the real gap.

Diagnostics: `results/newgen/scratch/{anatomy,linecover,farcands,scarce}.py` on `results/newgen/ifleet10`;
event catalogue `results/newgen/scratch/events.npz` (kind 0 = ecliptic-node crossings with |r-1| < 0.08 AU, 1669;
kind 1 = local minima of the distance to the 1 AU circle below 0.08 AU, 2666).

## 1. What the 10-craft fleet says
- 298 flybys, 145.8 km/s, **0.49 km/s per flyby**; per route 0.37 (route 12) to 0.68 (route 09).
- Out-of-plane impulses are 37-66 % of each route's Delta-v; the costliest 20 % of legs carry 50 % of the fleet Delta-v.
- Median flyby geometry: 0.037 AU off the 1 AU radius, 0.024 AU off the ecliptic. 22 of the 60 costliest legs target
  asteroids that have NO crossing within 0.03 AU of the 1 AU circle (asteroids 40, 163, 217, 296, 11, 31, 149, 150, ...):
  their cost is geometric. 111 of 298 asteroids have no such event; 20 have none within 0.05 AU (the "hard" set:
  1, 31, 39, 50, 66, 70, 74, 75, 96, 107, 118, 140, 149, 169, 179, 186, 269, 270, 280, 295).
- Rigid orbits (no thrust after launch) cover few targets: greedy cover with in-ecliptic craft at node crossings reaches
  28, 24, 17, 14, ... targets per craft (168 with 14 craft); inclined rigid craft 152. Phasing and plane changes are
  unavoidable, so per-flyby cost cannot drop by a large factor: re-choosing events for the costliest legs is worth
  ~10 km/s (~0.2 J) by a linear estimate.
- **Leaderboard arithmetic**: J = N + 2 + sum(x + x^2). At our 0.49 km/s per flyby and 298/N flybys per craft:

  | N | flybys/craft | Delta-v/craft | tank | J_i | J |
  |---|---|---|---|---|---|
  | 10 | 30 | 14.6 km/s | 873 kg | 1.23 | 14.3 (actual 14.37) |
  | 9 | 33 | 16.2 | 909 | 1.27 | 13.4 |
  | 8 | 37 | 18.2 | 957 | 1.32 | 12.6 |
  | 7 | 43 | 20.9 | 1024 | 1.40 | 11.8 (0.5 N cap binds) |

  The leaders' ~12.5 is 8 craft at OUR efficiency, not 10 craft at a third of our Delta-v. Each craft removed is worth
  ~0.9 J; route 01 already flies 37 flybys at 0.45 km/s each, so dense routes are not less efficient.
- Why dissolves fail: the dissolved route's scarce targets (few usable events) find no host that is phased right at one of
  those events. Hard/scarce targets per route (`scarce.py`): route 01 (37 flybys) and 03 have no hard target; 04 and 13
  have 4 each; 09 holds 39 (52 deg inclination, no event within 0.05 AU), 270, 118.
- Time coefficient: k grows 0.36 %/day (0.05 J-equivalent per day at J 14.4); a 9-craft result (-0.9 J) pays back within
  ~2.5 weeks of extra work, an 8-craft result (-1.8 J) within ~5 weeks.

## 2. Ideas (ranked by expected gain / effort)
A. ~~**Miss-tolerant dissolve + global candidates (9 craft, ~1-2 days).**~~ **DONE AND REFUTED.** Implemented as
   `insert_block` (adaptive backtracking homotopy), `w_block` (block insertion with ejection of the host's flybys in
   the window) and `dissolve_lns` (ruin and recreate, priced per net placed target). Block moves DO place targets the
   old operator could not (8, 39, 168, 270 went onto routes 03 and 06), but the dissolve of route 09 placed only 9 of
   28 targets for 1.241 J and stalled: 0.10-0.25 J per target against the 0.044 J the dissolve can afford. Original
   text kept below for the record.
   Original: Let a dissolve leave unplaceable targets as
   misses (+1 each, the true objective) instead of aborting; afterwards treat every missed target as a prize of 1.0 for
   insertion anywhere. Candidates come from the event catalogue (all events of the target within 0.05 AU, ranked by a
   linear phasing estimate 2 dphi v / (3 n T_lead) + offset terms), not only from close approaches of the current
   trajectory (< 0.3 AU at the same time), so the host may drift for years to reach the event. Add ejection chains:
   if a missed target cannot enter host H, evict H's cheapest-saving target and re-place it. Dissolve routes with no
   hard targets first (03, 01) or the cheapest tanks (12, 10); anneal allowance 1.5-2.0.
B. **Anchored construction for 8-9 craft (~2-3 days).** Reverse the failure mode: assign the 20 hard and ~40 scarce
   targets first to K anchor routes at their few events (small MIP; compatibility = a smooth phase line through the
   events in the linear model), then grow each route by cheapest insertion over the catalogue (grow_fleet.py logic with
   global candidates and the impulsive twin), then ALNS with A. The old from-scratch growth failed because the easy
   targets were taken first and 101 hard ones were left over.
C. **Column generation with an event-space pricing (~1 week).** Reuse the master of `ctoc14/colgen.py` (set cover,
   craft cap, duals); replace the Lambert beam search by a beam search over catalogue events in time order whose state is
   the impulsive route (re-optimised incrementally: appending a target changes only the last rows) with duals as prizes.
   Then set partitioning at N = 8, 9 and the ALNS. Largest potential (fresh sequences built for phasing), largest effort.
D. **Efficiency squeeze at fixed N.** Implemented as `run_ialns.py retime` (every OTHER pass of each asteroid from the
   event catalogue, ranked by the reduced model). NOTE the ranking's plane term was 7x too pessimistic in the first
   run; `phasemodel.K_Z = 8.8 km/s per AU` is the calibrated rate.
E. **Heavy-route care.** With 33-37 flybys tanks reach 950-1000 kg and the 0.5 N cap binds; keep `enforce_cap`
   (0.43 N x bin / tank) in every settle and expect more 0.05/0.025 d polish retries at export.

## 3. Workflow (revised)
1. DONE: block insertion, LNS dissolve and global retiming are in `tools/run_ialns.py`
   (`search --lns`, `retime`) and `ctoc14/globalopt.insert_block`; the reduced cost model is `ctoc14/phasemodel.py`.
2. Do NOT spend more compute on dissolve/relocate variants of the 10-craft fleet: they are over budget by ~4x.
3. The only remaining lever is a route GENERATOR that keeps the average per-flyby cost at 33-37 targets per craft.
   It must carry the eccentricity vector in its state (a coast arc has 4 free parameters and hits exactly 2 events),
   so the sound basis is the existing Lambert beam search (`ctoc14/search.py`) rescored with the reduced model, driven
   by column-generation duals (`ctoc14/colgen.py`) at N = 9. The phase-plane DP (`ctoc14/phasedp.py`) is NOT usable:
   it grants the eccentricity slack independently at every event and its 10 d / 1 deg grid lets the phase slip ~10 deg
   per step.
4. Time cost of continuing: the leaderboard coefficient grows 0.1/28 per day = 0.051 J-equivalent per day at J 14.37,
   so a redesign attempt must be expected to win more than ~0.05 J per day it takes.
