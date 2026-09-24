# Stage 9 — what the GTOC winners actually do, and the architecture we were missing (2026-09-21)

Target changed to **J < 13.0**. That is rank-3 territory (NUAA 12.97, THU 12.37, HIT 12.52 raw) and cannot come from
tuning: stages 6-8 established that every fleet we build lands on 10 craft near J 14.35, with a measured 3x premium
for any target forced onto a finished route.

## 1. The literature (structurally identical problems)

| competition | team | what they did |
|:--|:--|:--|
| GTOC9 "Kessler run" (123 debris, MANY missions, cost = fixed + variable per mission — **isomorphic to ours**) | NUDT (2nd) | Three levels: **top = partition debris into missions by ANT COLONY OPTIMISATION (DCB_ACO)** — each ant builds ALL missions in one constructive pass, pheromone learns across passes; middle = hybrid-encoding GA for sequence + transfer times; bottom = differential evolution for the precise trajectory. 123 debris in 12 missions. |
| GTOC9 | JPL (winner) | Branch-and-bound for long rendezvous chains exploiting natural nodal drift, **beam search for campaign synthesis**, ACO, GA — and "**databases of transfers between all bodies on a fine time grid ... an easy-to-compute yet accurate estimate of the transfer Delta-V**". |
| GTOC11 "Dyson ring" | Tsinghua + SISE | Pre-analysis to shrink the search space; **beam search over a pre-built database of 3-8 asteroid flyby FRAGMENTS** to make many single-mothership trajectories; **a genetic algorithm then selects the 10 motherships**; a separate database of optimal rendezvous times for all 83 453 asteroids by phase angle; greedy assignment. 388 asteroids flown by. |
| GTOC4 (44 NEA flybys) | Moscow State (winner) | Global sequence search with **bi-impulse** approximation first (48 flybys), low-thrust refinement second (44). |
| 2025, arXiv 2508.02904 | — | DP with the Markov state "only velocity and mass at the flyby epoch matter", global optimality for a FIXED sequence with an explicit error bound; improved the GTOC4 winner by 20.2 kg; neural-network fuel predictors + GPU. |

**The common architecture, and the one piece we never had:**
1. **Precompute a transfer database** so trajectory design becomes table lookup. Every winner does this. We re-solved
   Lambert inside every beam node, so one route cost 5-10 minutes and the fleet layer only ever saw ~500 routes.
2. **Generate at fragment/chain granularity**, not whole routes (Tsinghua's 3-8 asteroid fragments).
3. **Assemble the fleet with a learning constructive metaheuristic over the WHOLE partition** (ACO/GA), not a greedy
   pass and not a selection over independently built routes. This is what dissolves our 3x premium: an ant assigns
   every target exactly once BY CONSTRUCTION, so nothing is ever forced onto a finished route, and pheromone learns
   which targets to save for a later craft — precisely what greedy construction cannot do.

## 2. What was built (all new, stage 9)

| tool | what it does | cost |
|:--|:--|:--|
| `tools/encounter_db.py` | every epoch a target can be met near 1 AU: each ring pass (r in [0.8, 1.25] AU, abs(z) < 0.12) SAMPLED every 6 d, not just its abs(z) minimum | 26 694 nodes, seconds |
| `tools/transfer_db.py` | Lambert arcs between all node pairs (20-420 d), every revolution branch, departure AND arrival velocity stored | 9-13 M edges, ~4 min |
| `tools/graph_beam.py` | beam on the LINE graph (state = last edge, since a flyby fixes no velocity: junction = abs(v_dep(p->q) - v_arr(o->p)), both in the table) | **0.6 s per route** (was 5-10 min) |
| `tools/graph_fleet.py` | greedy and **ACO** fleet search; covered targets kept as zero-prize stepping stones | **0.5 s per 12-craft fleet** |

Validation of the graph against a real 10-craft fleet's 298 flybys: 99 % lie within 5 d of a database node, 70 %
within 2 d, none missing (the first, minima-only database held only 20 % within 2 d — hence the 6 d sampling).

## 3. Where it stands: 1000x faster, not yet better

| generator | best route |
|:--|:--|
| carrier DP + insertion (stage 7-8) | 40-48 flybys, 992-1143 kg, **0.032-0.034 J per flyby** |
| t10d's best route | 37 flybys, 918 kg, 0.0346 |
| **graph beam (stage 9)** | 24-32 flybys, 896-1064 kg, **0.042-0.055** |

**Diagnosed cause.** The graph prices every junction as one impulse, abs(v_dep - v_arr). Our best routes make their
money on exactly the legs where that is wrong: long, near-180-degree, phasing legs whose single-revolution Lambert
solution is degenerate and expensive, but which a continuous-thrust arc flies cheaply (`impulsive.LIN_FALLBACK`,
`ctoc14/linleg.py`). Evidence: capping the junction at 0.6 km/s collapses the route to 7 flybys, and at 0.35 km/s to
4 — in the pure-Lambert graph cheap junctions barely exist, while the real fleet averages 0.49 km/s per flyby.
Storing every revolution branch (9.1 M edges) did not fix it, because the problem is the impulsive model itself.

**The fix, and it is precomputable exactly like the rest:** the edge cost must be
`min(Lambert junction, LinLeg continuous-thrust cost)` for long legs. `ctoc14/linleg.py` already solves that leg model
in batch and it is the model the planner used in stages 5-6; it just has to be evaluated once per long edge and
stored next to the Lambert velocities. Until then the graph is a fast search over a cost model that misprices the
best legs.

## 4. Honest position against the clock

The banked 13.8019 needs raw J < 14.153 today, tightening ~0.05/day to 13.849 on 09-27 (t_end 2026-09-28 11:59,
28 d span, from our own k = 13.801851/14.366055 = 0.96073). Best in hand: 14.355 (exact set cover over 489 routes,
LP bound = MIP optimum). J < 13.0 needs ~1.4 J, i.e. the architecture above working end to end, not a tuning pass.
