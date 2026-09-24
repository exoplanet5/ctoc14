# Column generation + neighbourhood search campaign (2026-09-15/16)

Goal: a fleet-building method that reaches J < 20 (docs/design_colgen.md, docs/design_alns.md). Start: validated
J = 30.789 (results/CTOC14_Result_TEAM.txt, 16 craft, 294 covered).

## 1. Speedups (results identical to 1e-12, regression in results/colgen/scratch/regress_*.json)
| change | effect |
|---|---|
| `Ephemeris.ast_states_at` (vectorised asteroid states for (index, time) pairs) replaces 53k scalar calls per 20 expands in the encounter-time refinement | refinement cost 132 ms -> 17 ms per expand |
| Lambert: Stumpff functions evaluated once per iteration (bit-identical) | Lambert 2x faster |
| `kepler.propagate_batch` (one universal-variable propagation for many initial states) in `LinLeg.__init__` | linear-leg setup 11x faster |
| total | refined beam search 3.5x faster (114 s -> 33 s on the regression case); a camp5-style round ~4-5 min on 5 processes |

`search.Params` gained `prize` (per-target prizes in J units, also used in child pruning), `w_t`, `first_targets`,
`dv_max_t` (per-target leg caps; tested, see section 4). `jointsearch` can collect its beam (`P.collect_joint`).

## 2. Plan-level re-timing is an artefact (analysis/retime_probe.md)
Coordinate descent on flyby times cuts planned propellant 15-16 % (camp1 4499 -> 3758 kg). Exact conversion of the two
largest savers: camp1 sc5 planned -92 kg -> converted +24 kg and 2 flybys lost; camp2 sc3 planned -100 kg -> converted
-5 kg. The converter already optimises flyby times; do not re-time plans for fuel.

## 3. Column generation (ctoc14/colgen.py, tools/run_colgen.py)
Column store = 282k pool routes (+ generated), cost model tank fixed point m0 = 620 / (1 - 1.15 F / m0_search).
Master LP exact over the whole store (iterative pricing of all stored columns by the duals, 6 s). Pricing = refined beam
search with the duals as prizes (w_fuel 3.43 = marginal J per kg, w_t 0.9), every beam state is a column.

| stage | LP (planned J) | notes |
|---|---|---|
| pool only, N = 12 / 11 / 10 | 18.70 / 24.70 / ~39 | N >= 11 degenerate duals (mu 13, 48 targets at 1.0) |
| run1: 4 rounds N = 12 (8 beams, 2 min each) | 17.10 | cap not binding (sum y 11.03) |
| run1: 5 rounds N = 11 | 16.32 (sum y 10.66) | still not binding |
| run1: 7 rounds N = 10 | 21.49 | forcing 10 craft costs ~5 J; duals degenerate (43-52 at 1.0) |
| run2: 2 more rounds N = 11 (1000/1100 kg) | 16.06 (sum y 10.47) | LP optimum at ~10.5 craft |

Integer fleets from these columns:
| method | craft | covered | planned J |
|---|---|---|---|
| HiGHS MIP over the working set (60-120 s) | 11-12 | 243-250 | 65-72 (useless) |
| MIP on the LP support + 300 nearest columns, N = 10 / 11 | 10 / 11 | 256 / 267 | 58.0 / 48.2 |
| LP-guided dive without re-pricing, N = 10 / 11 | 10 / 11 | 258 / 266 | 56.1 / 49.4 |
| **priced dive** N = 11 (2 beams per level, look-ahead on top-3) | 11 | 278 | 37.2 |

The dive's look-ahead LP climbs ~5 per level after the second fixed route (16.7, 17.5, 23.4, 28.4, 33.2, ... 37.2):
fixed whole routes leave residuals that the pricing cannot cover with complementary routes. The integrality gap
(LP 16 vs integer 37) is the whole problem.

## 4. Neighbourhood moves with the planner (tools/run_lns.py)
State after flyby k reconstructed exactly (tools/exp_retime.py Chain, max |ddv| 1e-9 on pool and generated routes).
Moves: tail re-plan from flyby k (prize 1 on targets no other route covers), spawn a new craft over the uncovered
targets (launch 0-2000 d), dissolve the weakest route, **pair move** = joint re-plan of two routes from launch (shared
visited mask, `jointsearch`), **group move** = remove the k weakest routes and re-plan k-1 or k craft jointly.

| path | craft | covered | planned J |
|---|---|---|---|
| no-pricing dive N=10 (snap1) | 10 | 258 | 56.06 |
| + tail re-plans, 3 spawns (12, 5, 5 flybys) | 13 | 282 | 35.31 |
| + pair move of the two 5-flyby spawns -> 9 + 9 | 13 | 290 | 27.45 |
| + pair move, spawn (2 flybys), pair move | 14 | 293 | **25.42** |
| priced dive N=11 + tail/pair/spawn (round 1-2) | 12 | 280 | 35.1 (stalls: full routes, rarest leftovers) |

Per-target leg cap 2.0 km/s into targets priced >= 0.5 (N = 10 duals): worse pricing (62 vs 517 negative columns,
rc_min -1.05 vs -1.91), because at N = 10 almost every target is priced >= 0.5 and expensive legs shorten routes.

## 5. Waypoints (2026-09-16 01:30) — the key fix
Every planner leg must END at an asteroid (Lambert <= 400 d, linear <= 600 d). Excluding the targets other routes already
cover therefore removes the stepping stones between sparse uncovered targets: full re-plans with exclusion returned 2-5
flybys. With covered targets visitable at prize 0.001 (`run_lns.py --prize-covered 0.001`; dive pricing in
`run_colgen.py` now keeps covered targets visitable, `--dive-exclude` = old behaviour) the 5-flyby craft became 17 flybys
and the 6-flyby craft 13, planned J 25.42 -> 20.64 in 5 minutes. The same exclusion exists in `convert_tour`'s re-plan
after a failed window (`global_excluded` = other tours' targets) — only relevant when a window fails; not yet changed.

## 6. Results after waypoints (2026-09-16 02:25)
| fleet | craft | covered | planned J | notes |
|---|---|---|---|---|
| results/colgen/fleet_2049 (lns_wp1: full re-plans of the weak craft with waypoints) | 14 | **298** | **20.49** | first fleet covering every reachable target |
| results/colgen/fleet_2048 (+ tail re-plans at split 0.9) | 14 | 298 | 20.48 | being converted (tools/convert_pool.py) |
| waypoint dive N=11 (run3) | 11 | 281 | 34.23 | dive stays greedy: 16.5, 18.0, 24.1, 26.2, ... 34.2 |

Exact conversion of the 25.42 fleet (results/colgen/conv_a, tanks 600 + 1.3 F + 20): **14/14 valid, zero dropped
flybys**, actual 64-451 kg vs planned 57-358 kg at the 1000 kg search mass (0.98-1.17x at the conversion mass). The
cost model (s = 1.15, R = 20) predicted the fleet J within ~0.1.
Dissolve (remove the weakest craft, absorb its 8-10 unique targets by tail re-plans of nearby craft at split 0.33):
rejected so far (absorbs ~3 of 8).

## 7. First validated files of the new approach (2026-09-16 03:00-03:50)
- Conversion of every route of the 20.48 fleet: 21 routes (tools/convert_pool.py + conv_a), all valid, **no flyby dropped
  in any pass**.
- Fragment MILP over all 177 fragments (old + new): 14 craft, 298 covered, **J = 20.339** (results/CTOC14_Result_cg2
  selection; cg1 = 20.542 with part of the fragments, being validated).
- Right-sized tanks still carry 26-52 kg unused (pass 2 sizes from the heavier pass-1 burn): `tools/tighten_frags.py`
  re-converts at 600 + used + 8 kg; expected gain ~0.5 J (-> ~19.8).
- CG-LNS (tools/run_cglns.py) with LP duals of the residual: no gain, the residual LP is degenerate (current routes
  optimal, dual weight on 1-2 targets) -> pricing now uses unit prizes on the residual.

## 8. Below 20 (2026-09-16 04:40)
- `results/CTOC14_Result_cg1.txt` = TEAM: **J = 20.5415**, 14 craft, 298 covered, validator + validator2 VALID.
- Tank tightening round 1 (`tools/tighten_frags.py --margin 8`): 13/13 fragments valid with the identical flyby set;
  m0 down 18-44 kg per craft (sum J_i -0.45). Fragment MILP (dry): **J = 19.889** (`results/CTOC14_Result_cg3.txt`, being
  validated). The tightened craft still end with 11-24 kg unused -> round 2 (`--margin 6`) running.
- Dissolve with 1150 kg full re-plans of the absorbing routes: absorbs 5 of 8 targets of the weakest craft (22.48 vs
  20.48) -> rejected; the small craft hold late-window targets that time-limited big routes cannot reach.
- Second fleet: waypoint dive (11 craft, 34.2) + waypoint LNS -> 12 craft, 293 covered, planned 22.46 (spawn of a
  19-flyby craft over 16 leftovers).

## 9. Status 2026-09-16 06:40
- TEAM = `results/CTOC14_Result_cg4.txt`, **J = 19.7173** (both validators): cg3 fragments after a second tightening
  round (`--margin 6`, all 13 kept).
- Three independent fleets, all covering 298 with 13-14 craft: fleet_2048 (14 craft, planned 20.48, converted -> cg4),
  lns_d3 (waypoint dive + LNS, 13 craft, planned 19.55; being converted/tightened by `tools/finish_fleet.sh ... cg5`),
  lns_d2w (priced exclusion dive + waypoint LNS, 13 craft, planned 19.32, still improving).
- Plateau: the big routes are time-limited (14-15 yr), small craft hold late-window targets, single-route moves cannot
  remove a craft once everything is covered; dissolve has not been accepted yet (best: absorbs 5 of 8).
- New structural experiment: greedy campaign WITH waypoints (`results/colgen/greedy_wp`: seed = camp5's four best routes,
  141 targets, then one waypoint spawn per round). Hypothesis: the exclusion greedy decays (camp1: 29..8, 3) because it
  loses waypoints; with waypoints it may hold 25-30 per round -> 10-11 craft.

## 10. Four fleets (2026-09-16 07:20)
| fleet | construction | craft | covered | planned J |
|---|---|---|---|---|
| fleet_2048 | no-pricing dive + LNS (+ waypoints) | 14 | 298 | 20.48 (converted: **cg4 = 19.7173 validated**) |
| lns_d3 / fleet_d3c | waypoint dive N=11 + LNS | 13 | 298 | 19.546 |
| lns_d2w / fleet_d2w | priced exclusion dive + waypoint LNS | 13 | 298 | 19.295 |
| greedy_wp | camp5's 4 best routes + 9 waypoint spawns (40 min) | 13 | 298 | 19.409 |
Greedy with waypoints decays like the exclusion greedy (27, 27, 24, 21, 21, 19, 19, 17, 15 new per round), so waypoints
fix the repair moves, not the greedy construction. All four fleets land at 13-14 craft and planned 19.3-19.5: the
structure (big time-limited routes + small late-launch craft for the rare targets) is a strong local optimum.
Dissolve attempts (up to 1150 kg full re-plans of the absorbing routes) are all rejected.
Next: convert + tighten every fleet's routes and let the fragment MILP mix them (`tools/finish_fleet.sh`).

## 11. Multi-revolution legs: no gain (2026-09-16 07:30)
`ctoc14/lambert.py` gained a vectorised multi-revolution solver (`lambert_multirev`, `lambert_all`, `lambert_best`,
`tests/test_lambert_multirev.py`: propagation-exact for nrev 0-2, both branches, 33k-167k rows/s). Re-pricing the 115
Lambert junctions of six planned tours with nrev 0-2: **0 legs improve** (sum dv 62.7 km/s unchanged). Multi-revolution
solutions need tof above ~700 d (median tof_min), while our legs are 100-400 d. Not a lever for this problem.

## 12. Mixing fleets in the fragment MILP + recombination (2026-09-16 09:30)
`tools/finish_fleet.sh <fleetdir> <tag> <nproc>` = convert the fleet's new routes (`convert_pool.py`, keyed by a
(launch, flyby-sequence) signature so nothing is converted twice) -> `fleet_frags.py` picks each route's best fragment ->
`tighten_frags.py --rounds 2` -> selection-only fragment MILP over every valid fragment.
- cg5 (fleet_2048 + lns_d3 fragments, 211+ fragments): **13 craft, 297 covered, J = 19.551** (11 routes from lns_d3,
  2 from conv_a; asteroid 27 was the single flyby the converter dropped, and covering it costs more than the miss).
- The MILP selection is itself a planned fleet (`results/colgen/fleet_mix5`, planned J 19.532 with all 298): feeding it
  back into the waypoint LNS is a recombination step (path relinking), running now.
- Recombination (LNS restarted from the MILP mix, `results/colgen/lns_mix5`): **no improvement** (19.532 unchanged after a
  full round of waypoint re-plans, spawn and dissolve) — the mixed fleet is already a local optimum of these moves.
- `results/CTOC14_Result_cg5.txt` = TEAM: **J = 19.5511** (13 craft, 297 covered, both validators; the only avoidable miss
  is asteroid 27, dropped by the converter on one route).

## 13. J = 18.60 (2026-09-16 12:10)
`results/CTOC14_Result_cgX.txt` = TEAM: **J = 18.5989**, 13 craft, **298 covered** (only 131/144 missed), both validators.
Two levers on top of the fleet search:
1. **Dropped-flyby recovery**: one selected route had lost its LAST flyby (asteroid 27, 1.1 km/s over 140 d) because the
   converter ran out of propellant (used 478.8 of 483.3 kg). Re-converting that tour with an 80 kg larger tank keeps all
   24 flybys; after two tightening rounds the craft costs +0.05 J and removes a 1.0 miss. **Always check the registry's
   `n_plan` against the fragment's flyby count and re-convert the losers with a bigger tank.**
2. **Tightening to 6-12 kg of residual propellant** per craft (two rounds, `--margin 6..8`).
Fragment MILP over 302 fragments (four fleets, all tightened) picks 13 craft: 2 from conv_a, 10 from the waypoint-dive
fleet, 1 the recovered big-tank route.

## 14. Final of the session (2026-09-16 14:30): J = 18.4917
`results/CTOC14_Result_TEAM.txt` = `results/CTOC14_Result_cgF.txt`: **13 craft, 298 covered (only 131/144 missed),
J = 18.4917**, validator PASS + validator2 VALID (`results/validator2_cgF.log`). Tanks 700-1130 kg (J_i 1.09-1.52).
Fragment MILP over 356 fragments from four fleets, each route tightened 3-4 times (residual propellant 4-12 kg).
Session path: 30.789 -> 20.541 -> 19.889 -> 19.717 -> 19.551 -> 18.599 -> 18.492.

### Where the remaining gap to the leaders (raw ~13.1-13.6 with 9-10 craft) sits
Our 13 craft cover 9-38 targets each; the LP bound over 300k routes is ~16 planned J for 11 fractional craft, so even a
perfect integer selection from the current route pool would be ~17. The leaders' fleets need ~30 unique targets per craft
with almost no overlap; our long routes carry 5-11 repeat flybys (used as waypoints) and the rare late-window targets
need dedicated small craft. Closing that needs a route generator that produces longer *disjoint* routes, not better
allocation of the current ones: all four independent constructions (dive, priced dive, greedy with waypoints,
recombination) converge to 13 craft at planned 19.3-19.5.
