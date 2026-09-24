# Stage 5 (2026-09-19): the "orphan wall" was a search artefact -- skeleton-first fill for N=8

Supersedes the conclusion of `docs/stage4_n8_results.md` ("N=8 structurally blocked"). Written while the experiments
run; the measurements below are from `results/n8/` scratch runs (scripts in the session scratchpad, logic recorded here).

## 1. What t10d's own routes say (`t10d_anatomy`)
| | easy legs (231) | hard legs (67, prize>=1.5) |
|---|---|---|
| leg dv mean / median / p90 | 0.48 / 0.36 / 1.07 km/s | **0.51 / 0.44 / 1.18 km/s** |
| leg tof mean / p90 | 175 / 284 d | 181 / 300 d |
- The 67 "orphans" cost the SAME as easy targets when the route is planned around them. The 100-280 kg figure of
  stage 4 was the cost of INSERTING one into an already time-full route, not a property of the targets.
- Every t10d route is TIME-full (13.3-14.9 yr) at ~175 d/flyby; the free beam paces 116 d/flyby. Eight craft need
  298/(8*15 yr) = 147 d/flyby. Time efficiency of the sequencing, not fuel, is the short resource.
- Stripping each t10d route to its hard flybys and re-settling (`results/n8/skel10`): 6-7 hard targets cost only
  1.8-6.7 km/s per skeleton (tanks 629-714 kg), i.e. ~0.58 km/s per hard target. Hard skeletons are cheap.

## 2. Why every beam-based N=8 method failed
The beam has no way to be FORCED through a rare-event target: it only has prizes, and a prize distorts the whole
route (prize beam: 26 flybys/14 hard; even the unwindowed beam over "easy + 6 specific hard at prize 2" gives 29
flybys at 0.69 km/s -- attracting the beam to six targets costs 17 flybys).  Insertion growth (`grow_skel8`) is the
other fill mechanism and it stalls at 265 covered with 1180-1270 kg tanks. Neither mechanism was the leaders' 8-craft
route (they exist: leaderboard N=10 rows are 8 craft + 2 unreachable, raw 12.5-13.1).

Leg-model check (`legmodel2`, t10d route 01, departure = the twin's real flyby state): twin 16.6 km/s; linear model
15.8; multi-rev Lambert 20.3; single-rev Lambert 57.9 (one 368-d leg alone = 39 km/s). Coast misses <= 0.25 AU.
So the beam's leg models CAN price the good routes; single-rev Lambert is the only bad one (legs > ~300 d) and
`lin_drmax 0.15 AU` excludes ~6 of 36 legs (0.15-0.25 AU). The wall is the constrained SEARCH, not the physics.

## 3. New beam capability: mandatory time windows (`ctoc14/search.py`)
`Params.win_lo/win_hi` (300,): a windowed target is mandatory -- legs to it are only generated inside the window and
a state that passes `win_hi` without visiting it is dropped (deadline prune, `_alive`).  `P.lookahead` (default on
when windows are set): children are ranked by score + fuel of the cheapest Lambert leg to their next unvisited
mandatory target (window start/mid/end), A*-style.  `None` = old behaviour, all existing tools unchanged.

Hard-only chains (E1): under the default leg limits a hard-only greedy covers 35/67 in 8 routes (16,6,4,3,2,2,1,1);
with 900-d Lambert / 1000-d coast legs **60/67** (25,17,10,3,2,1,1,1) -- but at ~1 km/s per hard leg. Better skeletons
come from t10d itself (section 1).

T1 (windowed fill, t10d skeleton epochs +-45 d, beam 300, no lookahead) vs the t10d route through the same hard set:
| route | waypoints | windowed beam | t10d |
|---|---|---|---|
| 01 | 6 | 29 fb, 0.63 km/s/fb (planner) | 37 fb, 0.45 |
| 09 | 7 | **33 fb**, 0.65 | 28 fb, 0.68 |
| 13 | 12 | 23 fb, 0.60 | 30 fb, 0.51 |
Right ballpark; T2 (running) is the search-strength ladder: lookahead, beam 1000, drmax 0.25, windows +-90 d.

## 4. Plan (if T2 reaches ~35+ flybys per skeleton)
1. Skeletons: `results/n8/skel10` (10 hard chains from t10d) -> merge to 8 by inserting the hard targets of the two
   smallest chains into the other eight skeletons (`grow_skeletons.py --hard`, cheap: 7 pinned flybys in 15 yr, not 30).
2. Fill: windowed beam per skeleton over the easy targets (greedy over the 8 skeletons with exclusion, best-first by
   the skeleton with the fewest options), settle the twin (`settle_tour`), keep the fleet in `IFleet` format.
3. Leftover easy targets: `grow_skeletons.py --max-dtank 60` (first ~5 insertions per route are cheap), then
   `run_ialns.py search --relocate` passes, export, `validator2.py`.
Budget: J < 13 <=> mean tank <= 995 kg at 298 covered (37 flybys at <= 0.52 km/s each).

## 5. Results of the search-strength ladder and the segmented fill (2026-09-19, 01:00-02:00)
Windowed (one beam, all windows) on route 01 (6 wp): beam 300 -> 29 fb; +lookahead(250 d) 30; beam 1000 -> 32;
+drmax 0.25 +-90 d -> 34 (planner 937) but such tours often fail to seed the twin (single-junction Lambert cannot do
the long linear legs; `impulsive.LIN_FALLBACK` now seeds those legs with the linear thrust profile); dv_max 2.0/2.5
buys depth only with fuel (33 fb at 1029 kg).  Settled: 33 fb -> 940 kg twin (0.53 km/s/fb).  Route 09 (7 wp): 33 fb
-> **877 kg twin, 0.448 km/s/fb** (t10d's route 09: 28 fb at 976 kg).  With 9 waypoints the windowed beam drops to
25-28 fb: every deadline kills most of the beam (284 of 5781 states survive all windows), so easy targets between
waypoints are packed at ~320 d each.

**Segmented fill (`fill_segmented`, `--mode segmented`)**: waypoint-to-waypoint beams, each segment must END at the
next waypoint (easy targets get a window ending at the segment end, non-mandatory via `Params.mandatory`), K=30 states
kept per waypoint, a free tail after the last one.  Route 05 (9 wp, the balanced skeleton): **39 flybys, twin 1018 kg,
J_i 1.388 in 190 s** (windowed: 25 at 791).  Pareto 37:954 38:962 39:996 (planner) = exactly the J<13 line
(37 fb at <= 995 kg; 39 fb allow 1023 kg pro rata).  Route 09 failed at the merged-in target 217 (skeleton leg 1.9 km/s
over 479 d, beyond the 400-d Lambert / 0.15 AU linear legs) -> per-segment relaxed retry (700-d Lambert, 1000-d
coast, 0.30 AU, waypoint cap 4 km/s), else the waypoint is DROPPED to the leftovers.

Skeleton pipeline used: `results/n8/skel10` (t10d stripped) -> `grow_skeletons.py --hard-only` (9 hard targets of
routes 02/04 into the other 8: 7 of 9 at <= 1 kg) -> `tools/rebalance_skel.py --max-hard 9` (13->9 on route 13, all
moves <= +5 kg) = `results/n8/skel8b` (9/9/9/8/8/8/8/8 hard, 642-734 kg).  Full fill: `results/n8/fill2`.

## 6. Full fill and the leftover problem (02:00-02:50)
`results/n8/fill2` (segmented, portfolio +-45/+-90 d, lam 0.05): 05:36@998 09:36@1061 13:33@939 03:33@919 01:30@896
06:25@815 12:25@840 10:23@867 -> **241 covered, sum J_i 10.24 (J 12.24 if the 57 leftovers were free)**.  The first
routes fill from a rich pool (27 easy each), the last from a depleted one (74 targets -> 17 easy): the greedy tail
moved from the hard targets to the awkward easy ones.  Two fixes that did NOT work: (a) refill of one route over its
own easy + the leftovers (84-target pool: shallower, 27 vs 36, and the tours no longer settle); (b) dual prices on the
leftovers x1.8 and a fresh pass (`tools/dual_fill.py`): every route gets shallower and dearer (05: 33@919, 13: 29@1067,
09 nothing settled) -- pricing awkward targets drags each route through bad legs.  Insertion can absorb ~15-20 of 57
(first ~5 per route at <= 30 kg), so the fill itself must reach >= ~280.  Next: order test (are the late skeletons
intrinsically shallow?) and `--assign` (geometric pre-assignment of easy targets to skeletons, balanced, + pad).

## 7. Capacity, columns and column generation (03:00-05:00)
- Order test: skeletons 12/10/06 filled FIRST from the full pool reach 39/33/32 (25/23/25 when filled last) -> the
  late-route shortfall is pool depletion, not the skeletons.  Geometric pre-assignment (`--assign`, group + pad =
  ~52 targets per route) is too thin a pool (28 fb).  Wide segmented search (beam 1000, K 100) on skeleton 12: 40 vs 39
  -> ~40 flybys is a skeleton's ceiling; per-route capacity ~37-40, fleet ~305 vs 298 needed.  Packing decides.
- Column mode (`--columns`, 5 diverse columns per skeleton by exclusion at prize 0.3): 247 columns, ~105 distinct easy
  targets reachable per skeleton; only 2 easy targets (41, 64) reachable by NO column, 6 by one skeleton, most by
  3-7.  MIP (`tools/skel_select.py`, one column per skeleton, 1 J per miss) over those: 223 covered -- columns are
  diverse within a skeleton but all skeletons' first columns take the same popular targets.
- Column generation (`tools/skel_cg.py`: LP relaxation duals -> prizes of the next columns; the beam's cost - sum prize
  is the reduced cost): LP coverage 240.8 -> 247.5 after one round (+7/round, dual-priced columns are sparse, ~21 easy).
- Materialisation (`tools/skel_materialise.py`): all 8 chosen columns settle; twin/planner tank within +-3-7%.
- CG result (4 rounds, floor 0.05): LP coverage 240.8 -> 247.5 -> 250.7 -> 253.0 -> 254.3; final MIP over 460 columns
  **229 covered** (integer gap ~25).  Worse than the sequential fill (241).  Sparse dual-priced columns (~21 easy) are
  the reason; next attempt should price with a floor ~0.35 so columns stay deep, and add more columns per round.
- Wide segmented search (beam 1000, K 100): 40 fb but the tours do not settle (legs the twin cannot seed) -> keep
  beam 300 / K 30 for production.

## 8. Decision 05:50 (2026-09-20): bank N=9 first
N=8 needs <= 2 misses; tonight's best is 241 covered.  N=9 skeletons (`results/n8/skel9b`: t10d chains minus route 02,
8/8/8/8/8/8/7/6/6 hard, 637-686 kg) have ~9 x 28 = 252 easy capacity vs 231 -> `results/n9/chain.sh` (fill, insert,
export + validate `CTOC14_Result_n9a.txt`, relocate passes, export `n9b`).  Expected raw ~13.3-13.6 (shown ~13.1-13.3
vs 13.81 banked).  N=8 / J<13 remains open: better pricing in CG, deeper fills, or waypoint epochs free to move.

## 9. N=9 result and where this stands (06:45, 2026-09-20)
`results/n9/fill1`: 9 skeletons -> 234 covered (03:37@1058 05:37@983 09:34@936 13:28 10:26 12:26 01:21 06:19, 04 stayed
bare -- its waypoint set lets 1-3 beam states through each segment).  Insertion (`results/n9/grow1`): 263 covered after 5
rounds, then feasibility dries up (12 of 113 trials) -> ~35 unplaceable.  The sequential fill saturates at ~235-245 for
N=8 and N=9 alike: the late routes always face the same awkward remainder, and insertion adds ~25.
**t10d stays the submission** (validator VALID, raw 14.366).

What is established tonight (all measured, see s1-s8): the hard targets are not the obstacle; skeleton-first works;
a skeleton route caps at ~40 flybys; full N=8 coverage is a packing problem at the capacity edge (8 x ~29 easy vs 231);
sequential fills, refills, dual prices, geometric assignment and a first column-generation loop all stop at 223-254.
The leaders' 8-craft fleets prove it is solvable.

Candidate next steps (each a few hours, none tried yet):
1. CG with a prize FLOOR ~0.35 (deep columns biased to the duals) and 3 columns/skeleton/round, 8+ rounds; then the
   MIP with a time limit of 30 min.  The LP bound rose every round (241 -> 254); the integer gap (25) was the loss.
2. Two-route LNS on the 241 fleet: re-fill routes A and B jointly over their union + leftovers (pool ~120), accept on
   coverage - cost; cycle over pairs.  Pools of ~120 give ~30 flybys per route (measured), so pairs can trade.
3. Free the waypoint epochs: each hard target has several crossings; let the skeleton choose per route among 2-3 epochs
   (the merge/rebalance step only ever used t10d's).  This changes which easy targets line up and may lift the ~40 cap.

## Decision (2026-09-20 06:50): STOP
User chose to stop. t10d (N=10, raw 14.366, validator VALID, banked 13.81 shown) remains the final submission. All
solver jobs stopped. The stage-5 tools and findings above are the starting point if the N=8 packing problem is picked
up again (s9 options 1-3).
