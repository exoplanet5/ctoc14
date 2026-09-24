# Stage 4 result: N=8 is structurally blocked; the lever is a cheaper N=10, not fewer craft
Written 2026-09-18 late. Records what the branch-A/B/C runs actually showed, so N=8 is not re-attempted blindly.

## What was run (all at m0 1600, vinf 4.0, beam 300, the stage-3 levers)
| branch | method | result |
|---|---|---|
| A | fleet tree (`tools/fleet_tree.py`), gamma 1.0 | 3 routes, **73 covered**, then every depth-4 child pruned as coverage-infeasible |
| A | same, gamma 0.6 | depth-1 identical to gamma 1.0 (the node-value hard penalty, not the beam, drives route choice) |
| B | k-routes (`tools/kroutes.py`), hardness prizes | round 0 **177 covered** |
| B | k-routes, uniform base | round 0 **192**, round 1 **172** (declining, not climbing) |
| C | MIP (`tools/partition_mip.py`) over ~1600 distinct columns | N=8 **177**, N=9 **178** |
| (prior) | `greedy_cover.py` m0 1600 | **231** covered + 67 orphans -- the real disjoint-coverage ceiling |

## The finding: a hard ~231-235 ceiling on disjoint beam-route coverage
Every method that builds DISJOINT deep routes by beam tops out at ~231-235 of 298 covered. The other ~63-67 are the
structural orphans (far, inclined: median a 2.32 vs 1.78 AU, i 16.4 vs 11.6 deg). They do not chain into deep routes
(greedy: the 67 orphans chain into a single 6-flyby route) and are expensive to insert (100-280 kg each, tank
compounding). So covering all 298 needs ~10 craft -- which is exactly what t10d is.

Two mechanisms were newly measured tonight and both are dead ends for N=8:
1. **Hardness prizes trade depth for hard density.** A prize beam gives a 26-flyby/14-hard route; uniform gives
   46-flyby/4-hard. So a "hard-first" sequential tree (branch A) makes shallow fleets (8x26=208 << 282), and prize
   k-routes captures *less* of each group (177 < uniform's 192). Prizes cannot buy coverage.
2. **k-routes rounds do not climb capture.** Re-beaming reassigned groups loses coherence: 192 -> 172. This answers
   the stage-3 open question ("does capture climb from 75%?") -- no.

Root cause is the same wall the whole prior campaign hit (colgen "converges to 10.2 fractional craft"). The stage-3
mass-ceiling finding raised how DEEP one route can go (34 -> 46 flybys), but it did not change the ORTHOGONAL fact
that the 67 orphans cannot be packed into disjoint deep routes. N=8/9 at full coverage is not reachable this way.

## The achievable lever instead: a cheaper N=10 (stage-3 fallback 3)
t10d is N=10, raw 14.366, all 298 covered, banked 13.81 shown. Its routes were built at m0 1000 / vinf 2.0. Rebuilding
the SAME 10-route partition at m0 1600 / vinf 4.0 made every fresh route 6-8% cheaper at equal depth (stage-3 s6.3),
i.e. ~0.5-0.9 J on the fleet with NO new algorithm and no coverage risk. Raw 14.37 - ~0.7 = ~13.7 raw -> ~13.1 shown,
under the 13.81 break-even. This -- not N=8 -- is the realistic next improvement. An N=9 via the growth/LNS pipeline
(not beam) is a larger, riskier bet for maybe another ~0.3 J.

## Recommendation
Stop investing in an 8-craft partition; it is blocked by the orphan structure, confirmed by branch A, B and C tonight
and by the entire prior campaign. Keep t10d as the submission. If pursuing more: rerun the 10-craft pipeline at
m0 1600 / vinf 4.0 for a cheaper N=10. Assets kept: `results/n8/pool.jsonl` (~8600 columns), `results/n8/prizes_H.json`,
`tools/fleet_tree.py`.

## Decision (2026-09-18): STOP at N=10
Reviewed the evidence and chose to keep **t10d** (N=10, raw 14.366, all 298 covered, banked 13.81 shown, submitted
2026-09-17) as the final answer. N=8 abandoned as structurally blocked. No further compute. The cheaper-N=10 rerun
remains available in stage-3 fallback 3 if revisited later.
