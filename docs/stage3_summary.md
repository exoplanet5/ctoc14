# Stage 3 summary: the mass ceiling, the free launch energy, and the partition that is still missing
Written 2026-09-18 16:20. Supersedes `docs/stage2_instructions.md` sections 1-2 (the neighbourhood-beam pool and the
skeleton MIP were designed around a mass-capped generator and cannot reach eight craft).

## 0. Where we stand
| | |
|---|---|
| Submitted | `results/CTOC14_Result_t10d.txt`, raw J **14.366**, 10 craft, 298 covered, sent 2026-09-17, banked at **13.81 shown** |
| Best fleet on disk | still t10d (`results/newgen/ifleet_t10d`); nothing today beat it end to end |
| Deadline | 2026-09-28 12:00; last safe export/validate day 09-27 |
| Break-even for a resubmission | raw J < 14.26 (09-19), < 14.16 (09-21), < 14.00 (09-24), < 13.85 (09-27) |
| Running jobs | none |

Two real findings today, both free improvements to every beam we will ever run, and one obstacle that is now
precisely characterised. Details in `docs/globalopt_campaign.md` section 10.

## 1. Finding 1: the planner mass ceiling (the reason for the ten-craft equilibrium)
`beam_search(..., m0=X)` refuses a child that would take the craft below `M_DRY + m_margin`, so the planned
propellant is capped at `X - 640` and the resulting craft at `CostModel(0.70, 2).tank(X - 640, X)`.
**Every column ever generated in this project used m0 = 1000 -> tank <= 805 kg -> ~11.4 km/s -> ~34 flybys.**
A 37-flyby route needs ~18 km/s (930-970 kg) and was never representable. Column generation, the waypoint dive,
the skeleton MIP and greedy growth were all searching a space where no craft can carry more than ~34 targets,
so 10 x 34 ~ 300 was arithmetic, not a local optimum.

One beam over all 298 targets, prize 1 J per target, `w_fuel` from the cost model
(`results/newgen/scratch/m0_test.py`):

| m0 | tank ceiling | best route | planner tank | twin tank | J per flyby | end |
|---|---|---|---|---|---|---|
| 1000 | 805 | 34 flybys | 804 | 844 | 0.0343 | 14.40 yr |
| 1600 | 1038 | **46 flybys** | 1027 | 1067 | 0.0304 | 14.64 yr |
| 2000 | 1149 | 45 flybys | 1011 | - | 0.0307 | 14.63 yr |

m0 2000 is no better: at 46 flybys the route is **time**-limited (~116 d per flyby), not fuel-limited. Use m0 1400-1600.
The cost model needs no recalibration in the new regime: implied s over seven fresh routes = 0.689 (planner/twin
tank ratio 0.993).

## 2. Finding 2: the launch v_inf cap was throwing away free energy
`Params.vinf_cap` was 2.0 km/s in every run of the campaign; the contest allows |v_inf| <= 4 km/s and that energy
comes from the launcher, not the tank. Same beam, m0 1600 (`results/newgen/scratch/lever_test.py`):

| vinf_cap | flybys | twin tank | J_i | J per flyby | |v_inf| used |
|---|---|---|---|---|---|
| 2.0 | 46 | 1067 | 1.445 | 0.0314 | 1.23 |
| 3.0 | 46 | 1004 | 1.372 | 0.0298 | 2.11 |
| 4.0 | 46 | 1004 | 1.372 | **0.0298** | 2.11 |

-63 kg of tank, -0.073 J per craft, ~-0.55 J on a seven-craft fleet. 3.0 and 4.0 find the same optimum, so pass
**`--vinf 4.0`**: never worse, costs nothing. Note the twin optimiser already exploited the full 4 km/s (t10d's
routes launch at |v_inf| 0.8-3.85), so the gain is in the sequence the beam *chooses*, not in the final polish.

**Levers that do nothing:** raising the Lambert TOF cap 400 -> 700 d or the linear-model cap 600 -> 900 d returns the
identical route and doubles the beam time. Keep `tofs` 15-400 d, `lin_tofs` 20-600 d.

## 3. Finding 3 (tooling): the degenerate-tour screen
Some beam tours contain a Lambert leg that `from_tour` cannot reconstruct. The signature is unmistakable **before**
any optimisation: a junction impulse of 27.9 km/s and an initial tank of 3750-4130 kg, against 986-1464 kg for a
healthy tour. Settling one costs ~83 s and never converges; one k-routes group spent 1529 s losing 11 of 12
candidates this way. `ctoc14.impulsive.settle_tour` now takes `screen_tank=2000` and rejects them in ~0.1 s -- the
same group now takes **86 s** for the same answer. Impulses above `tcap` are *not* a problem: healthy tours settle
with 11-31 of them (`enforce_cap` absorbs them).

## 4. The obstacle: the fleet must be BALANCED, and our generators produce the opposite shape
J_i = 1 + x + x^2 is convex in the route length, so for a fixed 298 flybys the cheapest fleet is N *equal* routes,
and fewer craft always wins because the "+1" per craft dominates. Target arithmetic at our measured ~0.49 km/s per
flyby: **N = 10 -> 14.28 (= t10d), N = 8 -> ~12.4, N = 7 -> ~11.9, N = 6 -> ~11.1** (6 needs 50 flybys per craft,
above the time limit). The leader HHIT-SAT is raw 12.52, so seven balanced craft would win outright.

Every generator we have produces a decreasing profile instead:

* **Greedy cover** (`tools/greedy_cover.py`, m0 1600): 46, 41, 37, 33, 29, 25, 20 targets (tanks 1115...731),
  231 covered, sum cost 8.60 -- and the **67 left over chain into a 6-flyby route**. With `--vinf 4.0` and a
  diverse-candidate settle the same shape is cheaper (46@989, 43@954, 38@954, 33@859) but no flatter.
* **The leftovers are structural.** They are the far, inclined objects: median a 2.32 vs 1.78 AU, aphelion 3.81 vs
  2.75, i 16.4 vs 11.6 deg, e 0.631 vs 0.566. They can only be covered *woven into* a deep route. A uniform 1 J
  prize never buys them because every route is already fuel- AND time-capped (947 of 960 kg, 14.6 of 15 yr), so the
  beam spends its budget on the cheapest available targets.
* **Merging an existing fleet does not balance it.** `ifleet_t10d` is a genuine disjoint partition of all 298 into
  10 routes; merging whole routes down to 7 gives groups of 27-56, and a merged 55-target group is internally
  incompatible -- its beam returned 28 flybys, against 37 for a coherent 37-target group.
* **Balanced grouping helps but does not close.** Re-beaming coherent groups with the new levers beats t10d route
  for route (37 targets at 854/883/895 kg against t10d's 919), but **a beam confined to a 48-target group captures
  only ~36 of them (~75%)**. Seven routes therefore reach ~250 of 298, and each uncovered target costs a full 1 J.

That 75% capture rate is the single number standing between us and a 12-ish score.

## 5. What was built today (all tools take `--m0 --vinf --dvmax --tofmax --lintofmax`)
| tool | what it does |
|---|---|
| `tools/greedy_cover.py` | beam over the not-yet-covered targets, prize 1 J each, settles a diverse top-k and keeps the best twin; `--prizes`, `--max-depth`, `--scan/--try-best`, `--settle-timeout/--settle-budget` |
| `tools/dual_greedy.py` | subgradient loop around it: multiply the prize of every uncovered target by ~1.8 and re-run (Lagrangian dual of the set cover), each pass scored with the true objective |
| `tools/kroutes.py` | Lloyd's algorithm with beam routes as centroids: assign every target to its nearest route under a capacity cap (~298/N), re-beam each group with a pad, iterate. Round 0 assignment is balanced (48/48/43/42/41/39/37) |
| `tools/partition_ls.py` | re-plan one route against the targets nobody covers, prize 1 J on anything no other route holds, accept on the true twin dJ |
| `tools/partition_mip.py` | pick N columns from twin/planner pools (`--twin-files`, `--plan-files`, several `--N` at once) |
| `tools/neighbourhood_beams.py` | plain/prize beams on 50-100-target neighbourhoods around seed routes |
| `tools/grow_skeletons.py` | now has `--max-dtank` (cap a single insertion) |
| `ctoc14/impulsive.py` | `settle_tour(..., screen_tank=2000)` |
| `results/newgen/pool16.jsonl` | ~2,250 columns at m0 1600 (1,450 with >= 34 targets) for the MIP |

Useful negative: the MIP over the *old* (mass-capped) pool covers only 206/298 with 8 routes and 235/298 with 10 --
that pool cannot express the fleet we want and should not be used again.

## 6. What to try next, in order
1. **Finish the k-routes loop** (`tools/kroutes.py`, now ~86 s per group instead of 25 min). Round 0 only re-beams
   the seed assignment; the point is rounds 1+, where targets missed by their group are reassigned to whichever new
   route passes closest. Run N = 7 and N = 8, 5 rounds, and watch whether capture climbs from 75%.
2. **If capture stalls, accept a hybrid:** k-routes for the backbone, then `tools/grow_skeletons.py --max-dtank 60`
   and `run_ialns.py search --lns` to weave the leftovers in -- the historical path from ~30-flyby beams to full
   coverage, but now starting from 36-46-flyby routes instead of 28-34.
3. **The safe fallback is worth taking regardless:** rerun the *existing* 10-craft pipeline with `--m0 1600 --vinf 4.0`.
   Every route today came out 6-8% cheaper at equal depth, which is ~0.5-0.9 J on the fleet for a few hours of
   machine time and no new algorithm.
4. **Do not** rebuild pools at m0 1000, rerun colgen root pricing, the waypoint dive, or uncapped growth.

## 7. Operating notes
- One heavy job at a time: other applications hold ~7.7 GB of swap, and two 8-worker jobs drove the load to 25.
- Every settle loop needs a per-candidate timeout (150 s) and a per-route wall-clock budget; without them one
  degenerate tour costs 25 minutes.
- A failed settle must never kill a run (`LinAlgError` from the linearisation is routine); all four drivers now
  catch it and drop the candidate.
