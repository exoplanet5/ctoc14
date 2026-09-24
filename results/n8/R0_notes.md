# R0 execution notes (branch A prep) -- 2026-09-18

Status of `docs/stage4_n8_search_tree.md` section 2 (R0). **Prep only, no branch A run yet.**

## R0.1 hardness prizes -> results/n8/prizes_H.json
Script `results/newgen/scratch/prizes_H.py` (has `--cs`, `--orphan`).
- Deep-column frequency is heavily skewed: p25(f)=0; **80 reachable targets never appear in any >=34-target column**;
  median f = 0.034. The doc's coefficient cs=1.2 therefore flagged **156** hard (prize>=1.5), which auto-trips branch-A
  prune (iii) (`hard_left > 12*(8-d)`, budget 96 at the root, calibrated for ~87 hard). Every depth-1 child would prune.
- Chose **cs=0.45, orphan=0.6**: |prize>=1.5| = the 67 structural orphans (<=96, prune-safe), the 80-target f=0 tail
  lifted to ~1.45 (preferred, not hard-flagged). median 1.28, max 2.05, 129 targets at exactly 1.0.

## R0.2 tool changes (all compile; `plan_route` proven to settle real routes)
- `tools/greedy_cover.py` -- route-loop body factored into `plan_route(avail, PZ, a, P, grid, nproc, say, tag, cf)`
  returning the settled diverse top candidates `[(ip, tank, targets, value), ...]`. CLI behaviour preserved.
- `tools/fleet_tree.py` (NEW) -- branch A: partial-fleet beam, `Node` value + prune (i)-(iv) as specified, terminals
  scored on true `sumJi + misses`, `plan_route` per expansion, `result.json` with coverage/leftovers/hard-leftovers.
  Synthetic unit tests pass (hard-first ranking, prune, diverse, IFleet round-trip). Live smoke: beam->plan_route->
  settle->Node->save->result.json all execute.
- `tools/kroutes.py` -- `--prizes --pad-prize 0.3 --dual-up 1.6 --dual-cap 4`; uncovered targets assigned first each
  round; own-group at the dual-scaled prize, pad outsiders at `--pad-prize`; uncovered prizes multiplied each round.
- `tools/partition_ls.py` -- `--prizes`; defaults now `--m0 1600 --vinf 4.0`. Beam biased by prize, accept test still 1 J/target.

## CALIBRATION FINDING -> branch-A design corrected (RESOLVED)
Direct depth-1 probe over all 298 (beam 300, full grid, m0 1600, vinf 4.0):
- **hardness H: first route 26 flybys / 14-16 hard** (tank ~885, dv/fb 0.59-0.60).
- uniform 1 J: **46 flybys / only 4 hard** (tank ~1000, dv/fb 0.42-0.47).
Prizes trade total depth for hard density (the beam score is prize-dominated at w_fuel 0.52). This is not a bug to
tune away -- it is the mechanism. The doc's original prune (ii) `covered < 34*d` (a per-route depth floor) is
**incompatible** with it: prized route 1 is 26 flybys < 34, so branch A would prune every root child and die.

**Correction applied to `tools/fleet_tree.py`:**
1. Prune (ii) replaced by a COVERAGE-FEASIBILITY test `covered + 44*n_left < 282` (can the fleet still reach the
   floor?) and a hard-feasibility test `hard_left > 18*n_left + 4`. No per-route depth floor.
2. Node value: `nbar` capped at 44 and a `forced_miss = max(0, |R| - 44*n_left)` term, so the estimate does not
   explode while the tail is large; the hard-residual penalty uses the FULL (un-softened) hardness.
3. `--prize-gamma` (default 1.0) softens the beam prize `1 + gamma*(H-1)` if needed; `--cover-goal` (282).
4. `--dvpf-max` 0.60 -> **0.65**: hard-dense routes measure 0.586-0.599, so a 0.60 cap spuriously prunes viable hard
   routes (and a missed hard target costs a full 1 J).
Bugs found and fixed by the smokes: prize script `parents[2]`->`[3]`; a `sed` that deleted `--dvpf-max` with `--dead-left`.

**Self-balancing confirmed** (N=8 width-1 live probe): depth 1 = 1 route/25 covered/51 hard-left (accepted, was
impossible under the old prune); depth 2 = 2 routes/49 covered/37 hard-left. Hard targets absorbed first (67->51->37),
so later routes have only easy targets left and go deep -- exactly the intended shape. Branch A is READY.

## Branch A run (the next step -- NOT run yet; heavy, ~3-4 h, one job at a time)
`~/.venvs/astro313/bin/python tools/fleet_tree.py results/n8/A --n 8 --width 3 --branch 3 --prizes results/n8/prizes_H.json --prize-gamma 1.0 --cover-goal 282 --m0 1600 --vinf 4.0 --beam 300 --nproc 8 --scan 30 --try-best 3 --settle-timeout 150 --columns results/n8/pool.jsonl`
Then gate G1 on `results/n8/A/result.json`. If final coverage < 282, lower `--prize-gamma` to 0.6-0.7 (more depth,
fewer hard per route) and rerun.

## width-1 probe full result (readiness signal, not the real search)
N=8, width 1, branch 1, beam 200, scan 12, try-best 2, gamma 1.0: advanced depths 1-3 (covered 25/49/70,
hard-left 67->51->37->27, routes 25/24/21 flybys), then banked at depth 4 -- route 4's two settle attempts both
failed (undersampled; the real run uses width 3 / scan 30 / try-best 3, which has many more attempts AND two other
partial fleets to fall back on). Logic and self-balancing are sound; coverage is the open question.

**Coverage caveat for the real run:** at gamma 1.0 the routes are shallow (21-25 flybys). If that held across all 8
they would sum to ~200, far under 282 -- but (a) once the 67 hard are absorbed (~depth 5) the remaining easy routes
go deep (~44), and (b) width-3's node value rewards coverage, so it should pick deeper lines than this greedy probe.
Watch `result.json` coverage: **if < 282, rerun at `--prize-gamma 0.6`** (coverage is the unrecoverable gate -- every
miss is a full 1 J -- while a few hard leftovers are repairable in branch L, so err toward depth/coverage).
