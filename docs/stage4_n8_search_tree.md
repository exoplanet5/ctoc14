# Stage 4 instructions: the N = 8 search tree, target raw J < 13
Written 2026-09-18 evening. Instructions only -- nothing here has been run. Builds on `docs/stage3_summary.md`
(read its sections 1-4 first); `docs/stage2_instructions.md` sections 3-5 (grow / polish / export / pitfalls) still stand.
Fallback at every leaf: `results/CTOC14_Result_t10d.txt` (raw 14.366, already submitted, banked at 13.81 shown).

## 0. What "J < 13 with eight craft" means in numbers
The validator's J = sum J_i + 2 (targets 131 and 144 are unreachable for everybody; t10d = 12.366 + 2).
So **J < 13  <=>  sum J_i < 11  <=>  sum (x + x^2) < 3.0 over the eight craft, minus 1.0 for every extra miss.**

| fleet-mean dv per flyby | balanced 8 x 37.25 | mild taper 46,43,40,38,36,34,32,29 | steep taper 46,44,42,40,38,34,30,24 | mean tank |
|---|---|---|---|---|
| 0.40 km/s | 11.92 | 11.94 | 11.95 | 879 kg |
| 0.45 | 12.27 | 12.29 | 12.31 | 922 |
| 0.49 (= t10d's efficiency) | 12.57 | 12.60 | 12.63 | 958 |
| **0.53 (the J = 13 line)** | 12.89 | 12.94 | 12.97 | 995 |

Three consequences that shape the tree:
1. **Balance is worth only 0.03-0.06 J.** Stage 3 section 4 overstated it. A tapered profile is fine as long as the
   eight routes SUM to 298; do not spend machine time flattening a fleet. What killed the greedy cover was not its
   shape but its *tail*: 46+41+37+33+29+25+20 = 231, then 67 unchainable orphans.
2. **Coverage is everything.** One extra miss moves the J = 13 line from 0.53 to 0.41 km/s per flyby; two make it
   unreachable. The plan must end at 298 covered (at most 297).
3. **There is a real budget for weaving the hard targets in.** Routes built with the new levers cost 0.37-0.44 km/s
   per flyby (37 @ 854/883/895 kg, 36 @ 868, 46 @ 1004). Against the 995 kg line that leaves ~115 kg per craft,
   ~900 kg fleet-wide. Insertions cost 40-60 kg (easy target) to 100-280 kg (hard target, `grow_skel8.out`), so the
   budget buys **10-16 insertions, no more**. Hence the rule that drives everything below:

> **The backbone (8 beam-planned routes, disjoint) must cover >= 282, and the <= 16 it leaves must be EASY targets.**
> Hard targets (far, inclined: stage 3 section 4) go into the backbone by prize; easy ones are what insertion is for.
> Every generator so far did the opposite: it spent the beam on easy targets and left the hard ones for insertion.

## 1. The tree at a glance
```
R0  prep: hardness prizes H, settings, seeds                              (30 min, no heavy compute)
 |
 +-- A  fleet tree search: hard-first, depth-capped, beam over partial fleets   <- main line (3-4 h)
 |      gate G1: >= 282 disjoint, leftovers <= 16 and <= 4 of them hard
 |        +-- pass  -> L
 |        +-- 270-281 -> A' (one re-run with duals raised on the leftovers) -> else C
 |        +-- < 270 -> B
 +-- B  k-routes with dual prizes (Lloyd + subgradient), seeded by A's best fleet      (2-3 h)
 |      gate G1 again;   fail -> C
 +-- C  MIP over everything generated (pool16 + A + B columns), N = 8, then G1          (20 min)
 |      fail -> F9
 |
 L   leftover repair, cheapest operator first                                          (2-3 h)
 |    L1 partition_ls (re-plan one route, ejection allowed)  -> L2 capped insertion -> L3 block LNS
 |    gate G2: covered >= 297, impulsive sum J_i <= 10.90
 |        +-- pass -> P          +-- 10.90-11.3 -> P anyway, then decide on the submit rule
 |        +-- covered <= 295 -> back to A' / B with the stuck targets' prizes doubled (once), else F9
 P   polish (relocate / retime / swap), export, validate
 |    gate G3: validator2 VALID, raw J < 13.0   -> SUBMIT (any raw J < break-even of the day is worth submitting)
 |
 F9  same tree at N = 9 (J ~ 13.3-13.6: still beats t10d's break-even until 09-27)
 F10 last resort: existing 10-craft pipeline rerun with --m0 1600 --vinf 4.0 (stage 3 section 6.3, ~13.5-13.9)
```
Depth-first: run A to its gate before touching B or C. One heavy job at a time (8 workers).

## 2. R0 -- preparation (shared by all branches)
Fixed settings for every beam: `--m0 1600 --vinf 4.0 --dvmax 1.2 --tofmax 400 --lintofmax 600 --beam 300`,
launch grid `0,800,20`, `--settle-timeout 150`. Do not re-test any of these (stage 3 sections 1-2).

**R0.1 Hardness prizes** -> `results/n8/prizes_H.json` (300 floats). DONE via `results/newgen/scratch/prizes_H.py`:
- `f_t` = share of the >= 34-target columns in `results/newgen/pool16.jsonl` that contain target t;
- `o_t` = 1 if t is in the orphan list of `results/newgen/gc16/log.txt` (the 67), else 0;
- `prize_t = clip(1 + cs * (1 - f_t / median f)+ + 0.6 * o_t, 1.0, 3.0)`; unreachable = 0.

**Coefficient chosen `cs = 0.45` (not the 1.2 first written here).** The deep-column frequency is heavily skewed
(p25(f) = 0; **80 reachable targets never appear in any deep column**; median f = 0.034), so `cs = 1.2` flags 156
targets as hard (prize >= 1.5) -- and prune (ii) can only tolerate ~87, so branch A would prune every root child.
`cs = 0.45` gives: 67 hard (= the structural orphans, <= the budget), the 80-target f=0 tail lifted to ~1.45
(preferred, not hard-flagged), median 1.28, max 2.05, 129 targets at exactly 1.0. "Hard" = prize >= 1.5.

**Measured prize-depth tradeoff** (depth-1 beam over all 298, beam 300, m0 1600, vinf 4.0): full hardness -> route 1
is **26 flybys / 14-16 hard** (dv/flyby 0.59, near the cap); uniform 1 J -> **46 flybys / only 4 hard**. Prizes buy
hard density at the cost of depth. This is why branch A has no depth floor and uses `--prize-gamma` (see section 3).

**R0.2 Tool changes** (small, do them all before the first run):
| file | change |
|---|---|
| `tools/greedy_cover.py` | factor the body of the route loop into `plan_route(avail, prize, a) -> list of (ip, tank, targets)` returning the top `--try-best` *settled, diverse* candidates instead of only the best. No behaviour change for the CLI. |
| `tools/fleet_tree.py` (new, ~150 lines) | the tree search of section 3, built on `plan_route` |
| `tools/kroutes.py` | `--prizes file` (own-group targets use the file's prize), `--pad-prize 0.3` (padded outsiders are worth little: another route owns them), `--dual-up 1.6 --dual-cap 4` (after each round multiply the prize of every uncovered target), and make the assignment step hand **uncovered targets first** to their nearest route before the capacity fill |
| `tools/partition_ls.py` | `--prizes file`; default `--m0 1600 --vinf 4.0` |
All four must keep the try/except + SIGALRM guard around `settle_tour` and append every kept beam state to
`--columns results/n8/pool.jsonl` (branch C lives off that file).

## 3. Branch A -- fleet tree search (hard-first, depth-capped)
A beam search whose *states are partial fleets*. This is the "scheme search tree" proper.

- **Node at depth d** = d settled routes (disjoint target sets) + the remaining set R.
- **Expansion:** one prize beam over R with prizes H, `--max-depth D_d`, then `plan_route` settles up to `b = 3`
  diverse candidates (overlap < 80% between siblings -- tighter than the 92% used inside one route, the point is to
  branch on genuinely different regions). Each candidate is a child. ~3-4 min per expansion.
- **Depth cap D_d** (saves easy connectors for later routes -- the uncapped greedy burnt them on route 1):
  `D_d = min(44, ceil(|R| / (8 - d)) + 4)`  -> 42, 41, 41, ... and free (44) for the last two routes.
- **Node value** (lower is better), all in J units:
  `V = sum J_i(done) + n_left * [1 + c(nbar)] + forced_miss + 0.35 * sum_{t in R} (H_t - 1)`,
  `n_left = 8 - d`, `nbar = min(|R|, 44*n_left)/n_left`, `forced_miss = max(0, |R| - 44*n_left)` (targets no
  remaining route can reach), `c(n)` = x + x^2 at tank 601.5*exp(0.45 n / 39.2266). Capping `nbar` at 44 stops the
  estimate exploding while the tail is large. The penalty uses the FULL (un-softened) hardness `H`, so the tree keeps
  preferring fleets that have absorbed the hard targets even when the beam runs at a softened `--prize-gamma`. A
  TERMINAL node (8 routes, or R empty) is scored on the truth, `sum J_i(done) + |R|` (each leftover a 1 J miss).
- **Beam width w = 3** partial fleets per depth; **prune** a child if
  (i) it can no longer reach the coverage floor: `covered + 44*n_left < 282` (`--cover-goal`),
  (ii) it can no longer bring the hard leftovers to <= 4: `hard_left > 18*n_left + 4`,
  (iii) the route's tank > 1150 kg or per-flyby dv > 0.65 km/s (`--max-tank`, `--dvpf-max`; hard routes measure ~0.59-0.60).
  **There is deliberately NO per-route depth floor** (the earlier `covered < 34d` rule): prizes make the FIRST routes
  shallow and hard-dense (26 flybys, ~15 hard, measured) and the LATER routes deep and easy (~44) once the hard
  targets are gone, so a 26-flyby route 1 is correct, not a failure. A depth floor would kill branch A at the root.
- **Cost:** <= 3 x 3 expansions per depth x 8 depths = 72 route builds ~ 3.5-4 h. Cache by frozenset(R) (siblings
  often meet again). Save every node (`results/n8/A/d{d}_{k}/` as an IFleet) so a crash or a manual stop loses nothing,
  and so A' and B can start from any node.
- **Dead-end rule:** from depth 5 on, if the best node can no longer reach the coverage floor with margin
  (`covered + 44*n_left < 282 + 6`), stop -- the tail will not close; go to B seeded with that node.

Run (after R0.2):
`~/.venvs/astro313/bin/python tools/fleet_tree.py results/n8/A --n 8 --width 3 --branch 3 --prizes results/n8/prizes_H.json --prize-gamma 1.0 --cover-goal 282 --m0 1600 --vinf 4.0 --beam 300 --nproc 8 --scan 30 --try-best 3 --settle-timeout 150 --columns results/n8/pool.jsonl`

`--prize-gamma` softens the beam's prize gradient to `1 + gamma*(H-1)` (1.0 = full hardness). With the coverage-
feasibility prune, full hardness self-balances (shallow-hard early, deep-easy late), so start at 1.0; if branch A
under-covers (final < 282), lower to 0.6-0.7 to trade a little hard density for depth. See `results/n8/R0_notes.md`.

**Gate G1** (read `results/n8/A/result.json`): disjoint coverage >= 282, leftovers <= 16, hard leftovers <= 4,
every tank <= 1100, sum J_i + misses printed. **A'** (only if 270-281): multiply the prize of each leftover by 1.6
(cap 4), restart from the best depth-3 node rather than from scratch (~2 h).

## 4. Branch B -- k-routes with dual prizes
Use when A's tail does not close. Seed = A's best complete or depth-8 fleet (else `results/newgen/ifleet_t10d`,
which takes the 8 longest of its 10 routes).

`~/.venvs/astro313/bin/python tools/kroutes.py <seed> results/n8/B --n 8 --rounds 6 --cap 1.15 --pad 12 --pad-prize 0.3 --prizes results/n8/prizes_H.json --dual-up 1.6 --dual-cap 4 --m0 1600 --vinf 4.0 --beam 300 --nproc 8 --scan 24 --try-best 3 --columns results/n8/pool.jsonl`

~8 groups x 3 min = 25 min per round. Watch two numbers per round: distinct coverage (must climb: 250 -> 282) and
duplicates (pad targets flown twice; should fall once `--pad-prize` is in). **Stop early** if coverage does not rise
by >= 5 over two consecutive rounds. Then gate G1.

## 5. Branch C -- MIP recombination
Everything A and B generated is in `results/n8/pool.jsonl` (planner costs), plus `results/newgen/pool16.jsonl`.
`~/.venvs/astro313/bin/python tools/partition_mip.py results/n8/C --N 8,9 --twin-files "" --plan-files results/n8/pool.jsonl,results/newgen/pool16.jsonl --w-hard 1.0 --w-easy 1.0 --overlap-pen 0.05 --max-tank 1150 --min-targets 28 --time 900 --save-states`
(w = 1.0 on both classes: a miss costs 1 J whatever the target; `--hard-thr` only decides where overlap is forbidden.)
Twin-recost the 8 chosen columns (`results/newgen/scratch/recost_pool.py`), save as an IFleet, apply gate G1.
If the MIP's N = 8 coverage is < 270 the generated columns are still mutually incompatible -> F9.
Never feed it `recost35.jsonl` or `results/colgen/*` (mass-capped, stage 3 section 5).

## 6. L -- leftover repair (input: a G1 fleet, <= 16 leftovers)
Cheapest and most global operator first; re-measure coverage and sum J_i after each.
1. **L1 re-plan with ejection:** `tools/partition_ls.py <fleet> results/n8/L1 --rounds 6 --knn 30 --m0 1600 --vinf 4.0 --prizes results/n8/prizes_H.json --nproc 8`.
   It accepts on the true twin dJ with 1 J per newly covered target, so a route may drop two easy targets to pick up
   one hard one *only if* those easy ones are then placed by L2 -- run L1 and L2 alternately, twice.
2. **L2 capped insertion:** `tools/grow_skeletons.py results/n8/L1 results/n8/L2 --hard <comma-separated leftover ids> --max-dtank 60 --max-tank 1100 --nproc 8`,
   then once more with `--max-dtank 130` for whatever remains (a 130 kg insertion at x ~ 0.27 costs ~0.15 J; a miss
   costs 1.0 -- it is still a bargain, the cap only exists to stop compounding on one route).
3. **L3 block LNS:** `tools/run_ialns.py search results/n8/L2 results/n8/L3 --lns --bmax 5 --hosts 3 --lns-rounds 40 --nproc 8`.

**Gate G2:** covered >= 297 and impulsive sum J_i <= 10.90 (leaves 0.10 for the export loss of 0.02-0.06 and polish
noise). Stuck targets after L3: double their prize and go back to A'/B **once**; a second failure -> F9.

## 7. P -- polish, export, submit
1. `tools/run_ialns.py search results/n8/L3 results/n8/P --relocate --retime --swap --nproc 8` (expect -0.05 to -0.15 J).
2. Copy `results/newgen/export_t10d.sh` to `results/n8/export_n8.sh`, change the three paths, run it
   (`run_ialns.py export` -> `combine_submission.py` -> `validator2.py`). Heavy routes (> 1000 kg) need the 0.05 /
   0.025 d polish retries; an INVALID fragment stops the script.
3. **Gate G3 / submit rule:** validator2 VALID, and raw J below the day's break-even against the banked t10d
   (14.26 on 09-19, 14.16 on 09-21, 14.00 on 09-24, 13.85 on 09-27). J < 13 is the *goal*; anything VALID under the
   break-even is still worth submitting, then keep improving.

## 8. Schedule and stop rules
| day | work | decision by end of day |
|---|---|---|
| 09-19 | R0 + tool changes (morning), branch A (afternoon/evening) | G1 on A |
| 09-20 | A' or B, then C if needed | a G1 fleet exists, or switch to F9 |
| 09-21 | L1-L3 | G2 |
| 09-22 | P + export + validate + submit | first N = 8 submission |
| 09-23 to 09-26 | second pass through the tree from the submitted fleet (B seeded by it, duals on its expensive targets) | resubmit if better by more than the k drift (0.05 J per day) |
| 09-27 | final export / validation only, no new search | |
Hard stop: if no G1 fleet exists by the end of 09-21, run F9 (same commands with `--n 9`, G1 >= 286, G2 sum J_i <=
11.5) and in parallel time F10 on 09-22; whichever validates under the break-even gets submitted.

## 9. Bookkeeping
- Everything under `results/n8/`; one log per node; every branch appends to `results/n8/pool.jsonl`.
- After each gate, append one line to `results/n8/tree_log.md`: node, coverage, leftovers (hard/easy), sum J_i,
  mean km/s per flyby, wall time, decision. That file is the tree as actually searched.
- Do not: rebuild pools at m0 1000, raise TOF caps, re-tune v_inf, run two heavy jobs at once, or let insertion
  run without `--max-dtank`.
