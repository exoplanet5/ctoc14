# Stage 2 instructions: eight craft by neighbourhood beams (written 2026-09-17 22:45)

> **Sections 1-2 are SUPERSEDED (2026-09-18).** They were designed around a generator that could not build a craft
> heavier than 805 kg (the beam ran at m0 = 1000), so no pool built that way can express a 37-flyby route.
> See `docs/stage3_summary.md`; sections 3-5 (grow, polish, export, budget, pitfalls) still stand.

Goal: an 8-craft fleet with raw J <= 13.3 (shown <= 13.2 by 09-25, rank 4-6), from a partition of the 298 targets
into eight 37-target sequences. Keep `results/CTOC14_Result_t10d.txt` (J 14.366, ten craft) as the fallback and
**submit it now if not done** (shown 13.81 today, no downside: the board keeps the best score).

## 0. What is established (do not re-test)
| fact | where |
|---|---|
| Leaderboard N = craft + misses; the top five are 8 craft at 0.48-0.54 km/s per flyby, i.e. our efficiency | strategy_stage2.md section 1 |
| Planner-born 36-38-target routes cost 0.32-0.36 km/s per flyby in the twin (tank ~820 kg) | campaign doc 9.1, `results/newgen/recost35.jsonl` |
| The planner's cost at s = 0.70, reserve 2 equals the twin cost to 2 % | 9.1 |
| Hard targets (87 = in <= 5 dense columns) are cheap in Delta-v when a sequence is built for them (24-31-target prize columns with 9-17 hard: 780-800 kg) | 9.4 |
| Column generation (any craft price / prize scale) converges to 10.2 fractional craft; the dive fixes easy columns first | 9.4 |
| Insertion-based growth is cheap for ~5 targets per route, then 100-280 kg each (pinning) -> 10-craft equilibrium | 9.5, `results/newgen/ifleet_skel8.log` |
| A PLAIN beam confined to a subset of >= 46 compatible targets returns 33-36-flyby routes with 9 hard targets each; subsets of <= 37 return 6-23 | `results/newgen/pbeam8.out` |
| Twin construction from a planner tour: junction nodes 1 d after departure, retry ladder 1 / 0.25 / 0.05 d (`impulsive.from_tour`, `settle_tour`); 150 s per-column limit in the driver | 9.2, 9.3 |

The missing piece is therefore a POOL of deep, hard-rich, mutually compatible columns. Nothing else in the chain
(twin costing, master, growth, polish, export) needs new code.

## 1. Build the neighbourhood-beam pool (day 1 morning, ~2-3 h machine time)
Write `tools/neighbourhood_beams.py` from `results/newgen/scratch/partition_beam.py`:
1. Seeds = every route we trust: the 8 skeletons (`results/newgen/skel8`, `skel8b`), the 10 flown routes
   (`results/newgen/ifleet_t10d`), the prize-mode twin states with >= 9 hard targets (`results/colgen/n8a/twin_states`,
   select via `twin_run1.jsonl` + the hard set), and the two good partition routes (`results/newgen/pbeam8/tour_01/02.json`).
2. For each seed, the neighbourhood = the seed's own targets + the K nearest other targets by closest approach of the
   asteroid to the seed trajectory (`run_ialns.candidates`, dmax 1.5 AU, per_target 1, distance-ranked), with
   K chosen so that |subset| = 50, 60 and 70 (three subsets per seed). Neighbourhoods overlap: that is intended.
3. Run the plain beam inside each subset (exclude everything else): `Params(beam=400, w_fuel=1.0, m_margin=40,
   dv_max=1.2, lin_tofs=20..600 d, lin_drmax=0.15 AU, vinf_cap=2.0, tof_refine=True)`, m0 1000, launch grid 0-600 d
   step 20 (also try 0-1500 d for late seeds). Set `P.collect = []` so every state is kept, keep states with >= 28
   flybys, write them as `tour_json` records with a tag (`seed:subset_size`) to `results/newgen/nb_pool/columns.jsonl`.
   One beam = 10-70 s on 8 processes; ~30 seeds x 3 sizes = ~1 h.
4. Twin-recost every column with >= 30 targets: `results/newgen/scratch/recost_pool.py` generalised to take
   `--columns results/newgen/nb_pool/columns.jsonl` (it already writes the twin state inline). Expect ~5-15 s each,
   1 % failures. Run it ALONE (8 workers) -- the machine swaps when two heavy jobs run.
5. Report: number of columns with >= 33 targets and >= 8 hard targets, their tank distribution, and how many of the
   87 hard targets appear in at least one such column. **Gate A: >= 80 of 87 hard targets appear in deep columns and
   >= 200 such columns exist; otherwise widen K and add seeds (prize-mode columns, dive fixes) before going on.**

## 2. Partition (day 1 afternoon, minutes)
`results/newgen/scratch/skel_mip.py` with `--twin-files` pointing at the new recost file (add the option; it currently
reads `recost35.jsonl` + `n8a/twin*.jsonl`): N = 8, `--w-hard 0.6 --w-easy 0.03 --overlap-pen 0.05 --time 600`.
Read `skeletons.json`: the number of distinct targets covered and the uncovered hard list. **Gate B: >= 250 distinct
targets covered by the 8 columns with <= 6 uncovered hard.** If the MIP covers fewer, the pool lacks compatible deep
columns: go back to step 1 with seeds = the MIP's own 8 columns (their neighbourhoods are the right ones) and iterate
once or twice (this is column generation by hand, each loop ~1 h).

## 3. Grow, polish, export (day 2)
1. `tools/grow_skeletons.py skel_dir ifleet_dir --hard <uncovered hard ids> --nproc 8` with a NEW option
   `--max-dtank 60` (accept an insertion only if the tank grows <= 60 kg; today's run without it drove tanks to 1200 kg).
   Expect 30-45 leftovers; each leftover is either placed by a block move later or is a miss (1 J) -- if more than ~15
   remain, the partition (step 2) was not good enough.
2. `tools/run_ialns.py search ifleet_dir out_dir --nproc 8` (relocate / retime / swap passes; add `--lns` with
   `--victims` = none, it is the block-insertion operator that places leftovers into the routes with ejection).
3. `tools/run_ialns.py export` -> `combine_submission.py` -> `validator2.py` (see `results/newgen/export_t10d.sh`).
   Exact export costs +0.02-0.06 J over the impulsive J; heavy routes need the 0.05 / 0.025 d polish retries.
4. Submit if shown J < 13.81 x (0.9615 / k_now); the formula and table are in strategy_stage2.md section 1.

## 4. Budget and stop rule
- Deadline 2026-09-28 12:00; last safe export/validate day 09-27; every day costs 0.05 J-equivalent.
- Stop if Gate A fails after one widening, or Gate B fails after two hand iterations: then the pool cannot express an
  8-craft partition and the remaining option is 9 craft from the same pool (worth ~0.3 J shown at best).

## 5. Pitfalls learned today
- Never run two heavy jobs (8 + 5 workers): 7.9 GB of swap belongs to other apps; load went to 267 and a 20-s column
  took 700 s. `renice 10` the secondary job or run sequentially.
- A twin recost step in a loop must have a per-column time limit (150 s) and one pass per master solve; stragglers
  otherwise cost 10-25 min per solve.
- Prize-driven beams trade depth for hard targets (24-27 deep at prize scale 0.3-0.5); use PLAIN beams on restricted
  subsets for depth, prizes only to pick among states afterwards.
- `ColumnStore.add` resets a column's cost to a cheaper planner estimate when the same target set is re-generated;
  re-apply twin costs after every pricing call (done in `run_colgen.py twin_recost`).
- Greedy insertion must be capped per insertion (`--max-dtank`); "any insertion below 1 J" is the wrong criterion
  because the tanks compound (x + x^2).
