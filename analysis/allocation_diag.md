# Allocation diagnosis: where the 30.789 fleet wastes cost, and what the route pools can and cannot cover

Agent "diag", 2026-09-15. Data in `results/diag/`, scripts in the session scratchpad (`diag/team_alloc.py`,
`team_greedy_delete.py`, `pool_coverage.py`, `lp_prices.py`). Nothing outside `results/diag/` and this file was written.

## 0. Summary

1. **The 30.789 fleet is two fleets glued together.** SC5 and SC9–16 (nine patient craft from one greedy-exclusion
   campaign) are mutually disjoint and cover 243 targets for ΣJ_i = 14.61 (0.060 per target). The other 51 covered
   targets cost 10.18 in seven craft (SC6–8 from a second patient fleet, SC1–4 = old chaser fragments) plus 4 avoidable
   misses: **the last 55 reachable targets cost 16.2 of the 30.8**. 382 flyby events cover 294 targets: 88 duplicate slots,
   81 of them on the four chaser fragments SC1–4 (89 flybys, only 8 unique targets, ΣJ_i 5.67 = 0.71 per unique target).
2. **The marginal-cost criterion K > c deletes nothing** (every craft has K ≥ 2 > c ≤ 1.63) — the criterion is
   satisfied with margins of only 0.53–0.65 by SC1–4. Deleting all four together loses 8 targets and raises J to 33.1.
   The waste is not "a craft that should be deleted"; it is "51 targets that should have been inside the nine core tours".
3. **The pool is not missing any target.** All 298 reachable targets appear in the union pool (335k routes), every
   one of them inside some route of ≥ 28 flybys. Rarity is graded, not binary: the 40 rarest targets appear in 4–8k
   routes each (median target 14k). They are large-a / long-period (a > 2.3 AU, P 3–7 yr, Spearman(#routes, a) = −0.37)
   and high-inclination (i > 30°) objects; ten of them (3, 15, 64, 77, 96, 148, 193, 216, 251, 287, all a 2.2–3.7 AU with
   q ≈ 1) are **never reached inside the first 3 legs of any route** — they can only be caught by a craft that is
   already on a shaped orbit when the asteroid comes to perihelion.
4. **The 10-craft cover is blocked by combination, not absence.** Column generation converges (no negative
   reduced cost left), so the LP bound is exact for this pool: **LP(10 craft) = 39.50**, LP(11) = 25.18, LP(12) ≈ 21.1.
   The 10-craft LP already leaves ~22 targets fractionally uncovered: the pool's long routes (28–38 flybys) all draw on
   the same ~200 "easy" ring targets, and every route that contains a rare target is a specialised short/expensive
   route. _(Integer results: see section 3.)_

## 1. Per-craft allocation of `results/CTOC14_Result_TEAM.txt` (J = 30.789)

Source: `results/diag/team_alloc.json`, `team_greedy_delete.json`. "unique" = no other craft in the file flies by that
asteroid; "K−c" = unique − J_i (marginal value; the craft is worth keeping iff K > c, brief §41).

| SC | m0 | fuel used | J_i | flybys | unique | dup | K−c | kg/flyby | end (yr) | family |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 1 | 985 | 366 | 1.351 | 3 | 2 | 1 | +0.65 | 122 | 4.7 | chaser fragment |
| 2 | 1018 | 376 | 1.388 | 25 | 2 | 23 | +0.61 | 15.1 | 14.4 | chaser fragment |
| 3 | 1086 | 378 | 1.468 | 31 | 2 | 29 | +0.53 | 12.2 | 14.9 | chaser fragment |
| 4 | 1084 | 426 | 1.466 | 30 | 2 | 28 | +0.53 | 14.2 | 14.5 | chaser fragment |
| 6 | 1160 | 482 | 1.560 | 18 | 10 | 8 | +8.4 | 26.8 | 14.8 | patient fleet B |
| 7 | 1088 | 332 | 1.470 | 16 | 9 | 7 | +7.5 | 20.8 | 13.8 | patient fleet B |
| 8 | 1094 | 375 | 1.477 | 16 | 10 | 6 | +8.5 | 23.4 | 13.2 | patient fleet B |
| 5 | 1211 | 520 | 1.627 | 32 | 21 | 11 | +19.4 | 16.3 | 14.9 | patient core |
| 9 | 1215 | 475 | 1.632 | 31 | 16 | 15 | +14.4 | 15.3 | 13.7 | patient core |
| 10 | 1211 | 500 | 1.627 | 30 | 21 | 9 | +19.4 | 16.7 | 14.6 | patient core |
| 11 | 1207 | 533 | 1.621 | 28 | 25 | 3 | +23.4 | 19.0 | 14.8 | patient core |
| 12 | 1217 | 512 | 1.634 | 27 | 17 | 10 | +15.4 | 18.9 | 14.5 | patient core |
| 13 | 1210 | 478 | 1.626 | 26 | 20 | 6 | +18.4 | 18.4 | 15.0 | patient core |
| 14 | 1211 | 509 | 1.627 | 24 | 20 | 4 | +18.4 | 21.2 | 14.9 | patient core |
| 15 | 1196 | 578 | 1.607 | 23 | 21 | 2 | +19.4 | 25.1 | 14.8 | patient core |
| 16 | 1197 | 461 | 1.608 | 22 | 16 | 6 | +14.4 | 20.9 | 14.6 | patient core |

- 382 flyby events, 294 distinct targets → **88 wasted (duplicate) flyby slots**; 80 targets are covered more than
  once (73 twice, 6 three times, 1 four times).
- Pairwise overlap matrix (`team_greedy_delete.json`): the nine core craft {5, 9, 10, …, 16} have **zero** mutual
  overlap (they are one greedy-exclusion fleet). SC6–8 overlap each other by 3/1/3 targets. All other duplicates are the
  chaser fragments SC1–4 re-flying core targets (SC3 shares 8 with SC9, 5 with SC10; SC4 shares 6 with SC12, 5 with SC5).
- Marginal-cost criterion: sequential deletion of the craft with the smallest K−c stops immediately — SC3 has K = 2
  uncovered targets vs c = 1.468, so removing it raises J by 0.53. **No craft is deletable in isolation.** Removing
  SC1–4 together loses 8 targets {23, 65, 75, 86, 98, 133, 232, 277} (J → 33.12).
- Cost structure: core 9 craft = 243 targets / 14.61 (0.060 per target, 12–27 kg per flyby but with 1.3× tank
  margin, x ≈ 0.43); the tail = 51 targets / 10.18 + 4 misses (8, 111, 177, 248) = 16.2 for 55 targets (0.29 each).
  A leader at raw 13.6 with 9 craft pays 0.046 per target across *all* 298.
- Tank sizing waste: the core craft carry 1196–1217 kg tanks but used 461–578 kg; x = 0.43 → J_i 1.61–1.63 whereas
  the fuel actually burned would justify m0 ≈ 1080–1200 (J_i 1.46–1.61). Right-sizing after conversion would save
  ≈ 0.05–0.15 per craft (≈ 1 total) — real but second-order next to the tail cost.

## 2. Route-pool coverage (335k routes, 13 files)

Source: `results/diag/pool_coverage.{json,csv}`, `pool_rarity_summary.json`, `pool_coverage.log`.
Files: camp1–5 `pool.jsonl` (60.9k / 62.3k / 62.6k / 40.1k / 50.6k routes) + `results/phasing/pool/*.jsonl` (58.6k).
Every line is a beam state, so a route's prefixes are also lines; counts below are therefore "beam states containing
the target", inflated for early-leg targets. 91k routes have ≥ 20 flybys; longest 38 (camp5, ref, win400).

Per-target statistics: median target is in 14.4k routes, 6.4k of them with ≥ 20 flybys, 740 within the first 3 legs.
**Absent from every route: only 131 and 144.** Present but never in a ≥ 20-flyby route: none. Never in a ≥ 28-flyby
route: none. Min cost-per-target of a route containing each target: 0.038–0.044 for *every* reachable target (because
every target sits in some 35–38-flyby route with a ~1090 kg tank).

### 40 rarest targets (by number of routes containing them)

| id | routes | ≥20 | ≥28 | first-3 | a (AU) | e | i (°) | q | Q | P (yr) | U (Öpik) | why rare |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 261 | 4269 | 2368 | 1052 | 57 | 1.260 | 0.264 | 23.3 | 0.93 | 1.59 | 1.41 | 0.47 | Earth-like orbit (a 1.26, e 0.26): slow synodic drift, rarely phased |
| 53 | 4890 | 1690 | 196 | 14 | 2.315 | 0.605 | 29.5 | 0.91 | 3.72 | 3.52 | 0.68 | high-i + large-a |
| 123 | 4920 | 2192 | 180 | 96 | 2.085 | 0.593 | 4.5 | 0.85 | 3.32 | 3.01 | 0.45 | large-a, P ≈ 3 yr (few perihelia) |
| 187 | 5198 | 2827 | 345 | 154 | 2.637 | 0.607 | 13.4 | 1.04 | 4.24 | 4.28 | 0.33 | large-a, q > 1 |
| 118 | 5260 | 4290 | 2116 | 30 | 1.490 | 0.367 | 68.4 | 0.94 | 2.04 | 1.82 | 1.22 | i = 68° |
| 159 | 5390 | 3491 | 509 | 3 | 2.691 | 0.677 | 5.4 | 0.87 | 4.51 | 4.41 | 0.47 | large-a, P 4.4 yr |
| 236 | 5418 | 1432 | 82 | 23 | 2.482 | 0.610 | 75.7 | 0.97 | 4.00 | 3.91 | 1.41 | i = 76° + large-a |
| 15 | 5730 | 1665 | 144 | 0 | 2.724 | 0.642 | 13.8 | 0.98 | 4.47 | 4.50 | 0.42 | large-a, q ≈ 1, never early |
| 208 | 5826 | 3308 | 897 | 144 | 1.643 | 0.383 | 3.4 | 1.01 | 2.27 | 2.11 | 0.17 | q > 1 (never inside Earth's orbit) |
| 157 | 5865 | 3386 | 384 | 3 | 2.747 | 0.706 | 20.7 | 0.81 | 4.69 | 4.55 | 0.66 | large-a + i 21° |
| 148 | 6229 | 2286 | 233 | 0 | 2.513 | 0.632 | 9.9 | 0.93 | 4.10 | 3.98 | 0.43 | large-a, never early |
| 200 | 6234 | 2601 | 474 | 1 | 2.108 | 0.602 | 16.1 | 0.84 | 3.38 | 3.06 | 0.55 | P 3 yr + i 16° |
| 179 | 6249 | 2259 | 436 | 107 | 1.628 | 0.776 | 34.3 | 0.37 | 2.89 | 2.08 | 1.03 | high-i + high-e |
| 22 | 6277 | 1672 | 217 | 516 | 2.122 | 0.654 | 34.9 | 0.74 | 3.51 | 3.09 | 0.85 | high-i |
| 287 | 6471 | 2371 | 522 | 0 | 2.812 | 0.645 | 9.8 | 1.00 | 4.63 | 4.71 | 0.35 | large-a, q ≈ 1, never early |
| 199 | 6478 | 2394 | 577 | 4789 | 2.323 | 0.578 | 28.1 | 0.98 | 3.67 | 3.54 | 0.61 | high-i + large-a (but often early: reachable near launch only) |
| 88 | 6481 | 1768 | 379 | 922 | 2.600 | 0.634 | 17.1 | 0.95 | 4.25 | 4.19 | 0.48 | large-a |
| 50 | 6491 | 3953 | 640 | 31 | 1.866 | 0.508 | 10.6 | 0.92 | 2.81 | 2.55 | 0.39 | phasing (node 0.12 AU off ring) |
| 92 | 6501 | 3557 | 610 | 1 | 2.799 | 0.675 | 6.0 | 0.91 | 4.69 | 4.68 | 0.43 | large-a |
| 136 | 6544 | 4763 | 1528 | 279 | 1.105 | 0.076 | 24.1 | 1.02 | 1.19 | 1.16 | 0.43 | Earth-like a, i 24°: crosses ring only at nodes |
| 64 | 6597 | 1537 | 15 | 0 | 3.655 | 0.724 | 44.1 | 1.01 | 6.30 | 6.99 | 0.91 | a 3.66, i 44°, P 7 yr: 2 perihelia in 15 yr |
| 51 | 6829 | 3418 | 1060 | 62 | 2.036 | 0.493 | 10.0 | 1.03 | 3.04 | 2.91 | 0.25 | q > 1 |
| 160 | 6887 | 1733 | 725 | 40 | 2.320 | 0.561 | 1.9 | 1.02 | 3.62 | 3.53 | 0.22 | q > 1, node 0.47 AU off ring |
| 251 | 6962 | 3383 | 275 | 0 | 3.069 | 0.695 | 16.4 | 0.94 | 5.20 | 5.38 | 0.51 | a 3.07, P 5.4 yr, never early |
| 247 | 6965 | 4851 | 1672 | 127 | 1.744 | 0.433 | 8.7 | 0.99 | 2.50 | 2.30 | 0.27 | phasing |
| 230 | 6989 | 3842 | 1014 | 102 | 1.336 | 0.603 | 8.3 | 0.53 | 2.14 | 1.54 | 0.65 | phasing (node 0.15 AU off) |
| 40 | 7117 | 2011 | 308 | 458 | 2.275 | 0.566 | 9.3 | 0.99 | 3.56 | 3.43 | 0.33 | large-a |
| 280 | 7121 | 2751 | 493 | 226 | 1.490 | 0.377 | 36.8 | 0.93 | 2.05 | 1.82 | 0.72 | high-i |
| 116 | 7246 | 3525 | 761 | 135 | 2.541 | 0.668 | 6.2 | 0.84 | 4.24 | 4.05 | 0.50 | large-a |
| 34 | 7519 | 5224 | 2215 | 1 | 2.627 | 0.631 | 2.8 | 0.97 | 4.28 | 4.26 | 0.33 | large-a, node 0.49 AU off ring |
| 3 | 7536 | 4327 | 739 | 0 | 2.201 | 0.859 | 17.5 | 0.31 | 4.09 | 3.26 | 1.05 | e 0.86, fast (U 1.05), never early |
| 259 | 7664 | 4341 | 1204 | 6 | 2.608 | 0.726 | 8.4 | 0.71 | 4.50 | 4.21 | 0.65 | large-a |
| 207 | 7680 | 3250 | 822 | 254 | 1.546 | 0.472 | 13.5 | 0.82 | 2.28 | 1.92 | 0.47 | phasing |
| 277 | 7769 | 2813 | 919 | 4459 | 2.058 | 0.533 | 40.4 | 0.96 | 3.16 | 2.95 | 0.82 | i 40° (early-only) |
| 77 | 7781 | 3297 | 1095 | 0 | 2.592 | 0.618 | 13.2 | 0.99 | 4.19 | 4.17 | 0.39 | large-a, q ≈ 1, never early |
| 223 | 7787 | 2103 | 277 | 895 | 2.493 | 0.825 | 10.9 | 0.44 | 4.55 | 3.94 | 0.92 | high-e + large-a |
| 39 | 7808 | 2126 | 245 | 56 | 2.005 | 0.570 | 52.4 | 0.86 | 3.15 | 2.84 | 1.04 | i 52° |
| 254 | 7821 | 3295 | 283 | 47 | 2.671 | 0.702 | 5.9 | 0.80 | 4.55 | 4.36 | 0.56 | large-a |

(131 and 144 are absent: a = 8.6 AU retrograde and a = 17.8 AU.)

What makes a target rare — rank correlation of #routes with the elements over the 298 reachable targets:
a −0.37, Q −0.34, q −0.27, P −0.37, e −0.13, i −0.12, U(Öpik) −0.03, node-gap +0.05. **Semi-major axis / period is
the dominant factor** (a > 2.3 AU means 3–7-yr periods, i.e. only 2–5 perihelion passages in the mission, and the
craft must be phased to a passage that is itself rare); inclination is the second factor but only when i ≳ 30°
(118, 236, 39, 86, 232, 277, 22, 179, 280). A few near-Earth-like orbits (261, 136, 208, 51, 160 with q ≥ 1 or a ≈ 1.1–1.3)
are rare for the opposite reason: they never cross the 1-AU ring from inside, so they are met only by matching their
orbit. Encounter speed (U) is irrelevant — as expected for flyby-at-any-velocity routing.

Cheap-early reachability: 10 targets are never inside the first 3 legs of any route — 3, 15, 64, 77, 96, 148, 193,
216, 251, 287 — all a = 2.2–3.7 AU with q = 0.3–1.0 (perihelion at or inside the ring, aphelion at 4–6 AU). They can only
be intercepted near perihelion by a craft already shaped to the right longitude, which is why they end up at the end
of long tours or in dedicated mop-up craft. At the other extreme, 199 and 277 (high-i, a ≈ 2) appear in the first 3
legs 4.5–4.8k times but in ≥ 28-flyby routes only 577/919 times: they are cheap **only** right after launch, when the
launch v∞ can be spent on the plane change — a route that starts with them is committed early.

Targets the TEAM file misses avoidably: 8 (i 42°, 16k routes, only 336 with ≥ 28 flybys), 111 (i 30°, 14k routes),
177 (a 2.68, P 4.4 yr, 12k routes, 284 with ≥ 28), 248 (a 1.9, 27k routes). None is absent; all are in long routes.

## 3. Set-cover MILP over the union pool (`tools/select_routes_milp.py --colgen 8`)

Runs: `results/diag/milp_n{10,11,12}.log` (union of all 13 pool files, 334,998 routes, fuel margin 0.3, reserve 20 kg,
MILP time limit 600 s). Loading takes 6 s, column generation 6 s; RSS during column generation 1.29 GB (n10, `ps`).
Column generation converged for K = 10 (round 7: "no negative reduced costs") and reached 2 remaining negative columns
for K = 11, so the LP values are (essentially) exact lower bounds over the **whole** pool, not just the restricted master:

| max craft | LP bound (planned J) | columns kept | integer result |
|---:|---:|---:|---|
| 10 | **39.499** (exact) | 6992 | HiGHS still running at the 600-s limit when this report was forced to close; see `results/diag/milp_n10/` and `milp_n10.log` |
| 11 | **25.179** (2 negative columns left, so true LP ≤ 25.18) | 7713 | running, `results/diag/milp_n11/` |
| 12 | ≈ 21.1 (earlier camp1+2+3+tail union, `docs/phasing_campaign.md`); union run `results/diag/milp_n12.log` | | running |
| 12 (earlier, camp1+2+3+tail) | 21.1 | | 284 covered, planned J 33.86, gap 3.8 % (`results/phasing/milp_u123cg_n12`: 31,31,30,29,27,26,25,22,22,21,21,18 flybys) |

The LP-bound curve is the key number of this campaign: **going from 12 to 10 craft raises the planned lower bound
from ~21 to 39.5**, i.e. the 10-craft cap alone forces about 22 fractional misses even with fractional (idealised)
routes. The pool contains no ten routes — not even ten fractional combinations — that cover 298 targets. The
integrality gap at 12 craft (LP 21.1 vs integer 33.9) shows the LP cover of the tail is made of fractions of many
specialised short routes, which an integer fleet cannot afford.

Per-craft flyby counts, wall time and `/usr/bin/time -l` peak memory of the integer runs are appended to the logs when
they finish (each log ends with an `EXIT <code>` line); to re-run: the exact command lines are at the top of each log's
process (see `results/diag/milp_n10.log` header) — `tools/select_routes_milp.py results/diag/milp_nK <13 pool files>
--max-craft K --colgen 8`.

## 4. What blocks a 10-craft cover

- **Absent from every route: only 131 and 144.** Every other target is in ≥ 4,269 routes and in ≥ 15 routes of
  ≥ 28 flybys (target 64: 15 such routes; 236: 82; 53: 196; 123: 180; 15: 144). So no target is missing from the pool.
- **The block is combinatorial**: the exact LP bound with 10 craft is 39.5 (≈ 22 fractional misses), 25.2 with 11,
  ~21 with 12. The long routes (28–38 flybys) all draw from the same ~200 ring targets (a < 2.3 AU, i < 30°); the
  rare targets (large a / P 3–7 yr, i > 30°, or q ≥ 1 Earth-like orbits) sit in routes that are either short (mop-up
  chains, ≤ 20 flybys) or long-but-overlapping with the core. Ten near-disjoint 30-flyby routes would need the rare
  targets *distributed among* the long routes; the pool was generated by greedy exclusion rounds, so rare targets only
  enter the pool in the late rounds, where the remaining target set is sparse and the routes are short (camp1 rounds:
  29 31 31 30 28 26 25 22 18 8 3).
- The targets that the LP itself cannot cover at 10 craft are the ones with the highest duals; the 40 rarest in the
  table of section 2 are the candidates (261, 53, 123, 187, 118, 159, 236, 15, 208, 157, 148, 200, 179, 22, 287, 199,
  88, 50, 92, 136, 64, 51, 160, 251, 247, 230, 40, 280, 116, 34, 3, 259, 207, 277, 77, 223, 39, 254, and the TEAM misses
  8, 111, 177, 248). The scripted LP-pricing / best-swap analysis (`lp_prices.py [results/diag/milp_n10]`) that names
  them exactly, with the fractional coverage z_a and the cheapest swap into the 10-craft selection, could not be run
  before this report was forced to close; it needs ~3 min and one process.
- Ten targets are never reached in the first 3 legs of any route (3, 15, 64, 77, 96, 148, 193, 216, 251, 287): they
  must be planned as **mid-tour perihelion interceptions** of an already-shaped craft, not as mop-up.

## 5. Implication for the fleet-building design

The 30.789 fleet already contains a 9-craft disjoint core at 0.060 per target; everything above ~17 is the tail
(55 targets for 16.2). Neither the marginal-cost deletion rule nor a better set cover over the existing pools can fix
this: the pools do not contain ten mutually disjoint ~30-flyby routes, and the LP bound proves it (39.5 at 10 craft).
The pool must be regenerated so that the rare targets are *inside* the long routes: (i) build routes with a per-target
price (LP duals / crossing rarity) high enough that each long tour absorbs 3–5 rare targets, especially the large-a
perihelion targets at their 2030–2044 perihelion dates; (ii) generate the 10 tours jointly or in exclusion rounds that
*start* from the rare targets (seed each tour with 2–3 rare targets and fill with ring targets), instead of greedy
rounds that leave the rare targets for the end; (iii) keep the insertion machinery (`ctoc14/insertion.py`) as the
repair operator to move the ~50 tail targets into core tours at 5–10 kg each rather than 0.29 J each.
