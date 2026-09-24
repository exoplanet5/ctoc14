# CTOC14 Problem A — Ballistic reachability scan

Code (numpy/scipy/pandas/matplotlib only, no pykep): `analysis/ballistic/`
- `ephem.py` — contest Keplerian ephemeris (PDF eqs. 4-5: `n = sqrt(mu/(a AU)^3)` rad/s, `M = M0 + n (t - t_eph)`, Earth epoch
  MJD 60676.0, asteroids MJD 61200.0), vectorised Newton Kepler solver, two-body propagators (universal-variable and
  element-based, elliptic + hyperbolic), element extraction.
- `lambert.py` — vectorised Izzo-2015 Lambert solver: N = 0 (prograde/retrograde) and N >= 1 (both branches), Householder
  iterations with bracketing safeguards; ~2 us per N = 0 solve, ~8 us per N >= 1 pair (1 core).
- `test_lambert.py` — validation (log `out/test_lambert.log`).
- `scan_earth_to_ast.py` — Task 1 (Earth -> asteroid intercept grid); `refine_unreachable.py` — fine grid for the unreachable
  ones + 1 d x 1 d local refinement of every optimum.
- `scan_junction.py` — Task 2 (asteroid -> asteroid Lambert legs, flyby-junction delta-v); `REPORT_ONLY=1` re-runs the report
  and plots from the saved `junction_hist.npz` / `junction_per_epoch.csv` (report in `out/scan_junction_report.log`).
- `estimate_fleet.py` — Task 3 (flybys per spacecraft, fleet size / cost; greedy chain model + fixed-TOF budget table).
- `scan_ballistic_pairs.py` — bonus: zero-fuel DOUBLE flybys (Earth -> A -> B on one conic). Note: the v_inf grid npz is loaded
  eagerly at import (an open `NpzFile` is not fork-safe; workers read garbage otherwise).

All outputs are in `analysis/ballistic/out/` (CSV, NPZ, PNG, one `.log` per script with the full printed report).
Runs with `NPROC=8 ~/.venvs/astro313/bin/python analysis/ballistic/<script>.py`; wall times below are on this 10-core (32 GB) machine:
test_lambert 10 s, scan_earth_to_ast 264 s, refine_unreachable ~250 s, scan_junction 2024 s (7 procs, 1.4 GB each), scan_ballistic_pairs 218 s, estimate_fleet 5 s.

## 0. Validation (`out/test_lambert.log`)

- Earth state at t0 reproduces the PDF sample launch line to 0.54 m; asteroid 174 at the PDF sample flyby epoch (t = 115.84 d)
  is 9.5 m from the sample position; the sample launch v_inf evaluates to 4.000000 km/s. Earth |r| 0.98333-1.01850 AU,
  period 365.76 d. Asteroid element round trip: |da/a| < 2e-14.
- Lambert N = 0: 300 000 random Earth->asteroid geometries (TOF 5-3000 d, 49 % hyperbolic): 100 % converged, mean 3.1
  iterations, Kepler time-of-flight consistency 1.8e-15 median / 1.6e-7 worst, re-propagation position error median 1e-6 km,
  99.9 % < 8e-3 km, max 0.6 km (element propagator, arcs with |e-1| > 1e-6). Retrograde branch: h_z < 0 in 100 % of cases.
- Lambert N = 1, 2, 3 (200 000 cells, both branches): 100 % of existing solutions (T >= T_min) converged; (tof - t_kepler)/P = N to
  2e-9; re-propagation error median 5e-6 km, max 1.9 km; T_min verified against a 4001-point brute-force scan of T(x).

## 1. Task 1 — Earth -> asteroid direct ballistic intercept

Grid: launch date 0 <= t_L < 15 yr step **5 d** (1096 dates) x TOF **30-1500 d step 10 d** (148 values) = 162 208 cells per asteroid,
arrival required inside the window (139 564 valid cells). Per cell the launch v_inf = |v_arc(Earth) - v_Earth| is minimised over the arc
types: N = 0 prograde, N = 0 retrograde, and N = 1..5 revolutions (both branches, prograde; N-rev arcs only for TOF >= 250 N d).
Flyby = position matching only; the arrival relative speed is recorded but unconstrained. Wall time 264 s (8 processes).
Outputs: `earth_to_asteroid_summary.csv` (one row per asteroid, 37 columns incl. `min_vinf`, `min_vinf_refined`, `min_vinf_n0`,
`best_nrev`, `best_tL_day`, `best_tof_day`, `n_cells_le4`, `n_launch_dates_le4`, `min_vinf_tof_le365`, `best_flyby_vrel`, ...),
`earth_to_asteroid_vinf_grid.npz` (`vinf[300,1096,148]` float16 km/s, `vinf_n0`, `nrev_best`, `vinf_min_tof`),
`fig_min_vinf_hist.png`, `fig_launch_calendar.png`.

**Result: 298 of 300 asteroids can be intercepted by a zero-fuel (cost 1) spacecraft** (v_inf <= 4 km/s, arrival in the window);
296 of them with a single-revolution arc alone. The two exceptions are physically out of reach for anything in this window:

| id | a [AU] | e | i [deg] | q [AU] | |r| during 2030-2044 | min launch v_inf (fine grid 2 d x 4 d, TOF 20-5476 d, N = 0..8) |
|---|---|---|---|---|---|---|
| 131 | 8.61 | 0.883 | 134.0 (retrograde) | 1.009 | 10.7-16.2 AU | 13.7 km/s |
| 144 | 17.81 | 0.947 | 19.6 | 0.944 | 33.8-34.7 AU (at aphelion all window) | 15.6 km/s |

Neither ever comes inside 1.05 AU during the window, so they are 2 certain misses (cost 2). Everything else has q <= 1.045 AU
and passes through the Earth annulus at least once.

Distribution of the minimum launch v_inf over the grid (298 reachable; the 1 d x 1 d local refinement lowers the coarse minima by only
0.048 km/s on average, max 0.57 km/s, so the 5 d x 10 d grid is adequate):

| TOF allowed | reachable (<= 4 km/s) | <= 3 | <= 2 | <= 1 | min v_inf percentiles 10 / 25 / 50 / 75 / 90 [km/s] |
|---|---|---|---|---|---|
| <= 90 d | 139 | 111 | 74 | – | 1.06 / 2.02 / 4.34 / 9.91 / 16.05 |
| <= 180 d | 223 | 193 | 142 | – | 0.68 / 1.15 / 2.19 / 4.10 / 6.11 |
| <= 365 d | 289 | 274 | 254 | 153 | 0.45 / 0.67 / 0.98 / 1.47 / 2.73 |
| <= 730 d | 298 | 296 | 283 | – | 0.36 / 0.51 / 0.72 / 1.02 / 1.61 |
| <= 1500 d (full grid) | 298 | 298 | 298 | 255 | 0.27 / 0.42 / 0.61 / 0.82 / 1.14 |

Threshold counts on the full grid: min v_inf <= 0.5 km/s for 114 asteroids, <= 1.0 for 255, <= 1.5 for 294, <= 2.0 for 298
(single-revolution arcs only: 39 / 153 / 232 / 260, <= 4 km/s: 296).

Structure of the cheap solutions:
- The optimum is a multi-revolution arc for 264/298 asteroids (best N: 0 -> 36, 1 -> 39, 2 -> 48, 3 -> 161, 4 -> 16): the spacecraft
  leaves Earth with a small v_inf onto a near-Earth orbit (best arcs: a = 0.96-1.04 AU, e = 0.01-0.06, i = 0.1-1.5 deg, 10-90 %
  ranges) and simply waits, near the asteroid's node/MOID crossing, for the asteroid to come by. Median best TOF 1285 d
  (25/75 %: 810/1398 d). Of all 1.16 million feasible cells (2.4 % of the 48.7 million grid cells), 87 % are won by N >= 1 arcs;
  only 10.7 % of feasible cells have TOF <= 365 d (of which 96.6 % are N = 0).
- Retrograde arcs are never competitive (0 feasible cells).
- Feasible cells per reachable asteroid: median 2596 (min 159 for id 64, max 29 450); launch dates with at least one feasible
  TOF: median 740 of 1096 (10/90 %: 439/959). From a given launch date, a median of 222 asteroids (max 242) are reachable with
  some TOF <= 1500 d, and a median of 41 with TOF <= 365 d. Per launch year 2030..2044 the number of asteroids with >= 1 feasible
  launch date is 290, 292, 291, 292, 292, 291, 289, 293, 290, 291, 291, 273, 240, 163, 57 (the tail is the window end cutting the
  long arcs).
- Hardest reachable targets (min v_inf 1.4-1.9 km/s, all needing multi-rev arcs of 1230-1340 d): ids 236 (i = 75.7 deg), 64 (a = 3.66 AU,
  i = 44 deg), 110, 7, 219, 216. With TOF <= 365 d the worst reachable cases need 5.5-9.1 km/s (ids 64, 121, 110, 114, 236, 148 —
  a >= 2.5 AU, i.e. large aphelia, arriving fast).
- Flyby relative speed at the best cell: median 16.3 km/s (5.7-42.7); median over all feasible cells per asteroid: 9.7 / 16.1 / 29.9 km/s
  (10/50/90 %). Approaches are fast: the spacecraft cannot linger, every subsequent target needs a new arc.
- Min v_inf correlates weakly with inclination (r = 0.20) and eccentricity (0.22), not with q (0.06).

Interpretation: **every asteroid except 131 and 144 is a legitimate target for a cost-1 ballistic spacecraft**, but a ballistic
spacecraft has only 4 free parameters (launch date, v_inf vector), so it generically hits at most 2 asteroids (section 4).

## 2. Task 2 — asteroid -> asteroid legs and the flyby-junction delta-v

Definition. After flying by B on the incoming Lambert arc A -> B (arrival velocity v_in at B), the next ballistic leg B -> C
needs the outgoing Lambert velocity v_out at B; the **junction delta-v** is dv_B = |v_out - v_in|. The next leg is called
feasible if the best next hop satisfies min_C dv_B <= eta * a * TOF_out (low-thrust capability along the next leg;
a = 0.25 / 0.5 / 0.8 mm/s^2 = 2000 / 1000 / 625 kg at 0.5 N, eta = 0.7 / 0.8 / 0.9).

Grid (`scan_junction.py`, `out/scan_junction.log` + `out/scan_junction_report.log`): junction epochs t_B = 0..5460 d step **30 d** (183),
TOF_in and TOF_out = **60..600 d step 30 d plus 365 d** (20 values), **all 300 asteroids as B** and all 299 others as A and C,
arc types = single-revolution prograde for every TOF plus the two 1-revolution branches for TOF >= 250 d (a ~1 AU spacecraft can
do 1.x revolutions). Outgoing arcs are exhaustive (299 x 20 x 3 per (B, t_B)); incoming arcs are randomly subsampled to 200 per
(B, t_B) within each arc class (7.66 M "plausible" and 1.99 M "earthlike" samples in total; 6.6e8 single-rev + 3.9e8 one-rev
Lambert pairs). Arc classes (applied to both arcs, so that hyperbolic / Jupiter-crossing "solutions" are excluded):
`plausible` = bound, q > 0.4 AU, Q < 3.0 AU (24 % of valid arcs, 1373 of 5617 per (B, t_B)); `earthlike` = q > 0.6 AU, Q < 1.8 AU,
i < 10 deg (1 %, 55 per (B, t_B); this is what an Earth launch with v_inf <= 4 km/s can produce). Greedy-hop statistics are
accumulated only for t_B <= 15 yr - 600 d so that the window end is not mistaken for "no hop". Wall time 2024 s (7 processes,
sharing the machine with the pair scan of section 4). Outputs: `junction_hist.npz` (histograms, 0.1 km/s bins),
`junction_per_epoch.csv` (54 900 rows), `junction_per_asteroid.csv`, `fig_junction_dv_cdf.png`, `fig_junction_feasible_vs_tof.png`.

### 2.1 Junction delta-v distributions (class `plausible`, both arcs)

| TOF_out [d] | dv_B to a random C, pct 10/50/90 [km/s] | best next hop min_C dv_B, pct 10/50/90 [km/s] | P(feasible) eta=0.8, a=0.25 / 0.5 / 0.8 | eta=0.7, a=0.25 / 0.8 | eta=0.9, a=0.25 / 0.8 | mean # of C with dv_B <= 1 / 2 / 3 / 5 km/s |
|---|---|---|---|---|---|---|
| 60 | 8.8 / 21.0 / 40.0 | 4.4 / 10.2 / 22.0 | 0.2 % / 1.4 % / 5.0 % | 0.1 / 3.5 % | 0.3 / 6.8 % | 0.00 / 0.01 / 0.04 / 0.16 |
| 90 | 8.4 / 20.5 / 40.1 | 3.2 / 7.2 / 15.9 | 1.5 / 9.7 / 27.7 % | 1.0 / 21.1 % | 2.1 / 34.3 % | 0.00 / 0.03 / 0.10 / 0.39 |
| 120 | 8.2 / 20.2 / 39.9 | 2.5 / 5.6 / 12.1 | 6.3 / 31.4 / 62.2 % | 4.4 / 53.3 % | 8.5 / 69.5 % | 0.01 / 0.06 / 0.19 / 0.76 |
| 150 | 7.9 / 19.6 / 38.8 | 2.1 / 4.6 / 9.9 | 17.1 / 58.7 / 84.5 % | 12.5 / 78.7 % | 22.3 / 88.6 % | 0.01 / 0.10 / 0.32 / 1.26 |
| 180 | 7.6 / 18.9 / 37.9 | 1.9 / 4.0 / 8.4 | 34.0 / 78.8 / 94.1 % | 26.1 / 91.2 % | 41.7 / 95.9 % | 0.02 / 0.15 / 0.48 / 1.88 |
| 240 | 7.0 / 17.9 / 36.8 | 1.5 / 3.2 / 6.7 | 68.1 / 95.1 / 99.2 % | 59.3 / 98.5 % | 75.1 / 99.5 % | 0.04 / 0.29 / 0.89 / 3.42 |
| 300 | 6.7 / 17.4 / 36.5 | 1.3 / 2.8 / 5.8 | 86.6 / 98.9 / 99.9 % | 80.9 / 99.8 % | 90.6 / 99.9 % | 0.06 / 0.44 / 1.36 / 5.14 |
| 365 | 6.4 / 17.1 / 36.3 | 1.2 / 2.6 / 5.3 | 94.5 / 99.7 / 100 % | 91.4 / 99.9 % | 96.4 / 100 % | 0.08 / 0.60 / 1.84 / 6.84 |
| 450 | 6.3 / 17.0 / 36.1 | 1.1 / 2.4 / 5.1 | 97.9 / 99.9 / 100 % | 96.5 / 100 % | 98.7 / 100 % | 0.10 / 0.74 / 2.26 / 8.23 |
| 600 | 6.6 / 17.2 / 35.9 | 1.1 / 2.5 / 5.0 | 99.6 / 100 / 100 % | 99.2 / 100 % | 99.8 / 100 % | 0.10 / 0.70 / 2.12 / 7.80 |

Class `earthlike` (both arcs; the population an Earth-launched spacecraft starts in):

| TOF_out [d] | random C, pct 10/50/90 | best next hop, pct 10/50/90 | P(feasible) eta=0.8, a=0.25 / 0.5 / 0.8 | mean # of C with dv_B <= 1 / 2 / 3 / 5 km/s |
|---|---|---|---|---|
| 60 | 3.5 / 7.7 / 14.8 | 2.7 / 6.0 / 12.1 | 0.7 / 5.1 / 16.5 % | 0.01 / 0.05 / 0.14 / 0.47 |
| 120 | 3.3 / 7.4 / 14.5 | 1.8 / 3.8 / 7.3 | 14.8 / 57.0 / 86.9 % | 0.02 / 0.15 / 0.44 / 1.42 |
| 180 | 3.4 / 7.6 / 15.3 | 1.5 / 3.0 / 5.2 | 53.4 / 95.8 / 99.6 % | 0.04 / 0.28 / 0.80 / 2.67 |
| 240 | 3.3 / 7.4 / 14.9 | 1.3 / 2.5 / 4.4 | 87.9 / 99.6 / 100 % | 0.06 / 0.44 / 1.29 / 4.13 |
| 365 | 3.2 / 7.3 / 14.7 | 1.2 / 2.4 / 4.2 | 98.8 / 100 / 100 % | 0.08 / 0.57 / 1.60 / 4.93 |
| 600 | 3.3 / 7.2 / 14.1 | 1.2 / 2.5 / 4.4 | 99.9 / 100 / 100 % | 0.07 / 0.49 / 1.39 / 4.55 |

**Feasible junction fraction for the requested TOFs (best next hop, eta = 0.8):**

| next-leg TOF | a = 0.25 mm/s^2 (2000 kg) | a = 0.5 (1000 kg) | a = 0.8 (625 kg) | earthlike arcs, a = 0.25 / 0.5 / 0.8 |
|---|---|---|---|---|
| 60 d | 0.2 % | 1.4 % | 5.0 % | 0.7 / 5.1 / 16.5 % |
| 120 d | 6.3 % | 31.4 % | 62.2 % | 14.8 / 57.0 / 86.9 % |
| 180 d | 34.0 % | 78.8 % | 94.1 % | 53.4 / 95.8 / 99.6 % |
| 365 d | 94.5 % | 99.7 % | 100 % | 98.8 / 100 / 100 % |

Reading: a random pair of consecutive legs is hopeless (median 17-21 km/s junction, only 1 % of random next targets are within
reach of a 2000-kg spacecraft over 180 d), but the *best* next target is cheap: median 4.0 km/s at 180 d, 2.6 km/s at
365 d, 2.4-2.5 km/s asymptotically (longer legs do not help beyond ~400 d; the 1-rev arcs take over there). The catalogue is thin,
though: at any junction only ~0.5 (180 d) to ~1.8-2.1 (365-600 d) asteroids lie within 3 km/s, and ~0.15-0.7 within 2 km/s. The
distributions barely depend on the junction asteroid: per asteroid the feasible fraction (180 d, a = 0.25) spans 0.24-0.41
(5-95 %), 0.64-0.84 at a = 0.5, and 0.90-0.96 at 365 d; the hardest junction asteroids are the high-inclination ones (ids 10, 66,
260, 86, 236, 29, 147, 125, 213 with i = 39-76 deg and id 297 with a = 0.70 AU: 19-22 % feasible at 180 d / a = 0.25, still
83-97 % at 365 d). Mean feasible fraction (180 d, a = 0.5) by inclination: 0.81 (i < 10 deg, 120 asteroids), 0.79 (10-20, 75),
0.76 (20-30, 47), 0.69 (> 30, 58).

### 2.2 Greedy "next flyby as soon as possible" statistics

For each sampled incoming arc, the shortest grid TOF_out with a feasible best next hop, and the junction dv on that leg
(`report_greedy` in `scan_junction_report.log`). No incoming arc lacks a feasible hop within 600 d (0.0 % for every setting).

| class | eta | a [mm/s^2] | shortest feasible TOF: mean, pct 25/50/75 [d] | implied flybys / yr | junction dv on that leg: mean, pct 25/50/75 [km/s] | dv per year of flying [km/s/yr] |
|---|---|---|---|---|---|---|
| plausible | 0.8 | 0.25 | 215, 180 / 210 / 240 | 1.70 | 2.69, 1.95 / 2.55 / 3.35 | 4.6 |
| plausible | 0.8 | 0.50 | 157, 120 / 150 / 180 | 2.32 | 3.79, 2.65 / 3.65 / 4.75 | 8.8 |
| plausible | 0.8 | 0.80 | 129, 90 / 120 / 150 | 2.83 | 4.77, 3.25 / 4.55 / 6.05 | 13.5 |
| plausible | 0.7 | 0.25 | 229, 180 / 210 / 270 | 1.60 | 2.52 | 4.0 |
| plausible | 0.9 | 0.25 | 204, 150 / 210 / 240 | 1.79 | 2.85 | 5.1 |
| earthlike | 0.8 | 0.25 | 175, 150 / 180 / 210 | 2.08 | 2.06, 1.45 / 2.05 / 2.55 | 4.3 |
| earthlike | 0.8 | 0.50 | 127, 120 / 120 / 150 | 2.88 | 2.81, 2.05 / 2.75 / 3.55 | 8.1 |
| earthlike | 0.8 | 0.80 | 105, 90 / 90 / 120 | 3.47 | 3.47, 2.45 / 3.35 / 4.45 | 12.0 |

Hopping as fast as the thruster allows costs 4.6 km/s per year at 2000 kg (58 % of the 7.9 km/s/yr continuous-thrust
capability) and 13.5 km/s/yr at 625 kg: a light spacecraft flies faster but burns its tank in ~3.5 yr of flight either way.

## 3. Task 3 — flybys per spacecraft, fleet size

(`estimate_fleet.py`, `out/estimate_fleet.log`, `fleet_estimate.csv`, `fleet_budget_table.csv`, `fig_fleet_estimate.png`.)
First leg: from any launch date some asteroid is reachable ballistically (v_inf <= 4) with TOF **60 d** (median over launch dates;
10/90 %: 30/100 d); per asteroid the shortest feasible TOF is 100 d median (30/230 d). Each further leg is charged the junction
dv x LT_penalty (1.0 / 1.25 / 1.5 for the low-thrust-vs-impulse inefficiency), mass by the rocket equation (v_e = 39.227 km/s),
and the acceleration a = 0.5 N / m is updated (0.25 mm/s^2 at 2000 kg -> 0.83 at 600 kg; leg statistics interpolated in log a).

**Model A — greedy shortest hop** (expected-value chain, eta = 0.8; `plausible` class unless stated):

| fuel [kg] | J_i | flybys, LT penalty 1.0 / 1.25 / 1.5 | chain ends after [yr] (pen 1.25) | J_i per flyby (pen 1.25) | spacecraft for 298 asteroids | J_total = N_sc J_i + 2 misses |
|---|---|---|---|---|---|---|
| 0 | 1.00 | 1 (ballistic; 2 with a double flyby, section 4) | – | 1.0 (0.5) | 298 (149) | 300 (151) |
| 100 | 1.08 | 2 / 2 / 1 | 0.5 | 0.54 | 149 | 162 |
| 200 | 1.16 | 3 / 3 / 2 | 0.9 | 0.39 | 100 | 118 |
| 350 | 1.31 | 5 / 4 / 3 | 1.4 | 0.33 | 75 | 100 |
| 500 | 1.48 | 6 / 5 / 4 | 1.8 | 0.30 | 60 | 91 |
| 700 | 1.75 | 8 / 7 / 6 | 2.8 | 0.25 | 43 | 77 |
| 1000 | 2.22 | 11 / 9 / 7 | 3.8 | 0.25 | 34 | 78 |
| 1400 | 3.00 | 14 / 11 / 10 | 5.1 | 0.27 | 28 | 86 |
| 1400, earthlike arcs | 3.00 | 18 / 15 / 12 | 5.6 | 0.20 | 20 | 62 |

Every fuelled chain is **fuel-limited, never time-limited**: a full tank flown greedily is empty after 5-6 yr (11-15 flybys at
2.7-5.4 km/s per leg, legs shortening from 215 d to 140 d as the spacecraft lightens). Cost per flyby is flat at 0.19-0.27
from 500 kg upwards, i.e. J ~ 60-80 for the catalogue with 25-45 spacecraft.

**Model B — fixed leg TOF T with the typical best-next-hop dv(T)** (a planner that trades time for fuel):
n = 1 + min( 47.2 km/s / (1.25 dv(T)), (15 yr - 60 d) / T ).

| T [d] | dv(T) median / 25 %-ile [km/s] | eta a T at a = 0.25 [km/s] | dv-limited (pen 1.0 / 1.25) | time-limited | n (median dv) | n (25 %-ile dv) | J_i / flyby | N_sc (298) |
|---|---|---|---|---|---|---|---|---|
| 120 | 5.6 / 3.7 | 2.1 (infeasible for most) | 8.4 / 6.7 | 45 | 7 | 11 | 0.43 | 43 |
| 180 | 4.0 / 2.7 | 3.1 (34 % feasible) | 11.8 / 9.4 | 30 | 10 | 14 | 0.30 | 30 |
| 240 | 3.2 / 2.2 | 4.2 | 14.8 / 11.8 | 22.6 | 12 | 18 | 0.25 | 25 |
| 300 | 2.8 / 1.9 | 5.2 | 16.9 / 13.5 | 18.1 | 14 | 19 | 0.21 | 22 |
| 365 | 2.6 / 1.7 | 6.3 | 18.2 / 14.5 | 14.8 | **15** | 15 | 0.20 | 20 |
| 450 | 2.4 / 1.6 | 7.8 | 19.7 / 15.7 | 12.0 | 13 | 13 | 0.23 | 23 |
| 600 | 2.5 / 1.7 | 10.4 | 18.9 / 15.1 | 9.0 | 10 | 10 | 0.30 | 30 |

(`earthlike` arcs: 16 flybys at T = 240-300 d, 21 with 25 %-ile dv.) The optimum is a leg of ~300-365 d (~1 revolution) at
~2.5-3 km/s per junction, where the fuel (47 km/s) and the window (15 yr) run out together: **~15 flybys per full-tank spacecraft
(18-19 if the planner consistently finds 25 %-ile junctions), 1.0-1.3 flybys per year, J_i/flyby ~ 0.16-0.20**. The whole
catalogue (298 reachable) then needs **~16-20 full-tank spacecraft, J ~ 50-62 (+2 unavoidable misses)**, against 149 cost-1
double-flyby spacecraft (J = 151) or 298 single-flyby spacecraft (J = 300). Intermediate fuel loads are not better per flyby
(J_i/flyby minimum 0.19-0.25 at 500-1000 kg in model A) because the J_i = 1 + x + x^2 curve is only mildly convex while
the 1-per-launch overhead is large: prefer full or nearly full tanks.

Caveats (direction of bias): (i) impulsive junction at B with ballistic Lambert legs — a real low-thrust arc is not a Lambert
arc and can exploit the whole leg, so the true feasible set is larger (optimistic for the planner, i.e. our n is a lower bound in
this respect); the LT penalty 1.25 is a guess. (ii) **Catalogue depletion**: only ~1.8 asteroids per junction are within 3 km/s at
365 d (0.5 at 180 d). The statistics above are for an untouched catalogue; as it is consumed the branching factor drops in
proportion, so the last ~30 % of the asteroids will cost markedly more per flyby than the first 70 % — the fleet count is a lower
bound in this respect. (iii) High-inclination asteroids (i > 30 deg, 58 of them) are 15-20 % less connected than the rest; ids
131 and 144 are unreachable by anything. (iv) Statistics use the median incoming arc; a planner chooses its state, so real chains
should sit between the "median" and "25 %-ile" columns.

## 4. Zero-fuel double flybys (Earth -> A -> B on one conic)

(`scan_ballistic_pairs.py`, `out/scan_ballistic_pairs.log`, `ballistic_pairs.csv` (accepted pairs, 12 significant digits),
`ballistic_pairs_refined_all.csv`, `ballistic_pair_candidates.csv`, `ballistic_pairs_summary.txt`; 218 s with 8 processes.)
A ballistic spacecraft has 4 free parameters (t_L, v_inf vector) and each flyby removes 2 (3 position components minus the free
flyby time), so it can hit exactly 2 asteroids at isolated parameter points. Method: for every asteroid A take up to 500 random
feasible cells of the Task-1 grid (v_inf <= 4, best arc type), coast the conic from the flyby of A to the window end (2 d
step) and record the closest approach to every other asteroid B (candidate if < 0.06 AU): 1.59 M (cell, B) hits, clustered into
856 205 (A, B, N, month-of-t_B) clusters (31 624 with miss < 0.01 AU); the best cluster of each ordered pair (85 452) is then
refined by a vectorised Levenberg-Marquardt on (t_L, t_A, t_B) with the miss vector at B as residual (the arc is the Lambert arc
Earth(t_L) -> A(t_A) of the cluster's N/branch, so the flyby of A is exact by construction). Acceptance: miss at B <= 1000 km,
v_inf <= 4 km/s, 0 <= t_L, t_B <= 15 yr.

Result: **57 499 ordered pairs (38 244 unordered) are exact zero-fuel double flybys**; a further 11 071 converge to a double
flyby but need v_inf > 4 km/s, 16 666 clusters do not converge to <= 1000 km. Every one of the 298 reachable asteroids has
between 179 and 285 (median 260) distinct partners, so essentially any two reachable asteroids can be paired; the least connected
are ids 64, 236, 139, 119, 114, 11, 218, 216 (179-208 partners each). A greedy matching covers all **298 asteroids with 149
cost-1 spacecraft: J = 149 + 2 = 151, 0.5 per asteroid** — the trivial baseline any fuelled design must beat. Pair properties:
launch v_inf min 0.20, median 2.81 km/s (45 853 pairs <= 3.5 km/s, 33 302 <= 3.0, 12 265 <= 2.0); arc types N = 0: 8 157,
N = 1: 13 757, N = 2: 17 369, N = 3: 14 895, N = 4: 3 269, N = 5: 52; time Earth -> A median 916 d (10/90 %: 313/1382 d), A -> B
median 1538 d (346/3235 d); flyby speeds ~16 km/s median at both asteroids. Verification: 300 random accepted pairs re-propagated
independently (universal-variable propagator, and DOP853 at rtol 1e-13 for the worst cases): miss at B median 0.04 km, max 0.95 km.
(The first run of this script wrote the times with 8 significant digits, which alone produced up to 5 400 km of miss on
re-propagation — 1e-4 d at 20-30 km/s relative speed; the CSVs are now written with 12 digits.) Ballistic triple flybys are
generically impossible (6 constraints vs 4 parameters); the chance that a third asteroid passes within 1000 km of one of the
57 499 pair arcs by accident is ~(1000 km / 0.06 AU)^2 x 10 hits/arc x 57 499 arcs ~ 0.01, so none are expected.

## 5. Summary for the planner

- 298/300 asteroids are ballistically reachable (ids 131, 144 never come inside 10 AU during the window); min launch v_inf median
  0.61 km/s (all arc types), 0.98 km/s single-rev; 289 reachable with TOF <= 365 d.
- Zero-fuel double flybys are abundant (38 244 unordered pairs, every reachable asteroid has >= 179 partners): the cost-1 baseline
  is 0.5 per asteroid, J = 151.
- Junctions: best next target median 4.0 km/s at 180 d, 2.6 km/s at 365 d, floor ~2.4 km/s; feasible fraction at 2000 kg / eta 0.8:
  0.2 % (60 d), 6 % (120 d), 34 % (180 d), 95 % (365 d); at 625 kg: 5 / 62 / 94 / 100 %. Only ~0.5 (180 d) to ~2 (365 d) candidate
  targets per junction within 3 km/s: the catalogue is thin and will deplete.
- Full-tank spacecraft: 11-15 flybys in 5-6 yr if hopping as fast as possible (fuel-limited), ~15 (up to ~19) flybys over the full
  15 yr with ~1-yr legs at ~2.5-3 km/s each; J_i/flyby ~ 0.16-0.27; the catalogue needs ~16-20 full-tank spacecraft in the optimistic
  fixed-TOF model and ~25-45 in the greedy model, i.e. J ~ 50-80 plus the 2 certain misses. Use full tanks; do not expect legs
  shorter than ~150-200 d early in a heavy spacecraft's life.
