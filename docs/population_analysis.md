# CTOC14 Problem A — Target Population Analysis (300 NEAs)

Code (numpy/scipy/pandas/matplotlib only):
- `analysis/population/kepler.py` — contest-exact Keplerian propagation (eq. 4-5 of the PDF).
- `analysis/population/population_analysis.py` — the analysis; runs in ~16 s with `~/.venvs/astro313/bin/python`; log in `run.log`.
- `analysis/population/verify.py` — independent cross-checks (PDF sample lines, the lead's `ctoc14/kepler.py`, 0.1-d brute-force
  recounts, brute-force MOID); ~85 s; results in `verify.json` / `verify.log` (section 6 below).

Outputs in `analysis/population/` (all floats written with 12 significant digits):
- `population_stats.csv` — one row per asteroid, 60 columns (elements, q, Q, period, n [rad/s and deg/day], node distances and
  longitudes, annulus/band flags, r and z at t0, r_min/r_max in window with dates, perihelion count and MJD list, ascending/descending
  crossing counts and `MJD:r` lists, MOID with location, min Earth distance with date and relative speed, Earth-approach counts,
  band/accessible days, access-window count, Earth-relative speed range, `difficulty`, `note`).
- `nodal_crossings.csv` — 4293 rows, one per ecliptic-plane crossing in the window (id, node, MJD, date, t[s], r, longitude, x, y, z,
  Earth distance, v_rel to Earth, v_z, i, in_band).
- `earth_close_approaches.csv` — 331 rows, local minima of Earth distance < 0.3 AU (id, MJD, date, t[s], d, v_rel, r, z).
- `access_windows.csv` — 3220 rows, contiguous intervals with 0.8 <= r <= 1.6 AU and |z| <= 0.1 AU.
- `summary.json` — every aggregate number quoted below; `verify.json` — every verification number in section 6.
- PNGs: `scatter_a_e_a_i.png`, `hist_i_q_Q.png`, `timeline_nodal_crossings.png`, `timeline_earth_approaches_and_band_count.png`,
  `polar_nodal_crossings.png`.

## 0. Propagation model

- `n = sqrt(mu_sun/(a*AU)^3)` rad/s, `M(t) = M0 + n*(t - t_eph)`, t in seconds (MJD difference x 86400); asteroid epoch MJD 61200.0,
  Earth epoch MJD 60676.0 with the Table-3 elements; mu = 1.32712440018e11 km^3/s^2, AU = 149597870.7 km; frame HEI J2000.
  Vectorised Newton Kepler solver (Danby starter; max residual 8.9e-16 rad over the whole 300 x 10958 grid, e up to 0.947).
- Mission window: MJD 62502.0 (2030-01-01) to MJD 67980.75 (2044-12-31), 5478.75 d. Time-grid scans use 0.5-d steps (10 958 samples);
  event epochs (perihelia, nodes, Earth approaches) are solved analytically from the mean anomaly or refined by bounded 1-D minimisation,
  never read off the grid.

## 1. Orbital-element census

| quantity | min | median | max |
|---|---|---|---|
| a [AU] | 0.635 | 1.901 | 17.81 |
| e | 0.074 | 0.586 | 0.947 |
| i [deg] | 0.45 | 13.5 | 134.0 |
| q [AU] | 0.140 | 0.852 | 1.045 |
| Q [AU] | 0.959 | 3.03 | 34.7 |
| period [yr] | 0.506 | 2.62 | 75.2 |

- Orbit classes: 264 Apollo (a >= 1, q < 1.017), 19 Amor (q 1.017-1.3), 15 Aten (a < 1, Q > 0.983), 2 Atira (Q < 0.983: IDs 296, 297;
  Q = 0.973 / 0.959 AU). **All 300 have q <= 1.045 AU, 271 have q < 1 AU**, so every orbit crosses the 0.7-1.5 AU annulus and the
  0.8-1.6 AU band (`crosses_annulus_0p7_1p5 = 300`, `crosses_band_0p8_1p6 = 300`).
- Inclination: 42 with i <= 5 deg, **120 with i <= 10 deg**, 160 <= 15, 195 <= 20, 242 <= 30, 58 > 30, 23 > 45, **1 retrograde
  (ID 131, i = 134.0 deg)**. Histogram peak at 5-10 deg.
- Periods: 17 objects with P < 1 yr, 92 with P < 2 yr, 56 with P > 4 yr. Ten objects have a > 3 AU
  (58, 64, 129, 138, 139, 152, 204, 251: P = 5.3-7.0 yr, 2-3 perihelia in the window; plus the giants 131: 8.61 AU, 144: 17.8 AU).
  148 objects have Q <= 3 AU, 62 have Q <= 2 AU, 37 have Q <= 1.6 AU.
- 14 objects have q < 0.3 AU (IDs 4 = Phaethon q = 0.140, 49, 60, 95, 113, 125, 165, 188, 198, 241, 265, 278, 293, 296).
- **Geometric MOID to the contest Earth ellipse: 298/300 below 0.05 AU** (the PHA criterion), 118 below 0.02 AU, 52 below 0.01 AU;
  minimum 0.00011 AU (ID 76, 1997 XF11). The only two at/above 0.05 AU are **ID 117** (0.0501 AU; JPL SBDB 0.0485) and
  **ID 280** (0.0519 AU; SBDB 0.0513) — both still trivially reachable. (The unreachable 131 and 144 have small MOIDs, 0.026 and 0.0047 AU;
  their problem is phasing, see section 3.) Cross-check with JPL SBDB via the lead's `data/nea_identification.csv`: median |difference|
  0.00027 AU; the largest discrepancies (IDs 278, 206, 17, 94, 187; up to 0.019 AU) are objects whose SBDB elements differ from MEA.txt,
  so the MEA.txt-based value is the one that matters for the contest.
- Speed relative to Earth (minimum over the window, per asteroid): quantiles 10 % 2.05, 25 % 3.29, 50 % 5.6, 75 % 8.8 km/s — a quarter of
  the targets are at some point within ~3.3 km/s of Earth's velocity, relevant for cheap ballistic first legs.

## 2. Mission-window event statistics (2030-01-01 to 2044-12-31)

- **Perihelion passages: 2146 total**, median 6 per asteroid (10 % have <= 3, 25 % <= 4, 75 % <= 8, max 30 for ID 296). Only 131 and 144 have none.
- **Ecliptic-plane crossings: 4293 total** (median 12 per asteroid, max 59): 716 at r < 0.8 AU, **2731 inside 0.8-1.6 AU**, 846 at r > 1.6 AU.
  Orbits are fixed, so each asteroid crosses at constant r_asc / r_desc (`r_asc_node_au`, `r_desc_node_au`); every crossing epoch and its
  r is listed in `asc_crossings_mjd_r` / `desc_crossings_mjd_r` and in `nodal_crossings.csv`.
- Heliocentric distance in window: for 298 objects r_min = q (a perihelion occurs in the window; dates in `date_r_min`); 298 objects spend time in
  the band. **On average 91 asteroids are inside 0.8-1.6 AU at any instant** (30 % of asteroid-time; range 67-105), and 32 are simultaneously
  inside the band and within |z| <= 0.1 AU of the ecliptic (range 19-49).
- Earth close approaches (local minima of the Earth distance, refined to 1e-6 d): 331 events below 0.3 AU, **86 below 0.1 AU (73 distinct
  asteroids), 33 below 0.05 AU (32 asteroids), 3 below 0.02 AU**:
  ID 243 (267131 2000 EK26) 0.00542 AU = 811 000 km on 2041-02-02 (v_rel 16.7 km/s);
  ID 65 (5693 1993 EA) 0.00588 AU = 879 000 km on 2037-01-07 (18.8 km/s);
  ID 102 (161989 Cacus) 0.0145 AU on 2041-09-11 (15.0 km/s);
  next: 226 (0.0218 AU, 2036-02-25, 23.9 km/s), 174 (0.0231 AU, 2030-04-07, 12.0 km/s — the target of the PDF's sample ballistic first leg),
  20 (0.0265 AU, 2043-09-05, 8.0 km/s).
  Relative speed at the < 0.1 AU approaches: 6.5-31.8 km/s, median 13.3 km/s, 23 events below 10 km/s, none below 6 km/s — so approaches are
  fast and never inside the 4 km/s launch cone; a ballistic spacecraft must be launched to intercept, not to co-orbit.
  Per-year counts of < 0.1 AU approaches (2030..2044): 5, 3, 8, 4, 6, 3, 12, 6, 5, 4, 9, 4, 6, 7, 4.

## 3. Difficulty classification

Thresholds (`summary.json['thresholds']`): band 0.8-1.6 AU; "near ecliptic" |z| <= 0.10 AU; **unreachable** if retrograde or r_min(window) > 2.5 AU;
**easy** if i <= 10 deg and q <= 1.3 AU with accessible time > 0; **moderate** if 10 < i <= 20 deg with accessible time > 0; **hard** if i > 20 deg and
(>= 1 in-band nodal crossing or accessible time > 0); **very_hard** otherwise. An *access window* is a contiguous interval with 0.8 <= r <= 1.6 AU
and |z| <= 0.1 AU (a near-1-AU, near-ecliptic spacecraft could be co-located without a large plane change).

| class | N | median i [deg] | median q [AU] | in-band nodal crossings (total) | access windows (total / median length) | median accessible days per asteroid |
|---|---|---|---|---|---|---|
| easy | **120** | 5.9 | 0.80 | 878 | 1250 / 81.5 d | 715 |
| moderate | **74** | 14.0 | 0.88 | 693 | 763 / 42.5 d | 263 |
| hard | **104** | 32.4 | 0.88 | 1160 | 1207 / 23.5 d | 166 |
| very_hard | **0** | – | – | 0 | 0 | – |
| unreachable | **2** (131, 144) | – | – | 0 | 0 | 0 |

Because every target has q <= 1.045 AU and (apart from 131/144) reaches the band during the window, the classification collapses to inclination
bins: **298 of 300 targets are geometrically reachable** by a spacecraft that stays near 1 AU and near the ecliptic. What differs is the
*width* of the timing windows: hard-class access windows are ~24 d long versus ~82 d for easy ones, and hard-class targets can only be met near
their nodes.

- All 104 hard-class objects have >= 1 in-band nodal crossing in the window; 103 have >= 3, 69 have >= 6 (median 8, max 33). Their in-band
  nodal distances cluster around 1 AU (quartiles 0.97 / 1.01 / 1.05 AU).
- Per-asteroid access windows: 296 asteroids have >= 3, 232 have >= 5 (quartiles 5 / 8 / 15, max 34). Only four have <= 2: 131 and 144 (0),
  138 (2; i = 47 deg, a = 3.29 AU), 251 (2; i = 16.4 deg, a = 3.07 AU).
- Eight objects have no in-band nodal crossing; six of them (18, 127, 192, 211, 267, 292) are low-i (< 6 deg) with nodes just outside the band
  (0.73-0.79 AU or 1.6-3.5 AU) and 350-2300 accessible days, so they stay "easy"; the other two are 131 and 144.
- Very low-q objects (q < 0.3 AU, 14 objects) are still classed by inclination: they cross the band twice per orbit and are best met on the
  inbound/outbound legs near 1 AU, where their speed relative to a ~1 AU spacecraft is 20-40 km/s (irrelevant for a position-only flyby).

The 23 objects with i > 45 deg (the hardest reachable ones), sorted by i, with the number of in-band nodal crossings / access windows in the window:
236 (75.7 deg, 4/4), 167 (73.6, 3/3), 147 (70.6, 4/4), 118 (68.4, 8/8), 24 (64.0, 21/21), 112 (62.1, 27/27), 27 (60.7, 7/7), 213 (59.3, 10/10),
139 (57.8, 3/3), 57 (55.9, 9/9), 86 (55.9, 12/12), 232 (55.0, 4/4), 219 (54.5, 3/3), 63 (53.2, 6/6), 39 (52.4, 5/5), 29 (51.3, 12/12),
31 (48.4, 3/3), 41 (47.8, 3/3), 260 (47.7, 12/12), 138 (47.0, 2/2), 36 (46.8, 6/6), 226 (45.9, 28/28) — plus 131 (134 deg, 0/0).
For these, each access window *is* a nodal passage (window count = crossing count), typically 25-40 d wide.

### The two unreachable objects (verified with dates)

- **ID 131 (2025 VP): a = 8.605 AU, e = 0.8827, i = 134.0 deg (retrograde), q = 1.009 AU, Q = 16.2 AU, P = 25.24 yr.**
  Last perihelion 2025-12-20 (MJD 61029.4), next 2051-03-19 (MJD 70249.9). Heliocentric distance in the window: 10.71 AU (2030-01-01),
  15.33 (2035-01-01), 16.07 (2040-01-01), aphelion 16.20 AU in 2038-39, 13.31 AU on 2044-12-31; minimum Earth distance 10.08 AU.
  M(t0) = 57.5 deg, M(t_end) = 271.4 deg. Retrograde on top of that.
- **ID 144 (1999 XS35): a = 17.81 AU, e = 0.9470, i = 19.6 deg, q = 0.944 AU, Q = 34.7 AU, P = 75.2 yr.**
  Last perihelion 1999-11-01, next 2075-01-07. It sits at aphelion through the whole window: r = 33.82 (2030), 34.59 (2035), 34.58 (2040),
  33.78 AU (2044-12-31; aphelion 34.68 AU in 2037-38); minimum Earth distance 32.9 AU. M(t0) = 144.4 deg, M(t_end) = 216.3 deg.

Both have Earth-crossing orbits (MOID 0.026 / 0.0047 AU) but are on the far side of their orbits for all 15 years. Each costs exactly 1 unit
as a miss; **the practical target set is 298 objects and N_miss >= 2** (cost floor 2 before any spacecraft cost).

## 4. Target-opportunity density seen by a spacecraft near 1 AU

- **Ecliptic crossings inside 0.8-1.6 AU: 2731 over 15 yr = 182 per year** (per-year std 4.2), from 292 asteroids; per year (2030..2044):
  185, 187, 181, 173, 188, 185, 178, 186, 184, 178, 180, 180, 187, 177, 182.
  Of these, **1160 (77 per year) belong to i > 20 deg objects** (per year 87, 75, 74, 79, 72, 81, 79, 72, 82, 80, 68, 79, 76, 83, 73) — the events
  that *must* be used for the hard class. By class: easy 878, moderate 693, hard 1160.
- Ecliptic crossings within 0.9-1.1 AU of the Sun: **1783 (119 per year), from 227 asteroids**; 894 of them from 99 of the 104 hard-class objects —
  nearly every hard target offers several near-1-AU crossings.
- Access windows (r in band and |z| <= 0.1 AU): **3220 windows = 215 per year** (per year 252, 213, 215, 199, 221, 209, 220, 206, 218, 217, 204,
  213, 217, 211, 205), median length 38.5 d; 174 000 asteroid-days of accessibility in total (32 asteroids accessible at any instant).
- Where in the ecliptic: in-band nodal crossings are spread uniformly in ecliptic longitude (`polar_nodal_crossings.png`); no preferred quadrant.
- Kinematics at in-band nodes: |v_z| median 8.1 km/s (quartiles 4.3-12.8), speed relative to Earth at the crossing median 43 km/s (10 % below
  15 km/s). Only 50 in-band crossings occur within 0.1 AU of Earth (121 within 0.2 AU). These do not affect flyby feasibility (position match
  only) but rule out rendezvous-style hopping; flybys must be sequenced as a near-1-AU orbit threading a stream of ~180-215 timing windows per year.

## 5. Take-aways for mission design

1. 298 geometrically reachable targets; 131 and 144 are hopeless (N_miss >= 2).
2. The population is uniformly PHA-like: q ~ 0.85 AU, MOID < 0.052 AU for all 300. Every target orbit passes within 0.05 AU of Earth's orbit, so a
   spacecraft on a slightly modified Earth-like orbit (v_inf <= 4 km/s gives roughly 0.75-1.35 AU radial reach) intersects every target orbit
   somewhere; the whole problem is timing/phasing, not energy.
3. Timing density: ~182 in-band nodal crossings/yr + ~215 access windows/yr, uniform over the 15 years and over ecliptic longitude. A spacecraft
   living 15 yr near 1 AU sees ~3000 windows; sequencing them with 0.25-0.83 mm/s^2 of acceleration is the real optimisation problem.
4. Inclination is the only real difficulty axis: 120 targets with i <= 10 deg can be reached anywhere on long (median 82 d) windows; 104 targets
   with i > 20 deg need node-timed encounters (median 24 d windows, median 8 in-band crossings each over 15 yr; the 23 with i > 45 deg have 2-28).
5. Earth close approaches (< 0.1 AU: 86 events / 73 asteroids) are candidates for near-ballistic (cost ~1) spacecraft, but relative speeds are
   6.5-32 km/s, so such spacecraft must be launched to intercept; two flybys per ballistic spacecraft remains the generic limit.

## 6. Verification (`verify.py`, numbers in `verify.json`)

1. PDF sample lines (sec. 6.2): Earth at t = 0 reproduces line 1 to **0.5 m**; the sample launch v_inf evaluates to **4.000000 km/s**;
   asteroid 174 at t = 1.0008479152e7 s is **9.5 m** from the sample flyby position (line 2). Epochs and angle conventions confirmed.
2. Against the lead's independent `ctoc14/kepler.py` (62 epochs x 300 asteroids + Earth): max position difference **4.4e-5 km**, max velocity
   difference 6.5e-11 km/s.
3. Kepler-equation residual over the full grid: 8.9e-16 rad.
4. Brute-force recount on a 0.1-d grid (54 789 samples): perihelion counts and ecliptic-crossing counts agree for **all 300** asteroids
   (2146 and 4293 in total); grid r_min exceeds the analytic r_min by 0-3.4e-6 AU (as it must); grid Earth-distance minimum exceeds the refined
   one by 0-9.2e-6 AU with epoch differences <= 0.05 d.
5. Nodal table: at the 4293 analytic crossing epochs |z| <= 0.22 km, |r - r_node| <= 1e-9 AU; all ascending nodes have v_z > 0, all descending v_z < 0.
6. MOID: independent brute force (2400 x 1200 true-anomaly grid + Nelder-Mead from the 6 best cells) agrees with the stored values to 5e-14 AU for all 300.
7. Re-derived from the CSVs: class counts 120/74/104/0/2, 2731 in-band crossings from 292 asteroids, 86/33/3 Earth approaches, Atira = {296, 297},
   retrograde = {131}, r_min(window) > 2.5 AU = {131, 144}.

Changes relative to the first version of this analysis: the `n_rad_s` column was unreadable (written as 0.000000) and is now written in
scientific notation with a `n_deg_per_day` companion; all CSVs now carry 12 significant digits (was 6 decimals); the statement that the two
MOID > 0.05 AU objects were 131/144 was wrong — they are 117 and 280 (counts unchanged).

## 7. Re-run and second-implementation check (2026-09-05)

- `population_analysis.py` re-run from scratch (16.6 s) and `verify.py` re-run on the regenerated outputs (143 s): every number in
  sections 1-6 and in `summary.json` / `verify.json` is reproduced unchanged.
- `analysis/population/independent_check.py` (log: `independent_check.log`) is a second, independent implementation of the contest
  propagation (own Kepler solver, rotation via argument of latitude, velocities from radial/transverse components, 0.1-d brute-force
  grid) that shares no code with `kepler.py`, `population_analysis.py` or the lead's `ctoc14/`. It reproduces: PDF sample lines
  (Earth at t0 to 0.0005 km, launch v_inf 4.000000 km/s, asteroid 174 to 0.0095 km); **2146 perihelia, 4293 ecliptic crossings,
  2731 in-band crossings from 292 asteroids** with the identical per-year histogram (185, 187, 181, 173, 188, 185, 178, 186, 184, 178,
  180, 180, 187, 177, 182); q_max = 1.045 AU, 271 with q < 1, all 300 crossing 0.7-1.5 AU; inclination bins 120 (i <= 10 deg) / 75
  (10-20 deg, incl. 144) / 105 (> 20 deg, incl. 131) / 1 retrograde; r_min(window) > 2.5 AU only for 131 (10.71 AU) and 144 (33.78 AU),
  which are also the only two with zero perihelion passages; the six closest Earth approaches (243, 65, 102, 226, 174, 20) with the same
  distances and epochs; 73 / 32 / 3 asteroids with d_Earth,min < 0.1 / 0.05 / 0.02 AU; per-asteroid perihelion, crossing and in-band
  crossing counts identical to the CSV for all 300; |r_min| and |d_Earth,min| grid-vs-CSV differences <= 3.4e-6 and 9.2e-6 AU.
- PDF re-read (Tables 2-3, eq. 4-5): epochs MJD 60676.0 (Earth) / 61200.0 (asteroids), t0 = MJD 62502.0, T_max = 4.73364e8 s = 5478.75 d
  and the six Earth elements are exactly those used here.

## 8. Re-run 2026-09-06

- Full pipeline re-executed from scratch (`population_analysis.py` 18.4 s, `verify.py` 87 s, `independent_check.py` 3.5 s; all exit 0).
  The regenerated `summary.json` and `population_stats.csv` are byte-identical to the previous versions (`diff` empty); every line of
  `verify.log` and `independent_check.log` is reproduced (2146 perihelia, 4293 crossings, 2731 in-band from 292 asteroids, classes
  120/74/104/0/2, MOID >= 0.05 AU only for 117 and 280, Earth approaches 86/33/3, PDF sample lines to 0.5 m / 9.5 m / v_inf 4.000000 km/s).
- Additional stand-alone spot check (numpy from `MEA.txt` only, no project imports): Atira = {296, 297}; the 23 objects with i > 45 deg and
  the 10 with a > 3 AU are exactly the lists in sections 1 and 3; the six low-i objects whose two nodes both lie outside 0.8-1.6 AU are
  18 (3.54/0.76 AU), 127 (2.19/0.77), 192 (0.79/2.44), 211 (0.77/3.97), 267 (0.73/0.44), 292 (0.79/1.60); inclination bins 120/75/105/1;
  ID 131 r = 10.711/15.331/16.069/13.313 AU at 2030-01-01/2035-01-01/2040-01-01/2044-12-31, M = 57.5 -> 271.4 deg, perihelia MJD 61029.4
  (2025-12-20) and 70249.9 (2051-03-19); ID 144 r = 33.819/34.592/34.580/33.784 AU, M = 144.4 -> 216.3 deg, perihelia MJD 51483.3
  (1999-11-01) and 78944.9 (2075-01-07). Both remain unreachable for the entire window.
