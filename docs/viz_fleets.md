# Fleet files for the 3D viewer (FLEETS converter, 2026-09-24)

Deliverable of docs/viz_spec.md section 2 / section 4 "FLEETS": `results/viz/data/fleets/<key>.json` for the nine fleets
in `results/catalog/viz_selection.json`, plus `results/viz/data/index.json`. Scripts (both run from the project root with
`nice -n 5 ~/.venvs/astro313/bin/python ...`):

- `results/viz/data/make_fleets.py [key ...]` -- the converter (6.4 s for all nine). Reads the submission files and the
  twin `route_*.npz` directories, imports `ctoc14.kepler`, `ctoc14.impulsive`, `ctoc14.constants` (nothing under tools/ or
  ctoc14/ is modified). Writes only `results/viz/data/fleets/*.json` and `results/viz/data/index.json`.
- `results/viz/data/check_fleets.py [key ...]` -- the checker: recomputes every asteroid's position at every flyby epoch
  with `ctoc14.kepler.Body` from the elements in `results/viz/data/bodies.json` (falls back to MEA.txt if that file is
  missing), reports the max |error| in AU per fleet, and checks t strictly increasing, array lengths, thrust arcs inside
  [launch, end], flyby / covered / missed / J_raw consistency, file size <= 1.5 MB and the index.json entries.
  Exit 0 iff everything passes.

## Measured result (check_fleets.py, 2026-09-24 17:10)

```
  key       kind craft flybys  cov samples  size MB maxdt  max flyby err AU  worst (craft, ast, t)  status
 b300 submission    10    277  277   22310    0.987   3.4         4.992e-06  (2, 6, 1330.149297)  PASS
   v4 submission    12    305  298   25704    1.133   5.0         4.992e-06  (1, 6, 1330.149297)  PASS
 mix1 submission    18    491  295   31728    1.399   5.0         4.992e-06  (1, 6, 1330.149297)  PASS
 t10d submission    10    298  298   27361    1.193   2.0         5.402e-06  (3, 296, 5230.062489)  PASS
 s16E submission     9    298  298   24762    1.081   2.0         5.595e-06  (9, 142, 5381.526813)  PASS
 s16a submission     9    298  298   24762    1.081   2.0         5.124e-06  (2, 54, 5334.070994)  PASS
  p13       twin     8    246  246   22376    1.010   2.0         5.710e-06  (2, 270, 1990.677055)  PASS
  p19       twin     9    262  262   24664    1.114   2.0         4.945e-06  (7, 159, 3273.158605)  PASS
   g2       twin     8    224  224   22166    1.005   2.0         4.789e-06  (7, 263, 3337.621755)  PASS
9 fleets, 10.003 MB total; max flyby position error 5.710e-06 AU (854.2 km); RESULT PASS
```

The flyby error is dominated by the models' own flyby miss (validator tolerance 1000 km = 6.7e-6 AU for the submission
files, <= 93 km for the twins) plus the 6-significant-digit rounding of the sample position (<= 1e-6 AU); the spec
tolerance is 2e-4 AU. Every fleet's `covered`, `n_flybys`, `J_raw` and `sumJi` agree with `viz_selection.json` to 1e-6
(the converter flags any mismatch on its output line; none occurred). Twin tanks reproduce fleet.json exactly
(sum tank - 600 n = fleet.json `fuel` to 1e-12 kg).

File sizes (bytes): b300 986911, v4 1133487, mix1 1398655 (3 d grid), t10d 1193167, s16E 1081016, s16a 1081111,
p13 1009724, p19 1114093, g2 1004894; index.json 5753.

## Conversion rules

Common: heliocentric ecliptic J2000; positions in AU with 6 significant digits; times in days since MJD 62502.0, integer
grid epochs written as integers, other grid epochs to 1e-4 d, flyby epochs to 1e-6 d (0.09 s: the asteroid moves < 4 km,
so the epoch rounding cannot spoil the flyby check); mass in kg with 2 decimals. JSON is written compact
(`separators=(',', ':')`). Every flyby epoch is an element of `t` (exact float equality after the JSON round trip; the
viewer's `data.js` looks it up that way). If a fleet exceeds 1.5 MB the grid step is raised 2 -> 3 -> 4 -> 5 d; only mix1
(18 craft, 87 836 rows) needed 3 d.

### Submission files (b300, v4, mix1, t10d, s16E, s16a)
- Rows: `Line SC Event Time[s] x y z[km] vx vy vz[km/s] m[kg] Tx Ty Tz[N] AstID`; craft id = the file's SC number.
- Every row is a sample; all six files sample Event=1 at exactly 1 d, with no Event=2 rows (coasts are Event=1 rows
  with ~0 N). Kept samples: the Event 0 row (launch), the Event 4 row (end), every Event 3 row with ITS OWN position
  (the flyby sample is what the validator checked), every 1 d row within +-4 d of a flyby, the last row before any
  sparser stretch of the file, and a greedy 2 d grid (3 d for mix1) elsewhere. The only file gaps > 2 d are one final
  coast per file into a craft's last flyby (3.4 d in craft 5 of b300/v4/mix1, 5.0 d in the 384-row mop-up craft of
  v4/mix1); they are kept as-is (linear interpolation over 5 d at 1 AU is a 1e-3 AU chord error, invisible at scene scale).
- `thrust`: contiguous runs of Event=1 rows with |T| >= 0.05 N (10 % of Tmax); a flyby row between two 1 d samples
  does not break a run; each arc is closed by the next Event=1 row (or the end row). The threshold is recorded in the
  file (`thrust_def`). Why a threshold: with every row Event=1, "contiguous Event=1 runs" would be the whole mission.
  Fraction of samples >= 0.05 N: b300 76 %, v4 76 %, mix1 55 % (full-tank chaser era), t10d 14 %, s16E/s16a 16 %.
  Arcs are 20 d median for t10d/s16 and 108-137 d for the chaser fleets; 99 % of them contain at least one sample pair
  so the viewer can draw them.
- `m0` = mass on the Event 0 row (6 decimals, so sum J_i reproduces the catalog value exactly), `m` = the file's mass
  column at the kept samples, `fuel_kg` = m0 - m(end row), `launch_d`/`end_d` = Event 0 / Event 4 epochs.
- Extra per-craft fields: `n_rows_file`, `thrust_on_days`. Fleet-level extras: `n_flybys` (mix1 has 491 flyby rows for
  295 distinct targets: duplicates count in n_flybys, once in `covered`), `fuel_kg`, `submitted`, `shown`, `md5`, `note`,
  `grid_d`, `near_flyby_1d_window_d`, `thrust_def`, `mass_def`.

### Twin fleets (p13, p19, g2 -- results/s18/*/fleet/route_*.npz, run_ialns "ist" format)
- `ImpulsiveProblem(Ephemeris(), tL, vinf, ts, Ts, tf, asts)` (the same constructor as `tools/run_ialns.ipr`);
  `ip.integrate()` gives the state just before every impulse node (Ys) and at every flyby (Yf). Each Kepler arc is
  resampled with `ctoc14.kepler.propagate_twobody` from the node's post-impulse state (Ys velocity + Ts) at the global
  2 d grid, plus 1 d epochs within +-4 d of each flyby; node epochs and flyby epochs are samples themselves (flyby samples
  = the integrate() states at ip.tf). Consistency check inside the converter: propagating each arc to the next node
  reproduces integrate()'s state there to 0.0 km (`max_arc_err_km` per craft). The craft ends at its last flyby
  (`end_d` = max tf); the 2-3 impulse nodes after it carry |dv| < 1e-21 km/s and are dropped.
- `m0` = ip.tank() = 601.5 exp(dv / ve); the mass drops by the rocket equation at every impulse node (m ends at 601.5, the
  twin's margin); `fuel_kg` = ip.tank() - 600 as specified (sum = fleet.json `fuel`). Extra per-craft fields: `route`
  (npz name), `vinf_kms`, `dv_kms`, `impulses` = [[t_d, |dv| km/s], ...] for every impulse >= 1 m/s.
- `thrust`: one arc per impulse node with |dv| >= 1e-3 km/s (1 m/s; the 20 d bins the optimiser left at 1e-4 km/s
  and below are numerical noise), centred on the node with half-width max(2 d, 10 d x |dv| / dv_cap) capped at 10 d,
  where dv_cap = 0.43 N x 20 d / m is the continuous-thrust capability of one bin (the quantity `ImpulsiveProblem.tcap`
  caps): an impulse at the cap is drawn as its whole 20 d bin, which is what it physically is. The spec's "+-1 d" is the
  floor, raised to the 2 d grid step because the viewer draws an arc only over sample pairs that lie inside it
  (with +-1 d no pair would). Result: 1448 / 1634 / 1562 arcs (p13 / p19 / g2), median 4 d, max 19 d, all visible.

### index.json
`{"generated_by", "frame", "selection", "fleets": [...]}` ordered by date then selection order; each entry carries key,
file, title, n_craft, covered, n_flybys, J_raw, sumJi, fuel_kg, shown, submitted, date, kind, source, priority, grid_d,
size_bytes, note. A partial run (`make_fleets.py t10d`) rewrites only that fleet and keeps the other entries.

## Known limits
- Submission-file samples are the file's own low-thrust integration rows; between kept samples the viewer interpolates
  linearly (2 d chord at <= 1.5 AU: < 2e-4 AU sagitta).
- The twin trajectories are the impulsive model, not a flown low-thrust conversion: velocity jumps at the nodes, and the
  thrust arcs are the equivalent-burn picture described above, not actual thrust samples.
- `m` for twins is piecewise constant between nodes (mass only changes at impulses).
- mix1 is on a 3 d grid (still 1 d within +-4 d of every flyby).
