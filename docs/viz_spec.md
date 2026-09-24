# Visualisation spec: 3D interactive solar system with the CTOC14 fleets (2026-09-24)

Goal: one interactive web page (Three.js) showing the Sun, Earth + Moon, all 300 catalogue asteroids on their orbits,
and, for a chosen fleet result, every spacecraft's trajectory animated over the 15-year mission, with flybys flashing
as they happen and a running count. At least five fleet results, representing different fleet configurations
(craft count, coverage, era of the project), selectable from a menu.

All work goes under `results/viz/` (page + data), `results/catalog/` (sorted results), `docs/viz_*.md`.
Nothing under results/ is deleted or modified by this work. Python: `~/.venvs/astro313/bin/python`.

## 1. Frames, units, epochs (the same for every file)
- Heliocentric Ecliptic J2000 (the contest frame). Positions in **AU**, times in **days since mission start**
  t0 = 2030-01-01 00:00 UTC = MJD 62502.0 (t0 = 0 d, mission end = 5478.75 d).
- Sources: MEA.txt (asteroid Keplerian elements, epoch MJD 61200.0 = 2026-06-09; columns ID a[AU] e i Ω ω M[deg]),
  Earth elements from docs/problem_summary.md (epoch MJD 60676.0), μ☉ = 1.32712440018e11 km³/s², AU = 149597870.7 km.
  ctoc14/kepler.py (load_mea, propagate_twobody) and ctoc14/constants.py already implement these.
- Submission files results/CTOC14_Result_*.txt: 15 columns `Line SC Event Time[s] x y z[km] vx vy vz[km/s] m[kg]
  Tx Ty Tz[N] AstID`; Event 0 launch, 1 thrust sample, 2 coast node, 3 flyby, 4 end.
- Twin fleets (not exported to submission format): `route_*.npz` in run_ialns "ist" format (tL, vinf, ts, Ts, tf,
  asts); load with `tools/run_ialns.ipr(st)` -> `ctoc14.impulsive.ImpulsiveProblem` (ip.integrate() gives the state
  at every impulse node; Kepler arcs between nodes -> sample with ctoc14.kepler.propagate_twobody). Flyby epochs ip.tf,
  targets ip.asts; launch ip.tL, v_inf ip.vinf; tank ip.tank().

## 2. Data files (JSON, written by the converters)
### `results/viz/data/bodies.json`
```
{ "epoch_mjd_t0": 62502.0, "units": {"length": "AU", "time": "day", "angle": "deg"},
  "mu_sun_km3s2": 1.32712440018e11, "au_km": 149597870.7,
  "earth": {"a": 1.0009175020, "e": 0.017566762041, "i": 0.002976847126, "Om": 189.953211282428, "om": 273.196254000254,
            "M0": 357.4135031077, "epoch_mjd": 60676.0},
  "moon":  {"note": "illustrative geocentric mean orbit, not part of the contest model",
            "a_km": 384400, "e": 0.0549, "i": 5.145, "Om": 125.08, "om": 318.15, "M0": 135.27, "epoch_mjd": 51544.5,
            "period_d": 27.321582, "Om_rate_deg_per_day": -0.05295, "om_rate_deg_per_day": 0.11140},
  "asteroids": [ {"id": 1, "a": ..., "e": ..., "i": ..., "Om": ..., "om": ..., "M0": ..., "epoch_mjd": 61200.0}, ... 300 ],
  "unreachable": [131, 144] }
```
Angles in degrees. The page propagates Kepler orbits itself (solve Kepler's equation), so no asteroid ephemeris tables.

### `results/viz/data/index.json`
```
{ "fleets": [ {"key": "t10d", "file": "fleets/t10d.json", "title": "...", "n_craft": 10, "covered": 298, "J_raw": 14.366,
               "shown": 13.8019, "submitted": "2026-09-17", "kind": "submission|twin", "note": "..."}, ... ] }
```
Ordered by date. At least 5 entries, chosen by the CATALOG agent (section 4).

### `results/viz/data/fleets/<key>.json`
```
{ "key": "t10d", "title": "t10d: 10 craft, 298 targets, raw J 14.366 (submitted 09-17)",
  "source": "results/CTOC14_Result_t10d.txt", "kind": "submission", "n_craft": 10, "covered": 298,
  "missed": [131, 144], "J_raw": 14.366055, "sumJi": 12.366055, "date": "2026-09-17",
  "craft": [ { "id": 1, "m0": 948.7, "launch_d": 41.0, "end_d": 5438.0, "n_flybys": 42, "fuel_kg": 348.7,
               "t": [41.0, 43.0, ...],                      # sample times, days, strictly increasing, 2 d step (1 d near flybys), ~2000-3000 points
               "xyz": [[x,y,z], ...],                       # AU, same length as t, 6 significant digits
               "thrust": [[t_on, t_off], ...],              # days; contiguous Event=1 runs (submission) or impulse nodes ±1 d (twin)
               "flybys": [{"t": 86.1, "ast": 57}, ...] }, ... ] }
```
Every flyby epoch must appear exactly in `t` (insert the sample), and the sample position there must match the asteroid's
position from bodies.json to better than 2e-4 AU (the CHECKER verifies this). Size target: <= 1.5 MB per fleet file.

## 3. The page: `results/viz/index.html` (+ `results/viz/js/*.js`, `results/viz/css/*.css`)
- Three.js and OrbitControls loaded ONLY from cdn.jsdelivr.net/npm/ or cdnjs.cloudflare.com (pinned versions); everything
  else inline or local. Must run from a plain static server (`python -m http.server` in results/viz) and offline except
  for the CDN. `<title>` = "CTOC14 Fleet Viewer". Dark theme, colour tokens on `:root`, works at phone width.
- Scene: Sun (point light + glow), ecliptic grid/reference circles at 0.5/1/1.5/2 AU, Earth orbit line + Earth sphere,
  Moon on its geocentric orbit (toggle "Earth–Moon zoom": camera follows Earth at ~0.01 AU scale so the Moon is visible;
  the Moon's orbit is exaggerated only in that mode if needed, say so in the legend), 300 asteroid orbit lines (faint,
  reachable = grey, the two unreachable = red dashed) and asteroid points moving in time; targets already flown by the
  current fleet turn green, remaining ones stay grey.
- Fleet layer: one colour per craft (distinct palette, up to 12 craft), trajectory drawn up to the current time
  (trail) with the future part faint; craft marker; thrust arcs highlighted (brighter/thicker) on the trail; flyby
  event: flash ring at the asteroid + the target turns green + a counter increments; a timeline at the bottom with all
  flyby ticks per craft.
- Controls: fleet selector (from index.json), play/pause, speed (1–200 days per second), time slider, current date
  (UTC) readout, per-craft show/hide, "follow craft N" camera, Earth–Moon zoom, orbit lines on/off, labels on/off.
- HUD: fleet title, craft count, covered/298 so far, total fuel so far (interpolate m along t if available or sum
  fuel per craft at launch), raw J, and a per-craft table (flybys so far / total, fuel).
- Performance: use BufferGeometry, update positions per frame from Kepler propagation (300 asteroids + Moon) and
  trajectory interpolation (binary search in t); target 60 fps on a laptop.

## 4. Agents and deliverables
- CATALOG: inventory every results/CTOC14_Result_*.txt (with results/validator2_*.log) and every twin fleet directory
  (results/s16/best/fleet, results/s18/*/fleet, results/s15b/best_closed/fleet, results/newgen/*fleet*, ...); for each
  record craft count, covered, misses, raw J, sum J_i, validity (from the validator log or tools/validator2.py if
  missing; twin fleets: from fleet.json), date, provenance (stage), md5. Write results/catalog/fleets.json and
  docs/results_catalog.md (sorted by date, with a "keep/derived/superseded" tag per file), and results/catalog/README.md
  explaining the layout. Then pick >= 6 fleets for the viewer that differ in configuration: the earliest valid one,
  the one with most craft, t10d (10 craft, first banked), s15a or s16E (9 craft, submitted), s16a (best valid),
  one 8-craft RHFA twin fleet (results/s18/p13_claims_ttl4, 246/298) and the 9-craft RHFA twin (p19, 262/298);
  write the selection with titles/notes to results/catalog/viz_selection.json.
- BODIES: write results/viz/data/bodies.json (section 2) from MEA.txt + docs/problem_summary.md; add a Python
  check that the page's Kepler propagation formula (documented in docs/viz_kepler.md with the exact equations and
  angle conventions the JS must use) reproduces ctoc14.kepler.propagate/ephemeris positions to < 1e-6 AU for 10
  asteroids at 5 epochs.
- FLEETS (after CATALOG): write results/viz/data/fleets/<key>.json for every selected fleet and index.json;
  submission files: parse, thin to the 2 d grid + flyby epochs, thrust arcs from Event=1 runs; twin fleets: propagate
  Kepler arcs between impulse nodes. Report sizes and the max flyby position error.
- VIEWER: build the page (section 3) against a self-made fixture fleet (2 craft, 5 asteroids, 200 d) so it does not
  wait for real data; then, when data lands, switch to real files.
- CHECKER: verify every fleet file (flyby positions vs asteroid positions from bodies.json < 2e-4 AU; t strictly
  increasing; sizes), serve the page locally, open it in Chrome (claude-in-chrome MCP tools via ToolSearch) or a
  headless check if Chrome is unavailable, screenshot each fleet at 3 epochs, list every console error, and fix what is
  broken (any file). Write results/viz/README.md (how to run) and docs/viz_report.md (what was built, screenshots list,
  known limits).
