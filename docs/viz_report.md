# CTOC14 fleet visualisation: checker / integration report (2026-09-24)

Deliverable of `docs/viz_spec.md` section 4 "CHECKER". Everything below is measured; paths are relative to
`/Users/mickey/solarsystem/ctoc14`. Run instructions and controls: `results/viz/README.md`.

## 1. What was built

One static web page, `results/viz/index.html` (+ `css/viewer.css`, `js/*.js`), served with `python -m http.server`
from `results/viz/`. Three.js 0.170.0 and OrbitControls come from cdn.jsdelivr.net (pinned, the only external
resource). The scene has the Sun, Earth on its contest orbit, the Moon on an illustrative geocentric mean orbit
(Earth-Moon zoom mode), 300 asteroid orbit lines and asteroid points propagated by Kepler's equation in JS
(`js/kepler.js`, a transcription of `ctoc14/kepler.py` documented in `docs/viz_kepler.md`), and a fleet layer with one
colour per craft: full path faint, flown trail bright, thrust arcs whitened, craft marker, flash ring at each flyby,
flown targets turning green, HUD counters and a per-craft timeline. Data files: `results/viz/data/bodies.json` (Earth,
Moon, 300 asteroids from MEA.txt), `results/viz/data/index.json` and nine `results/viz/data/fleets/<key>.json`
(6 submission files thinned to a 2-3 d grid with every flyby epoch as an exact sample, 3 impulsive twin fleets
resampled on Kepler arcs between the impulse nodes). Converters and checkers live next to the data
(`make_bodies.py`, `make_fleets.py`, `check_bodies.py`, `check_fleets.py`). Agent docs: `docs/viz_kepler.md`,
`docs/viz_fleets.md`, `docs/viz_viewer.md`; results inventory: `results/catalog/` (CATALOG agent).

## 2. The nine fleets (results/viz/data/index.json, oldest first)

| key | kind | craft | covered | flyby events | raw J | sum J_i | fuel kg | date | submitted / shown | file MB | what it represents |
|---|---|---|---|---|---|---|---|---|---|---|---|
| b300 | submission `results/CTOC14_Result_b300.txt` | 10 | 277 | 277 | 49.196882 | 26.196882 | 11002.3 | 09-06 | no | 0.987 | earliest validator2-VALID file: sequential full-tank tours (2000 kg craft), 23 misses, the "chaser" era |
| v4 | submission `results/CTOC14_Result_v4.txt` | 12 | 298 | 305 | 32.259473 | 30.259473 | 12513.0 | 09-07 | 09-07 19:35, shown 29.876511 | 1.133 | first 298-cover and first leaderboard entry; full-tank chasers + mop-up craft |
| mix1 | submission `results/CTOC14_Result_mix1.txt` | 18 | 295 | 491 | 43.843002 | 38.843002 | 15557.0 | 09-10 | no | 1.399 (3 d grid) | most craft among VALID files: greedy mix of swarm + chaser fragments |
| t10d | submission `results/CTOC14_Result_t10d.txt` | 10 | 298 | 298 | 14.366055 | 12.366055 | 2741.8 | 09-17 | 09-17 12:04, shown 13.801851 | 1.193 | impulsive fleet search (relocate + dissolve); first banked score, rank 8 |
| s16E | submission `results/CTOC14_Result_s16E.txt` | 9 | 298 | 298 | 13.746635 | 11.746635 | 3076.4 | 09-23 | 09-23 17:05, shown 13.511581 | 1.081 | closed 9-craft fleet after iterated MILP moves + swaps + polish; the file on the board (rank 10) |
| s16a | submission `results/CTOC14_Result_s16a.txt` | 9 | 298 | 298 | 13.739947 | 11.739947 | 3070.5 | 09-23 | no (withheld) | 1.081 | same target sets as s16E, tighter conversions; best raw J of the project |
| p13 | twin `results/s18/p13_claims_ttl4/fleet` | 8 | 246 | 246 | 64.862766 | 10.862766 | 3123.2 | 09-24 | no | 1.010 | stage 18 rolling-horizon fleet auction with plan claims (ttl 4); best 8-craft RHFA fleet |
| p19 | twin `results/s18/p19_claims_N9/fleet` | 9 | 262 | 262 | 50.292377 | 12.292377 | 3574.6 | 09-24 | no | 1.114 | RHFA 9-craft diagnostic: coverage scales with N at ~29 flybys per craft |
| g2 | twin `results/s18/g2_plan/fleet` | 8 | 224 | 224 | 86.940599 | 10.940599 | 3189.3 | 09-24 | no | 1.005 | RHFA first full run with the planbeam proposer (no claims): the lane-coherence loss that claims fixed |

Total 10.003 MB of fleet data; index.json 5753 bytes; bodies.json 48902 bytes.

## 3. Data checks (all PASS)

- `nice -n 5 ~/.venvs/astro313/bin/python results/viz/data/check_fleets.py` -> `RESULT PASS max_flyby_err_au=5.709694e-06`,
  exit 0. Per fleet max flyby position error (sample position vs asteroid position recomputed with
  `ctoc14.kepler.Body` from bodies.json elements), tolerance 2e-4 AU: b300/v4/mix1 4.992e-06, t10d 5.402e-06,
  s16E 5.595e-06, s16a 5.124e-06, p13 5.710e-06 (craft 2, ast 270, t 1990.677055 d = 854 km), p19 4.945e-06,
  g2 4.789e-06. The residual is the models' own flyby miss (validator tolerance 1000 km) plus 6-digit rounding.
- `nice -n 5 ~/.venvs/astro313/bin/python results/viz/data/check_bodies.py` -> `RESULT PASS max_check_err_au=2.041120e-15`,
  exit 0: the documented JS formula vs `ctoc14.kepler.Ephemeris`, 10 asteroids x 5 epochs max 2.04e-15 AU, Earth 0.0,
  all 300 asteroids x 110 epochs max 1.07e-14 AU.
- Independent re-parse by the checker (own script, no shared code): index.json lists 9 fleets (>= 5 required), every
  fleet file parses, every `t` strictly increasing (0 violations in 225,733 samples), every flyby epoch is an exact
  element of `t` (0 missing out of 2699 flyby events), `xyz` and `m` lengths match `t`, every thrust arc lies inside
  [launch, end], the distinct-target count equals `covered` in the file and in index.json for all 9, every craft ends
  <= 5478.75 d, every file <= 1.5 MB (largest mix1 1.399 MB).
- In the browser, at mission end every fleet's HUD equals index.json: covered 277/298/295/298/298/298/246/262/224,
  fuel 11002.3/12513.0/15557.0/2741.8/3076.4/3070.5/3123.2/3574.6/3189.3 kg, raw J and sum J_i to 4 decimals; the
  number of green asteroid points equals `covered` for each fleet.

## 4. Browser check (Chrome, claude-in-chrome, 2026-09-24 17:08-17:18)

Served `results/viz` with `python -m http.server 51203` (free port picked at random; server stopped and port freed
afterwards; the tab was closed). Console after load: `[viewer] loaded REAL data (results/viz/data/): 9 fleet(s),
300 asteroid(s)` and `[viewer] fleet b300: 10 craft, 277 flyby events, 22310 samples, tmax 5443.37 d`; **zero console
errors or warnings** after load, after switching through all nine fleets (50-56 ms each), during playback, slider
jumps, Earth-Moon zoom, follow-craft, per-craft hide, orbit/label toggles and the screenshot runs (checked with
`read_console_messages` errors-only at the end: none).

Measured while the tab was visible (the tab has to be foreground: `requestAnimationFrame` is suspended in hidden
tabs and the page then says so): 120 fps (the display's refresh rate) with the real b300 fleet playing at 200 d/s;
the clock advanced 600.0 d in 3.0 s wall time; 31 flyby flashes fired on that stretch (t 1000 -> 1600 d). Time
slider: t = 3000 d -> `2038-03-20 00:00 UTC`, covered 151/298 for b300. Earth-Moon zoom: camera 0.011 AU from Earth,
near plane 2e-5, Moon at 395,854-401,670 km geocentric (mean orbit), legend note shown. Follow craft 3: OrbitControls
target coincides with the craft position (distance 0.000000 AU), marker visible. Hide craft 1: trail node hidden, table
row greyed. Orbit lines off: group hidden. The "N fps" printed in the corner of the screenshots below (1-8 fps) is
the readout during the screenshot tool's own blocking capture, not the running frame rate.

### Screenshots (`results/viz/screenshots/`, 1232x960 JPEG, all nine fleets at three epochs + three extra views)

| fleet | t = 600 d (2031-08-24) | t = 2800 d (2037-09-01) | mission end |
|---|---|---|---|
| b300 | `b300_t0600.jpg` (25/298) | `b300_t2800.jpg` (139/298) | `b300_t5443_end.jpg` (277/298, 11002.3 kg) |
| v4 | `v4_t0600.jpg` (26/298) | `v4_t2800.jpg` (153/298) | `v4_t5475_end.jpg` (298/298, 12513.0 kg) |
| mix1 | `mix1_t0600.jpg` (41/298, 18 craft) | `mix1_t2800.jpg` (202/298) | `mix1_t5475_end.jpg` (295/298, 15557.0 kg) |
| t10d | `t10d_t0600.jpg` (14/298) | `t10d_t2800.jpg` (142/298) | `t10d_t5457_end.jpg` (298/298, 2741.8 kg) |
| s16E | `s16E_t0600.jpg` (21/298) | `s16E_t2800.jpg` (147/298) | `s16E_t5475_end.jpg` (298/298, 3076.4 kg) |
| s16a | `s16a_t0600.jpg` (21/298) | `s16a_t2800.jpg` (147/298) | `s16a_t5475_end.jpg` (298/298, 3070.5 kg) |
| p13 | `p13_t0600.jpg` (27/298) | `p13_t2800.jpg` (144/298) | `p13_t5469_end.jpg` (246/298, 3123.2 kg) |
| p19 | `p19_t0600.jpg` (30/298) | `p19_t2800.jpg` (157/298) | `p19_t5453_end.jpg` (262/298, 3574.6 kg) |
| g2 | `g2_t0600.jpg` (28/298) | `g2_t2800.jpg` (130/298) | `g2_t5472_end.jpg` (224/298, 3189.3 kg) |

Extra: `b300_t1500_earth_moon_zoom.jpg` (Earth-Moon zoom, Moon on its orbit ring), `t10d_t1210_follow_craft1.jpg`
(camera following craft 1), `t10d_t3300_top_view.jpg` (top view, 10 trails at 172/298).

## 5. Fixes made by the checker (results/viz only)

1. `js/ui.js` per-craft fuel column mixed two meanings: launched craft showed fuel used so far, unlaunched craft
   showed their total loaded fuel (b300 at t = 0 read "0.0 / 1357.7 / 0.0 / 1382.3 ..."). Now every row reads
   "used / total" kg (0 before launch); header renamed in `index.html`.
2. `js/ui.js` total fuel read "-0.0 kg" at t = 0 (m0 to 6 decimals minus a mass column rounded to 0.01 kg); used
   fuel is clamped at 0.
3. `js/ui.js` twin fleets ended 1.5 kg per craft below their catalog fuel (HUD 3111.2 vs index 3123.2 kg for p13,
   3561.1 vs 3574.6 for p19, 3177.3 vs 3189.3 for g2): the twin mass arrays end at the 601.5 kg margin while
   `fuel_kg = tank - 600` (docs/viz_fleets.md). A finished craft (t >= end) now reports its accounted `fuel_kg`, so
   the end-of-mission HUD equals the catalog; during flight m(t) is still used. The "fuel used (m(t))" label no
   longer flips to "fuel loaded" once all craft are finished.
4. `js/scene.js` Earth-Moon zoom toggled while paused left the Moon at the Earth's centre (inside the enlarged Earth
   sphere, geocentric distance 0 km) until the clock next changed, because the Moon is only propagated in zoom mode
   inside `update(t)`. The scene now remembers the last time and re-propagates when zoom is switched on (measured
   395,854 km right after the toggle).
5. `js/ui.js` timeline height cap raised 180 -> 230 px on desktop so 18 craft (mix1) keep 9 px rows instead of 6.7 px.

No data file needed a change. Nothing under `results/` other than `results/viz/` was written; `tools/`, `ctoc14/`
and `docs/viz_spec.md` were not modified. Files written by the checker: `results/viz/README.md`,
`results/viz/screenshots/*.jpg` (29 files, 6.0 MB), `docs/viz_report.md`, plus the edits above in
`results/viz/index.html`, `results/viz/js/ui.js`, `results/viz/js/scene.js`.

## 6. Known limits

- WebGL draws 1 px lines: thrust arcs are distinguished by a brighter (whitened) colour and by the thick bars in the
  timeline, not by line width in the 3D view.
- Flyby flashes fire only when playback crosses the event with a step < 60 d; slider jumps recolour targets and
  update counters without a flash.
- The browser caches the ES modules: after editing `js/*.js` hard-reload (Cmd+Shift+R).
- `requestAnimationFrame` is suspended in background tabs: playback only advances while the tab is visible.
- Submission-file trajectories are linear interpolations of the file's own rows on a 2 d grid (3 d for mix1; 1 d
  within +-4 d of a flyby); twin trajectories are the impulsive model (velocity jumps at the nodes) and their thrust
  arcs are the equivalent-burn picture described in docs/viz_fleets.md, not flown low-thrust samples; twin `m` is
  piecewise constant between nodes.
- The Moon is an illustrative mean orbit (not part of the contest model), only shown in Earth-Moon zoom; Earth and
  Moon spheres are enlarged for visibility (the legend says so).
- Asteroid ID labels are 300 DOM nodes (about +4 ms per frame when on); off by default.
- Screenshots were taken at a 1400x857 CSS viewport (device pixel ratio 2, canvas 2800x2182); phone layout (<= 760 px)
  was verified by the VIEWER agent at 386 px, not re-checked here.
