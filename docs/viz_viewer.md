# CTOC14 Fleet Viewer: the page (docs/viz_spec.md section 3), built 2026-09-24

## Files
- `results/viz/index.html` (title "CTOC14 Fleet Viewer"), `results/viz/css/viewer.css` (dark theme, colour tokens on `:root`,
  phone layout below 760 px), `results/viz/js/`:
  - `kepler.js` Kepler propagation (below), Moon mean orbit, MJD -> UTC.
  - `data.js` data loading with fixture fallback, fleet preparation (typed arrays, first-detection table), binary-search interpolation.
  - `scene.js` Sun (point light + additive glow sprite), ecliptic circles at 0.5/1/1.5/2 AU + vernal-equinox line, Earth orbit +
    sphere, Moon on its geocentric orbit (Earth-Moon zoom), asteroid orbit lines (one `LineSegments` for reachable, red dashed
    `LineLoop`s for unreachable), asteroid `Points` with per-vertex colour (grey / green flown / red unreachable), OrbitControls.
  - `fleet.js` per-craft layer: faint full path, bright flown trail (`drawRange`), current partial segment, thrust arcs
    (`LineSegments` of samples inside `thrust` intervals, whitened craft colour, drawn up to the current time; faint copy for the
    future), craft marker sprite, flyby flash rings (camera-facing, size proportional to camera distance, 1.4 s).
  - `ui.js` panel, HUD, per-craft table, timeline canvas (one row per craft: launch-end bar, thick thrust bars, white flyby ticks,
    cursor; click/drag scrubs), HTML label overlay.
  - `main.js` clock, fleet switching, flyby-event detection, render loop, `window.__viz` handle for automated checks.
- Fixture (never collides with the real files): `results/viz/data/fixture/{bodies.json,index.json,fleets/demo.json}` written by
  `results/viz/data/fixture/make_fixture.py` (2 craft, 200 d, 3 flybys on Lambert arcs sampled with
  `ctoc14.kepler.propagate_twobody`; asteroids 174, 162, 181 flown, 131 unreachable, 5 extra; max flyby position error
  1.08e-6 AU from the 6-significant-digit rounding). `results/viz/data/fixture/test300/` = the same demo fleet on a snapshot of
  the real 300-asteroid `bodies.json` (performance test, `?data=data/fixture/test300/`).

## External resources (pinned, nothing else)
- `https://cdn.jsdelivr.net/npm/three@0.170.0/build/three.module.js` (import map name `three`)
- `https://cdn.jsdelivr.net/npm/three@0.170.0/examples/jsm/` (import map prefix `three/addons/`, used for
  `controls/OrbitControls.js` only). 0.170.0 is the last single-file module build (0.171+ splits off three.core.js).

## Data selection
`data/index.json` + `data/bodies.json` are tried first; if either is missing the page falls back to `data/fixture/` and logs
`[viewer] loaded FIXTURE data ...` (and shows a status banner). `?data=<dir>/` forces a directory. Fleet `file` paths resolve
relative to the loaded index. Optional per-craft `m` array (same length as `t`) switches the HUD fuel from "loaded at launch" to
"used, m0 - m(t)".

## Kepler propagation in JS
`docs/viz_kepler.md` did not exist when the viewer was started (it landed at 16:46, while kepler.js was being written), so
`kepler.js` was written to mirror `ctoc14/kepler.py` `Body.state`; it follows viz_kepler.md steps 1-5 literally (the doc calls
the +-1 rad Newton clip optional; kepler.js keeps it, tolerance 1e-13, and its Moon uses the equivalent per-day mean motion).
Worked examples of viz_kepler.md section 4 reproduced by kepler.js: asteroid 1 at t = 1000 d max component error 4.8e-13 AU
(doc prints 11 digits), Earth at t = 0 and 1000 d <= 4.4e-10 AU (doc prints 9 digits), velocity identical to 9 decimals.
Formulas: `M = M0 + n (t + t_off)` with `n = sqrt(mu / a_km^3)` rad/s and `t_off = (62502 - epoch_mjd) * 86400` s; Newton on
`E - e sin E = M` (M reduced to [0, 2pi), start `M + e sin M` for e < 0.8 else pi, step clipped to +-1, tol 1e-13, <= 60 it);
`x_p = a (cos E - e)`, `y_p = a sqrt(1-e^2) sin E`; `r = R3(-Om) R1(-i) R3(-om) [x_p, y_p, 0]`, km -> AU.
Measured: node cross-check against `ctoc14.kepler.Ephemeris` for asteroids 1, 3, 4, 12, 57, 100, 131, 200, 250, 300 + Earth at
t = 0, 365.25, 1234.5, 3000, 5478.75 d: max |JS - python| = 1.52e-14 AU (55 positions).
Moon: illustrative geocentric mean orbit from the `moon` block (period 27.321582 d, Om/om linear rates), true distance; Earth and
Moon spheres are enlarged (Earth 0.012 AU normal / 3e-4 AU in zoom, Moon 1e-4 AU); the legend says so in zoom mode. The Moon is
hidden outside zoom mode because it would sit inside the enlarged Earth marker.

## Controls / HUD
Fleet selector (index.json order), play/pause (space), speed 1-200 d/s (log slider, default 30), +-10 d, reset, time slider and
timeline scrub, UTC date + t readout, per-craft show/hide (table checkbox), follow craft (menu or "f" button), Earth-Moon zoom,
orbit lines on/off, labels on/off + asteroid IDs, home/top view. HUD: title/kind/source/note, craft count, launched, covered so
far / reachable (= asteroids - unreachable), fuel, raw J, sum J_i, per-craft flybys so far / total and fuel.

## Performance (measured in Chrome 2026-09-24, MacBook, main-thread time per frame incl. render call, tab hidden so driven by hand)
- fixture (5 asteroids): 4.09 ms/frame.
- 300 asteroids (test300; 300 Kepler solves + 107,280 orbit-line vertices + fleet + HUD + timeline): 4.23 ms; with 300 asteroid
  ID DOM labels 7.89 ms; Earth-Moon zoom 2.89 ms. All well inside 16.7 ms (60 fps).
- Everything is BufferGeometry; asteroid positions and colours are two dynamic attributes updated per frame; the trail uses
  `drawRange` (no per-frame geometry rebuild); craft position by binary search + linear interpolation; timeline statics are cached
  in an offscreen canvas and only the cursor is redrawn. Note: `requestAnimationFrame` is suspended in hidden tabs, so the fps
  readout shows "paused (tab hidden)" there.

## Verified in Chrome (claude-in-chrome), no console errors
Load + fixture fallback log, play/scrub, flown targets turning green, flyby flash at t = 88 d, thrust arcs, follow craft (controls
target == craft position), Earth-Moon zoom (camera 0.011 AU from Earth, Moon at 0.0025 AU), orbit toggle, asteroid IDs, per-craft
hide, top view, phone width 386 px (panel collapses to a menu sheet, no horizontal scroll). A canvas-size feedback loop in the
timeline (flex item growing with its own dpr-scaled intrinsic size) was found and fixed (wrapper + absolute canvas + size cap).

## Known limits
- Trail line width is 1 px (WebGL); thrust arcs are distinguished by a brighter colour, not thickness.
- Flyby flashes fire only when playback crosses the event with a step < 60 d; jumps (slider, +-10 d over an event) recolour the
  target without a flash.
- Labels for 300 asteroid IDs are DOM nodes (cost above); off by default.
