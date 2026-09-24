# CTOC14 Fleet Viewer

Interactive 3D solar system (Three.js) showing the Sun, Earth + Moon, all 300 catalogue asteroids on their orbits and,
for a chosen CTOC14 fleet result, every spacecraft's trajectory animated over the 15-year mission (2030-01-01 to
2044-12-29), with flybys flashing as they happen and running counters. Nine fleet results are selectable
(10/12/18-craft "chaser" era files, the 10-craft t10d, the 9-craft s16E/s16a, and three 8/9-craft RHFA twin fleets).

## Run

The page is static; it only needs a local web server (ES modules and `fetch` do not work from `file://`) and internet
access to cdn.jsdelivr.net for Three.js 0.170.0 (the only external resource).

```
cd /Users/mickey/solarsystem/ctoc14/results/viz
~/.venvs/astro313/bin/python -m http.server 8765
```

then open <http://localhost:8765/> in Chrome (any modern browser with WebGL works). Any free port is fine.
The page loads `data/index.json` + `data/bodies.json` and logs `[viewer] loaded REAL data ...` in the console; if those
files are missing it falls back to the small fixture under `data/fixture/` and shows a banner. `?data=<dir>/` forces
another data directory (e.g. `?data=data/fixture/`).

After editing anything under `js/`, hard-reload (Cmd+Shift+R): the browser caches the ES modules.

## Controls

| Control | Where | Effect |
|---|---|---|
| Fleet | panel select | switch fleet (index.json order, oldest first); note under the select shows kind/source/provenance |
| play / pause | button or space | animate the clock at the chosen speed |
| speed | log slider, 1-200 days per second (default 30) | playback speed |
| -10 d / +10 d | buttons or arrow keys | step the clock |
| reset | button | pause and go to t = 0 |
| time slider | top of the timeline strip | scrub anywhere in the mission |
| timeline | bottom canvas, one row per craft | launch-end bar, thick segments = thrust arcs, white ticks = flybys, cursor = now; click/drag to scrub |
| asteroid orbit lines | checkbox | 300 faint orbit lines on/off (reachable grey, unreachable 131/144 red dashed) |
| labels / asteroid IDs | checkboxes | Sun/Earth/Moon/craft labels; asteroid ID labels (300 DOM nodes, off by default) |
| Earth-Moon zoom | checkbox | camera parks 0.011 AU from Earth and follows it; the Moon (true mean-orbit distance) and its orbit appear; Earth/Moon spheres are enlarged, as the legend says |
| follow | select or the "f" button in the craft table | camera target locked to that craft |
| home view / top view | buttons | reset the free camera |
| on (checkbox per craft) | craft table | show/hide a craft's trail, marker and timeline row |
| mouse | scene | drag rotate, wheel zoom, right-drag pan (OrbitControls) |

## HUD

covered so far / 298 reachable, craft count, launched so far, fuel used so far (m0 - m(t) interpolated along each
craft's mass array; a finished craft counts its accounted `fuel_kg`), raw J and sum J_i of the fleet, and per craft:
flybys so far / total and fuel used / total in kg. Flown targets turn green in the scene; a ring flash marks each
flyby during playback.

## Files

- `index.html`, `css/viewer.css`, `js/{main,data,scene,fleet,ui,kepler}.js` - the page (see `docs/viz_viewer.md`).
- `data/bodies.json` - Earth, Moon and 300 asteroid elements (`docs/viz_kepler.md` documents the propagation).
- `data/index.json`, `data/fleets/<key>.json` - the nine fleets (`docs/viz_fleets.md` documents the conversion).
- `data/make_bodies.py`, `data/make_fleets.py` - regenerate the data; `data/check_bodies.py`, `data/check_fleets.py` -
  verify it (run from the project root with `nice -n 5 ~/.venvs/astro313/bin/python ...`, both exit 0 when everything passes).
- `data/fixture/` - the 2-craft development fixture (never collides with the real files).
- `screenshots/` - each fleet at t = 600 d, 2800 d and mission end, plus Earth-Moon zoom, follow-craft and top views.
- `docs/viz_report.md` - what was built, the fleet list with numbers, checks, fixes and known limits.
