# CTOC14-A Solver Design Brief

## 0. What exists already (validated against the PDF example)
Package `ctoc14/` (import with project root on sys.path):
- `constants.py` — all official constants, `cost_sc(m0)`.
- `kepler.py` — `Ephemeris` (`earth_state(t)`, `ast_state(k,t)` k = ID-1, `all_ast_states(t)`), `Body`, `solve_kepler`,
  `propagate_twobody(r0,v0,dt)` (universal variables, vectorised over dt, ~1e-6 km over 15 yr). t = seconds since t0.
- `lambert.py` — `lambert(r1, r2, tof, prograde=True)` vectorised single-rev (12 µs/solve), NaN where unsolved.
- `thrust.py` — exact validator thrust model `thrust_interp(ts, Ts, t)` (sliding-window cubic Lagrange), `arc_max_thrust`.
- `validator.py` — replica of the official checker: `validate(path) -> Report` (errors, max errors, flybys, J). `integrate_segment`.
- `submission.py` — `SCTrajectory` (launch → coast / thrust_arc(ts,Ts,flybys) / flyby / end) + `write_submission`.
Tests in `tests/` reproduce PDF example lines to < 10 m and mass to 1e-7 kg; `results/test_example_submission.txt` passes.

Analyses (read their docs first): `docs/population_analysis.md` (298 reachable, 131 & 144 impossible; ~180 in-band nodal
crossings/yr), `docs/ballistic_reach.md` / `analysis/ballistic/out/earth_to_asteroid_summary.csv` (296 targets reachable by a
zero-fuel spacecraft, median min v_inf 0.98 km/s; junction delta-v statistics), `docs/methods_survey.md` (GTOC4 lessons).

## 1. Physics of the problem in one paragraph
Every target passes within 0.05 AU of Earth's orbit; the problem is phasing, not energy. A flyby needs only position
coincidence (1000 km) — arrival velocity is free — so high-inclination targets are met at their ecliptic nodes and
fast-moving targets cost nothing extra. The spacecraft (0.25–0.83 mm/s², 47 km/s total, 3.5 yr of thrust) threads a stream of
~200 timing windows per year on a near-1-AU orbit. Position-only targeting has a leverage penalty: thrust applied late in a
leg barely moves the flyby position, so a junction impulse Δv is deliverable within the following leg of duration Δt only if
Δv ≲ ½·a·Δt (constant thrust) and costs up to 2Δv of propellant if spread uniformly; thrusting near the junction (end of the
previous leg + start of the next) makes it nearly impulsive (cost ≈ Δv). Hence legs are optimised in a 2-leg sliding window.

## 2. Cost logic (drives every decision)
- Missing a target costs 1. A zero-fuel spacecraft costs 1 and generically reaches ≤ 2 targets (4 DOF: launch time + v∞ vector
  vs 2 net constraints per flyby) ⇒ pairing leftovers is J-neutral at worst, gains 1 per pair.
- Full tank costs 3 and buys 47 km/s / 3.5 yr of thrust; if it yields n flybys the cost per target is 3/n (n ≈ 40–60 expected).
  Marginal fuel cost dJ/dx = 1 + 2x ≤ 3 per 1400 kg ⇒ the last 140 kg (≈2.9 km/s) costs 0.29 ⇒ worth it if it buys ≥ 1 flyby.
- Therefore: maximise flybys per full-tank spacecraft; add spacecraft while (flybys gained) > (cost added); mop up leftovers
  with cheap spacecraft (x ≈ 0–0.1, cost 1.0–1.11) reaching 2–4 targets each; accept 131 and 144 as misses (floor J ≥ 2 + ΣJ_i).

## 3. Pipeline
```
reach (Earth→target ballistic table)  ─┐
                                        ├─► search.py : beam search, one spacecraft, Lambert-junction model ─► Tour JSON
fleet.py : sequential allocation (remove visited targets, repeat), cost accounting, stop rule, mop-up via pairs.py
lowthrust.py : convert each Tour to continuous thrust (exact validator thrust model), 2-leg sliding window, 1000 km miss ─► SCTrajectory
submission.py ─► CTOC14_Result_TeamID.txt ─► validator.py (must PASS) ─► J
```

## 4. Module contracts
### 4.1 Tour JSON (output of search, input of lowthrust)
```json
{"m0": 2000.0, "t_launch": 1.2e7, "vinf": [vx,vy,vz], "legs": [
   {"ast": 174, "t_flyby": 2.2e7, "dv_est": 0.0, "tof": 1.0e7},
   {"ast": 57,  "t_flyby": 3.1e7, "dv_est": 0.63, "tof": 0.9e7}, ...],
 "dv_total_est": 31.2, "fuel_est": 1120.0, "n_flybys": 47, "t_end": 4.6e8}
```
`dv_est` = Lambert junction Δv (km/s) at the START of that leg (0 for the first leg if within v∞ ≤ 4).

### 4.2 search.py — `beam_search(eph, excluded:set[int], m0, t_launch_grid, params) -> list[Tour]`
- Root: for each launch date in `t_launch_grid` (default: every 10 d over the first 3 yr, extendable), each target ∉ excluded
  ∪ {131,144}, each TOF in the grid (15–400 d, ~40 values): Lambert Earth→target; feasible if |v1 − v_E| ≤ 4 + η·a·TOF
  (the excess over 4 km/s must be paid by thrust: Δv_est = max(0, |v1−v_E| − 4)).
- Expansion of a state (t, r, v, m, visited): vectorised Lambert to all unvisited targets × TOF grid; Δv = |v1 − v|;
  feasible if Δv ≤ η·a(m)·TOF (η default 0.6, calibrated by lowthrust.py); fuel Δm = m·(1 − exp(−κΔv/VE)) with κ = 1 + Δv/(2·a·TOF)
  (leverage penalty); child velocity = Lambert arrival velocity v2; keep the best `n_per_target` TOFs per target.
- Constraints: t_flyby ≤ T_MISSION; m ≥ 600 (+ margin); optional min leg TOF 15 d (need ≥ 8640 s spacing for samples).
- Beam scoring (expose weights): score = n_flybys − w_t·(t/T_MISSION) − w_dv·(Δv_total/Δv_max); beam width W (default 300);
  also keep a diversity rule (no more than K children of the same parent). Return the best tours (by n_flybys, then Δv).
- Performance: vectorise (298 targets × 40 TOF = 12k Lambert per state ≈ 0.15 s); use multiprocessing across beam states.
  Budget: one full-tank tour in ≤ 30 min on 10 cores.

### 4.3 lowthrust.py — `convert_tour(eph, tour, params) -> SCTrajectory` (+ per-leg diagnostics)
- Thrust parametrised by samples every h (default 1 d, ≥ 8640 s) on the validator's cubic-Lagrange model; cap |T_k| ≤ 0.45 N at
  samples and verify the interpolant max ≤ 0.5 N (`arc_max_thrust`), shrink if violated.
- Per 2-leg window: variables = thrust samples over legs k, k+1 (+ flyby times within ±δ); constraints = miss distance at
  flybys k and k+1 (target 0, accept ≤ 900 km); objective = propellant. Linearise: sensitivities ∂r(t_f)/∂T_j by variational
  integration or finite differences on a vectorised fixed-step RK integrator; solve the small least-norm / SLSQP problem; re-integrate
  exactly (validator integrator) and iterate until miss < 900 km. Freeze leg k, slide. Coast where thrust is negligible.
- Launch: v∞ direction/magnitude from the tour (|v∞| ≤ 4 − 1e-6 km/s), position = Earth exactly.
- Must produce rows via `SCTrajectory` so the file is self-consistent; run `validate` on the result.
- Calibration deliverable: for ~50 random junctions report achieved propellant vs Lambert Δv (κ) and the max feasible Δv/(a·TOF) (η).

### 4.4 pairs.py — cheap spacecraft for leftovers
`find_chains(eph, targets:set[int], max_fuel_kg, params) -> list[Tour]`: Earth→a(→b→c) chains with v∞ ≤ 4 and tiny Δv (grid +
root-finding on (t_L, v∞) for exact 2-flyby ballistic solutions; small thrust for a 3rd/4th). Greedy set cover of leftovers.

### 4.5 fleet.py — allocation and cost report
Sequential: tour_1 = best search over all targets; exclude; repeat with full tank while (flybys − 3) > 0 and it beats
alternatives; optimise the last tank size x; mop-up with pairs.py; anything left = miss. Produce `results/fleet_report.md` with
per-spacecraft cost, flybys, cost per target, and the marginal-cost table answering "which targets are not worth visiting".

### 4.6 tools/run_solver.py — end-to-end CLI producing `results/CTOC14_Result_<team>.txt`, validator report, J.

## 5. Verification requirements
- Every produced file must PASS `ctoc14.validator` (and an independently written second validator).
- Unit tests: kepler vs PDF; lambert round trips; thrust rule vs PDF masses; search feasibility model vs lowthrust results.
- Report honest numbers: flybys per spacecraft, propellant used, J. No claims without a validated file.
