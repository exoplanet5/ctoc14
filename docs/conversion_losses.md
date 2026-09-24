# Low-thrust conversion losses: why legs are lost and what the fallbacks recover

Scope: `ctoc14/lowthrust.py::convert_tour` only (planner untouched). Test tours: `results/fleet_b300/tour_sc3.json` (44 planned),
`tour_sc8.json` (19), `tour_sc7.json` (25), converted with the tour's own m0 = 2000 kg, `verbose=False`, 3 conversions in parallel
(fork). Every output was validated with `ctoc14.validator` (PASS). Outputs and logs: `results/review2/` (`<tag>_scN.txt`,
`<tag>_scN_info.json`, `<tag>_scN.log`, `<tag>_summary.json`); runner: `tools/convert_review2.py`.

## 1. Why legs fail (from the baseline logs)

The window solver (2-leg sliding window, Newton/IRLS on thrust samples + flyby times, 12 iterations, then re-plan → drop)
loses legs through four mechanisms, all visible in `results/fleet_b300/convert_sc*.log` and `results/review2/base_sc*.log`:

1. **Slow-but-converging windows cut off at `max_iter = 12`.** Several dropped windows were still converging geometrically
   when the iteration limit hit: SC7 ast 125 at 340 km (tolerance 300 km), SC3 ast 141 at 1314 km, ast 98 at 2953 km, ast 123 at
   9078 km, SC7 ast 216 at 3593 km (baseline). Cost: a flyby, plus a re-plan that usually reshuffles the rest of the tour.
2. **Saturation (max |T| = 0.45 N through the whole window, miss stalls at 1e6–1e8 km).** The junction Δv needed after the
   frozen flyby is larger than the two free legs can deliver (SC3 ast 106 after a 35-day leg; ast 256/117/183; SC8 ast 90, 41, 77;
   SC7 ast 15). The thrust that would have prepared the junction lies in the previous, already frozen, leg.
3. **Cascading after a failure.** Once a window fails, `Ts[free]` is reset and the rest of the tour is re-planned by a beam search
   from the actual state with `eta = 0.5`. From an off-plan state the Lambert-junction model is very pessimistic (the actual
   post-flyby velocity differs from the Lambert arrival velocity by several km/s; with the planner's default of 2 TOFs per target
   the arrival velocities that make the *next* leg feasible are discarded), so a re-plan typically keeps 1–2 of 6–8 remaining
   targets ("replanned from leg 37: 1 legs (was 7)" on SC3, "2 legs (was 6)" on SC8, "6 legs (was 8)" on SC7).
4. **Propellant exhaustion at the end.** The plans use 1348–1398 of the 1400 kg; the conversion's per-leg cost is on aggregate
   close to the plan (cumulative actual/plan 0.9–1.3, launch leg excluded) but individual legs cost 3–8× plan (SC3 ast 117: 84 kg
   vs 11.5 planned; SC8 ast 192: 174 vs 42; SC7 ast 12: 88 vs 11), so the last planned legs meet an empty tank: SC3 froze its 38th
   leg at m = 613 kg and the beam re-plan from there could afford one more target; SC8 reached leg 13 with 61 kg.

A latent bug was also found and fixed: `info.setdefault('replans_at', {})[i] = info['replans_at'].get(i, 0) + 1` evaluates the
right-hand side first and raised `KeyError` on the first re-plan whose best sequence did not start with the same two targets
(the whole conversion then failed; `tools/convert_fleet.py` would have reported the tour as `ok=False`).

## 2. What was implemented (all keyword options of `convert_tour`)

| option | mechanism | default |
|---|---|---|
| `adaptive_iter` (+ `max_iter_ext=36`) | keep iterating past `max_iter` while the largest miss still halves every two iterations | see §3 |
| `widen_on_fail` (+ `max_widen=8`) | (a) when a 2-leg window fails, un-freeze the previous leg and re-solve legs (i−1, i, i+1) as one 3-leg window (the previous leg's thrust is free again, so the tight junction gets the capacity of a whole extra leg); if that fails too, restore the state at the original failure and continue with re-plan/drop | see §3 |
| `smart_replan` | (b) before the full beam search: re-time the remaining sequence *in its original order* from the actual state (Lambert model, planner's η = 0.6, **all** feasible TOFs per target kept, beam stratified by position in the list); accept it if complete; else compare the order-preserving version with ≤ 2 skips against the full beam search and take the one with more targets; any candidate whose first two (target, time) pairs match a window that already failed at this leg (±10 d) is rejected | see §3 |
| `fuel_guard` (+ `fuel_guard_factor='auto'`) | (c) at each window: planned propellant of the remaining legs × calibration factor (actual/planned ratio of the flown legs, launch leg excluded, clipped to [1, 1.3], active after 300 kg planned) vs. remaining tank; if short by > 10 kg, drop the fewest legs — the most expensive planned legs first, without overshooting the deficit, never the imminent leg | see §3 |

With all four flags `False` the function reproduces the original behaviour exactly (checked: identical frozen-leg lines against
the pre-refactor logs).

Two designs were tried for (c) and rejected: a fixed-time Lambert chain re-evaluated from the actual state (a 2.4-day shift of
one flyby time changed the estimate for SC7 from 1385 to 1827 kg — single-revolution Lambert junctions are hypersensitive to
small time changes and mostly measure model artefacts), and greedy removal on that chain (removing a middle leg turns the next
junction into a multi-revolution transfer that single-rev Lambert calls infeasible, so it always stripped the tail). The
order-preserving re-timer with skips is kept for (b) but is not used by the guard.

## 3. Results

| configuration (options on) | sc3 (44 planned) | sc7 (25) | sc8 (19) | total |
|---|---|---|---|---|
| base (all off) | 38 | 22 | 15 | 75 |
| all four on | 44 | 20 | 10 | 74 |
| all except fuel_guard | 44 | 25 | 12 | 81 |
| all except smart_replan | 44 | 25 | 15 | **84** |
| all except widen_on_fail | 38 | 25 | 15 | 78 |
| all except adaptive_iter | 44 | – | 15 | – |

(converted flyby counts, every file validated PASS; from results/review2/*_info.json)

## 4. Recommendation

Defaults set by the lead (2026-09-07): `adaptive_iter=True`, `widen_on_fail=True`, `smart_replan=False`, `fuel_guard=False`.
The 3-leg widening removes the sc3 losses entirely and adaptive iteration removes the sc7 losses; the smart re-planner and the
propellant guard lose targets on sc8 and are kept as opt-in options. The re-plan counter was additionally changed to count per
target instead of per leg index (a per-index counter disabled re-planning after a few drops at the same index and caused cascades,
e.g. plan_E1 SC7 lost 19 of 25).
