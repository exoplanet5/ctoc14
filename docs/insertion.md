# Inserting uncovered targets into existing planned tours (`ctoc14/insertion.py`)

## 1. What it does
`ctoc14/insertion.py` takes planned tours (JSON from `search.tour_json`: `t_launch`, `vinf`, `m0`, `legs[(ast, t_flyby, dv_est, tof)]`)
and tries to insert uncovered targets into them with the **same Lambert-junction model as the planner** (`search.roots` /
`search.expand`): junction Δv at the start of leg i is `dv[0] = max(0, |v1[0] − vE| − 4)` and `dv[i] = |v1[i] − v2[i−1]|`
(Lambert departure / arrival velocities), feasible iff `dv ≤ η·a·tof` (η = 0.6, a = Tmax/m with m the mass at the start of the leg),
propellant `m(1 − exp(−κ dv / v_e))`, κ = 1 + dv/(2 a tof). Tours are re-derived from `(t_launch, m0, asts, times)` alone
(`build_tour`), so the JSON `dv_est`, `tof`, `fuel_est`, `vinf` are only used for the unit check.

**Candidates** for target j in tour T: for every gap k = −1 (launch → first flyby), 0 … n−2 (between flybys k and k+1) and n−1
(after the last flyby), every intermediate time t_j on a 1 d grid (legs ≥ 15 d), and every shift δ ∈ {−10, −7.5, …, +10} d of the
flyby time of k+1. The detour changes legs A = k→j, B = j→k+1 and, when δ ≠ 0, C = k+1→k+2 (shortened by δ, Lambert arc recomputed)
and the junction at k+2 (its arrival velocity changed); all other legs keep their Lambert arcs but are re-accumulated with the lower
mass. Every candidate is checked for (i) the feasibility rule on the modified legs (unchanged legs cannot get worse: the mass only
decreases downstream), (ii) timeline `t ≤ 15 yr − 5 d`, (iii) budget `fuel_new ≤ m0 − 600 − margin` (margin default 60 kg), optional
(iv) `u_max`: max dv/(η a tof) over the modified legs. All (gap, t_j, δ) candidates of one (tour, target) pair — about 40 000 for a
49-leg tour — are evaluated in one vectorised pass (3 Lambert calls + one vectorised mass chain), 0.4 s single-core.

**Greedy driver** (`insert_targets`): evaluate every (tour, target); apply the globally cheapest feasible insertion (extra propellant);
rebuild that tour from scratch (`apply_insertion` → `build_tour`, asserting the rebuilt fuel equals the vectorised estimate); re-evaluate
only that tour; repeat until nothing fits. Targets already in a tour's plan are skipped for that tour; 131/144 are skipped
(`try_unreachable=False`). Optional second phase `grow_m0`: for tours with m0 < 2000 kg, find the smallest m0' that makes the cheapest
leg-feasible candidate fit (bisection on the mass chain; leg feasibility only gets worse with m0) and rank by ΔJ = J(m0') − J(m0).

**Output**: every tour (modified or not) as `outdir/<name>.json` in the planner format (plus `inserted: [...]`), so the directory is a
drop-in input for `tools/convert_fleet.py`; `outdir/insertion_<tag>.json` with the history and, for every remaining target and tour,
the cheapest leg-feasible candidate ignoring the budget (extra kg needed vs budget available).

```
P=~/.venvs/astro313/bin/python
$P tools/run_insertion.py results/fleet_b300 results/insertion/run --targets 15,137,138,216 --margin 60 --nproc 4
$P tools/run_insertion.py results/fleet_b300 results/insertion/run2 --targets-json results/fleet_b300/tail3_targets.json --full-tank-only --margin 0
   options: --tours tour_sc1,tour_sc4  --eta 0.6 --tstep 1 --delta-max 10 --delta-step 2.5 --tof-min 15 --tof-max-last 400
            --no-launch-gap --u-max 1.0 --grow-m0 --try-unreachable --from-report results/CTOC14_Result_b300_report.json
$P tests/test_insertion.py
```

## 2. Verification
- Reconstruction vs stored plan (all 10 `results/fleet_b300/tour_sc*.json`, the 5 `tail3_m*.json`, `extra_sc11/12.json`):
  max |dv − dv_est| = **2.7e-10 km/s** (requirement 1e-3), fuel to 1e-8 kg, v∞ vector exact, all legs feasible (min slack +0.014 km/s on
  sc2, +0.077 on sc1 — the plans sit close to the η rule).
- Vectorised candidate evaluation vs rebuilt tour: agree to < 1e-6 kg on random candidates including shifted ones (`tests/test_insertion.py`).
- Planner formulas reproduced exactly (`mass_chain` vs `search.expand`).

## 3. Results

### 3.1 Remaining propellant of the planned tours (the binding constraint)
| tour | m0 | n | fuel_est | left = m0−600−fuel | budget @ margin 60 | budget @ margin 0 |
|---|---|---|---|---|---|---|
| sc1 | 2000 | 49 | 1383.7 | 16.3 | −43.7 | +16.3 |
| sc2 | 2000 | 48 | 1396.2 | 3.8 | −56.2 | +3.8 |
| sc3 | 2000 | 44 | 1347.6 | 52.4 | −7.6 | +52.4 |
| sc4 | 2000 | 39 | 1364.8 | 35.2 | −24.8 | +35.2 |
| sc5 | 2000 | 36 | 1373.1 | 26.9 | −33.1 | +26.9 |
| sc6 | 2000 | 29 | 1377.0 | 23.0 | −37.0 | +23.0 |
| sc7 | 2000 | 25 | 1397.7 | 2.3 | −57.7 | +2.3 |
| sc8 | 2000 | 19 | 1363.5 | 36.5 | −23.5 | +36.5 |
| sc9 | 1000 | 6 | 342.1 | 57.9 | −2.1 | +57.9 |
| sc10 | 850 | 2 | 159.6 | 90.4 | +30.4 | +90.4 |

The fleet_b300 plans were made with the planner's default 2 kg margin, so every full tank ends at 602–653 kg: **with the 60 kg margin
the budget rule is violated before any insertion** and nothing can be inserted into sc1–sc9. Everything below is therefore reported
for margin 60 (the rule as specified) and margin 0 (use the plan's leftover only), plus the unconstrained "what would it cost" table.

### 3.2 Test 1 — uncovered list [15, 131, 137, 138, 144, 216] against all 10 tours (`results/insertion/test1_*`)
- margin 60: **0 inserted** (131, 144 skipped as unreachable; 137/138 are already in sc10's plan; 15/216 already in sc7's plan).
- margin 0: **1 inserted — #137 into sc4** between 300 and 274 at t = 4357 d, flyby of 274 shifted −5 d: extra **30.3 kg**
  (junction Δv 0.81 / 0.24 / 1.79 km/s, utilisation 0.87 of the η rule), sc4 left with 5.0 kg → `results/insertion/test1_m0/tour_sc4.json` (40 legs).
- #138 and #216: **no feasible detour in any gap of any tour** (not even ignoring the budget). #15: only by appending it after sc8's last
  flyby (77 → 15 at 5262 d, 9.15 km/s over 282 d, 159.7 kg).
- `--grow-m0` (sc9, sc10): nothing — raising m0 lowers the acceleration, so the tight detours become infeasible before the tank is big
  enough (137 into sc9 would need ~200 kg at m0 = 1000 with u = 0.79).

### 3.3 Test 2 — the 22 `tail3_targets.json` targets against the 8 full-tank tours (`results/insertion/test2_*`)
- margin 60: **0 inserted**. margin 0: **1 inserted — #141 appended after sc8's last flyby** (77 → 141 at 5203 d, 1.41 km/s over 223 d,
  **23.5 kg**, u = 0.16), sc8 left with 12.9 kg → `results/insertion/test2_m0/tour_sc8.json` (20 legs).
- Sensitivity (margin 0): shift range ±20 d, η = 0.7, or a 600 d append leg all give the same single insertion. With η = 0.7 two more
  detours become leg-feasible (41 into sc6 for 118 kg, 160 into sc6 for 230 kg) but do not fit any budget.
- Cheapest leg-feasible candidate per target ignoring the budget (margin-60 run, i.e. original tours; "—" = no feasible detour in any
  tour; tours where the target is already planned are excluded):

| target | cheapest option | extra kg | 2nd option | extra kg |
|---|---|---|---|---|
| 15 | sc8 append after 77 at 5262 d, Δv 9.15 | 159.7 | — | |
| 21 | sc8 append after 77 at 5331 d, Δv 3.22 | 53.4 | sc8 between 90/… at 3620 d (δ −2.5), 4.97/1.30/3.64 | 100.4 |
| 41 | — | | | |
| 50 | — | | | |
| 76 | sc8 append after 77 at 5213 d, Δv 4.16 | 72.0 | — | |
| 77 | — | | | |
| 84 | sc8 append after 77 at 5343 d, Δv 3.98 | 66.1 | — | |
| 96 | — | | | |
| 116 | sc8 after 233 at 2865 d (δ +10), 2.39/0.61/5.39 | 73.4 | — | |
| 138 | — | | | |
| 141 | sc8 append after 77 at 5203 d, Δv 1.41 | **23.5** | — | |
| 160 | — | | | |
| 166 | — | | | |
| 183 | sc8 append after 77 at 5315 d, Δv 6.14 | 103.5 | — | |
| 204 | — | | | |
| 212 | sc8 append after 77 at 5317 d, Δv 8.43 | 143.0 | — | |
| 216 | — | | | |
| 229 | sc5 after 14 at 1238 d, 1.26/1.74/0.96 | **27.3** | sc7 after 171 at 1793 d (δ −5), 2.43/0.37/2.06 | 37.3 |
| 240 | sc8 append after 77 at 5358 d, Δv 7.15 | 119.4 | — | |
| 268 | — | | | |
| 272 | sc7 after 171 at 1645 d, 1.26/1.20/3.07 | **31.7** | sc6 append after 52 at 5443 d, Δv 3.84 | 63.7 |
| 285 | sc2 append after 85 at 5449 d, Δv 0.88 | **13.7** | sc6 append after 52 at 5441 d, Δv 4.46 | 74.6 |

Reading: (i) only 9 of the 22 have any feasible detour at all, and only 4 of those cost < 40 kg (285: 13.7, 141: 23.5, 229: 27.3,
272: 31.7); (ii) sc8 is the only tour whose plan ends early (4980 d), which is why almost every option is an "append after 77" with a
long final leg — cheap in kg because the leg is long, but 3–9 km/s of thrusting on an almost empty spacecraft; (iii) 13 targets
(41, 50, 77, 96, 138, 160, 166, 204, 268 and, outside the tours that already plan them, 15, 84, 212, 216, 240) have no feasible
k→j→k+1 detour anywhere with η = 0.6 — these are rare-window targets that need dedicated chains or a re-plan, not a local insertion.

### 3.4 Converted fragments as base tours (`--from-report`, `results/insertion/test3_*`) — not usable with this model
Using the converted flyby sequences/times (`results/CTOC14_Result_b300_report.json`) as bases would expose the propellant the
conversion actually left unused (sc3 107 kg, sc8 56 kg, sc1 20 kg). But the single-rev Lambert chain through the *converted* flyby
points is a different trajectory from the thrusted one: the planner model prices it at 1390 kg (sc2, +32 kg vs actual) up to 1988 kg
(sc7, +640 kg) and finds infeasible legs in 9 of the 10 fragments (min slack down to −13 km/s). The greedy then "finds" negative extra
propellant by re-timing those legs (e.g. −482 kg for 138 in sc10) — model artefacts. The CLI therefore excludes infeasible converted
bases by default (`--allow-infeasible-base` to force) and none of the 22 targets fits into sc2 (the only feasible base, budget −50 kg).
Conclusion: insertions into converted trajectories must be priced with the low-thrust machinery (e.g. `convert_tour(..., allowed_extra=)`
re-planning from the actual state), not with the Lambert model.

## 4. What this means for the fleet (data for the lead)
1. **The planned full-tank tours have no room**: the planner spends the tank to the 2 kg margin, so local insertion can only use the
   2–52 kg it happened to leave. Under the specified rule (60 kg reserve) the answer for both test sets is **0 insertions**; with the
   reserve dropped, 137 → sc4 (30 kg) and 141 → sc8 (24 kg) are the only fits, and both are worth +1 target each if they survive conversion
   (`results/insertion/test1_m0/tour_sc4.json`, `results/insertion/test2_m0/tour_sc8.json`; sc4's candidate uses 0.87 of the η rule —
   marginal for the converter; sc8's is a long, gentle 223 d leg).
2. To make insertion useful the room must be created upstream: plan the full tanks with `--mmargin 60–100` and then fill with the
   cheapest uncovered targets (285 13.7 kg, 141 23.5 kg, 229 27.3 kg, 272 31.7 kg are the candidates at hand), or implement a *replace*
   move (swap a planned target that is also reachable by another tour for an uncovered one) — not part of this module.
3. Insertion cannot rescue the 13 rare-window targets listed in 3.3; those need the tail/mop-up chains or target-weighted re-planning.

## 5. Files
- `ctoc14/insertion.py` — `InsParams`, `mass_chain`, `TourModel`, `build_tour`, `tour_from_json`, `tour_from_converted`,
  `check_reconstruction`, `evaluate_insertions`, `best_insertion`, `apply_insertion`, `min_m0_for`, `insert_targets`, `load_tours`.
- `tools/run_insertion.py` — CLI (options above); `tests/test_insertion.py` — reconstruction / consistency / greedy checks.
- `results/insertion/<run>/` — updated tour JSONs + `insertion_<run>.json` + `<run>.log` for: `test1_m60`, `test1_m0`, `test1_grow`,
  `test2_m60`, `test2_m0`, `test2_m0_d20`, `test2_m0_eta07`, `test2_m0_last600`, `test3_m60`, `test3_m0`, `test1c_m60` (test3/test1c =
  converted bases, superseded by the default exclusion described in 3.4).
