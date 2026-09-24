# Review of `docs/problem_summary.md` against `CTOC14_problem.pdf`

Reviewer: rules-extraction agent, 2026-09-05. All 10 PDF pages read (images + `pdftotext`).
Companion files: `docs/rules_spec.json` (machine-readable rules), `docs/thrust_interpolation_rule.md`
(validator thrust model, index-base resolution, overshoot data).

## 1. Verdict

`problem_summary.md` is accurate. Every constant, epoch, constraint, scoring formula and format rule in
it matches the PDF, and every derived number checks out to the stated precision. No hard errors were
found. The items below are (a) small wording/precision corrections, (b) rules the summary omits, and
(c) additional derived facts the lead should have.

## 2. Line-by-line verification of the summary

| summary claim | PDF | status |
|---|---|---|
| HEI J2000 frame; units km, km/s, s, kg, N | §3.1, Table 1 | OK |
| EOM, `m' = -T/(Isp g0)`, Isp 4000 s, g0 9.80665, Tmax 0.5 N | eq (1)-(3), Table 2 | OK |
| mu = 1.32712440018e11 km³/s², AU = 149597870.7 km | Table 2 | OK |
| fixed Keplerian Earth/asteroids, `M = M0 + n(t - t_eph)`, `n = sqrt(mu/(a AU)^3)` rad/s | §3.3, eq (4)-(5) | OK (t - t_eph in seconds) |
| Earth elements, epoch MJD 60676.0 = 2025-01-01 | Table 3 | OK, all six values verified digit by digit; date verified |
| MEA.txt columns, epoch MJD 61200.0 = 2026-06-09 | §8, Table 2 | OK; date verified; IDs 1..300 consecutive verified |
| no gravity assists, no Earth return | §4.4 | OK |
| t0 = MJD 62502.0 = 2030-01-01; window 0 ≤ t ≤ 4.73364e8 s = 5478.75 d | eq (6), (14) | OK (15 × 365.25 d × 86400 s = 473 364 000 s exactly) |
| launch: r = r_Earth (1 km), v_inf ≤ 4 km/s (0.01 km/s) | eq (8)-(9), §7 | OK |
| m_dry 600, m(t) ≥ 600, m0 ≤ 2000 | eq (10)-(12) | OK |
| flyby d ≤ 1000 km at declared Event=3 time; first detection only; repeats unpenalised | eq (13), §4.3, §7 | OK |
| J_i = 1 + x + x², x = (m0-600)/1400; 1 empty, 3 full | eq (15) | OK |
| J = ΣJ_i + N_miss | eq (16) | OK |
| k = 1 − 0.1 (t_end − t_submit)/(t_end − t_start) ∈ [0.9, 1]; J_final = kJ; earlier better | eq (17)-(18), §5.4 | OK |
| 15 columns, order, event codes, first Event=0 / last Event=4, strictly increasing Time | Table 1, §6.1-1 | OK |
| arc definition, Event=3 inside arc allowed, arcs separated by 0/2/4 | §6.1-2 | OK |
| cubic sliding-window Lagrange, window j−1..j+2 clipped to [0, n−4], n<4 → all samples; zero outside [t1,tn] | §6.1-3 | OK, and the centred reading is now **confirmed** by the PDF example (see §4 below) |
| ≥ 8640 s spacing; |T(t)| ≤ 0.5 N at all times incl. between samples | §6.1-4, §7 | OK |
| Event 0/2/3/4 thrust (0,0,0); mass by `-|T(t)|/(Isp g0)`, 0.01 kg | §6.1-6/7 | OK |
| validator integrates line to line; 1 km / 1 m/s / 0.01 kg | §7 | OK |
| "Output ≥ 10 significant digits" | §6.1-8 | **wording**: the PDF *recommends* (建议) ≥ 10 digits; it is not a format requirement. See §3.3 — use `%.17g`. |
| Isp·g0 = 39.227 km/s; full-tank Δv = 47.2 km/s | derived | OK (39.2266 km/s; 47.228 km/s) |
| ṁ(Tmax) = 1.2746e-5 kg/s; 1400 kg lasts 1.098e8 s = 1271 d = 3.48 yr | derived | OK (1.27465e-5 kg/s; 1.0983e8 s; 1271.2 d; 3.480 yr) |
| accel 0.25 → 0.83 mm/s²; ≈ 7.9 km/s per year at 2000 kg | derived | OK (7.889 km/s/yr at constant 2000 kg). **Refinement**: with mass depletion a full-tank spacecraft gains 8.8 km/s in its first year of continuous thrust, and 26.3 km/s/yr at dry mass. |
| ballistic spacecraft: 4 DOF, 2 net constraints per flyby → at most 2 asteroids generically | derived | OK as a generic-counting argument (see §3.5 for caveats) |

## 3. Corrections and refinements (do not edit the summary; apply in the lead's copy)

### 3.1 Precision / wording
* "Output ≥ 10 significant digits" → "≥ 10 significant digits *recommended* by the PDF; use full
  double precision (`%.17g`)". Reason: at 10 significant digits a 30 km/s velocity is rounded by up to
  5e-9 km/s, which over a 1e8 s coast segment alone is ~0.5 km of the 1 km budget (`%.10e` = 11
  digits gives 0.05 km; `%.17g` makes it negligible).
* "≈ 7.9 km/s per year" is the constant-mass figure; 8.8 km/s in year 1 with depletion, 26 km/s/yr at
  dry mass.
* "Adjacent Event=1 samples ≥ 8640 s apart": the rule applies *within an arc* and *ignores Event=3
  lines in between* (§6.1-4). There is no minimum spacing for Event=2/3/4 lines (only strict increase).

### 3.2 Rules present in the PDF but absent from the summary
* File name `CTOC14_Result_TeamID.txt`, UTF-8; fields separated by spaces or tabs; blank lines and
  `#` comment lines are ignored (§6, §6.1-8, §7).
* File sorted by `SC_ID` ascending, `SC_ID` consecutive from 1; `Line` is a positive integer (§6, Table 1).
* `Asteroid_ID` must be 1–300 on Event=3 lines and 0 otherwise (Table 1, §7).
* `m0,i` used in the cost is read from the **Event=0 line's mass column** (§7).
* Event=3 lines "record the state at the flyby time and do not change the thrust state"; the thrust at
  a flyby time inside an arc is the arc interpolant (§6.1-3/5), yet the Event=3 line's thrust columns
  must be `(0,0,0)` (§6.1-6, §7).
* Event=2 lines are "pure trajectory nodes with zero thrust, usable to split long integration segments"
  (§6.1-5) — and they *separate* arcs (§6.1-2).
* Validator's format check requires **every line** (including Event=2 and Event=4) to satisfy
  0 ≤ t ≤ Tmax (§7 bullet 1), which is stricter than §5.1 ("flybys outside the window are not counted").
* Time coefficient: computed from the **last valid submission**; submissions after t_end are not
  scored; t_start/t_end are published on the website (§5.3).
* The validator checks each segment starting from the *submitted* state of the previous line (§7
  bullet 2) — no error accumulation; each Event=3 line is individually required to satisfy d ≤ 1000 km
  (§7 bullet 5), so a redundant or marginal flyby declaration can fail the file.
* Problem statement §2 item 4: partial coverage is explicitly allowed ("未全覆盖的方案仍可提交").

### 3.3 Derived facts worth adding
* **Fuel is cheap at low loads.** `m_fuel = 600 (e^{Δv/39.2266} − 1)`; J_i − 1 = x + x²:
  1 km/s → 15.5 kg, +0.011; 5 km/s → 81.6 kg, +0.062; 10 km/s → 174 kg, +0.140; 20 km/s → 399 kg,
  +0.366; 30 km/s → 689 kg, +0.735; 47.2 km/s → 1400 kg, +2.0. The binding limits for a lightly
  fuelled spacecraft are thrust-time (74 d of full thrust for 5 km/s at 600–680 kg) and the calendar,
  not the fuel cost. Consequently a "ballistic + small correction" spacecraft (a few 100 m/s, ΔJ ~1e-3)
  can generically pick up a third target that a pure ballistic double misses by ≲ 1e6 km.
* Cost-per-asteroid benchmark: a ballistic double costs 0.5 per asteroid; an x = 0.5 spacecraft
  (30 km/s) must cover ≥ 4 to beat that (0.43); a full tank must cover ≥ 7.
* Since an empty spacecraft costs exactly one missed asteroid, a spacecraft covering a single asteroid
  is never strictly beneficial; every spacecraft must cover ≥ 2.
* Earth elements give a period of 365.76 d (a = 1.00092 AU): do **not** substitute a JPL ephemeris;
  the PDF example launch position agrees with the Table-3 Keplerian Earth to 0.5 m.
* MEA.txt population (verified): a 0.64–17.8 AU, e 0.07–0.95, i 0.45–134° (retrograde objects exist),
  q 0.14–1.05 AU, Q 0.96–34.7 AU; 206 have e > 0.5, 105 have i > 20°, 131 have a > 2 AU. Full coverage
  is not realistic; the cost function makes skipping an asteroid cost exactly 1.

### 3.4 Index-base ambiguity of the Lagrange window — resolved
The PDF mixes 1-based sample labels (t_1…t_n) with a 0-based clip range [0, n−4]. Reconstructing the
§6.2 example arc (|T_k| = 0.5 sin(πk/40), 9°/sample rotation about +z, elevation 17.1887°) and
integrating line 4 → line 5 reproduces the printed mass (1999.8273580 kg) to 2e-10 kg only with the
0-based centred window `l = clip(j−1, 0, n−4)`, window {j−1, j, j+1, j+2}; the forward-window reading
is off by 1.5e-5 kg and 57 m. The summary's "window j−1..j+2" is therefore correct with j 0-based.
Full pseudo-code in `docs/thrust_interpolation_rule.md`.

### 3.5 Caveats on the "ballistic ≤ 2 asteroids" argument
The DOF count (4 launch DOF; each flyby = 3 position equations − 1 free time = 2) is right and gives a
discrete set of ballistic doubles; a ballistic triple can only occur by a ~(1000 km / 1 AU)² ≈ 5e-11
coincidence. But the count says nothing about *existence* under |v_inf| ≤ 4 km/s and the 15-year
window: most catalogue orbits (e up to 0.95, i up to 134°) are unreachable ballistically, so the
number of feasible ballistic doubles is an empirical question. Conversely, the 1000 km ball plus a few
100 m/s of (almost free) Δv relaxes the count considerably (§3.3).

## 4. Traps and ambiguities a solver must respect

1. **Absolute window.** All lines, including the final Event=4 and any Event=2, must have
   0 ≤ Time ≤ 473 364 000 s; a single late line is a format violation for the whole file. Launching
   later only shortens the mission. t = 0 is a legal launch time (used in the example).
2. **Flybys count only when declared.** Each counted asteroid needs an explicit Event=3 line at the
   flyby time with `Asteroid_ID` 1–300 and d ≤ 1000 km computed from that line's position; undeclared
   close passes are ignored. Repeats are harmless for the score but every Event=3 line must pass the
   distance test — never declare a marginal or duplicate flyby you have not verified.
3. **v_inf tolerance (0.01 km/s) and launch position tolerance (1 km) are for rounding**, not design
   margin. The example uses exactly 4.000000 km/s. Design to ≤ 4.000 km/s.
4. **Time coefficient.** k = 0.9 → 1.0 linearly from t_start to t_end and is taken from the *last*
   submission: resubmit only if J_new·k_new < J_old·k_old. The 10 % swing equals ~10 % of J (tens of
   cost units for J ~ 100–200), so early submission of a decent solution followed by disciplined
   resubmission matters.
5. **Thrust model fidelity.** Submitted states must be produced with the validator's exact model
   (0-based centred cubic Lagrange per component, mass from the interpolated vector norm, zero thrust
   outside the arc); a 1 % thrust-model discrepancy at 0.5 N / 2000 kg gives ~9 km per 86 400 s
   segment vs the 1 km tolerance.
6. **Interpolation overshoot.** Steps inside an arc overshoot by 6.4 %, a two-sample plateau by
   12.5 %, ramp-to-plateau junctions by 1–3 %, sudden 90° turns at full thrust by 6.6 %. Safe
   primitives: constant-magnitude arcs with instant on/off at the arc ends (0 % overshoot) and
   concave bumps whose peak is a sample (0 %). Always run the exact per-interval polynomial check and
   cap samples at e.g. 0.4999 N (also guards against `|T| = 0.5 + 1e-10` after decimal rounding).
7. **Sample spacing.** Event=1 lines inside an arc must be ≥ 8640 s apart; thrust cannot be shaped
   finer than 0.1 d. Event=3 lines are exempt (the example has an Event=3 line 10 s before an arc).
8. **Event=2 splits arcs.** Never put an Event=2 between samples of one arc. Do use Event=2 nodes to
   split long coasts (≲ 100 d) so print rounding and the validator's integrator error stay ≪ 1 km.
9. **Zero-fuel spacecraft cannot thrust at all** (m must stay ≥ 600 kg), and m0 on the Event=0 line
   is what the cost uses; a spacecraft with m0 = 600 must have no non-zero Event=1 line.
10. **Event=3 inside an arc**: thrust columns (0,0,0) but the mass must include the burn to that
    instant; the segment before and after it uses the arc interpolant.
11. **Units**: Time in seconds from t0 (not MJD), km, km/s, kg, N; T/m is m/s² → divide by 1000.
12. **No gravity assists**: automatically true in the two-body model (planets have no gravity in the
    EOM); do not add any perturbation or Earth/asteroid gravity.
13. **Segment-wise checking** (each segment restarted from the submitted state with 1 km / 1 m/s /
    0.01 kg slack) could in principle be abused as free Δv; treat the slack as numerical margin only —
    the organisers can inspect and reject such files. Keep own integration errors ≤ 1e-3 km.
14. **Sorting/IDs**: sort by SC_ID then Time, SC_IDs consecutive from 1, positive line numbers (use
    global 1..N), Asteroid_ID 0 on non-flyby lines, first line Event=0 with mass = m0, last Event=4.
15. **n < 4 arcs**: n = 1 gives thrust only at an instant (no effect, no mass change); n = 2 linear,
    n = 3 quadratic. Avoid n = 1.
16. **Mission window vs Event=4**: yes, Event=4 must be ≤ Tmax (format check); a spacecraft may end
    immediately after its last flyby (no return, no coast required).
17. **Catalogue realism**: several targets have Q up to 35 AU or retrograde orbits; plan on skipping
    them (cost 1 each) rather than chasing full coverage.
