# viz_kepler.md: the Kepler propagation the viewer's JavaScript must implement (2026-09-24)

Contract for `results/viz/js/*.js`: how to turn one element block of `results/viz/data/bodies.json` into a heliocentric
ecliptic J2000 position in AU at a mission time `t` in days since t0. It is a literal transcription of
`ctoc14/kepler.py` (`Body.__init__` + `Body.state`), the contest ephemeris, so a JS implementation of exactly these
steps agrees with the Python solver and with the flyby positions in `results/viz/data/fleets/*.json`.

Verified: `results/viz/data/check_bodies.py` evaluates these steps in plain Python (own Newton solver, no code shared
with ctoc14.kepler) against `ctoc14.kepler.Ephemeris`: 10 asteroids (1, 2, 4, 57, 100, 131, 144, 200, 250, 300) at
5 epochs (0, 1000, 2500, 4000, 5478.75 d): **max error 2.04e-15 AU** (tolerance 1e-6 AU); Earth at the same epochs
0.0 AU; all 300 asteroids on a 50 d grid over the whole mission (110 epochs): max 1.07e-14 AU. Run from the project
root: `nice -n 5 ~/.venvs/astro313/bin/python results/viz/data/check_bodies.py` (exit 0 = PASS).

## 1. Inputs (all from bodies.json)
| symbol | JSON field | unit | meaning |
|---|---|---|---|
| μ | `mu_sun_km3s2` = 1.32712440018e11 | km³/s² | solar gravitational parameter |
| AU | `au_km` = 149597870.7 | km | astronomical unit |
| t0 | `epoch_mjd_t0` = 62502.0 | MJD | mission start 2030-01-01 00:00 UTC; viewer time `t` is days after this |
| a | `a` | AU | semi-major axis |
| e | `e` | – | eccentricity (all bodies elliptical, max 0.947 for asteroid 144) |
| i | `i` | deg | inclination (asteroid 131 has i = 134°, retrograde; the formulas need no special case) |
| Ω | `Om` | deg | longitude of the ascending node |
| ω | `om` | deg | argument of perihelion |
| M0 | `M0` | deg | mean anomaly at the element epoch |
| t_e | `epoch_mjd` | MJD | element epoch: 61200.0 for every asteroid, 60676.0 for Earth |

Blocks: `earth` (one block), `asteroids[k]` (300 blocks, `id` 1..300 in order, so `asteroids[id-1]`), `moon`
(illustrative geocentric orbit, section 6), `unreachable` = [131, 144] (still drawn, red dashed per viz_spec section 3).

## 2. The five steps (per body, per frame)
Angles in **radians** from here on: multiply every degree field by π/180. Precompute steps 1 and 5 once per body.

**Step 1: mean motion, with a in km.**
```
a_km = a * AU
n    = sqrt(μ / a_km³)                       [rad/s]
```
(For a per-day form use n_d = 86400 n [rad/day]; keep the seconds form for exact agreement with the solver.)

**Step 2: mean anomaly at viewer time t [days since t0].**
```
Δt   = (t + (t0 − t_e)) * 86400              [s since the element epoch]   (t0 − t_e = 1302 d asteroids, 1826 d Earth)
M    = M0 + n Δt                             [rad, unreduced]
M    = M mod 2π, mapped into [0, 2π)         (JS: ((M % TWO_PI) + TWO_PI) % TWO_PI)
```

**Step 3: Kepler's equation E − e sin E = M by Newton's method.**
```
E    = (e < 0.8) ? M + e sin M : π           initial guess
repeat (max 50 times):
    ΔE = (E − e sin E − M) / (1 − e cos E)
    E  = E − ΔE
    stop when |ΔE| < 1e-12
```
f'(E) = 1 − e cos E ≥ 1 − e > 0.05 for every catalogue body, so Newton converges monotonically from these starts;
in practice 3–5 iterations. (ctoc14/kepler.py additionally clips |ΔE| ≤ 1 rad; with these starting values the clip
never triggers for the catalogue, and the check above confirms identical results without it. Add it if you like.)

**Step 4: perifocal position (x toward perihelion, y along motion at perihelion, z = 0).**
```
x_p  = a_km (cos E − e)
y_p  = a_km sqrt(1 − e²) sin E
r    = a_km (1 − e cos E)                    (radius, only needed for the velocity)
```

**Step 5: rotate perifocal → heliocentric ecliptic J2000 and convert to AU.**
R = R3(−Ω) · R1(−i) · R3(−ω), the classical rotation with cΩ = cos Ω, sΩ = sin Ω, ci = cos i, si = sin i,
cω = cos ω, sω = sin ω (exactly the matrix built in `Body.__init__` of ctoc14/kepler.py):
```
      | cΩ cω − sΩ sω ci     −cΩ sω − sΩ cω ci      sΩ si |
R  =  | sΩ cω + cΩ sω ci     −sΩ sω + cΩ cω ci     −cΩ si |
      | sω si                 cω si                 ci    |

x = (R[0][0] x_p + R[0][1] y_p) / AU
y = (R[1][0] x_p + R[1][1] y_p) / AU
z = (R[2][0] x_p + R[2][1] y_p) / AU
```
The third column of R multiplies z_p = 0 and is never used for the position. Result (x, y, z) in AU, heliocentric
ecliptic J2000 = the frame of the submission files and of `fleets/*.json` (`xyz`). No axis swap for Three.js is part
of this contract: if the page uses y-up it must apply one fixed swap to *everything* (bodies and trajectories).

## 3. Velocity (optional, same convention; ctoc14/kepler.py `Body.state`)
```
fac  = sqrt(μ a_km) / r
vx_p = −fac sin E
vy_p =  fac sqrt(1 − e²) cos E             [km/s]
v    = R · (vx_p, vy_p, 0)                 [km/s, ecliptic J2000]
```
Only needed if the page wants velocity vectors; the fleet trajectories are sampled positions and need none.

## 4. Worked numeric example: asteroid 1 at t = 1000 d
Elements (bodies.json `asteroids[0]`): a = 1.7689528 AU, e = 0.4232224, i = 22.14242°, Ω = 336.03924°,
ω = 27.87452°, M0 = 252.06689°, t_e = 61200.0.

| step | quantity | value |
|---|---|---|
| 1 | a_km | 264 631 572.248803 km |
| 1 | n | 8.462395747545927e-08 rad/s (= 7.311509925879681e-03 rad/day, period 859.355 d) |
| 2 | Δt | (1000 + 1302) · 86400 = 198 892 800.0 s |
| 2 | M0 | 4.399397165762 rad |
| 2 | M unreduced | 4.399397165762 + 8.4623957475e-08 · 198 892 800 = 21.230493015137 rad |
| 2 | M mod 2π | 2.380937093599 rad (136.417646749°) |
| 3 | E start (e < 0.8) | 2.380937093599 + 0.4232224 sin(2.380937093599) = 2.672705122408 |
| 3 | Newton 1 | ΔE = 7.297e-02 → E = 2.599737400263290 |
| 3 | Newton 2 | ΔE = 3.914e-04 → E = 2.599345976137045 |
| 3 | Newton 3 | ΔE = 1.227e-08 → E = 2.599345963862454 |
| 3 | Newton 4 | ΔE = 3.3e-16 → stop, **E = 2.599345963862454 rad** |
| 4 | x_p, y_p | −338 668 556.414253 km, 123 732 521.082095 km (= −2.263859471, 0.827100817 AU) |
| 4 | r | 360 563 625.282758 km = 2.410218966 AU |
| 5 | R | row 0: [ 0.983664841163, −0.094730559943, −0.153067309612] |
|   |   | row 1: [ 0.036745410222,  0.938092814141, −0.344429451245] |
|   |   | row 2: [ 0.176219338004,  0.333178620368,  0.926249832304] |
| 5 | r (km) | (−344 857 602.757556, 103 628 073.868074, −18 454 918.145304) |
| 5 | **r (AU)** | **(−2.305230690, 0.692710888, −0.123363508)** |
| 3 | v (km/s) | (−7.135097195, −12.277275312, −5.744460236) |

Reference `ctoc14.kepler.Ephemeris().ast_state(0, 1000*86400)`: r = (−2.30523069041, 0.692710888084,
−0.123363508177) AU, v = (−7.135097194791, −12.277275312477, −5.744460235796) km/s; difference 0.0 km in every
component. A JS implementation that prints these numbers for asteroid 1 at t = 1000 is correct.

Second reference value (Earth block, t = 0 d, Δt = 1826 · 86400 = 157 766 400 s): n = 1.988246716529e-07 rad/s,
M mod 2π = 6.189968562284, E = 6.188304310055, **r = (−0.130351646, 0.974736357, −0.000051052) AU**
(= `Ephemeris().earth_state(0)`). Earth at t = 1000 d: (1.004593546, 0.001658423, 0.000008937) AU.

## 5. Reference JavaScript (drop-in; matches steps 1–5 exactly)
```js
const D2R = Math.PI / 180, TWO_PI = 2 * Math.PI;
function makeBody(el, B) {            // el: bodies.json element block, B: bodies.json root
  const a = el.a * B.au_km, e = el.e, i = el.i * D2R, Om = el.Om * D2R, om = el.om * D2R;
  const cO = Math.cos(Om), sO = Math.sin(Om), ci = Math.cos(i), si = Math.sin(i), cw = Math.cos(om), sw = Math.sin(om);
  return { a, e, n: Math.sqrt(B.mu_sun_km3s2 / (a * a * a)), M0: el.M0 * D2R,
           off: (B.epoch_mjd_t0 - el.epoch_mjd) * 86400, se: Math.sqrt(1 - e * e), au: B.au_km,
           R00: cO * cw - sO * sw * ci, R01: -cO * sw - sO * cw * ci,
           R10: sO * cw + cO * sw * ci, R11: -sO * sw + cO * cw * ci,
           R20: sw * si,                R21: cw * si };
}
function solveKepler(M, e) {
  M = ((M % TWO_PI) + TWO_PI) % TWO_PI;
  let E = e < 0.8 ? M + e * Math.sin(M) : Math.PI;
  for (let k = 0; k < 50; k++) {
    const dE = (E - e * Math.sin(E) - M) / (1 - e * Math.cos(E));
    E -= dE;
    if (Math.abs(dE) < 1e-12) break;
  }
  return E;
}
function positionAU(b, tDays, out) {   // out: [x, y, z] in AU, heliocentric ecliptic J2000
  const E = solveKepler(b.M0 + b.n * ((tDays + 0) * 86400 + b.off), b.e);
  const xp = b.a * (Math.cos(E) - b.e), yp = b.a * b.se * Math.sin(E);
  out[0] = (b.R00 * xp + b.R01 * yp) / b.au;
  out[1] = (b.R10 * xp + b.R11 * yp) / b.au;
  out[2] = (b.R20 * xp + b.R21 * yp) / b.au;
  return out;
}
```
Orbit lines: sample E uniformly on [0, 2π) (e.g. 256 points), apply steps 4–5 directly; no Kepler solve needed.

## 6. Moon (illustrative only, not part of the contest model)
`bodies.json.moon` is a mean geocentric orbit: a_km = 384400, e = 0.0549, i = 5.145°, Ω = 125.08°, ω = 318.15°,
M0 = 135.27° at MJD 51544.5 (J2000), period 27.321582 d, with the node regressing at −0.05295°/day and the perigee
advancing at +0.11140°/day. Propagate with the same steps but:
- n = 2π / (27.321582 · 86400) rad/s (do NOT use μ☉; the orbit is geocentric),
- Δt = (t + (62502.0 − 51544.5)) · 86400 s, and use Ω(t) = Ω + Ω̇ · Δt_days, ω(t) = ω + ω̇ · Δt_days (degrees) when
  building R, so R is rebuilt per frame for the Moon only,
- the result is a geocentric offset in km; divide by AU and add the Earth position from section 2. The inclination
  is referred to the ecliptic, so no extra frame rotation. Magnitude 0.00257 AU: invisible at solar-system scale,
  which is why viz_spec section 3 has the Earth–Moon zoom mode (any exaggeration must be stated in the legend).

## 7. Files
- `results/viz/data/bodies.json` (48 902 bytes): written by `results/viz/data/make_bodies.py` from MEA.txt (via
  `ctoc14.kepler.load_mea`) and `ctoc14.constants` (Earth elements, μ, AU, epochs); values are the MEA.txt decimals
  unchanged (7 decimals), Earth elements as in docs/problem_summary.md.
- `results/viz/data/check_bodies.py`: the verification described at the top (prints a table of reference positions
  and per-point errors; `RESULT PASS max_check_err_au=2.041120e-15`).
