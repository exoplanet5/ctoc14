# CTOC14 validator thrust model — precise specification

Source: `CTOC14_problem.pdf` §6.1 rules 2–7 and §7 (校验说明). Verified numerically against the §6.2
example (see "Evidence" at the end). Constants: `Isp = 4000 s`, `g0 = 9.80665 m/s²`, `T_max = 0.5 N`,
`Isp·g0 = 39226.6 m/s`.

## 1. What is an arc

* Lines of one spacecraft are sorted by `Time` (strictly increasing).
* A **thrust arc** is a maximal run of consecutive `Event=1` lines. `Event=3` (flyby) lines that fall
  between two `Event=1` lines of the run **do not break the arc** and are **not** samples.
  `Event=0`, `Event=2`, `Event=4` lines terminate an arc (an `Event=2` placed between two `Event=1` lines
  splits them into two arcs — never do this inside an intended arc).
* The arc's samples are `(t_k, T_k)`, `k = 0..n-1` (0-based), `t_k` = the line's `Time`,
  `T_k = (Tx, Ty, Tz)` in newtons. Constraint: `t_{k+1} - t_k >= 8640 s` for every k (Event=3 lines in
  between are ignored for this spacing check), and `|T_k| <= 0.5`.
* An `Event=3` line is "inside" the arc iff `t_0 < t_flyby < t_{n-1}`. An `Event=3` line just before
  `t_0` or just after `t_{n-1}` is in coast (the example puts a flyby 10 s before `t_0`; that is legal).

## 2. Thrust vector at an arbitrary time

```
function THRUST(arc, t):                        # arc.t[0..n-1], arc.T[0..n-1] (vectors, N)
    n = arc.n
    if t < arc.t[0] or t > arc.t[n-1]:          # outside the arc: instantaneous on/off
        return (0, 0, 0)
    if n < 4:                                   # use all samples, degree n-1
        idx = [0 .. n-1]                        #   n=1: constant on the single instant t_0
    else:                                       #   n=2: linear, n=3: quadratic
        j = largest k in [0, n-2] with arc.t[k] <= t     # interval index, 0-based
                                                # (any t == t_k gives the same value from either
                                                #  neighbouring window because the polynomial
                                                #  interpolates the node exactly)
        l = clamp(j - 1, 0, n - 4)              # window start, 0-based
        idx = [l, l+1, l+2, l+3]
    # Lagrange interpolation, applied to EACH CARTESIAN COMPONENT separately
    T = (0, 0, 0)
    for i in idx:
        L_i = 1
        for k in idx, k != i:
            L_i *= (t - arc.t[k]) / (arc.t[i] - arc.t[k])
        T += L_i * arc.T[i]
    return T                                    # vector; NOT re-normalised, NOT clipped
```

Windows in words (n ≥ 4, 0-based):

| interval `[t_j, t_{j+1}]` | window samples            | character                   |
|---------------------------|---------------------------|-----------------------------|
| j = 0                     | 0,1,2,3                   | edge (interval at window start) |
| 1 ≤ j ≤ n−3               | j−1, j, j+1, j+2          | centred (two nodes each side)   |
| j = n−2                   | n−4, n−3, n−2, n−1        | edge (interval at window end)   |

The interpolant is a cubic polynomial per component on every interval; it is continuous at the nodes
but its derivative is not (different windows on the two sides), except at nodes interior to a shared
clipped window.

## 3. Mass flow and dynamics used by the validator

```
for each pair of consecutive lines (A, B) of one spacecraft:      # every line, all event types
    state  = (r_A, v_A, m_A)   taken from line A (NOT from the previous integration)
    arc    = the thrust arc containing (t_A, t_B), or none
    integrate from t_A to t_B:
        T      = THRUST(arc, t)  if arc else (0,0,0)
        r''    = -mu * r/|r|^3  +  (T / m) / 1000        # T/m in m/s^2 -> km/s^2
        m'     = -|T| / (Isp*g0)                          # norm of the INTERPOLATED VECTOR
        check |T| <= 0.5 at every evaluation ("全程检查") -> violation if exceeded
    require |r(t_B) - r_B| <= 1 km, |v(t_B) - v_B| <= 1 m/s, |m(t_B) - m_B| <= 0.01 kg
```

Notes:

* Because every sample is itself a line, each validator segment lies inside a single interpolation
  interval (or is split by an Event=3 line inside it). The thrust is therefore a smooth polynomial on
  every segment; the only discontinuities (arc start/end) sit exactly at segment boundaries.
* The mass on an `Event=3` line inside an arc must reflect the burn up to that instant; its thrust
  columns must nevertheless be `(0,0,0)` (they are placeholders, the validator uses the interpolant).
* `m' = -|T(t)|/(Isp g0)` uses the norm of the interpolated vector. For a rotating direction this is
  slightly *smaller* than the interpolated magnitude (chord effect); do not integrate mass with a
  magnitude-only model.
* The check is segment-wise from the *submitted* state, so errors do not accumulate along the file.

## 4. Index-base ambiguity and its resolution

The PDF writes the samples as `t_1 < … < t_n` (1-based) but clips the window start to `[0, n−4]`
(only meaningful 0-based). Two readings:

* **Reading A (0-based, centred)** — as coded above: `j` is 0-based, window `{j−1, j, j+1, j+2}`
  clipped to the arc. This is what a numpy-style implementation
  (`j = searchsorted(t, tt, 'right') - 1; l = clip(j-1, 0, n-4)`) produces, and the lower clip bound
  `0` is needed precisely for `j = 0`.
* **Reading B (1-based j, forward window)** — `j` 1-based, `l = j−1` as a 0-based offset, window
  `t_j … t_{j+3}` (the interval's left node plus three nodes to the right), clipped at the end.

Both give identical results on the first interval of an arc; they differ from the second interval on.

**Resolution: Reading A is confirmed by the PDF example.** The §6.2 arc was reconstructed exactly
(41 samples, 86 400 s apart, `|T_k| = 0.5·sin(πk/40)`, direction rotating 9°/sample about +z at
elevation 17.188734°, azimuth 40.107046° + 9°·k; the reconstruction reproduces the printed samples
k=1,2 to 2e-13 N). Integrating line 4 → line 5 (the *second* interval, where A and B differ):

| reading | mass error vs printed 1999.8273580 kg | position error |
|---------|----------------------------------------|----------------|
| A       | 2e-10 kg                               | 3 m            |
| B       | 1.5e-5 kg                              | 57 m           |

Only A matches the printed 11 significant digits (line 4, the first interval, is reproduced by both
to 2e-9 kg / 3 m). Both readings would pass the 0.01 kg / 1 km tolerances on this gentle profile, so a
solver using Reading A is safe even in the unlikely event the validator uses B, provided profiles are
smooth; for abrupt profiles the difference between readings can reach several km per segment, which
is one more reason to keep sample-to-sample changes small.

## 5. Overshoot facts (cubic Lagrange, uniform spacing, confirmed windows)

Lebesgue constants: centred interval 1.25, edge intervals 1.631. Hence the only *a-priori* bound is
`|T(t)| <= 1.25·max|T_k|` (centred) / `1.631·max|T_k|` (edge) — far too conservative to design with;
use the exact check of §6 instead. Measured peak `|T(t)|` for magnitude patterns at 0.5 N with fixed
direction:

| sample pattern (N)                                   | peak |T(t)|      |
|------------------------------------------------------|------------------|
| constant `[.5,.5,.5,…]`, instant on/off at arc ends  | 0.5000 (exact)   |
| half-sine bump `0.5·sin(πk/K)` (peak is a sample)    | 0.5000 (exact)   |
| hard step inside arc `[…,0,0,.5,.5,…]`               | 0.5321 (+6.4 %)  |
| hard step at arc start `[0,.5,.5,…]`                 | 0.5321 (+6.4 %)  |
| hard step-off `[…,.5,.5,0,0,…]`                      | 0.5321 (+6.4 %)  |
| two-sample plateau `[0,.5,.5,0]`                     | 0.5625 (+12.5 %) |
| linear ramp, 1/2/3/5 intermediate samples → plateau  | +3.2 / +2.1 / +1.6 / +1.1 % |
| cosine ramp over 2/3/4 intervals → plateau           | +3.2 / +1.6 / +0.9 % |
| constant 0.5 N, direction rotating θ per sample      | centred intervals: ≤ 0.5 (1−O(θ⁴)); edge intervals: +1.9e-5 N (10°), +2.9e-4 N (20°), +1.4e-3 N (30°), +5.8e-3 N = +1.2 % (45°) |
| constant 0.5 N, sudden turn of 60/90/120/180°        | +3.4 / +6.6 / +9.8 / +12.8 % |

Take-aways: (i) any plateau reached *through a ramp* overshoots at the ramp/plateau junction;
(ii) a constant-magnitude arc that switches on/off instantly at the arc ends does **not** overshoot
(the constant is reproduced exactly, and the rule already makes thrust zero outside `[t_0, t_{n-1}]`);
(iii) a bump whose maximum is a sample and which is concave around it does not overshoot;
(iv) direction changes at full magnitude overshoot only in the two edge intervals of an arc and only
noticeably for > ~15° per sample.

## 6. Safe design choice (recommended)

1. **Generate the submitted states with exactly the model of §2–3** (Reading A, componentwise cubic,
   `m'` from the interpolated vector norm, thrust zero outside the arc). Any other thrust
   parametrisation (piecewise-constant, spherical-angle interpolation, magnitude interpolation)
   deviates by ~10 km per 86 400 s segment at 0.5 N / 2000 kg for a 1 % model difference — i.e. fails
   the 1 km check. Integrate each segment to ≤ 1e-3 km so the whole 1 km tolerance is left for the
   validator's own integrator.
2. **Prefer constant-magnitude arcs with instant on/off** (bang-bang is what fuel/time-optimal
   low-thrust control gives anyway): all samples `|T_k| = T_c`, direction sampled from the smooth
   optimal steering law at spacing ≥ 8640 s chosen so that the turn per sample is ≤ ~10°.
   Do not encode a switch by ramp samples inside an arc; end the arc (next line `Event=2`) and start
   a new arc instead. For partial-throttle needs use separate constant-magnitude arcs.
3. **Cap the sample magnitude below 0.5 N**, e.g. `T_c = 0.4999 N` (0.02 % loss), so that (a) the tiny
   rotation-induced edge-interval overshoot (≤ 2e-5 N at 10°/sample) and (b) decimal rounding of the
   printed components (`|T|` of a rounded triple can exceed 0.5 by ~1e-10) can never trip the hard
   `<= 0.5` test.
4. **Prove it before writing the file** — exact per-interval maximisation, which is cheap because the
   interpolant is polynomial:

```
for each arc, for each interval j:
    build the three cubic polynomials p_x, p_y, p_z on [t_j, t_{j+1}] (window of §2, in a local
        variable τ = (t - t_j)/(t_{j+1} - t_j) for conditioning)
    S(τ) = p_x² + p_y² + p_z²                # degree 6
    candidates = {0, 1} ∪ { real roots of S'(τ) (degree 5) inside (0,1) }
    peak_j = max over candidates of sqrt(S)
    require peak_j <= 0.5 - margin           # margin ≥ 1e-6 N recommended
also: every sample |T_k| <= 0.5 - margin; every Event=0/2/3/4 line thrust == (0,0,0) exactly;
      t_{k+1} - t_k >= 8640 s within each arc (ignoring Event=3 lines);
      re-parse the written text file and repeat both checks on the parsed (rounded) numbers.
```

   If a peak exceeds the cap, scale all samples of that arc by `(0.5 - margin)/peak` (the interpolant
   is linear in the samples, so the scaled arc is provably compliant), then re-integrate and re-target.
5. Print all floats with `%.17g` (10 significant digits, the PDF's minimum suggestion, leaves only
   ~0.5 km of slack on 1e8-second coasts through velocity rounding alone).
6. Optional double check with an independent brute-force grid (e.g. 1000 points per interval); it must
   agree with the exact peak to ~1e-6 N.

## 7. Evidence (numbers reproduced from the PDF example with this model)

* Line 1 position vs Table-3 Earth at t=0: 0.5 m; launch `v_inf` = 4.000000 km/s.
* Line 1 → line 2, two-body coast of 1.0008e7 s: 3.2 m, 5.6e-10 km/s.
* Line 2 (`Event=3`, asteroid 174 at t = 1.0008479152e7 s): distance 9 m, relative speed 11.78 km/s.
* Line 3 → 4 (first interval): 2.9 m, 5e-10 km/s, 2e-9 kg (both readings).
* Line 4 → 5 (second interval): Reading A 2.8 m / 5e-10 km/s / 2e-10 kg; Reading B 57 m / 1.3e-6 km/s / 1.5e-5 kg.
* Exact peak of the example arc's interpolant: 0.500000000 N (peak sample 0.5 N at k=20, sine bump).

Scripts used (scratch, not part of the repo): `verify.py`, `overshoot.py` in the session scratchpad.
