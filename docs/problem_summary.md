# CTOC14 Problem A (题目甲) — Problem Summary

**第十四届全国空间轨道设计竞赛 (CTOC14) 题目甲: 近地小行星防御多航天器遍历探测轨道设计与优化**
Setter: 太空系统运行与控制全国重点实验室 (State Key Laboratory of Space System Operation and Control).

## Goal
Design N ≥ 1 continuous-low-thrust spacecraft launched from Earth that **fly by** (pass within 1000 km of) as many
as possible of 300 catalogued near-Earth asteroids (MEA.txt) inside a fixed 15-year window, minimizing total cost J.

## Decision variables
- Number of spacecraft N ≥ 1.
- Per spacecraft i: launch epoch t_launch,i ≥ t0; initial mass m0,i = 600 + m_fuel,i ≤ 2000 kg;
  hyperbolic excess vector v∞ (|v∞| ≤ 4 km/s, any direction); thrust history T(t) (0 ≤ |T| ≤ 0.5 N);
  asteroid visit sequence.

## Frame, dynamics, ephemerides
- Heliocentric Ecliptic Inertial J2000 (HEI). Units in submission: km, km/s, s, kg, N.
- r'' = −μ r/r³ + (T/m) u ; ṁ = −T/(Isp g0). Isp = 4000 s, g0 = 9.80665 m/s², Tmax = 0.5 N.
- μ☉ = 1.32712440018e11 km³/s², AU = 149597870.7 km.
- Earth and asteroids move on **fixed Keplerian orbits** (two-body). M(t) = M0 + n (t − t_eph), n = sqrt(μ/(a·AU)³) [rad/s].
- Earth elements (epoch MJD 60676.0 = 2025-01-01 00:00 UTC): a=1.0009175020 AU, e=0.017566762041, i=0.002976847126°,
  Ω=189.953211282428°, ω=273.196254000254°, M0=357.4135031077°.
- Asteroid elements: MEA.txt, columns ID a[AU] e i Ω ω M [deg], common epoch MJD 61200.0 (2026-06-09 00:00 UTC).
- No gravity assists allowed (planet or asteroid). No Earth return required.

## Constraints
- Mission start t0 = MJD 62502.0 (2030-01-01 00:00 UTC). All events (launch, flybys, end) in 0 ≤ t ≤ 15 yr = 4.73364e8 s = 5478.75 d.
- Launch: r_sc = r_Earth(t_launch) (tolerance 1 km), |v_sc − v_Earth| ≤ 4 km/s (tolerance 0.01 km/s).
- Mass: m_dry = 600 kg; m(t) ≥ 600 kg always; m0 ≤ 2000 kg.
- Flyby detection: d ≤ 1000 km between spacecraft and asteroid at the declared Event=3 time. Only first detection of each asteroid counts. Repeats not penalized.

## Cost (lower is better)
- J_i = 1 + x_i + x_i²,  x_i = (m0,i − 600)/1400 ∈ [0,1]   → empty spacecraft costs 1, full tank costs 3.
- J = Σ_i J_i + N_miss,  N_miss = 300 − N_covered (each missed asteroid costs exactly 1 = one empty launch).
- Time coefficient k = 1 − 0.1 (t_end − t_submit)/(t_end − t_start) ∈ [0.9, 1.0]; J_final = k·J. Earlier submission is better.

## Submission file (CTOC14_Result_TeamID.txt), 15 columns per line
Line SC_ID Event Time[s] x y z[km] vx vy vz[km/s] m[kg] Tx Ty Tz[N] Asteroid_ID
- Event: 0 launch, 1 thrust sample, 2 coast node, 3 asteroid flyby, 4 end. Per spacecraft: first line Event=0, last Event=4, Time strictly increasing.
- A run of consecutive Event=1 lines is one thrust arc; Event=3 lines may be inserted inside an arc without breaking it; arcs are separated by Event 0/2/4.
- Thrust between samples: **3rd-order sliding-window Lagrange interpolation** (window of 4 samples: j-1..j+2, clipped to [0, n-4]); if n<4 use all n samples ((n-1)-th order). Thrust is zero outside [t1, tn] (instant on/off at arc ends).
- Adjacent Event=1 samples ≥ 8640 s (0.1 d) apart. |T(t)| ≤ 0.5 N at ALL times incl. between samples (interpolation overshoot is a violation).
- Event 0/2/3/4 lines carry (0,0,0) thrust. Mass follows ṁ = −|T(t)|/(Isp g0) (tolerance 0.01 kg).
- Validator integrates each line's state to the next line with the interpolated thrust: position error ≤ 1 km, velocity ≤ 1 m/s, mass ≤ 0.01 kg. Output ≥ 10 significant digits.

## Derived capability numbers
- Exhaust velocity Isp·g0 = 39.227 km/s. Full tank Δv = 39.227·ln(2000/600) = 47.2 km/s.
- Mass flow at Tmax = 1.2746e-5 kg/s → full tank (1400 kg) lasts 1.098e8 s = 1271 d = 3.48 yr of continuous full thrust (out of 15 yr).
- Acceleration 0.25 mm/s² (2000 kg) → 0.83 mm/s² (600 kg). ≈ 7.9 km/s per year of continuous thrust at 2000 kg.
- Ballistic spacecraft (cost 1): launch time + v∞ vector = 4 DOF; each flyby = 3 position constraints − 1 free time = 2 net → at most 2 asteroids generically.
