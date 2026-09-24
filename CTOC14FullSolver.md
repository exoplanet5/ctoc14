You are an expert in spacecraft trajectory design, astrodynamics, optimal control, combinatorial optimization, and scientific computing.

Your task is to help me systematically solve **Problem A of the 14th China Trajectory Optimization Competition (CTOC14)**:

**“Multi-Spacecraft Tour Design and Optimization for Near-Earth Asteroid Defense Reconnaissance.”**

I will provide two source files:

1. `CTOC14_problem.pdf`
2. `MEA.txt`

These two files are the authoritative definition of the problem. Always follow their exact dynamics, constraints, parameters, scoring rules, ephemeris definitions, event definitions, and submission format. Do not silently replace any problem definition with external assumptions.

The final goal is not merely to produce a feasible baseline solution, but to build a complete, modular, reproducible solver architecture with the potential to achieve a competitive CTOC14 score.

---

# 1. Optimization objective

Within the fixed 15-year mission window starting from 2030-01-01 00:00:00 UTC, design a set of continuously thrusting spacecraft to perform flyby reconnaissance of as many as possible of the 300 given near-Earth asteroids, while minimizing

\[
J=\sum_{i=1}^{N}J_i+N_{\rm miss},
\]

where

\[
J_i=1+x_i+x_i^2,
\]

and

\[
x_i=\frac{m_{0,i}-600}{1400}.
\]

Here:

- \(N\) is the number of spacecraft,
- \(m_{0,i}\) is the initial mass of spacecraft \(i\),
- \(N_{\rm miss}\) is the number of asteroids not successfully visited.

The official final score also includes the submission-time coefficient specified in the problem statement, but trajectory optimization should initially minimize the raw mission cost \(J\).

---

# 2. Exact problem model

## 2.1 Reference frame

All position and velocity states are expressed in the:

**Heliocentric Ecliptic Inertial frame, HEI**

with:

- origin at the Sun barycentric center used by the problem,
- \(x-y\) plane equal to the J2000.0 ecliptic,
- \(x\)-axis pointing toward the J2000.0 vernal equinox,
- \(z\)-axis completing a right-handed system.

Use:

- position in km,
- velocity in km/s,
- internal integration time in seconds.

---

# 2.2 Spacecraft dynamics

The spacecraft follows heliocentric two-body motion with continuous thrust:

\[
\ddot{\mathbf r}
=
-\frac{\mu_\odot}{r^3}\mathbf r
+
\frac{\mathbf T}{m}.
\]

The thrust magnitude must satisfy

\[
0\leq \|\mathbf T\|\leq0.5\ {\rm N}.
\]

The mass flow is

\[
\dot m
=
-\frac{\|\mathbf T\|}{I_{sp}g_0}.
\]

Use the exact constants from the problem statement, including:

\[
\mu_\odot
=
1.32712440018\times10^{11}\ {\rm km^3/s^2},
\]

\[
I_{sp}=4000\ {\rm s},
\]

\[
g_0=9.80665\ {\rm m/s^2},
\]

\[
T_{\max}=0.5\ {\rm N}.
\]

Be extremely careful with unit conversion between N, kg, m/s², and km/s².

---

# 2.3 Spacecraft mass constraints

Dry mass:

\[
m_{\rm dry}=600\ {\rm kg}.
\]

Initial mass:

\[
600\leq m_0\leq2000\ {\rm kg}.
\]

Fuel mass:

\[
m_{\rm fuel}=m_0-600.
\]

At all times:

\[
m(t)\geq600\ {\rm kg}.
\]

---

# 2.4 Earth departure condition

Mission epoch:

2030-01-01 00:00:00 UTC,

corresponding to

\[
{\rm MJD}=62502.0.
\]

Each spacecraft may have its own launch time:

\[
t_{\rm launch}\geq0,
\]

provided all events remain within the 15-year mission window.

At launch:

\[
\mathbf r_{\rm SC}
=
\mathbf r_\oplus.
\]

The Earth-relative hyperbolic excess velocity satisfies:

\[
\left\|
\mathbf v_{\rm SC}
-
\mathbf v_\oplus
\right\|
\leq4\ {\rm km/s}.
\]

The direction is unrestricted.

This \(v_\infty\) does not enter the spacecraft fuel cost directly and should therefore be treated as an extremely valuable free initial velocity control variable.

The solver should exploit this aggressively.

---

# 2.5 Earth and asteroid ephemerides

All celestial bodies use the heliocentric two-body propagation model defined in the problem.

Earth orbital elements are given in the problem PDF at epoch:

\[
{\rm MJD}=60676.0.
\]

The 300 asteroid orbital elements are given in:

`MEA.txt`

with common epoch:

\[
{\rm MJD}=61200.0.
\]

Implement propagation from classical orbital elements:

\[
(a,e,i,\Omega,\omega,M_0)
\]

to Cartesian states:

\[
(\mathbf r,\mathbf v)
\]

at arbitrary epochs.

Implement a numerically robust Kepler equation solver.

Do not substitute JPL Horizons or any external ephemeris for the competition ephemeris.

---

# 2.6 Asteroid reconnaissance condition

The mission requires only a flyby.

Asteroid \(j\) is considered successfully visited if

\[
\left\|
\mathbf r_{\rm SC}(t)
-
\mathbf r_j(t)
\right\|
\leq1000\ {\rm km}.
\]

Critically, there is no velocity-matching constraint.

Do not impose

\[
\mathbf v_{\rm SC}
=
\mathbf v_j.
\]

This is a moving-target flyby problem, not a rendezvous problem.

Asteroids do not provide gravity assists.

The spacecraft state remains continuous through a flyby.

Repeated visits to the same asteroid do not provide additional credit.

---

# 2.7 Mission duration

All events must satisfy

\[
0\leq t\leq5478.75\ {\rm days}.
\]

This is a fixed absolute 15-year mission window.

Flybys after the end of the window do not count.

---

# 3. Interpretation of the scoring function

Each spacecraft costs

\[
J_i=1+x_i+x_i^2.
\]

Therefore:

- launching a dry spacecraft costs 1,
- launching a fully fueled 2000 kg spacecraft costs 3,
- missing one asteroid also costs 1.

The optimization objective is therefore not simply minimum fuel.

The real objective is:

**maximize useful asteroid coverage per spacecraft cost.**

For a route \(R\) covering \(K_R\) previously uncovered asteroids with estimated normalized fuel \(x_R\), define the route utility

\[
U(R)
=
K_R
-
(1+x_R+x_R^2).
\]

Use this idea in route generation, pruning, and route selection.

Do not prescribe the spacecraft count \(N\) in advance.

The global optimization should determine:

- how many spacecraft to launch,
- which asteroids belong to which spacecraft,
- which difficult asteroids should potentially be left uncovered.

---

# 4. Overall solver architecture

Do not formulate the complete problem immediately as one enormous nonlinear optimization involving:

- all 300 asteroid assignments,
- all launch times,
- all encounter times,
- all continuous thrust controls.

That formulation is too large, highly nonconvex, and likely to fail.

Instead use a hierarchical architecture:

\[
\boxed{
\text{Ephemeris}
\rightarrow
\text{Ballistic Atlas}
\rightarrow
\text{State-Expanded Graph}
\rightarrow
\text{Route Search}
\rightarrow
\text{Set Cover}
\rightarrow
\text{Low-Thrust Refinement}
\rightarrow
\text{Submission Validation}
}
\]

---

# 5. Phase 0 — core astrodynamics

First implement a clean astrodynamics foundation.

Suggested structure:

```text
ctoc14/
    constants.py
    ephemeris.py
    kepler.py
    lambert.py
    dynamics.py
    propulsion.py
    data.py
    tests/
```

Implement at least:

```python
state_from_elements(elements, mjd)
earth_state(mjd)
asteroid_state(asteroid_id, mjd)
propagate_two_body(...)
```

Requirements:

- double precision,
- consistent units,
- vectorized propagation where possible,
- efficient batch evaluation for all 300 asteroids.

---

# 6. Mandatory regression test

Use the spacecraft trajectory example contained in the official PDF.

The sample spacecraft departs Earth on 2030-01-01 and later performs a ballistic flyby of asteroid 174.

Reconstruct the corresponding Earth-to-174 Lambert transfer.

Verify:

- Earth state,
- asteroid 174 state,
- epoch conversion,
- HEI orientation,
- Lambert solution,
- Earth-relative launch \(v_\infty\).

The reconstructed solution should reproduce approximately

\[
v_\infty\simeq4\ {\rm km/s}.
\]

This must become an automated regression test.

If this test fails, stop and fix the underlying dynamics before proceeding.

---

# 7. Phase 1 — Earth-to-asteroid ballistic atlas

Initially ignore electric propulsion.

Search ballistic transfers of the form

\[
Earth(t_L)\rightarrow Asteroid_j(t_A).
\]

Decision variables:

\[
t_L,\quad t_A.
\]

Constraints:

\[
t_A>t_L,
\]

and

\[
v_\infty
=
\left\|
\mathbf v_{\rm Lambert}(t_L)
-
\mathbf v_\oplus(t_L)
\right\|
\leq4\ {\rm km/s}.
\]

Consider reasonable Lambert branches:

- short-way,
- long-way,
- optionally multi-revolution later.

The first implementation may use only zero-revolution Lambert solutions.

---

# 8. Ballistic-atlas search strategy

Use a coarse search first.

Example launch-time spacing:

30–90 days.

Example time-of-flight range:

50–1500 days.

For each asteroid, identify promising regions.

Then locally refine the best windows by optimizing

\[
\min_{t_L,t_A}v_\infty.
\]

Do not store only one best window per asteroid.

Store approximately 50–200 distinct promising arrival states per asteroid if computationally feasible.

Define:

```python
class ArrivalState:
    asteroid_id
    launch_time
    arrival_time
    r_arr
    v_arr_sc
    v_arr_ast
    vinf
```

The critical state variable is

\[
\mathbf v_{\rm SC}(t_A),
\]

not the asteroid velocity.

---

# 9. State-expanded graph

The graph must not simply use asteroid identity as a node.

For the same asteroid \(A\), different flyby epochs and spacecraft arrival velocities can produce completely different future opportunities.

Represent nodes as:

\[
(A,t_A,\mathbf v_{\rm SC}).
\]

Thus transitions are of the form:

\[
(A,t,\mathbf v)
\rightarrow
(B,t',\mathbf v').
\]

This is a state-expanded temporal graph.

---

# 10. Phase 2 — asteroid-to-asteroid impulsive surrogate

Given a spacecraft state at asteroid \(A\),

\[
S_A=
(A,t_A,\mathbf r_A,\mathbf v_{\rm SC,A}),
\]

enumerate future targets \(B\) and future epochs \(t_B>t_A\).

Solve a Lambert transfer from

\[
\mathbf r_A(t_A)
\]

to

\[
\mathbf r_B(t_B).
\]

Let the Lambert departure velocity be

\[
\mathbf v_L.
\]

Define an impulsive surrogate cost

\[
\Delta v_{\rm proxy}
=
\left\|
\mathbf v_L-\mathbf v_{\rm SC,A}
\right\|.
\]

Use the Lambert arrival velocity as the next spacecraft state:

\[
\mathbf v_{\rm SC,B}.
\]

This generates candidate transitions

\[
S_A\rightarrow S_B.
\]

The quantity \(\Delta v_{\rm proxy}\) is not the true optimal continuous-thrust cost.

It is only a fast heuristic used for:

- ranking,
- pruning,
- route generation,
- initial guesses.

Clearly label it as a surrogate.

---

# 11. Graph-edge pruning

Do not construct a dense all-to-all graph.

For each state, retain only a limited number of promising successor states.

Use a score based on quantities such as:

\[
\Delta v_{\rm proxy},
\]

flight time,

remaining mission duration,

orbital geometry,

future target density,

and number of reachable successors.

For example:

\[
C_{AB}
=
w_1\Delta v_{\rm proxy}
+
w_2\Delta t
-
w_3N_{\rm future}.
\]

Retain only roughly 20–100 best outgoing transitions per state.

---

# 12. Continuous-thrust feasibility surrogate

Use

\[
a_T=\frac{T_{\max}}{m}
\]

to estimate whether an impulsive transition is plausible under finite thrust.

A simple first-order estimate is

\[
t_{\rm thrust}
\approx
\frac{\Delta v_{\rm proxy}}{T/m}.
\]

Fuel can also be approximated using

\[
m_f
=
m_i
\exp
\left(
-\frac{\Delta v}{I_{sp}g_0}
\right).
\]

These approximations should only be used for:

- pruning,
- route ranking,
- initial fuel estimates.

The final trajectory must be optimized under the true continuous-thrust dynamics.

---

# 13. Phase 3 — route generation

Search routes of the form

\[
Earth
\rightarrow
A_1
\rightarrow
A_2
\rightarrow
\cdots
\rightarrow
A_K.
\]

Implement at least one of:

1. Beam Search
2. GRASP
3. Adaptive Large Neighborhood Search, ALNS

Start with Beam Search.

Each search state should contain approximately:

```python
class RouteState:
    time
    position
    velocity
    estimated_mass
    visited_mask
    asteroid_sequence
    encounter_times
    launch_time
    vinf_vector
    estimated_cost
```

---

# 14. Beam-search objective

Do not rank routes only by accumulated \(\Delta v\).

Possible score:

\[
F
=
-N_{\rm visited}
+
\lambda_1J_{\rm est}
+
\lambda_2\Delta v_{\rm remaining}
+
\lambda_3\frac{t}{T_{\max}}
-
\lambda_4N_{\rm future}.
\]

Alternatively maximize

\[
F
=
N_{\rm visited}
-
\lambda J_{\rm spacecraft}
+
\eta N_{\rm future}.
\]

Experiment with different weights.

The objective is to discover high-target-density routes.

---

# 15. Explicit search for ballistic chains

Implement a dedicated search for nearly ballistic routes such as

\[
Earth\rightarrow A\rightarrow B\rightarrow C.
\]

Because Earth departure \(v_\infty\leq4\ {\rm km/s}\) is free, a dry 600 kg spacecraft that visits two or more asteroids with no propulsion can be highly valuable.

Implement:

```python
search_ballistic_chains()
```

The search can parameterize departure by:

\[
t_L,\quad
\mathbf v_\infty
\]

or use Lambert-derived initial states.

---

# 16. Optional Monte Carlo ballistic screening

Also implement a Monte Carlo mode.

Randomly generate:

\[
t_L
\]

and

\[
\mathbf v_\infty,
\qquad
\|\mathbf v_\infty\|\leq4\ {\rm km/s}.
\]

Propagate the resulting heliocentric ballistic spacecraft orbit through the entire mission.

Precompute asteroid positions at time samples.

Use spatial search methods such as:

- KD-tree,
- Ball tree,
- spatial hashing.

Search for near misses such as

\[
d<10^5\ {\rm km}
\]

or even

\[
d<10^6\ {\rm km}.
\]

Promising near misses should then be locally optimized until the miss distance satisfies

\[
d\leq1000\ {\rm km}.
\]

---

# 17. Phase 4 — route pool

Build a large database of candidate routes.

For example:

```python
class Route:
    id
    launch_time
    vinf
    asteroid_sequence
    encounter_times
    estimated_fuel
    estimated_initial_mass
    estimated_cost
    asteroid_mask
```

Aim for at least thousands of routes.

If computationally feasible, generate

\[
10^4-10^6
\]

candidate routes before dominance pruning.

---

# 18. Route dominance pruning

If two routes \(R_1,R_2\) satisfy

\[
S_{R_1}\supseteq S_{R_2}
\]

and

\[
c_{R_1}\leq c_{R_2},
\]

then \(R_2\) is dominated and can be discarded.

For identical asteroid sets, retain a small Pareto set emphasizing:

- minimum fuel,
- earliest completion,
- best future extensibility.

---

# 19. Phase 5 — global route selection

Formulate route selection as a weighted set-cover or prize-collecting set-cover problem.

Binary variable

\[
y_r\in\{0,1\}
\]

indicates whether route \(r\) is selected.

Binary variable

\[
z_j\in\{0,1\}
\]

indicates whether asteroid \(j\) is missed.

Solve:

\[
\min
\sum_r c_r y_r
+
\sum_j z_j.
\]

Subject to

\[
\sum_{r:j\in S_r}y_r+z_j\geq1.
\]

Potential solvers:

- `scipy.optimize.milp`
- HiGHS
- OR-Tools
- Gurobi if available

Prefer an implementation that does not require commercial software.

---

# 20. Route overlap is allowed

Multiple spacecraft may visit the same asteroid.

Repeated visits do not provide additional credit, but they are not prohibited.

Therefore route selection is not an exact-cover problem.

Do not impose

\[
\sum_r y_r=1
\]

for each asteroid.

---

# 21. Hard-target analysis

Before route generation, analyze all 300 asteroid orbits statistically.

Identify targets with:

- high inclination,
- retrograde motion,
- very large semimajor axis,
- large heliocentric distance throughout 2030–2045,
- poor accessibility from Earth under the free \(v_\infty\) constraint.

Pay special attention to asteroid IDs:

- 131
- 144

Their orbital elements in `MEA.txt` make them obvious dynamical outliers relative to most of the target set.

Do not allow a very small number of extreme targets to distort the search strategy for the main asteroid population.

It is acceptable to maintain:

```text
main_target_pool
hard_target_pool
```

and design dedicated missions for hard targets.

---

# 22. Strategy for asteroid 131

Do not incorrectly require the spacecraft to match its retrograde velocity.

This remains a position intercept problem.

Investigate:

- early departure,
- substantial out-of-ecliptic motion,
- long-duration low thrust,
- intercept before the target moves even farther outward.

Solve an Earth-to-131 minimum-initial-mass problem with variables including:

\[
t_L,
\]

\[
\mathbf v_\infty,
\]

and continuous thrust history.

---

# 23. Strategy for asteroid 144

First propagate its heliocentric distance over the full 2030–2045 mission interval.

Then investigate a fast heliocentric escape/intercept trajectory:

\[
Earth\rightarrow144.
\]

Again, velocity matching is not required.

Determine whether interception is possible within 15 years.

If feasible, solve for

\[
m_{0,\min}.
\]

Then compare its spacecraft cost

\[
1+x+x^2
\]

with the cost of simply leaving the asteroid uncovered:

\[
N_{\rm miss}=1.
\]

A dedicated spacecraft costing more than 1 and visiting only one new asteroid is generally not advantageous.

A hard-target mission becomes interesting if it can also visit additional asteroids or otherwise improve the global score.

---

# 24. Phase 6 — true low-thrust refinement

Once an asteroid sequence is fixed,

\[
Earth
\rightarrow
A_1
\rightarrow
\cdots
\rightarrow
A_K,
\]

solve the true continuous-thrust trajectory optimization problem.

Preferred approaches:

- direct multiple shooting,
- direct collocation.

Start with direct multiple shooting.

---

# 25. Low-thrust decision variables

At minimum include:

\[
t_{\rm launch},
\]

\[
\mathbf v_\infty,
\]

encounter epochs

\[
t_1,\ldots,t_K,
\]

and a discretized thrust-control history.

For example, divide each leg into \(N_c\) control intervals with thrust vector

\[
\mathbf T_k.
\]

Impose

\[
\|\mathbf T_k\|\leq0.5\ {\rm N}.
\]

---

# 26. Flyby constraints

Do not impose rendezvous.

Incorrect:

\[
\mathbf r_{\rm SC}
=
\mathbf r_{\rm ast}
\]

and simultaneously

\[
\mathbf v_{\rm SC}
=
\mathbf v_{\rm ast}.
\]

Correct:

\[
\|
\mathbf r_{\rm SC}(t_k)
-
\mathbf r_{A_k}(t_k)
\|
\leq1000\ {\rm km}.
\]

For initial numerical stabilization, it is acceptable to use

\[
\mathbf r_{\rm SC}
=
\mathbf r_{A_k}
\]

as an equality condition.

After obtaining a feasible solution, release this to the actual 1000 km inequality.

---

# 27. Low-thrust objective

For a fixed asteroid route, minimize initial mass

\[
m_0
\]

or equivalently fuel mass.

Because

\[
J_i=1+x+x^2
\]

is monotonic in \(m_0\), minimizing fuel also minimizes that route's spacecraft score.

Always report the actual \(J_i\) afterward.

---

# 28. Initial guess generation

Low-thrust NLPs are highly sensitive to initial guesses.

Use the Phase-2 Lambert route to generate an initial control guess.

Convert each impulsive correction

\[
\Delta\mathbf v
\]

into a finite-duration thrust arc.

Approximate the burn duration as

\[
\tau
\sim
\frac{m\Delta v}{T}.
\]

Initialize thrust direction with

\[
\hat{\mathbf T}
=
\frac{\Delta\mathbf v}{\|\Delta\mathbf v\|}.
\]

Then allow the nonlinear optimizer to refine the entire trajectory.

---

# 29. Numerical tools

Use Python as the primary implementation language.

Useful packages include:

- NumPy
- SciPy
- Numba
- joblib
- multiprocessing
- pandas
- matplotlib

Potential optimization tools:

- `scipy.optimize`
- `scipy.optimize.least_squares`
- `scipy.optimize.minimize`
- `scipy.integrate.solve_ivp`

If available, strongly consider:

- CasADi
- IPOPT through CasADi
- cyipopt

for the low-thrust nonlinear optimization.

The baseline solver should not depend on commercial software.

---

# 30. Performance

Large-scale Lambert screening will likely be a computational bottleneck.

Use:

- vectorization,
- parallel processing,
- cached ephemerides,
- precomputed asteroid-state grids,
- Numba where appropriate.

Implement concepts such as:

```python
EphemerisCache
LambertCache
```

Avoid repeated propagation of identical states.

---

# 31. Internal time representation

Use a single internal time convention.

Recommended:

```python
t_sec
```

defined as seconds since MJD 62502.0.

Implement:

```python
t_to_mjd(t_sec)
mjd_to_t(mjd)
t_to_datetime(t_sec)
```

Do not mix MJD, days, and seconds inside core dynamics without explicit conversion.

---

# 32. Official submission thrust representation

The competition does not directly use the optimizer's internal control discretization.

The submitted Event=1 thrust samples are reconstructed using the official moving-window cubic Lagrange interpolation rule.

Therefore implement the exact competition interpolation separately:

```python
submission_thrust_interpolator()
```

It must reproduce the official definition exactly.

---

# 33. Thrust-sample constraints

Within one thrust arc, adjacent Event=1 samples must satisfy

\[
\Delta t\geq8640\ {\rm s}.
\]

The magnitude constraint must hold not only at samples:

\[
\|\mathbf T_k\|\leq0.5\ {\rm N},
\]

but also everywhere between samples after interpolation:

\[
\|\mathbf T(t)\|\leq0.5\ {\rm N}.
\]

Cubic Lagrange interpolation may overshoot.

Therefore do not automatically place thrust nodes exactly at 0.5 N.

Either use margin, for example

\[
0.47-0.49\ {\rm N},
\]

or explicitly maximize

\[
\|\mathbf T(t)\|
\]

inside each interpolation interval.

---

# 34. Independent submission validator

Implement an independent validator:

```python
validate_submission(path)
```

that reproduces the official checks as closely as possible.

Validate:

- file structure,
- SC_ID ordering,
- event sequence,
- strictly increasing time,
- launch state,
- \(v_\infty\),
- spacecraft mass,
- thrust limits,
- interpolated thrust,
- asteroid flybys,
- mission window,
- event fields.

From every submitted state row, numerically integrate the official dynamics to the next row using the official reconstructed thrust history.

Internally aim for tighter consistency than the official tolerances, for example:

position:

\[
<0.1\ {\rm km},
\]

velocity:

\[
<0.1\ {\rm m/s},
\]

mass:

\[
<0.001\ {\rm kg}.
\]

The official acceptance limits must be read from the PDF and used as the hard validation thresholds.

---

# 35. Output precision

Use sufficient floating-point precision in the final submission.

For example:

```python
"%.17g"
```

or

```python
"%.12e"
```

Do not sacrifice numerical consistency to reduce output file size.

---

# 36. Recommended software structure

A possible final project structure is:

```text
ctoc14/
│
├── data/
│   └── MEA.txt
│
├── constants.py
├── time_utils.py
├── kepler.py
├── ephemeris.py
├── lambert.py
├── dynamics.py
├── propulsion.py
│
├── ballistic_atlas.py
├── asteroid_graph.py
├── beam_search.py
├── monte_carlo.py
├── route_pool.py
├── set_cover.py
│
├── low_thrust/
│   ├── transcription.py
│   ├── shooting.py
│   ├── controls.py
│   └── refine.py
│
├── submission/
│   ├── interpolate.py
│   ├── export.py
│   └── validator.py
│
├── analysis/
│   ├── target_statistics.py
│   ├── hard_targets.py
│   ├── route_statistics.py
│   └── score_analysis.py
│
├── cache/
├── results/
├── tests/
│
└── main.py
```

---

# 37. First implementation milestones

Do not attempt to implement the entire solver in one step.

Proceed through the following milestones.

## Milestone A — target statistics

Read `MEA.txt`.

For all 300 asteroids compute:

- \(a\),
- \(e\),
- \(i\),
- perihelion distance \(q\),
- aphelion distance \(Q\),
- orbital period.

Generate plots including:

- \(a-e\),
- \(a-i\),
- \(q-i\),
- heliocentric-distance envelope over 2030–2045.

Automatically identify dynamical outliers.

---

## Milestone B — competition ephemeris

Implement the exact two-body ephemeris specified by CTOC14.

Validate Earth and asteroid states.

Use the official Earth-to-asteroid-174 example as a regression test.

---

## Milestone C — Lambert ballistic atlas

Search

\[
Earth\rightarrow\text{all 300 asteroids}
\]

for ballistic flyby windows.

For every stored solution output:

```text
asteroid_id
launch_time
arrival_time
time_of_flight
vinf
arrival_spacecraft_velocity
```

Store multiple distinct windows for each asteroid.

---

## Milestone D — accessibility statistics

Determine how many asteroids can be reached ballistically with

\[
v_\infty\leq4\ {\rm km/s}.
\]

Analyze targets by:

- minimum \(v_\infty\),
- earliest reachable epoch,
- number of launch windows,
- inclination,
- semimajor axis.

Identify hard targets.

---

## Milestone E — first multi-target routes

Construct asteroid-to-asteroid Lambert surrogate transitions.

Initially test on only the easiest 10–50 targets.

Implement Beam Search and find:

- 2-target routes,
- 3-target routes,
- 5-target routes,
- longer routes if possible.

---

# 38. Every stage must produce numerical evidence

Do not only write code.

After each module:

1. run it,
2. report actual numerical results,
3. quantify numerical errors,
4. inspect physical plausibility,
5. save results,
6. only then move to the next stage.

Examples:

```text
Earth -> asteroid 174 regression:
position error = ...
velocity error = ...
vinf = ...
```

and

```text
Ballistic atlas:
targets reachable = ...
Lambert solves attempted = ...
successful Lambert solutions = ...
CPU time = ...
```

---

# 39. Preserve failure information

Do not hide failures.

Log cases including:

- Lambert failure,
- optimizer failure,
- invalid time ordering,
- negative remaining fuel,
- thrust interpolation overshoot,
- shooting divergence,
- inability to satisfy a flyby constraint.

Use structured logging.

Failures are useful information for refining the search.

---

# 40. Key scientific questions after the baseline

Once the baseline is working, investigate:

### A

Can purely ballistic trajectories naturally encounter multiple asteroids?

### B

Can small thrust corrections connect one ballistic chain to another?

### C

Which asteroids behave as high-degree hubs in the temporal transfer graph?

### D

Which asteroids are isolated dynamical outliers?

### E

How many additional asteroids can typically be covered per extra 100 kg of fuel?

### F

Does the best mission favor:

- many low-cost spacecraft,
- an intermediate number of moderately fueled spacecraft,
- or a small number of heavily fueled long-tour spacecraft?

### G

Is full coverage actually optimal?

Since leaving one asteroid uncovered costs only

\[
1,
\]

some isolated asteroids may rationally be abandoned if they require an expensive dedicated mission.

---

# 41. Marginal-cost criterion

For a new spacecraft,

\[
c=1+x+x^2.
\]

If it covers only \(K\) previously uncovered targets, then it is strictly better than missing those \(K\) targets only if

\[
K>c.
\]

Use this criterion explicitly during:

- route generation,
- hard-target analysis,
- local search,
- spacecraft creation and deletion.

---

# 42. Local improvement

Once a complete mission solution exists, implement local-search operators such as:

- route merge,
- route split,
- asteroid insertion,
- asteroid removal,
- asteroid swap,
- subsequence replacement,
- launch-time shift,
- encounter-time shift.

Example:

SC1:

\[
A\rightarrow B\rightarrow C
\]

SC2:

\[
D\rightarrow E
\]

Try alternatives such as:

\[
A\rightarrow B\rightarrow D\rightarrow E
\]

and determine whether \(C\) can be inserted into another route.

Treat this as a vehicle-routing-style local optimization problem.

---

# 43. ALNS operators

At a later stage, implement Adaptive Large Neighborhood Search.

Destroy operators:

- random asteroid removal,
- expensive-asteroid removal,
- entire-route removal,
- hard-cluster removal.

Repair operators:

- cheapest insertion,
- regret-2 insertion,
- regret-\(k\) insertion,
- new-spacecraft creation.

Evaluate candidate solutions using the true score \(J\), or the most accurate available surrogate when full low-thrust refinement is too expensive.

---

# 44. Full coverage is not mandatory

Always retain

\[
N_{\rm miss}
\]

as an explicit optimization term.

Do not force all 300 targets to be visited.

If feasible, maintain both:

- best full-coverage solution,
- best unconstrained minimum-\(J\) solution.

---

# 45. Maintain an internal leaderboard

Every time a better mission is found, record:

```text
run_id
timestamp
J
N_spacecraft
N_covered
N_missed
total_fuel
max_fuel
mean_targets_per_spacecraft
max_targets_per_spacecraft
```

Also save:

```text
solution.json
routes.csv
submission.txt
```

Do not overwrite previous best solutions.

---

# 46. Reproducibility

Every stochastic method, including:

- Monte Carlo,
- random restarts,
- GRASP,
- ALNS,

must record its random seed.

Store configuration in something like:

```yaml
config.yaml
```

All important search parameters should be reproducible.

---

# 47. Visualization

Generate at least:

1. orbital-element distribution of all 300 asteroids,
2. ballistic accessibility histogram,
3. asteroid temporal transfer graph,
4. route-length distribution,
5. fuel versus number of visited asteroids,
6. spacecraft mission timelines,
7. 3D heliocentric trajectory plots for representative spacecraft,
8. encounter chronology,
9. overall target coverage,
10. best-score evolution versus computation time or iteration.

---

# 48. Distinguish exact quantities from surrogates

Always explicitly distinguish:

### Exact CTOC14 quantities

The quantities defined by the official dynamics and scoring rules.

### Surrogate quantities

Examples:

- Lambert \(\Delta v\),
- impulsive approximations,
- simple thrust-time approximations,
- approximate fuel estimates.

Never report surrogate quantities as if they were exact mission costs.

---

# 49. Central modeling insight

Always remember:

This is a

**moving-target flyby routing problem**

not a

**multiple-rendezvous problem**.

At an asteroid encounter, the spacecraft velocity

\[
\mathbf v_{\rm SC}
\]

may be completely different from the asteroid velocity.

This distinction is fundamental to obtaining efficient solutions.

---

# 50. Working methodology

Work iteratively:

\[
\text{implement}
\rightarrow
\text{run}
\rightarrow
\text{validate}
\rightarrow
\text{analyze}
\rightarrow
\text{improve}.
\]

Do not generate thousands of lines of untested code at once.

For each phase:

1. explain what is being implemented,
2. run the implementation,
3. report quantitative results,
4. verify the underlying physics,
5. save intermediate outputs,
6. continue to the next phase.

If numerical experiments invalidate an earlier assumption, revise the solution strategy instead of forcing the original design.

---

# 51. Immediate task

Start with Phase 0 and Phase 1.

Perform the following tasks now:

1. Read `CTOC14_problem.pdf`.
2. Read `MEA.txt`.
3. Extract and cross-check all official constants and constraints.
4. Analyze the orbital statistics of all 300 asteroids.
5. Implement the exact competition two-body ephemeris.
6. Reproduce the asteroid-174 example from the PDF as an ephemeris and Lambert regression test.
7. Build the first Earth-to-asteroid ballistic atlas.
8. Compute the minimum ballistic \(v_\infty\) found for all 300 targets.
9. Determine how many targets admit a ballistic transfer satisfying
   \[
   v_\infty\leq4\ {\rm km/s}.
   \]
10. Identify hard targets and dynamical outliers.
11. Based on these results, propose concrete Phase-2 graph-search parameters.

Do not proceed directly to full low-thrust optimal control before validating these steps.

---

# 52. Execution permission

You have permission to:

- create and modify project files,
- run Python code,
- run numerical experiments,
- install reasonable open-source Python dependencies if necessary,
- benchmark alternative algorithms,
- parallelize expensive searches,
- cache intermediate results,
- revise the solver architecture when numerical evidence suggests a better approach.

Continue autonomously through Phase 1 until the ballistic atlas and the Earth-to-174 regression test are fully validated.

The goal is not merely to produce code, but to produce a numerically trustworthy, competition-oriented CTOC14 solver.