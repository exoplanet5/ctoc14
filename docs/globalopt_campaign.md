# Whole-trajectory optimisation + impulsive fleet search (2026-09-16/17)

Start: TEAM = `results/CTOC14_Result_cgF.txt`, J = 18.4917 (13 craft, 298 covered). Goal: J <= 16.
Result so far: **J = 14.366055, 10 craft, 298 covered** (`results/CTOC14_Result_t10d.txt`, validator PASS + validator2 VALID).

| step | file | craft | J (validated) |
|---|---|---|---|
| start | CTOC14_Result_cgF.txt | 13 | 18.4917 |
| whole-trajectory min-fuel SCP, same flyby sequences | CTOC14_Result_go1.txt | 13 | 16.8231 |
| remove repeated (waypoint) flybys | CTOC14_Result_dd2.txt | 13 | 16.7200 |
| impulsive fleet search: relocate + dissolve route 11 | CTOC14_Result_d11.txt | 12 | 15.8097 |
| + dissolve route 08 + relocate passes | CTOC14_Result_r15.txt | 11 | **15.0374** |

## 1. Why the old pipeline wasted propellant
The converter (`lowthrust.convert_tour`) is a feasibility Newton on 2-leg windows that stops once every flyby is within 300 km.
It never minimises propellant, and a correction can only start one or two legs before the flyby. The whole-trajectory
problem

    min sum_s |T_s|   s.t.  r(t_j) = r_ast_j(t_j) for every flyby j (flyby times free), |T_s| <= Tcap, |v_inf| <= 4 km/s

over all thrust samples of the mission cuts propellant by 25-55 % with the flyby sequences unchanged (SC7 493 -> 218 kg, SC9 417
-> 212 kg, SC11 116 -> 55 kg). Tanks shrink from 720-1101 kg to 656-900 kg.

## 2. The solver (`ctoc14/globalopt.py`)
- IRLS closed form for the L1 objective (weights 1/|T_s| + proximal rho), rows = 2 per flyby: only the B-plane components
  (perpendicular to the relative velocity) are constrained; the flyby-time shift follows from the along-velocity component.
  Jacobi-scaled normal equations.
- Sensitivities must be the exact discrete ones of the RK4 model: batched central differences per sample interval on the
  7-state and on the 4 Lagrange-window samples, chained with cumulative 7x7 products (`GlobalProblem.chain`, `rows_at`).
  Chained Kepler STMs with finite-difference noise give ~2 % errors, which the optimiser exploits (steps rejected forever).
- Pitfalls fixed on the way:
  1. a trace-scaled regulariser next to the huge sensitivities of free flyby times and v_inf silently violated the other
     constraints: eliminate the times by B-plane projection;
  2. in thrust space a large step changes how much propellant is burned early and therefore every later acceleration
     a = T/m: optimise **mass-free** (samples = accelerations x m0, massless RK4), convert back with
     m0 = (600 + margin) exp(dv / ve) (`to_massless`, `to_thrust`);
  3. **tank tightening is exact by scaling every sample by m0_new/m0_old** (m(t) is proportional to m0 at a fixed acceleration
     history); changing m0 alone moved flybys by 1 AU;
  4. inserting a target 0.1-0.3 AU off the trajectory needs an aim-point homotopy (`insert_homotopy`: the target point
     moves from the current relative position to the asteroid in stages of <= 1.5e6 km with a short restoration after each);
     plain L1 or L2 steps stall;
  5. the RK4 model at 0.25 d drifts 1000-8000 km from the validator (DOP853) on impulse-derived thrust profiles: final polish
     at 0.1 d, retry at 0.05 / 0.025 d until the validator accepts.

## 3. Fast impulsive twin (`ctoc14/impulsive.py`)
Impulses every 20 d + Kepler arcs, same interface as `GlobalProblem`, so `optimise` / `restore` / `insert_homotopy` /
`remove_flybys` run unchanged. It reproduces exact tanks within 0.2-1 % (e.g. 824.8 vs 826.7 kg) and insertion costs within
~0.01-0.02 J, at ~1.5 s per full route optimisation and ~10 s per insertion trial (exact model: ~10 min / ~7 min).
`to_exact` turns an impulsive solution into a continuous-thrust mass-free problem (constant acceleration per bin); the exact
SCP then restores the flybys (6e6 km -> 26 km in 9 s) and the exact tank comes out within 0-3 kg of the impulsive tank.

## 4. Fleet search (`tools/run_ialns.py`)
- `import`: impulsive twins of an exact fleet (`tools/run_gfleet.py` craft states).
- `relocate` pass: removal saving of every target (route re-optimised without it), then insertion trials of the most
  expensive targets into other routes at their close approaches (< 0.3 AU); apply non-conflicting moves with
  saving(A) > insertion cost(B).
- `dissolve V`: insert all unique targets of V into the other routes in rounds (one placement per host per round, cheapest
  first; candidate radius 0.3 AU widening to 0.6 AU; 3 failures = miss); accept if the fleet J drops. `--anneal x`: keep a
  trial up to x worse (no misses) and run relocate passes on it.
- `export`: validated continuous-thrust fragments (`frags/frag_<route>.txt`) -> `tools/combine_submission.py`.

Observed: removal of any single target saves at most ~0.03 J on a relocated fleet, but far insertions cost 0.1-0.2 J;
dissolves pay only for routes whose targets pass within ~0.05 AU of other routes. After a dissolve the greedy placements
are poor and relocate passes recover a lot (11-craft fleet 15.554 -> 14.940 impulsive).
Blocked dissolves (11 -> 10 craft): route 06 (+1.27 J vs 1.20, asteroid 34 unplaceable), route 09 (+1.58 vs 1.21),
route 07 (+1.31 vs 1.16).

## 5. Structural facts (useful for further work)
- Every craft flies an Earth-like orbit (a 0.93-1.10 AU, e 0.02-0.16, i 0.1-3.5 deg) and every flyby happens at r ~ 1.00 AU
  (the PHAs' Earth-orbit crossings). Craft phases relative to Earth drift over +-180 deg during the mission.
- Per flyby the optimised fleet uses ~0.33 km/s (leaders' J implies ~0.14 km/s): a route is a nearly rigid phase curve, so a
  target costs ~0.028 km/s per deg/yr of slope change; only targets within ~1-2 deg (0.02-0.03 AU) of a route are cheap.
- Pool LP bound with the recalibrated cost (s = 0.70, reserve 2 kg): 13.79 at 10.26 fractional craft; integer selections from
  the old pool are poor (J >= 19 even with a 0.15 miss penalty), so the pool is not the lever.

## 6. Reproduce
```
P=~/.venvs/astro313/bin/python
$P tools/run_globalopt.py results/CTOC14_Result_cgF.txt results/newgen/go1 --iters1 120 --iters2 60 --nproc 9
$P tools/run_gfleet.py init results/newgen/fleet_go1 --from-sub results/CTOC14_Result_cgF.txt --from-ckpt results/newgen/go1
$P tools/run_gfleet.py dedupe results/newgen/fleet_go1 results/newgen/fleet_dd2 --iters 60 --nproc 9 --keep-big
$P tools/run_ialns.py import results/newgen/fleet_dd2 results/newgen/ifleet0 --nproc 8
$P tools/run_ialns.py search results/newgen/ifleet0 results/newgen/ifleet1 --nproc 8 --relocate --rounds 6
$P tools/run_ialns.py search results/newgen/ifleet1 results/newgen/ifleet2 --nproc 8 --rounds 8          # dissolves
$P tools/run_ialns.py search results/newgen/ifleet2 results/newgen/ifleet3 --nproc 8 --rounds 6 --relocate --dissolve-k 4
$P tools/run_ialns.py export results/newgen/ifleet3 results/newgen/efleet --nproc 2
$P tools/combine_submission.py results/CTOC14_Result_X.txt $(ls results/newgen/efleet/frags/frag_*.txt | sort)
$P tools/validator2.py results/CTOC14_Result_X.txt
```

## 7. Ten craft (2026-09-17 04:20-09:20)
| step | file | craft | J (validated) |
|---|---|---|---|
| relocate passes on 11 craft (impulsive 14.940) | CTOC14_Result_r14.txt | 11 | 14.9501 |
| annealed dissolve of route 07 (`--anneal 0.7`: trial 15.559 -> 2 relocate passes -> 14.788) | CTOC14_Result_t10a.txt | 10 | 14.8454 |
| relocate passes (impulsive 14.523) | CTOC14_Result_t10b.txt | 10 | 14.4959 |
| relocate passes to saturation (impulsive 14.378) | CTOC14_Result_t10c.txt | 10 | **14.3950** |

- Annealing is essential for dissolves: greedy placements are poor, relocation recovers 0.6-0.8 J afterwards; with
  +0.35 allowance every 11 -> 10 attempt failed, with +0.7 route 07 succeeded.
- 10 -> 9 attempts (routes 06, 03) fail clearly (+2.2 J placements plus unplaceable targets).
- Impulsive model optimism on heavy routes: 20-d impulses above what 0.45 N delivers on a ~1000 kg craft; fixed with a
  mass-dependent cap (`ImpulsiveProblem.tcap`, `enforce_cap` homotopy); after rebalancing all tanks are < 1000 kg and exports
  agree within 0-3 kg.
- Negative results: launch-date shifts (+-15-30 d) gain 0.009 J in total; retiming a target to another crossing of its own
  route gains nothing; regrowing a route from scratch on its own targets inserts only 17 of 28; simultaneous from-scratch
  growth of 11 routes covers 197 targets cheaply (sum J_i 12.09) then stalls on 101 hard targets.
- Swap moves (`--swap`: exchange expensive targets between two routes) and retiming within a route find no improvement on
  the saturated 10-craft fleet (one 0.006 J relocation in a full pass).
- 10 -> 9 with `--anneal 0.9`: route 09 (the most expensive, J_i 1.341) leaves 4 targets unplaceable (8, 39, 168, 251) after
  +1.80 J of placements; route 06 leaves 3 (7, 171, 208) after +1.61 J. Insertion-based moves cannot reach 9 craft: a
  different route generator (sequences built for cheap phasing from the start) is needed for the next step.
- Final export of the saturated fleet (impulsive 14.357): `results/CTOC14_Result_t10d.txt`, exact J 14.366055, validator PASS + validator2 VALID = TEAM.

## 8. Block moves, the reduced cost model, and why 10 -> 9 fails (2026-09-17 15:00-17:00)

### 8.1 Adaptive and block insertion (`globalopt.insert_block`, `run_ialns.w_block`)
Diagnosis of the failed 10 -> 9 dissolves: all 23 single-target insertion trials of route 09's hard targets
(8, 39, 168, 251) into the other nine routes fail the same way -- the fixed-stage aim-point homotopy restores one or
two stages, then every restore step is rejected, the fuel stops growing and the miss grows linearly with the remaining
stages (`results/newgen/scratch/diag_ins.py`). The mechanism is not a missing candidate but the bump cost: displacing a
route that is pinned by ~30 flybys by dr over a free window T needs dv ~ 4 dr / T, which the per-bin cap cannot deliver
over the ~170 d gaps of a dense route.

Two new operators follow:
- `insert_block(gp, items, ...)`: ADAPTIVE joint aim-point homotopy. All aim points move together with s in [0, 1];
  a step whose restoration fails is REVERTED (thrust, flyby times, v_inf) and retried with half the step, so the
  homotopy backtracks instead of diverging. It handles one or many targets.
- `w_block`: insert a time-contiguous block of targets into a host after ejecting the host's own flybys inside the
  block window (+- margin). Block insertion WITHOUT ejection essentially always fails; with ejection it succeeds at
  0.013-0.11 J per target (`results/newgen/scratch/t_block.py`). The host pays one excursion to the victim's phase line
  instead of one bump per target.
- `dissolve_lns`: ruin and recreate. Remove the victim, then place the unassigned targets in blocks, greedy on
  cost per NET placed target against a price lam that starts at (victim cost)/(targets) and is raised only when the
  round has no move under it. Ejected targets return to the pool (tabu `--tabu` rounds).

### 8.2 Result: redistribution cannot reach 9 craft
`run_ialns.py search ... --lns --victims 09 --anneal 1.0` places 9 of route 09's 28 targets for 1.241 J and then
stalls: 19 unplaceable at any price below 0.25 J each. Block moves DO place the targets the old operator could not
(8, 39, 168, 270 went onto routes 03 and 06), so the capability gap is closed; the economics are not:

| quantity | value |
|---|---|
| saved by dissolving a route | ~1.24 J (one craft) |
| needed per target (28 targets) | 0.044 J |
| achieved per target by block moves | 0.10-0.25 J |

A marginal target costs 2-5x what the dissolve can pay because the hosts' phase lines are 0.2-1.5 AU from the victim's
targets at the times those targets are reachable. 10 -> 9 needs routes DESIGNED together, not redistribution.

### 8.3 Reduced cost model of a route (`ctoc14/phasemodel.py`) -- validated to 2 %
Every flyby happens near r = 1 AU, so a craft is a slowly moving point in element space and its propellant is the path
length in the Delta-v metric:

    dv = K_TOTAL * [ K_DRIFT * TV(drift rate) + K_PLANE * TV(inclination vector) ],
    K_DRIFT = 0.028 km/s per deg/yr, K_PLANE = 29.8 km/s per rad, K_TOTAL = 0.82

drift = n(a) - n_E [deg/yr]. Measured on the 10 optimised routes: tangential part corr 0.994 (ratio 1.01), normal part
corr 0.997 (ratio 1.01), total corr 0.994 with ratio 1.013 +- 0.019 (`results/newgen/scratch/phasemodel[23].py`).
Charging the eccentricity vector separately DOUBLE COUNTS: one tangential impulse moves a and e by the same relative
amount (d(a)/a = 2 dv / v = de), which is why the fleet's "ecc part" equals its "drift part" to three digits.

Consequences: a route is a phase curve whose slope changes cost 0.028 km/s per deg/yr; the fleet's 145.8 km/s is
TV(drift) 270-490 deg/yr plus TV(i_vec) 0.17-0.37 rad per route. A flyby sits a median 4.1 deg (p90 10.7, max 16.8)
from the craft's MEAN longitude -- the eccentricity slack 2e sin(L - peri).

### 8.4 Negative: the phase-plane DP is not a valid generator (`ctoc14/phasedp.py`)
A max-plus DP over (time, phase bin, drift bin) with one event per ring pass prices routes in 2 s and is the natural
pricing problem for column generation. Its reachability test (the craft's eccentricity vector solves both the radius
and the longitude match, |(e_r, e_t)| <= 0.16) is however far too permissive: it grants the slack independently at
every event, while a craft carries ONE eccentricity vector. The DP therefore claims 53 asteroids for 0.8 km/s and
163 for 5.2 km/s, against 24-37 per route at 11-19 km/s in reality; capping the phase slack at +-1 deg still leaves it
2x optimistic, because the radius match is relaxed the same way. A usable generator must carry the eccentricity vector
in the state (a coast arc has 4 free parameters and can hit exactly 2 events, which is what the Lambert beam search in
`ctoc14/search.py` already models).

A second, independent defect is discretisation: with a 10 d step the craft's phase moves ~10 deg between grid points
(against a 1 deg phase bin), and the per-step drift is rounded to a whole number of phase bins, so the drift error
accumulates over 548 steps. Even the self-consistent regime (radius within 3 % of the craft's a, phase slack +- 2 deg,
|z| <= 0.005 AU at the asteroid's node crossing) then claims 34 asteroids for 0.8 km/s and 63 for 17.1 km/s, against
the best real route's 37 for 16.6 km/s. A finer grid (1 d, 0.25 deg) needs ~2 GB for the back-pointers, so the DP would
have to be re-cast (continuous phase offset per state, or a beam instead of a dense grid) before it is worth anything.

### 8.5 The break-even that reframes the fleet-size question
J = N + 2 + sum(x + x^2) rewards removing a craft by ~0.9 J, but the remaining craft must absorb the work. With the
fleet's total Delta-v D spread evenly over N craft (tank = 601.5 exp(D/N/ve)):

| N | tank if D stays 145.8 | J | break-even D | allowed growth |
|---|---|---|---|---|
| 10 (now) | 873 kg | 14.32 (actual 14.357) | - | - |
| 9 | 909 kg | 13.43 | 180.1 km/s | +34.3 (23.5 %) |
| 8 | 957 kg | 12.56 | 201.9 km/s | +56.1 (38.5 %) |
| 7 | 1024 kg | 11.75 | 213.2 km/s | +67.4 (46.2 %) |

Measured total Delta-v of this campaign's fleets: 13 craft 117.4, 12 craft 122.8, 11 craft 129.8, 10 craft 145.8 km/s
-- increments of +5.4, +7.0, +16.0 km/s per craft removed, i.e. accelerating. 10 -> 9 has a 34.3 km/s budget, but the
LNS dissolve of route 09 spent 1.24 J (~5.5 km/s per target) on 9 targets alone: redistribution runs ~4x over budget
because a MARGINAL target costs 1.6-9 km/s while the fleet AVERAGE is 0.49 km/s per flyby.

So the leaders' ~12.5 is not "8 craft at our efficiency after consolidation": at 0.49 km/s per flyby and 298 flybys,
8 craft gives 12.56 -- they hold the average per-flyby cost while flying ~37 targets per craft, whereas our per-flyby
cost degrades every time we consolidate (0.394 at 13 craft -> 0.489 at 10). The lever is therefore NOT another
consolidation move but route sequences that keep the average cost when a craft carries 33-37 targets.

### 8.6 Plane cost: already near the achievable optimum
The fleet spends 62.6 km/s (43 % of 145.8) on out-of-plane Delta-v and flies by at a median |z| of 0.024 AU, which
suggests moving flybys onto the asteroids' ecliptic-node crossings. It is a mirage:
- The flybys are ALREADY within a median 6 days of a node crossing (258 of 298 within 30 d); an asteroid near
  perihelion leaves the ecliptic at ~7 km/s, so 0.024 AU IS 6 days. Hitting a node needs the craft's phase to be right
  within ~1 day, and the SCP already moves flyby times freely -- these times are its optimum trade.
- The naive round-trip estimate 2 K_PLANE |z| / r = 59.6 km/s per AU overstates the real marginal rate SEVENFOLD: the
  fleet's 62.6 km/s over 298 flybys of median |z| 0.024 AU is 8.8 km/s per AU (`K_Z` in phasemodel.py). Each flyby only
  pins the inclination vector to a LINE, and one smooth i_vec path serves the whole route.
- The shortest path touching those lines in order (convex, `results/newgen/scratch/planeopt.py`) is 50.9 km/s against
  the actual 62.6, i.e. a 19 % bound -- and NOT extractable: 400 extra SCP iterations move the tanks by 0.2-0.6 kg, so
  the routes are converged. The gap is the price of sharing each impulse between the in-plane and out-of-plane jobs
  (the L1 objective minimises |T|, not the components).

### 8.6b Global retiming is exhausted too
`run_ialns.py retime` tries every OTHER ring pass of each asteroid (event catalogue, 3 candidates per target ranked by
the reduced model): 842 trials over the 298 targets, 20 min on 8 cores. It found TWO improving moves worth 0.0024 J in
total (asteroid 120 on route 06, asteroid 65 on route 01; impulsive J 14.3574 -> 14.3550, saved in
`results/newgen/rt1`). That is below the impulsive-to-exact export noise (0-3 kg per route, ~0.003 J each), so it was
not exported and TEAM stays at `CTOC14_Result_t10d.txt`, J 14.366055. Note the ranking's plane term was 7x too
pessimistic in this run (fixed as `phasemodel.K_Z`), but with gains this small a re-ranked rerun cannot matter.

### 8.7 What a redesign has to achieve (the target to aim at)
The break-even above is a total-Delta-v budget; the useful form is the per-flyby cost, since that is what a route
generator controls. Per-flyby cost of this campaign's fleets: 0.394 (13 craft), 0.412 (12), 0.435 (11), 0.489 (10) km/s
-- each consolidation degraded it by 4.6 %, 5.6 %, 12.4 %. With 298 flybys spread over N craft:

| degradation of 0.489 km/s per flyby | 9 craft: J (gain vs 14.357) | 8 craft: J (gain) |
|---|---|---|
| 0 % | 13.42 (+0.93) | 12.19 (+2.17) |
| 10 % | 13.79 (+0.57) | 12.96 (+1.40) |
| 20 % | 14.18 (+0.18) | 13.40 (+0.96) |
| 30 % | 14.60 (-0.24) | 13.88 (+0.48) |

So a 9-craft fleet pays off while its routes stay within +23.5 % of today's per-flyby cost, and an 8-craft fleet within
+38.5 %. Given that consolidation from 11 to 10 cost +12.4 %, ONE more step is plausible and worth ~0.4-0.6 J, and
8 craft is worth ~1 J -- but only from a generator that builds 9 or 8 dense routes directly. Consolidating the existing
fleet cannot do it: the LNS dissolve pays ~5.4 km/s per relocated target (11x the average), which would add ~150 km/s
against a 34 km/s budget.

### 8.8 The old column pool cannot build 9 craft (and is 3.2 J worse than the fleet we fly)
LP relaxation of the set-covering master over the existing pool of 304 589 columns with the recalibrated cost
(s = 0.70, reserve 2 kg, `results/newgen/scratch/lp_recal.py`):

| craft cap N | LP value | fractional craft used |
|---|---|---|
| 11-13 (cap not binding) | 13.787 | 10.26 |
| 10 | 17.562 | 10.00 |
| 9 | 36.980 | 9.00 (coverage collapses) |

Two conclusions. (1) With 9 columns the pool cannot cover the targets at all -- its routes are too short, because the
beam search that produced them costs a leg with a Lambert-junction estimate and never builds 30+ target routes.
(2) At N = 10 the pool's LP bound (17.56) is 3.2 J WORSE than the 10-craft fleet actually flown (14.357), which was
built by consolidation plus whole-trajectory SCP. So the pricer, not the master, is the weak part: any new
column-generation run must cost its columns with the impulsive twin (~1.5 s per route) instead of the old estimate.

The pool is not short of DENSE columns -- 12 001 columns hold >= 30 targets and 1 665 hold >= 35 (max 38), with a best
ratio of 32.7 targets per unit cost against 20.8-29.0 for the flown routes (on the optimistic s = 0.70 costing). They
OVERLAP: dense columns are built from the ~197 easy targets, so no 9 of them cover the 101 hard ones. Forcing a route
to carry its share of the hard targets is exactly what the duals are for, so the next attempt is a column-generation
run whose pricer (a) is costed by the impulsive twin and (b) is driven by duals until the easy targets are worthless.

## 9. Stage 2: eight craft by column generation with twin-costed columns (2026-09-17 evening)
Plan and rationale: `docs/strategy_stage2.md` (leaderboard N = craft + misses; the top five are 8-craft fleets at our
per-flyby cost, so the gap is the partition, not route efficiency).

### 9.1 Gate F: planner-born dense columns ARE eight-craft material
`results/newgen/scratch/recost_pool.py` rebuilds a planner tour as a Lambert chain and settles it with the impulsive
twin (`impulsive.from_tour` + `run_ialns.settle`). 830 of the pool's 1 665 columns with >= 35 targets (ordered by
planner cost per target; 722 camp5, 107 pool) before the machine ran out of memory (`results/newgen/recost35.jsonl`,
twin states included):

| quantity | p10 | p50 | p90 |
|---|---|---|---|
| Delta-v per flyby (km/s) | 0.316 | 0.336 | 0.355 |
| tank (kg) at 35-38 targets | 806 | 821 | 836 |
| twin cost / planner cost (s = 0.70, reserve 2) | 1.004 | 1.016 | 1.028 |

So (a) 36-target routes cost 0.34 km/s per flyby in the twin -- cheaper than our best flown route (01: 37 at 0.45);
(b) the planner's cost estimate at s = 0.70 is unbiased to 2 %, so the existing beam pricer can drive column generation
and the twin only needs to correct the LP support. The catch is the same as in section 8.8: the dense columns are
built from the easy targets (46 targets appear in none of the 829, 89 in at most five), so the duals must force the
hard targets into new columns.

### 9.2 Twin construction lesson (`impulsive.from_tour`, `impulsive.settle_tour`)
Junction impulses placed 1 s after each flyby make the restoration reject EVERY step: the SCP shifts flyby times, and a
flyby that moves past its own departure impulse invalidates the linearisation (the impulse then acts before the flyby
and the miss jumps by Delta-v x shift). With the impulse node 1 d after departure the same columns settle in 4-11 s
instead of 25-35 s and the 4 % failure rate disappears; some tours then start 0.5-1 AU off (an impulse of 3 km/s one
day late), so `settle_tour` retries with lags 1 d, 0.25 d, 0.05 d (829/830 succeeded; 4 needed the 0.25 d retry).

### 9.3 Driver: `run_colgen.py --twin-recost W --twin-files ...`
After every master solve the LP support (y > 0.02) is recosted by the twin (W workers), the twin tanks replace the
planner costs (and are re-applied after each pricing call, because `ColumnStore.add` would otherwise reset a column to
a cheaper planner estimate), twin states go to `outdir/twin_states/col_<j>.npz` for a later `run_ialns import`, and
the MIP selections are recosted before acceptance. Run `results/colgen/n8a`: all pools + previous columns (304 589),
homotopy N = 9 (3 rounds) then N = 8 (12 rounds), 2 beams of 300 per round, m0 900/1000, dive with 2 beams per level.

Machine note: the run competes with ~8 GB of swap in use by other applications (load average 200+ when the recost ran
8 workers next to it); run one heavy job at a time.

### 9.4 Root pricing stalls at ~10.2 fractional craft; hard-target-rich columns are CHEAP
- Hard cap N = 9: craft dual 20.8, 112 targets at 1.0 -> a column prices out only with ~22 hard targets; the beam
  saturates after one round (rc_min -1, deepest 24). Replaced by a fixed cost per craft (`RMP(craft_cost)`,
  `--craft-cost`) under a loose cap: duals moderate (median 0.11-0.18), 900-1700 negative columns per beam.
- Beam depth vs prize: plain (scale 0) 36 deep but nothing prices out; scale 0.15 / 0.3 / 0.5: rc_min -0.90 / -1.33 /
  -0.58 at depth 31-33 / 31-33 / 24-26 (`results/newgen/scratch/prize_scale_test.py`); `--prize-scale 0.3` chosen.
- Even so the LP moves 0.03 per round (13.88 -> 13.85, 10.24 -> 10.22 craft; craft cost 2 -> 4 changes nothing): the
  count is capacity-limited by column overlap, not by price.
- The twin costs of the prize-mode columns overturn section 8's marginal-cost picture: routes of 24-31 targets holding
  9-17 of the 87 hard targets (targets in <= 5 dense columns) settle at 780-800 kg (0.37-0.41 km/s per flyby), e.g.
  column 283800: 29 targets, 17 hard, 803 kg. Hard targets cost the AVERAGE when the sequence is built for them; the
  5.4 km/s marginal price was the pinned-route artefact. The remaining problem is purely the partition (overlap), so
  the run switched to the waypoint dive (fix a column, re-price the residual with covered targets at prize 0) at
  craft cost 6, prize scale 0.3, 2 beams per level.

### 9.5 Skeleton growth (H) also equilibrates at the 10-craft level
`results/newgen/scratch/skel_mip.py` picks 8 twin-costed skeletons (disjoint on the hard targets): 70/87 hard, 195
distinct targets after dedupe, tanks 692-806 kg. `tools/grow_skeletons.py` then inserts the rest by cheapest tank
increase (one accepted insertion per route per round). Single-route test (`grow_skel.py`, column 283800: 29 targets,
17 hard, 803 kg): +8 easy targets at 0.66 km/s each -> 37 targets, 919 kg, J_i 1.280, 0.449 km/s per flyby -- the
flown route 01's efficiency WITH 17 hard targets. Fleet growth (`results/newgen/ifleet_skel8`): the first ~5 insertions
per route cost 0-55 kg, then 60-285 kg each (pinning), and 8 routes ended at 1184-1273 kg with 265 covered and 33
unplaced (J 48.2). Stopping the growth at a cost threshold leaves ~50 targets for 2 more craft = the 10-craft
equilibrium again (~14.3). Conclusion: the leaders' advantage is TIME efficiency (2.5 flybys per craft-year with hard
targets included); hard targets cost time (long legs), not Delta-v, and only a sequence built for 37 targets from the
start reaches it. Last test of the day: partition first (nearest skeleton), plain beam inside each subset
(`partition_beam.py`).

## 10. The planner mass ceiling (2026-09-18) -- why everything equilibrated at ten craft

### 10.1 The bug in the search settings, not in the search
Every column ever generated -- the phasing pools, the colgen pools, `recost35.jsonl`, the partition beams -- was produced
by `beam_search(..., m0=1000)`. The beam refuses a child that would take the craft below `M_DRY + m_margin`, so the
planned propellant is capped at `m0 - 640 = 360 kg`, and `CostModel(s=0.70, reserve=2).tank(360, 1000) = 805 kg`:
**no column in any pool can describe a craft heavier than ~805 kg**, i.e. more than 11.4 km/s of real dv.
A 37-flyby route at our measured 0.47-0.49 km/s per flyby needs 17-18 km/s, i.e. a 930-970 kg craft. It was never
representable. That is the whole explanation of the "10-craft equilibrium": column generation, the waypoint dive, the
skeleton MIP and greedy growth were all searching a space in which no craft can carry more than ~34 targets.

The planner mass is a *ceiling*, not a commitment: `m0` sets both the modelled acceleration (0.5 N / m0) and the
propellant budget, and the cost model converts the planned fuel fraction into a tank. Reachable tank by m0:

| m0 (planner) | fuel cap | tank ceiling | J_i ceiling |
|---|---|---|---|
| 1000 | 360 | 805 | 1.168 |
| 1200 | 560 | 894 | 1.254 |
| 1400 | 760 | 971 | 1.335 |
| 1600 | 960 | 1038 | 1.411 |
| 2000 | 1360 | 1149 | 1.546 |

### 10.2 Measurement (`results/newgen/scratch/m0_test.py`, `results/newgen/m0_test.out`)
One beam over all 298 targets, prize 1 J per target, `w_fuel` from `CostModel.w_fuel(0.5*(m0-640), m0)`,
beam 300, launch grid 0-800 d:

| m0 | best route | planner fuel | planner tank | J_i | J per flyby | end |
|---|---|---|---|---|---|---|
| 1000 | 34 flybys | 359 (at the cap) | 804 | 1.167 | 0.0343 | 14.40 yr |
| 1600 | **46 flybys** | 947 | 1027 | 1.399 | **0.0304** | 14.64 yr |

At m0 1000 the route is fuel-capped AND time-capped; at m0 1600 it buys shorter transfers and 12 more flybys for
0.23 J. The whole final beam is at the same depth (p10 = p50 = p90), so this is the regime, not one lucky state.

### 10.3 Consequence for the fleet
With 298 targets at ~0.47 km/s per flyby, J(N) = N (1 + x + x^2) + 2 with tank = 601.5 exp(0.47*298/N/VE):
N = 8 -> 12.36, N = 7 -> 11.57, N = 6 -> 10.87. The binding constraint is now mission TIME (15 yr, ~116 d per flyby
at 46 flybys), not the fuel. Search plan: `tools/greedy_cover.py` (beam on the not-yet-covered targets, prize 1 J
each, twin recost of every accepted route), then `tools/partition_ls.py` (re-plan one route against the uncovered
set) and the usual growth / polish / export chain.

### 10.4 The second free lever: the launch v_inf cap (2026-09-18)
`Params.vinf_cap` was 2.0 km/s in every run of the campaign, but the contest allows |v_inf| <= 4 km/s and that
energy is FREE (it comes from the launcher, not the tank). One beam over all 298 targets, m0 1600, prize 1 J per
target, `results/newgen/scratch/lever_test.py`:

| vinf_cap | flybys | planner tank | twin tank | J_i | J per flyby | |v_inf| used |
|---|---|---|---|---|---|---|
| 2.0 | 46 | 1027 | 1067 | 1.445 | 0.0314 | 1.23 |
| 3.0 | 46 | 956 | 1004 | 1.372 | 0.0298 | 2.11 |
| 4.0 | 46 | 956 | 1004 | 1.372 | **0.0298** | 2.11 |

The cap was binding: with it open the route buys 2.11 km/s of free launch energy, keeps 46 flybys and saves 63 kg
of tank (-0.073 J per craft, ~-0.55 J on a 7-8 craft fleet). 3.0 and 4.0 give the same route here (the optimum sits
at 2.11 km/s), so use 4.0: it is never worse and costs nothing. Use `--vinf 4.0` everywhere
(`greedy_cover.py`, `neighbourhood_beams.py`, `partition_ls.py`, `dual_greedy.py` all take `--vinf/--dvmax/--tofmax/--lintofmax`).

### 10.5 Why a greedy cover cannot finish the job
Greedy at m0 1600 (prize 1 J per target, `results/newgen/gc16`): 46, 41, 37, 33, 29, 25, 20 targets (tanks
1115...731), 231 covered, sum cost 8.60 -- and then the **67 remaining targets chain into a 6-flyby route**. They
are the far, inclined objects (median a 2.32 vs 1.78 AU, aphelion 3.81 vs 2.75, i 16.4 vs 11.6 deg) and can only be
covered woven INTO a deep route. With a uniform 1 J prize the beam is indifferent between targets and every route is
already fuel- AND time-capped (947 of 960 kg, 14.6 of 15 yr), so it always leaves the orphans.
Fix: `tools/dual_greedy.py`, a subgradient loop that multiplies the prize of every uncovered target by ~1.8 and
re-runs the greedy (the Lagrangian dual of the set cover), scoring every pass with the true objective.
The cost model is unchanged in this regime: implied s over 7 fresh routes = 0.689 (planner/twin tank ratio 0.993).

### 10.6 Levers that do NOT matter
Same beam (m0 1600, vinf 4.0, prize 1 J per target): raising the Lambert TOF cap `Params.tofs` from 400 d to 700 d,
or the linear-model cap `lin_tofs` from 600 d to 900 d, returns the IDENTICAL route (46 flybys, twin 1004 kg) and
doubles the beam time (229 s vs 103 s). Keep `tofs` 15-400 d and `lin_tofs` 20-600 d.
