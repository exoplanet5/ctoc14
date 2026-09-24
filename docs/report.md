# CTOC14 Problem A — Analysis, Solver and Optimisation Direction

## 1. What the contest asks
第十四届全国空间轨道设计竞赛 (CTOC14, 2026) 题目甲, set by 太空系统运行与控制全国重点实验室: *near-Earth-asteroid defence
multi-spacecraft traversal flyby design*. N ≥ 1 low-thrust spacecraft (600 kg dry, ≤ 2000 kg wet, Isp 4000 s, 0.5 N) launch from
Earth (v∞ ≤ 4 km/s, any direction, launch on/after 2030-01-01) and must pass within 1000 km of as many of 300 catalogued NEAs
as possible before 2044-12-31 (15-year absolute window). Sun-centred two-body dynamics, Keplerian Earth and asteroids, no
gravity assists. Cost J = Σ_i (1 + x_i + x_i²) + N_miss, x_i = fuel_i/1400 kg; a missed asteroid costs exactly one empty launch.
The final score is k·J with k rising linearly from 0.9 (contest start) to 1.0 (deadline), so early submission is rewarded.
Contest facts (docs/contest_background.md): organised by 中国力学学会, hosted by NUDT; 甲题 released 2026-08-31 12:00 Beijing time (t_start),
deadline 2026-09-28 12:00 (t_end) → k = 0.9 + (days elapsed)/280 (0.921 on 09-06, 0.95 on 09-14, 0.975 on 09-21); k is taken from the
LAST valid submission; at most 30 online validations per team per day; site https://ctoc14.educoder.net. Delaying a submission pays only
if J drops by more than ~0.36 % per day.
Full digest: `docs/problem_summary.md`; machine-readable rules: `docs/rules_spec.json`; trap list: `docs/rules_review.md`.

## 2. Are the 300 targets real? Yes — the 300 brightest potentially hazardous asteroids
JPL SBDB nearest-neighbour matching, JPL Horizons osculating elements at the contest epoch (JD 2461200.5) and the MPC NEA
catalogue all identify every row of MEA.txt. The file is a verbatim extract of MPCORB `NEA.txt` (epoch K2669 = 2026-06-09) for
**PHAs with H ≤ 18.49** (299 flagged today plus 2004 RW164), sorted by brightness: #1 (3122) Florence, #2 Cuno, #4 Phaethon,
#9 Geographos, #12 Toutatis, #131 = 2025 VP (retrograde, a = 8.6 AU), #144 = 1999 XS35 (a = 17.8 AU). 286 objects agree with
Horizons to < 6e-4; the 14 others are short-arc orbits where JPL and MPC solutions differ (the contest uses the MPC values).
Details: `docs/nea_identification.md`, table `data/nea_catalog_ids.csv`.

## 3. What the population looks like (docs/population_analysis.md)
- q ≤ 1.045 AU for all, MOID < 0.05 AU for 298: every target orbit passes close to Earth's orbit, so **energy is never the
  issue; phasing is**. A flyby needs only position coincidence, so high-inclination targets are met at their ecliptic nodes.
- 120 targets with i ≤ 10°, 104 with i > 20°; a spacecraft near 1 AU sees ~180 in-band nodal crossings per year, uniform in time
  and longitude.
- Only two targets are impossible in 2030–2044: **#131** (r = 10.7–16 AU, next perihelion 2051) and **#144** (at aphelion,
  r ≈ 34 AU). Floor: N_miss ≥ 2.
- Ballistic reach from Earth (Lambert scan, `docs/ballistic_reach.md`): 296 targets are interceptable by a zero-fuel spacecraft
  (median minimum v∞ 0.98 km/s); a zero-fuel spacecraft can generically hit at most 2 targets (4 launch DOF vs 2 net constraints
  per flyby).

## 4. Capability numbers that drive the design
| quantity | value |
|---|---|
| exhaust velocity | 39.23 km/s |
| full tank Δv | 47.2 km/s (1400 kg) |
| continuous full thrust with a full tank | 3.48 yr of the 15-yr window |
| acceleration | 0.25 → 0.83 mm/s² (2000 → 600 kg) |
| marginal cost of fuel | dJ/dx = 1 + 2x ≤ 3 per 1400 kg (last 140 kg ≈ 2.9 km/s costs 0.29) |
| leverage penalty of position-only targeting | junction Δv deliverable in the next leg if Δv ≲ ½·a·Δt; ≈ impulsive if thrust straddles the flyby |

## 5. Solver (package `ctoc14/`, design in docs/solver_design.md)
1. **Exact models** — Keplerian ephemeris, universal-variable propagation, the validator's sliding-window cubic-Lagrange thrust
   model, a replica of the official checker and a writer that integrates every row with the same model. All reproduce the PDF
   example to < 10 m and < 1e-7 kg (`tests/`).
2. **Sequence planner** (`search.py`) — depth-synchronous beam search over (target, time-of-flight) with vectorised Lambert
   legs; a junction Δv is accepted if Δv ≤ η·a·TOF (η = 0.6) with propellant m(1 − e^{−κΔv/v_e}), κ = 1 + Δv/(2aTOF); coast-distance
   prefilter makes one expansion ≈ 10 ms; score = absolute time + w·fuel. One full-tank spacecraft: **49 planned flybys**
   (beam 300, 85 s on 8 cores).
3. **Fleet allocation** (`tools/run_fleet_search.py`) — sequential: plan the best full-tank tour, remove its targets, repeat;
   tail spacecraft are planned with launch anywhere in the window and small tanks (`tools/run_tail_search.py`).
4. **Low-thrust conversion** (`lowthrust.py`) — thrust samples every day on the exact validator model; two-leg sliding window;
   Newton iteration on (thrust samples, flyby times) with 3-D miss constraints, finite-difference sensitivities on a vectorised
   RK4 integrator (agrees with the exact integrator to metres for smooth profiles), L1-reweighted least-norm steps, 0.45 N cap
   (interpolation overshoot margin), exact re-integration at every frozen flyby, a fuel guard, and **re-planning from the actual
   state** (beam search from state) whenever a leg cannot be closed. Output rows are produced by the validator's integrator.
5. **Fleet conversion** (`tools/convert_fleet.py`) — converts all tours in parallel, re-sizes each tank to the propellant actually
   used (pass 2), merges, validates, and reports J.

## 6. Results (validated with the replica checker, results/CTOC14_Result_b300.txt, 42 838 rows)
**First complete solution: 10 spacecraft, 277 of 300 targets covered, J = 49.197** (Σ J_i = 26.197, N_miss = 23).
Max dynamics errors 7e-8 km / 1e-6 m/s / 2e-8 kg; max interpolated thrust 0.484 N.

| spacecraft | m0 [kg] | flybys (validated) | propellant [kg] | J_i | J_i per target |
|---|---|---|---|---|---|
| SC1 | 2000 | 46 | 1380 | 3.000 | 0.065 |
| SC2 | 2000 | 48 | 1358 | 3.000 | 0.062 |
| SC3 | 1888 | 38 | 1181 | 2.766 | 0.073 |
| SC4 | 2000 | 39 | 1382 | 3.000 | 0.077 |
| SC5 | 1966 | 36 | 1294 | 2.928 | 0.081 |
| SC6 | 1973 | 29 | 1305 | 2.942 | 0.101 |
| SC7 | 2000 | 22 | 1349 | 3.000 | 0.136 |
| SC8 | 2000 | 15 | 1344 | 3.000 | 0.200 |
| SC9 | 985 | 3 | 366 | 1.351 | 0.450 |
| SC10 | 850 | 1 | 44 | 1.210 | 1.210 |

Missed: [15, 21, 41, 50, 76, 77, 84, 96, 116, 131, 138, 144, 160, 166, 183, 204, 212, 216, 229, 240, 268, 272, 285] (131 and 144 impossible; the others were lost in the low-thrust conversion of the depleted-catalogue tours).
Baselines: do nothing J = 300; all-ballistic doubles J ≈ 151; planning-model estimate of this fleet ≈ 30 (conversion losses cost ~19).
A mop-up pass over the 22 reachable leftovers found one full-tank chain covering 18 of them (planned) and cheap chains for the rest;
Section 9 tracks the merged result.

## 7. Should every NEA be flown by? Cost logic and answer
- A miss costs 1. A spacecraft is worth adding only if (targets gained) > (1 + x + x²), and a target is worth adding to an
  existing tour if its marginal propellant/time cost is below 1 unit — which is essentially always true inside a full-tank tour
  (0.06–0.15 per target).
- **Full-tank tours**: every target that fits must be visited. Nothing in the catalogue except #131/#144 is "too expensive".
- **Leftovers** (targets the sequential tours could not absorb or that were lost in low-thrust conversion) obey a simple table:

| leftover situation | best action | cost | vs. missing them |
|---|---|---|---|
| 18 leftovers chainable by one full tank (this fleet) | add spacecraft (x = 1) | 3.0 | saves 15 |
| 3 leftovers chainable by an 850 kg craft | add spacecraft (x = 0.18) | 1.21 | saves 1.8 |
| 2 leftovers, ballistic double | add zero-fuel spacecraft | 1.0 | saves 1 |
| 1 isolated leftover (#138 here) | **miss it** | 1.0 | a dedicated craft costs ≥ 1.0 (J-neutral at best) |
| #131 (retrograde, r > 10 AU) and #144 (r ≈ 34 AU) | **miss** | 1.0 each | unreachable in 2030–2044 |

- So the answer is: **fly by everything reachable except isolated singles; 297 of 300 is the practical target; 2 misses are
  forced by physics and 1 by economics.** Fuel is never the reason to skip a target; the only real currency is spacecraft
  count (3 units per full tank) and conversion losses (each lost planned flyby = +1).
- Corollary for strategy: the dominant lever is flybys per full-tank spacecraft (each spacecraft saved is worth 3 units, more than
  any mop-up decision), followed by minimising conversion losses, then mop-up chains, then tank right-sizing (~0.05–0.2 each).

## 8. Optimisation direction (ranked by expected gain)
1. **Raise flybys per full-tank spacecraft** (each spacecraft removed = −3): larger beam, launch-date grid, fuel/time weight
   sweeps, target-value weighting (visit rare/hard targets early so they do not end up as expensive leftovers), local re-optimisation
   of each tour with the other tours' targets excluded, and 2-opt style exchanges between tours.
2. **Reduce conversion losses**: legs dropped in conversion are lost flybys; calibrate η/κ per leg length from conversion
   statistics, use 3-leg windows for tight junctions, allow the planner to use the actual arrival velocity.
3. **Right-size tanks and reserve margins** (small but free): x from 1.0 to 0.95 saves 0.15 per spacecraft.
4. **Mop-up**: pair/triple chains among leftovers with launch anywhere in the window; accept the isolated ones as misses.
5. **Submit early and iterate**: k = 0.9 → 1.0 is worth ~10 % of J, i.e. as much as one spacecraft.

## 9. Final validated solution (results/CTOC14_Result_v4.txt = results/CTOC14_Result_TEAM.txt, 12 spacecraft)
**Covered 298 of 300, N_miss = 2 (only the unreachable #131 and #144), Σ J_i = 30.259, J = 32.259** — verified by the
independent `tools/validator2.py` (VALID) and the in-house replica checker.

| spacecraft | m0 [kg] | distinct targets | J_i | J_i per target |
|---|---|---|---|---|
| SC1 | 2000 | 48 | 3.000 | 0.062 |
| SC2 | 2000 | 46 | 3.000 | 0.065 |
| SC3 | 1959 | 42 | 2.914 | 0.069 |
| SC4 | 2000 | 39 | 3.000 | 0.077 |
| SC5 | 1966 | 36 | 2.928 | 0.081 |
| SC6 | 1973 | 29 | 2.942 | 0.101 |
| SC7 | 1874 | 24 | 2.739 | 0.114 |
| SC8 | 2000 | 15 | 3.000 | 0.200 |
| SC9 | 2000 | 17 | 3.000 | 0.176 |
| SC10 | 985 | 3 | 1.351 | 0.450 |
| SC11 | 813 | 2 | 1.175 | 0.588 |
| SC12 | 850 | 4 | 1.210 | 0.303 |

How it was obtained (2026-09-07): the 12 tours of the first fleet were re-converted with the improved converter (3-leg window fallback,
adaptive iteration, expensive-leg drop, window clamp); old and new fragments were then selected greedily by marginal gain
(distinct new targets minus J_i), which gave 11 spacecraft / 294 covered / J = 35.05; the four remaining reachable targets
(12, 98, 117, 268) were chained by one 850 kg spacecraft (200 kg used, J_i = 1.21) → 298 covered.
Score trajectory: do-nothing 300 → all-ballistic doubles ≈ 151 → first fleet 49.2 → with mop-up chains 35.2 → **32.26**.
The remaining levers are all about spacecraft count (12 × ~2.5 average): fewer, fuller tours (planning with rarity weighting
converted poorly, see docs/conversion_losses.md), or merging the small mop-up craft.

## 10. Update 2026-09-10: patient-swarm campaign and fragment MILP

A second strategy was built and tested end-to-end (details in docs/j20_research.md §6): many small-tank "patient" spacecraft
planned jointly (`ctoc14/jointsearch.py`) so that each target is caught by whichever craft is cheaply phased. Planned J was 24.7
(12 craft × 1200 kg covering 281 + 3 tail craft covering 17), but exact conversion lost it: four craft needed 1250–1420 kg tanks
(planned fuel under-estimated by 10–40 % on legs with dv > 1.2 km/s), and the 23-target tail cost 10–12 by any method.
The pure swarm file validates at J = 34.98 (291 covered, 16 craft).
A set-cover MILP over all 81 converted fragments of both families (`tools/select_frags_milp.py`, duplicates allowed) gives the
new best file **results/CTOC14_Result_TEAM.txt = results/CTOC14_Result_milp1.txt: 12 spacecraft, 298 covered, J = 32.111**
(replica checker PASS, validator2 VALID) — all chaser-family fragments, SC3 swapped for its lighter 1888 kg variant.

## 11. Update 2026-09-16: column generation, waypoint neighbourhood search, tank tightening (J 30.79 → 19.55)

**Result**: `results/CTOC14_Result_TEAM.txt` = `results/CTOC14_Result_cg5.txt`, **J = 19.5511**, 13 spacecraft, 297 of 300
asteroids (misses: 27, 131, 144), both validators PASS/VALID (`results/validator2_cg5.log`). Previous best 30.789.

### 11.1 What changed conceptually
1. **Waypoints.** Every leg of the planner must END at an asteroid (Lambert ≤ 400 d, linearised ≤ 600 d). All earlier
   campaigns *excluded* the targets other craft already covered, which removes the stepping stones between the sparse
   remaining targets: a re-planned route then collapses to 2-5 flybys. Letting a route fly by covered asteroids again
   (prize 0.001 instead of exclusion; a repeat flyby costs nothing in J) turned 5- and 6-flyby craft into 17 and 13
   flybys and took a fleet from planned J 25.4 to 20.6 in five minutes. This single change is the difference between
   ~30 and ~19.
2. **Fleet-level neighbourhood search with the planner itself** (`tools/run_lns.py`): re-plan a route's tail from the
   exact state after flyby k (or the whole route), spawn a craft over the uncovered targets, re-plan two routes jointly,
   dissolve a craft; every candidate is priced with the fleet cost J = Σ cost_sc(tank_i) + misses, accepted only if J drops.
3. **Column generation** (`ctoc14/colgen.py`, `tools/run_colgen.py`) gives the route pool, the duals that price targets
   in J units, and the dives that build the initial fleets; the LP over 300k routes bounds the planned J at ~16 for
   11 craft. The integrality gap (LP 16 vs integer 19.5) is what remains.
4. **Tank tightening** (`tools/tighten_frags.py`): the standard right-sizing (600 + 1.03 × pass-1 burn + 10) leaves
   26-52 kg unused because the lighter craft burns less; re-converting at 600 + actual burn + 8 kg, twice, removes
   0.6 J over the fleet with the identical flyby sets.
5. **Fragment MILP across fleets**: four independently built fleets (13-14 craft, planned 19.3-19.5) were converted into
   a shared fragment store keyed by (launch epoch, flyby sequence); the set-cover MILP mixes their routes.

### 11.2 Pipeline (reproducible)
```bash
P=~/.venvs/astro313/bin/python
$P tools/run_colgen.py results/colgen/runX --pools 'results/phasing/camp*/pool.jsonl' 'results/phasing/pool/*.jsonl' \
   --homotopy 12,11,10 --rounds 4,4,8 --beams 2 --beam 300 --nproc 5 --dive --dive-beams 2     # route pool + dive fleet
$P tools/run_lns.py results/colgen/runX/dive results/colgen/lnsX --rounds 4 --splits=-1,0.33 \
   --prize-covered 0.001 --beam 250 --nproc 4 --spawn                                          # waypoint LNS
tools/finish_fleet.sh results/colgen/lnsX cgN 4        # convert new routes, tighten twice, fragment MILP over all
$P tools/combine_submission.py results/CTOC14_Result_cgN.txt $(cat results/CTOC14_Result_cgN.txt.selection)
$P tools/validator2.py results/CTOC14_Result_cgN.txt
```

### 11.3 What did not work (measured, not assumed)
- Plan-level re-timing of flyby times: −16 % planned propellant, 0 to −2 flybys after conversion (`analysis/retime_probe.md`).
- Multi-revolution Lambert legs: 0 of 115 planned junctions improve (legs are 100-400 d, multi-rev needs > 700 d).
- Relaxed leg cap (2.0 km/s) into high-priced targets: fewer and worse columns.
- Dissolve / group moves (also with 1150 kg tanks): the small craft hold late-window targets that the time-limited big
  routes cannot absorb; all attempts rejected.
- Recombination (LNS restarted from the MILP mix): no improvement.

## 12. Update 2026-09-17: whole-trajectory optimisation and impulsive fleet search (J 18.49 → 14.37)
Full record: docs/globalopt_campaign.md. The converter only closed flybys; a whole-trajectory L1 minimum-propellant SCP over
all thrust samples, flyby times and the launch v_inf direction (ctoc14/globalopt.py) cut propellant 25-55 % with the same flyby
sequences (J 16.823). A fast impulsive twin (ctoc14/impulsive.py, exact to ~0.5 %) made a fleet search affordable
(tools/run_ialns.py): relocating targets between craft and dissolving whole craft by aim-point-homotopy insertions reduced the
fleet from 13 to 11 craft. Validated: results/CTOC14_Result_t10d.txt, J = 14.366055 (10 craft, 298 covered); annealed dissolves (accept a worse trial, recover by relocation) took the fleet from 11 to 10 craft.

## 13. Update 2026-09-17 (pm): why the 10-craft fleet is a floor for marginal search
Full record: docs/globalopt_campaign.md sections 8.1-8.8. Three new operators were built and tested, a reduced cost model
was validated, and the fleet-size question was settled quantitatively.

**Built.** `globalopt.insert_block` (adaptive aim-point homotopy that BACKTRACKS instead of diverging, one or many
targets at once); `run_ialns.w_block` (insert a time-block of targets into a host after ejecting the host's flybys in
that window); `run_ialns.dissolve_lns` (ruin and recreate, priced per net placed target); `run_ialns.py retime`
(try every other pass of each asteroid from an event catalogue). Block insertion places targets the old operator could
not -- 8, 39, 168 and 270 went onto routes 03 and 06 after 23 of 23 single-target trials had failed.

**Settled.** Removing a craft is worth ~0.9 J but the remaining craft absorb the work: 9 craft pays only while the
fleet's total Delta-v stays below 180.1 km/s (+23.5 % of 145.8), 8 craft below 201.9 (+38.5 %). Redistribution costs
~5.4 km/s per relocated target -- 11x the 0.49 km/s fleet average -- so dissolving a 28-target route would add ~150 km/s
against a 34 km/s budget. The LNS dissolve of route 09 confirmed it: 9 of 28 targets placed for 1.241 J, then stalled.

**Reduced cost model** (`ctoc14/phasemodel.py`), accurate to 2 % on all ten routes:
dv = 0.82 [0.028 km/s per deg/yr of drift-rate variation + 29.8 km/s per rad of inclination-vector variation].
Charging the eccentricity vector separately double counts (one tangential impulse moves a and e by the same relative
amount). Useful corollary: the marginal plane cost is 8.8 km/s per AU of flyby |z|, not the 59.6 a round-trip estimate
suggests, so moving flybys onto the asteroids' node crossings is a mirage -- they already sit a median 6 days from one.

**Negatives.** The phase-plane DP (`ctoc14/phasedp.py`) is not a valid generator: it grants the eccentricity slack
independently at every event and its 10 d grid lets the phase slip ~10 deg per step. The routes are converged (400
extra SCP iterations move the tanks by 0.2-0.6 kg). Asteroids 131 and 144 are provably unreachable in the window
(10.7-16.2 AU and 33.8-34.7 AU), so N_miss = 2 is a floor. The old 304 589-column pool cannot cover the targets with
9 columns at all, and at N = 10 its LP bound (17.56) is 3.2 J worse than the fleet we fly.
