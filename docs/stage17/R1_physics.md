# Stage 17 R1: PHYSICS lens (low-thrust route physics and geometry)

2026-09-23 22:00 CST, agent PHYSICS. MEASURED numbers come from scripts/outputs in `results/s17/physics/`;
everything else is an estimate or recall.

**Verdict.** Our leftovers are not a geometric family but a TEMPORAL one: long-period targets with ~4 node
passages in 15 yr ("pins"). Pins cost nothing extra inside a deep route. Tails are expensive because the first
routes take the "fillers" (multi-passage targets that make any route deep). The lever is how pins and fillers are
distributed across routes; no new orbit shape helps.

## §1 Diagnosis: what actually binds

**Definitions.** A node *event* = local minimum of a target's distance to the 1-AU ecliptic circle (torus
distance) < 0.05 AU inside the 15-yr window: 2318 events over the 298 targets (`events.py` -> `events.json`).
Classes by event count: **pins** <= 4 (95 targets, median a 2.3-2.6 AU); **medium** 5-9 (123); **fillers** >= 10
(80, median a 0.97-1.26 AU).

### D1. One geometry only; tails are not a different orbit family (MEASURED, `summ.txt`)
| s16a | deep r1-r4 | tails r5-r9 |
|:--|:--|:--|
| r at flyby, IQR (AU) | 0.957-1.040 | 0.969-1.064 |
| median |ecliptic latitude| | 1.2° | 1.5° |
| craft osculating a / e / i | 1.00 / 0.07 / 1.9° | 1.01 / 0.08 / 2.5° |
| target dnode (node offset from 1 AU) | 0.042 AU | 0.044 AU |
| v_rel at flyby | 16.2 km/s | 16.7 km/s |

- Every target is a PHA, so its Earth-MOID node is the only cheap access point.
- The pins' OTHER node lies at a median 2.34 AU; only 3% fall in [0.75, 1.30] AU (`rarity.txt`).
- **So:** inclination, eccentricity and sma families, second-node and perihelion/aphelion sweeps are dead by data.
  Any "specialist geometry" still flies the same 1-AU torus.

### D2. "Unwanted" = temporally rare, not geometrically far (MEASURED, `rarity.txt`, `lagan.txt`)
- Spearman ρ with the number of deep pool routes (>= 36 fb) carrying a target: event count +0.55 (p 8e-25),
  a -0.52, i -0.23, dnode +0.11 (p 0.06, n.s.).
- Least-wanted 99 vs most-wanted 99: events (< 0.1 AU) 4 vs 12; a 2.30 vs 1.59 AU; i 16.2° vs 9.1°.
  83% of the least-wanted sit in the s16a tails.
- The stage-14 isolated 8 have event counts 3, 3, 4, 2, 3, 3, 7, 3: 7 of 8 are pins (`history_classes.txt`).
- **Why "99 wanted by no route":** a pin exists at only ~4 points of the (t, lag) plane, one per ~4.2 yr; a free
  route passes one only by chance. Resonant phasing does not rescue pins (one flyby suffices; a resonant craft just
  repeats the same mismatch). Fillers ARE the near-1:1 resonant ones: a ~ 1 AU, a node pass every year at a slowly
  drifting lag, so any craft can take one whenever convenient.

### D3. A pin costs nothing extra inside a deep route (MEASURED, `pool_legs.txt`, `rare_deep.txt`)
59,655 legs of the 2490 honest pool routes:
- Deep routes (>= 38 fb): leg into a pin median 0.404 km/s / 124 d; into a filler 0.408 km/s / 104 d.
- 459 distinct deep routes: kg/fb = 0.010 x n_pins + 9.94 (zero slope). 47 of them carry >= 12 pins at
  9.6-11.2 kg/fb; mean 8 pins, max 14.
- Shallow routes (< 30 fb): every leg costs 0.47-0.59 km/s over 190-232 d, whatever the target class.
- **Leg cost is a property of the ROUTE, not of the target.**

### D4. Depth is bought with fillers; tails are filler-starved (MEASURED, `fillers.txt`, `history_classes.txt`)
- Over 459 deep routes: corr(fillers, depth) = 0.56, corr(fillers, kg/fb) = 0.04.
- s16a fillers per route: r1 21, r2 15, r3 11, r4 15 | r5 6, r6 6, r7 4, r8 0, r9 2. The deep four hold 62 of
  80 fillers; the tails hold 58 of 95 pins.
- Tails: legs 188 d vs 121 d (median), 0.52 vs 0.31 km/s per leg = 0.64 vs 0.43 km/s per flyby. Targets within
  0.035 AU of the track: deep 45-52, tails 27-38 (law 7 seen from the other side).
- **Disjointness is mostly fillers:** they are 57% of the targets shared by two deep pool routes (27% of the
  catalogue); a filler sits in 101 deep routes on average, a medium in 57, a pin in 39. Law 2 is largely filler hoarding.

### D5. Past failures re-read in these classes (MEASURED on old files, `history_classes.txt`)
- **Stage-5 fill2 (241):** the first two routes filled took 17 and 15 fillers, the last four 4-10 each. The 57
  leftovers were 24 pins + 27 medium (median 6 events) + 6 fillers; the prize-defined "hard" skeletons omitted 24 pins.
- **Stage-8 J9grow:** 41 of 76 leftovers are pins (base rate 32%).
- Both failures = fillers consumed early + pins never designed in.

### D6. Checks of laws 6 and 7
- **Law 6 CONFIRMED from a new angle** (MEASURED, `lagcal.json`): a phase-slope cost on the TRUE-longitude lag of
  the s16a flybys prices routes at 27-73 km/s vs twin 16-21 (2-4x over; the equation of centre, e 0.07 = +-8°).
  On the MEAN-longitude lag (osculating orbit) it gives 8-15 km/s = 49-78% of the twin; the rest is e/i-vector
  steering. Any (t, lag) proxy must carry the craft's e-vector, i.e. a carrier.
- **Law 7 CONFIRMED and refined:** the density that matters is FILLER density; filler count is the best single
  predictor of depth.
- **Element-space (Edelbaum-like) routes:** the only element space that matters is (mean lag, e-vector, i-vector)
  about the 1-AU circle. Secular drift there is the stage-10 drifting carrier: it rotates the track, adds no windows.
  **A target's window count is fixed by its period and the 15-yr window; no element drift changes it.**

### D7. What binds for 8 craft
- **Fuel does not.** At the 09-27 bar 3858 kg = m0 1082 = 23.0 km/s per craft = 0.62 km/s/fb at 37.25 fb
  (J < 13: m0 1024, 20.9 km/s, 0.56 km/s/fb). The tails' 0.64 km/s/fb is nearly affordable.
- **Pace does:** 8 craft need a leg every <= 147 d; the tails fly 188 d. Pace comes from fillers: 80 = 10 per craft.
- **The pool cannot be recombined into this shape** (T0, §3): 8 deep routes with pairwise-disjoint pin sets cover
  at most 62/95 pins.
- **Per-craft recipe:** ~12 pins + 15 medium + 10 fillers, 37 fb, <= 147 d/leg, <= 11 kg/fb. Each ingredient exists
  singly (47 deep routes carry >= 12 pins; 33 deep routes of 38 fb carry only 8-10 fillers at 10.56 kg/fb, best
  8.52). The complementary combination does not.

## §2 Proposals

### P1 (top). "Pins and grout": joint pin lanes first, then a class-quota concurrent fill
**(a) Precedent.**
- CAS NSSC 2022 (Hao, Zheng, Li; 深空探测学报 9(4)): flyby epoch = the PHA's ecliptic-crossing time; beam over
  crossing events; 18-21 PHAs per 10 yr from an ecliptic orbit (abstract, verified in `docs/methods_survey.md`).
- GTOC4 nodal-intercept teams (Beihang, PoliTo, Glasgow, CAS): node passages as the discrete event set (survey, verified).
- GTOC5 ACT/PoliTo: phase-free Edelbaum ranking first, phasing second ("orbit first, phase second").
- GTOC12 ACT & Friends: phasing-indicator pre-selection + max-weight independent set over non-overlapping ship sets.
- No GTOC team, to my knowledge, classified targets by window count [recall, unverified].

**(b) Mechanism.**
1. **Pin lanes:** 8 fixed-carrier phase lanes chosen jointly to hold ~90 of 95 pins, each at one of its ~4 windows
   (carrier-DP events restricted to pins + column generation + exact max coverage; `rare_lanes_b.py`, built).
   Realise each lane as a twin skeleton: `carrier_dp.seed_ip` -> `globalopt.insert_block` -> `run_ialns.settle` ->
   `s15b_regrid`.
2. **Concurrent fill with class quotas**, pin-richest skeleton first, segmented mode. Hard per-route caps: fillers
   <= 10, medium <= 16. Caps count members, never name them, so each route still chooses from the WHOLE remaining
   class (pool stays ~2.5x the route; no pre-assignment). Fillers go last, as grout: each has >= 10 windows, so
   one lies near any track at some epoch.
3. **Close and polish:** s15b_close, s16_iter, s15b_relocate, normal export.

**(c) Laws.**
- Law 3 (sequential starvation): measured cause = early routes hoard fillers (D4, D5); the quota removes it at the generator.
- Law 1 (forcing premium): pins are designed into the lane before the fill (stage-5 skeletons ~0.58 km/s per hard
  target; D3: pins cost the same as fillers inside deep routes). Fillers are never forced.
- Law 7: not escaped but used; the quota spreads the density resource evenly.
- **Risks (honest):**
  - Structurally this is skeleton-then-fill (241, dead). New parts: the pin/filler classification, pins by design
    (87/95 in T1b vs 71/95 flown in total by stage-5 fill2), and the quota.
  - The quota may just move starvation into the medium class (stage-5 leftovers: 27 medium).
  - T1b model dv 8.3 km/s per ~12-pin lane (~0.7 km/s per pin) = ~40% of a craft's budget before any fill.

**(d) Falsification test** (<= 60 min, 2 cores; thresholds fixed now). A/B on stage-5 fill2: same skeletons
(results/n8/skel8b), same segmented mode, portfolio and order; the ONLY change is a filler cap of 10 per route.
Needs `tools/s17_quota_fill.py` (copy of skel_fill.py) + a search.py copy with a class counter in the beam state
(1.5 h build). Baseline: 241 covered; last four routes 25/23/25/25.
- **PASS:** covered >= 262 AND last four average >= 29 fb AND fleet planner kg/fb up by <= 1.0.
- **KILL:** covered <= 248, OR the first two routes both drop below 32 fb.
- Step 2 (only on PASS): swap in T1b's 8 pin lanes as skeletons (87 pins instead of 71).

**(e) Payoff.** Balanced 8 x 37 at 10.5 kg/fb -> raw 12.86; at 12 kg/fb -> 13.37; both beat the 09-27 bar
(13.56). P(8-craft submittable fleet by 09-27 06:00) = **4%**; P(useful 9-craft rebalance >= 200 kg as a
by-product) = **6%**.

**(f) Effort.** Build ~12-16 h (quota fill 3-4, lane->skeleton 2, glue 2, closer/polish reuse 2); compute ~1.5
days on 6-8 procs. Reuses rare_lanes_b.py, carrier_dp, globalopt.insert_block, skel_fill, s15b_regrid,
s15b_close, s16_iter, s15b_relocate, the export chain.

### P2. Remaining-window (earliest-deadline) prize inside the existing generators
- **(a) Precedent:** the nodal-event teams above (CAS 2022, GTOC4 Beihang/PoliTo); JPL GTOC9 chains on nodal
  drift, where a target is reachable only while its RAAN aligns (survey, stage-9 notes). EDF is textbook scheduling;
  I know of no GTOC paper naming it [recall, unverified].
- **(b) Mechanism:** a target's prize at epoch τ depends on its REMAINING windows after τ:
  p = 1 + β·[last or second-to-last window], no premium otherwise. A route takes a pin only as it is about to
  vanish and leaves fillers (many future windows) to later routes. Needs a time-dependent prize in a copied beam
  (search.py / s14_twinbeam).
- **(c) Laws:** differs from static scarcity prizes (dead: they pull routes toward pins at every epoch); a deadline
  prize is zero except near a pin's last windows. Honest: still a prize; law 1's depth-for-hard-density trade may recur.
- **(d) Test (<= 60 min):** four sequential free routes from the full pool, search.py beam 300 (fill2 portfolio),
  baseline vs β = 1. **PASS:** pins covered +8 or more AND total flybys down <= 4 AND kg/fb up <= 0.5.
  **KILL:** pins +3 or fewer, OR total flybys down >= 8.
- **(e) Payoff:** ~2% alone for 8 craft; an add-on raising P1 and the LNS generator yield.
- **(f) Effort:** 2-3 h build (prize callback in a copied beam), 1 h compute.

### P3. Orbit-state leg surrogate (GTOC12-style) for honest event-plane pricing
- **(a) Precedent:** GTOC12 ACT & Friends NN (2x50 nodes; relative state + Lambert vectors + m0), cumulative
  error < 10-15 kg after > 15 hops vs 240 kg Lambert-only (survey, verified); Zhang et al. 2025, arXiv 2508.02911,
  0.78% dv error (survey, verified).
- **(b) Mechanism:** gradient-boosted regression (no GPU) of twin leg dv on orbit-state features: Δt, Δ mean-lag
  and slope change, Δ e-vector, Δ i-vector, node gap, v_rel, mass; trained on the 59,655 honest legs. Replaces
  the carrier model's 1.5-2x calibration fudge in the P1 lanes and the Lambert/LinLeg pricing of law 5.
- **(c) Laws:** addresses law 5 and uses law 6. Does NOT escape law 1: all training legs were CHOSEN legs (cf. the
  R² 0.01 of facility location on assigned targets).
- **(d) Test (<= 60 min):** integrate the pool (2490 routes x ~0.2 s), 80/20 split by route. **PASS:** held-out
  R² >= 0.70, median |err| <= 0.08 km/s, AND R² >= 0.5 on the s16a tail legs (out of distribution).
  **KILL:** held-out R² < 0.5.
- **(e) Payoff:** enabling only, ~1% alone; improves P1 lanes and the other agents' screening speed.
- **(f) Effort:** 3-4 h, reusing `lagcal.py` features and `run_ialns.ipr`.

**Killed by my data (do not spend cores on these):** inclination/eccentricity/sma families, second-node or
aphelion sweeps (D1); resonance-tuned specialists (D2); a pure pin-specialist craft (one lane holds <= 16 pins
(T1), and the same carrier flies only 20-29 fb when given all targets (T1 sanity): a specialist is a tail by
construction; pins must be spread ~12 per route).

## §3 Measurement actually run (~17 min compute, 1-2 cores)
Thresholds were written to `results/s17/physics/THRESHOLDS.txt` BEFORE each pre-registered run.

**T0 — can the existing pool supply complementary pin skeletons?** (`rare_milp.py` -> `rare_milp.json`). Exact
MILP (HiGHS) over the 459 distinct honest routes with >= 36 fb: choose 8 or 9 routes with pairwise-DISJOINT pin
sets, maximising pin coverage, then the union.

| variant | pins covered | union covered |
|:--|--:|--:|
| 8 routes, pin-disjoint | 62/95 | 213-220/298 |
| 9 routes, pin-disjoint | 65/95 | 225-227/298 |
| 8 routes, no disjointness | 74/95 | 249/298 |

**Verdict: KILL** (62 <= 75) of the shortcut "select pin-complementary deep routes from the pool": deep pool routes
all carry the same convenient pins, so new lanes are required.

**T1 — pin-lane capacity in the stage-7 carrier model** (`rare_lanes.py` -> `rare_lanes.log/.json`; carrier_dp
imported, not modified). Events of the 95 pins only; eps 0.03, kgap 10, kphase 1.5, lam 0.5; 40,000 launch
carriers, top 3000 DP-valued. One lane holds <= 16 pins (median 10). Greedy 8 lanes: 81 pins; after the
pre-registered re-pricing **79/95**, model dv 2.9-11.1 km/s per lane (mean 6.2).
**Verdict: INCONCLUSIVE** (between 70 and 85; mean dv > 5 km/s). Sanity: the same carriers given ALL targets fly
20-29 fb at model 10.6-18.6 km/s, consistent with stage 7.

**T1b — exploratory, NOT pre-registered** (`rare_lanes_b.py` -> `rare_lanes_b.json`). 6 rounds of column
generation up-weighting rarely-carried pins (18,000 lanes), then exact max coverage (HiGHS hit its 300 s limit, so
these are incumbents).

| | 8 lanes | 9 lanes |
|:--|:--|:--|
| pins covered | **87/95** | 92/95 |
| pins per lane | 9-13 | 10-14 |
| model dv per lane | 5.5-11.9 km/s (mean 8.3) | 3.9-12.7 km/s |
| pins missed | 36, 64, 90, 114, 132, 151, 192, 280 | 140, 169, 217 |

Reading: pin coverage by 8 designed lanes is physically available (87, vs 62 from the pool and 71 flown in total
by stage-5 fill2). Model price ~0.7 km/s per pin; stage-5 real skeletons measured 0.58.

**Supporting class re-analysis** (`history_classes.py` -> `history_classes.txt`) on s16a, stage-5 fill2 and
stage-8 J9grow gives D4/D5: filler hoarding + pins left out explain both historical 8-craft failures.

**Meaning for P1:** the pin half is feasible in the model (87/95 with 8 lanes). The filler/medium half — does a
quota fill keep late routes deep? — is untested; that is P1's 60-min A/B test, which I would run first in round 2
if cores are granted.

**File index (results/s17/physics/):** `events.json` (3041 torus events < 0.15 AU), `fleet_flybys.json`,
`fleet_lags.json` (s16a per-flyby geometry), `pool_popularity.npz`, `track_dmin.json`, `summ.txt`, `lagan.txt`,
`rarity.txt`, `pool_legs.txt`, `rare_deep.txt`, `fillers.txt`, `history_classes.txt`, `rare_milp.json`,
`rare_lanes*.{log,json}`, `lagcal.json`, `THRESHOLDS.txt`.
