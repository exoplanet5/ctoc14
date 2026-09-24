# Phasing campaign (2026-09-14) — applying CTOC14FullSolver.md to beat J = 32.1

## 1. Where the leaders are (leaderboard/20260913_0600)
Displayed J = k·J_raw with k = 0.9 + 0.1·(t_submit − t_start)/28 d (our 32.259 of 09-07 shows as 29.877, k = 0.926).
Top: HIT-SAT 12.34 (N = 10, raw ≈ 13.1), NUAA 12.44 (N = 10), THU-LAD 12.83 (N = 9, raw ≈ 13.6). Since ΣJ_i ≥ N, these fleets
miss ≤ 3 targets and average J_i ≈ 1.11–1.28: **~30 targets per craft on 150–330 kg of propellant (5–10 kg per flyby)**, 2 flybys
per craft-year for 15 years. Our TEAM file (raw 32.1) is the chaser regime (12 craft, 28 kg/flyby).

## 2. Experiments (results/phasing/)
| test | result |
|---|---|
| ballistic screening, 300 random (t_L, v∞) orbits (`tools/exp_ballistic_screen.py`, §16 of the brief) | mean 0.5 targets within 0.01 AU of the craft, 7 within 0.03 AU; ring-hugging orbits (v∞≈0) have 45 targets within 0.01 AU of the *path* vs 20–35 for |v∞| ≥ 2 km/s. Pure ballistic chains are rare → not the lever. |
| single craft 700 kg, w_fuel 5, linear legs ≤ 1000 d, beam 200 (`probe_m700`) | 18 flybys / 95 kg / 13.8 yr = 5.3 kg per flyby but only 1.3 per year |
| exact conversion of that tour (720 kg tank) | 16 flybys, 117 kg (1.4× plan); the 430-d leg with 0.35 AU coast miss cost 27 kg vs 7 planned → linearisation invalid beyond ~0.2 AU (use lin_drmax 0.15 AU) |
| joint 10×900 kg, w_fuel 5 (`K10_m900`) | 154 covered, every craft time-limited at 14–15 yr with half the fuel left |
| joint 10×1000 kg, w_fuel 2 / 1, dv_max 1.0 / 1.2 (`K10_m1000_A/B`) | 181 / 200 covered, ~20 flybys per craft, 16 kg per flyby, fuel- and time-limited |
| single craft 1000 kg, w_fuel 1–2, beam 300–500, v∞ cap 4 / 0.7 / 1.5 (`single_m1000_*`) | **31 / 31 / 34 flybys on 353–360 kg in 14–15 yr (10.5–11.6 kg per flyby)** — the per-craft target of the leaders is reachable; the v∞ cap changes little |
| two-leg (pre-shaping) linear optimisation vs sequential single legs (11 legs of the probe) | joint is only 17 % cheaper → not the missing lever |
| linear leg model vs analytic / exact SLSQP on synthetic phasing offsets | exact min-fuel (SLSQP, exact dynamics) 0.130 km/s for a 0.02 AU along-track shift over 200 d; LinLeg IRLS 0.179, best-single-segment 0.128; 1-D analytic 0.052 is not attainable because 3-D position targeting must cancel the radial side-effect. The model is 10–30 % pessimistic, not 3×. |
| cheap-encounter density from random ring states (LinLeg, 40–250 d) | 0.5 targets reachable ≤ 0.5 km/s, 0.7 ≤ 0.8 km/s per state — cheap encounters are rare; a good sequence must be *searched* (wide beam), the joint beam of 150 states shared by 10 craft (~15 per craft) is far too narrow |

Conclusion: the leaders' regime is reproduced per craft (31–34 flybys, ~360 kg, m0 ≈ 1000 kg, J_i ≈ 1.37–1.47) by the hybrid
Lambert + linearised-leg planner with w_fuel 1–2, dv_max 1.0–1.2, wide beams; the open problem is **allocation**: 10 such craft
must cover 298 targets with little overlap. The joint beam search cannot (200 covered); the route-pool + set-cover architecture
(§17–§19 of the brief) with diversified searches is the next step (`tools/run_pool_campaign.py`, `tools/select_routes_milp.py`).

## 3. Tools added
- `ctoc14/search.py`: `Params.vinf_cap` (launch |v∞| cap), `Params.collect` (route-pool collection of every beam state).
- `tools/run_search.py`: `--pool-out`, `--vinf-cap`, `--launch-min-days`; `tools/run_joint.py`: `--lin-tofmax/--lin-step/--lin-drmax/--vinf-cap`;
  `tools/run_tail_search.py`: argparse with the same options and `--pool-out`.
- `tools/select_routes_milp.py`: set cover over planned routes (tour JSONs, directories, jsonl pools), dominance pruning,
  right-sized tanks 600 + fuel·(1+margin) + reserve, optional fleet-size cap.
- `tools/run_pool_campaign.py`: greedy rounds (exclude covered, search, dump pool) + final MILP.
- `tools/exp_ballistic_screen.py`: Monte-Carlo ballistic near-miss screening.

## 4. Route-pool campaigns (2026-09-14 evening)
`tools/run_pool_campaign.py`: greedy rounds (exclude the targets covered so far, beam-search one craft, dump every beam state
to the pool), then `select_routes_milp.py` over the pool.
| campaign | settings | rounds (new targets per round) | greedy total | own-pool MILP (12 craft) |
|---|---|---|---|---|
| camp1 | 1000 kg, w_fuel 1.5, dv_max 1.2, no rarity | 29 31 31 30 28 26 25 22 18 8 3 | 251 (47 hard leftovers) | 278, planned J 39.4 |
| camp2 | 1100 kg, w_fuel 1.0, dv_max 1.5, rarity 1 on camp1's 47 leftovers, w_rare 0.3 | 23 22 25 28 31 32 30 27 26 19 4 | 267, **all 47 hard targets taken** | 287, planned J 32.1 |
| camp3 | as camp2, rarity = LP duals of the pool (`tools/price_targets.py`), w_rare 0.5 | 22–24 every round, no decay | 272 | – |
| union camp1+2+3+chaser tail, column generation (`--colgen 8`) | LP bound 21.1 | | | 284, planned J 33.9 (gap 3.8 %) |

Exact conversion (`convert_fleet.py`, tanks = 600 + 1.3·fuel_est + 20): camp1 fleet 12/12 valid, 10 lossless (SC4 dropped 7,
SC12 dropped 2), fuel 355–482 kg (plan 340–360); camp2 fleet 12/12 valid, one flyby dropped, fuel 460–580 kg (plan 430–450).
The fragment set-cover MILP over all 110 converted fragments (old chaser family + new) still returns the chaser solution
(`results/CTOC14_Result_milp2.txt` = 32.111): the patient fleets miss 12–14 hard targets and cost ~1.6 per craft after conversion.

Encounter-time resolution: refining the flyby times of a planned 32-leg tour at 0.25 d (coordinate descent, ±6 d) lowers the
Lambert-chain Δv from 24.6 to 17.5 km/s (−29 %) with shifts of 1–3 d. Implemented in `search.expand` (`Params.tof_refine`,
CLI `--tof-refine`): every plausible candidate is refined ±4 d at 1 d then ±0.75 d at 0.25 d before the feasibility filter.
Single craft 1000 kg, w_fuel 1.5: **37 flybys on 360 kg** (34 without refinement), expansion 4× slower.

Efficiency notes: `select_routes_milp.py` now stores routes as packed bitmasks (156k routes in 3 s), dedupes with numpy, and
prices columns by LP duals (column generation) instead of the O(n²) dominance pass; the greedy campaigns are CPU-bound
(5–10 processes, < 0.3 GB each); do not run two fleet conversions at once (the log-watcher kills are harness-side, the jobs survive).
