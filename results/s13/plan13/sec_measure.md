All four probes ran on the best **settled** 8-route pass that exists, FREE's N6a:
- 263 covered, sum J_i 10.2068;
- depths 48, 39, 42, 36, 34, 26, 23, 15, at 1048, 1039, 993, 930, 881, 836, 848, 699 kg;
- source `scratchpad/passes/N6a`.

It leaves 35 leftovers. Their distance to the nearest host (`N6a_leftovers.json`) is:

| nearest host within | 0.03 AU | 0.05 AU | 0.08 AU | 0.12 AU | 0.20 AU | 0.35 AU |
|---|--:|--:|--:|--:|--:|--:|
| leftovers | 4 | 7 | 13 | 19 | 30 | 34 |

### M-A. Planner segment re-plan as a closer (the "SR / re-sequencing" operator): 0 of 35

**Setup** (`sr_fill.py`, `run_sr_N6a.sh`, `sr_N6a/NN/log.txt`):
- every host is re-planned segment by segment with all of its own flybys pinned as waypoints (+-45 d);
- the 35 leftovers are available at prize 1;
- beam 300, npt 6, max_children 150, timing variants kept at each waypoint;
- the shared `skel_fill` is used through a wrapper with one fix: segment 0 now collects root states.

**Result**
- **No host inserted a single leftover.**
- Four hosts re-fly exactly (01 48/48 @ 1048 kg, 04 36/36 @ 931, 05 34/34 @ 881, 06 26/26 @ 834).
- 02 re-flies 38/39, 03 35/42, 07 21/23 @ 890 kg; 08 does not settle.

With PAR M4 (a free re-plan of an at-capacity route absorbs 0 of 15), this closes the question: a beam-born host is
at the planner's own capacity, so only the twin can add targets to it.

### M-B. A 9th "sweeper" route over the leftovers

`sweeper_probe.py`; strict pool; launch grid 0-1500 d.

| pool | dv 1.2 / 0.15 AU / npt 6 | dv 2.0 / 0.25 AU / npt 6 |
|---|---|---|
| all 35 leftovers | **5 flybys @ 634.5 kg** (settled) | planner front to 21 @ 924 kg; the three deepest (19-21) **do not settle**; SWEEP_ALL_PLACEHOLDER |
| the 18 left after closing (M-C) | n/a | SWEEP_DEAD_PLACEHOLDER |

The leftovers of a good pass are anti-chainable. They cannot seed a route of their own.

### M-C. Sequential twin closing (the quantity FREE's Gate 2 and every judge called decisive)

**Setup** (`close_probe.py`, `close_N6a/closing.json`, `close_N6a.out`). Rounds of:
1. `RI.candidates`: 3 per target, radius 0.08, widened to 0.15 and then 0.30 AU on failure;
2. `RI.w_insert` trials;
3. cheapest insertion per host, **one per host per round**, re-priced every round, accepted if dJ < 0.25;
4. a 2x-longer homotopy (stages 20, iters 10) on the 4th try.

**Result**
- **17 of 35 placed for +1.7256 J**; covered 263 -> 280, sum J_i 10.2068 -> 11.9325.
- Sorted dJ: 0.002, 0.0036, 0.013, 0.019, 0.032, 0.037, 0.049, 0.056, 0.084, 0.103, 0.119, 0.145, 0.171, 0.201,
  0.203, 0.243, 0.245.
- **The cheapest 8 average 0.026 J; the next 9 average 0.17 J.**
- **17 targets could not be placed at all**, even with the long homotopy: 3, 15, 33, 48, 53, 58, 114, 118, 132,
  137, 139, 158, 196, 219, 236, 243, 247. One more (17) priced above 0.25.
- Several unplaceable targets are only 0.06-0.09 AU from a *deep* host (33 -> 02, 53 -> 01, 118 -> 01/05, 114 ->
  03). The deep hosts (39-48 flybys) accepted almost nothing cheap.
- The shallow tail (08, 15 flybys, 699 kg) accepted far targets cheaply: 110 at 0.25 AU for 0.019 J, 216 at 0.17 AU
  for 0.049 J. It also accepted 62 at 0.119.

**Relocate/swap afterwards:** RELOCATE_PLACEHOLDER

So the closer's price curve is:
- about 8 leftovers at ~0.026 J each;
- about 9 at ~0.17 J;
- then half the leftovers unplaceable (1 J each as misses). Block moves with ejection on the t10 fleet placed such
  targets at 0.10-0.25 J, or stalled (`docs/globalopt_campaign.md` section 8.2).

FREE's "~13.5 with point insertion" assumed ~0.04 J for all 35. On N6a the measured closer lands at
**RESULT_J_PLACEHOLDER**.

### M-D. Tail portfolio on N6a's head (the capacity test)

**Setup** (`tail_probe.py`, `tail_N6a_k5/log.txt`):
- keep N6a routes 1-5 (199 covered);
- re-plan routes 6-8 over the 99 left with six portfolio members:

  | member | settings |
  |---|---|
  | T1 | dv 1.2, npt 6 |
  | T2 | dv 1.6, drmax 0.20 |
  | T3 | dv 2.0, drmax 0.25, npt 2 |
  | T4 | T1 + launch grid 0-1500 d |
  | T5 | T2 + w_t 0.5 |
  | T6 | dv 2.0, npt 6, grid 0-1500 d |

- for each route, take the deepest candidate that settles with twin tank <= min(1.15 x planner, 1050 kg).

**Result:** TAIL_PLACEHOLDER

### What the four probes mean

J is decided by **settled pass capacity**, not by closing. The closer can finish about 8 leftovers at ~0.026 J
each; every leftover beyond that costs 0.1-1 J. With sum J_i(pass) near 10.2-10.8 that gives:

| outcome | settled 8-route pass must reach | at sum J_i | leftovers left for the closer |
|---|--:|--:|--:|
| J < 13 | ~290 | <= 10.8 | <= 8, all within reach of a shallow host |
| beat t10d's banked 13.80 shown (raw < 13.90 on 09-26) | ~283 | <= 10.6 | <= 15, none unplaceable |

Today's best settled pass is 263. Planner-level passes reach 272-278, but with unsettled tails.
