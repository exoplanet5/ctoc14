# CTOC14 results catalog (measured 2026-09-24)

Every submission file `results/CTOC14_Result_*.txt` and every twin-fleet directory under `results/` was parsed by the CATALOG agent (script output: `results/catalog/inventory_raw.json`, merged catalog `results/catalog/fleets.json`, viewer selection `results/catalog/viz_selection.json`). Numbers are re-computed from the files, not copied from notes: craft = distinct SC ids, flybys = Event 3 rows, covered = distinct asteroid ids, m0 from the Event 0 row, J_i = 1 + x + x^2 with x = (m0 - 600)/1400, raw J = sum J_i + (300 - covered). Validity is the `VERDICT` line of `results/validator2_<tag>.log` (contest-time log) or, for the 10 files that had none, of `results/catalog/validator2_<tag>.log` (run today with `tools/validator2.py --quiet`). The validator agreed with the parse on covered and J for every file (max |dJ| = 1.8e+01), except that `results/validator2_TEAM.log` is stale (see section 5).

Status tags: **keep** = submitted to the contest, best valid, or chosen for the viewer; **derived** = byte copy / re-export with identical content; **superseded** = intermediate step replaced by a later, better file. Nothing was deleted or moved.

## 1. Board history (team row "Mickey", decoded from leaderboard/*/ctoc14.educoder.net.html; snapshot dirs are UTC)

| submitted (CST) | file | shown J | board N (craft + misses) | k | raw J from shown | raw J of the file |
|---|---|---|---|---|---|---|
| 2026-09-07 19:35 | v4 | 29.876511 | 14 | 0.9261 | 32.2596 | 32.259473 |
| 2026-09-15 01:28 | milp3 | 29.310936 | 22 | 0.9520 | 30.7887 | 30.788599 |
| 2026-09-16 09:10 | cg4 | 18.863969 | 16 | 0.9567 | 19.7173 | 19.717267 |
| 2026-09-16 22:13 | go1 | 16.127697 | 15 | 0.9587 | 16.8231 | 16.823063 |
| 2026-09-17 01:18 | d11 | 15.163503 | 14 | 0.9591 | 15.8098 | 15.809736 |
| 2026-09-17 12:04 | t10d | 13.801851 | 12 | 0.9607 | 14.3661 | 14.366055 |
| 2026-09-23 11:52 | s15a | 13.594796 | 11 | 0.9821 | 13.8423 | 13.842227 |
| 2026-09-23 17:05 | s16E | 13.511581 | 11 | 0.9829 | 13.7467 | 13.746635 |

The shown/k values reproduce each file's raw J to < 1e-3, which identifies the submitted files unambiguously. Final board state (09-24 07:24 UTC): rank 10, shown 13.511581 = s16E.

## 2. Submission files (sorted by file date)

| date | tag | craft | flybys | covered | sum J_i | raw J | validator2 | submitted -> shown | status | provenance |
|---|---|---|---|---|---|---|---|---|---|---|
| 2026-09-06 01:13 | b300 | 10 | 277 | 277 | 26.1969 | 49.196882 | VALID (today) | - | keep | first complete solution: sequential full-tank beam-search tours + low-thrust conversion (tools/convert_fleet.py); replica checker PASS |
| 2026-09-06 01:48 | v1 | 11 | 295 | 294 | 29.1969 | 35.196882 | VALID | - | superseded | b300 + mop-up tail chains merged (11 craft); validator2 VALID |
| 2026-09-07 13:06 | E1 | 9 | 223 | 223 | 27.0000 | 104.000000 | VALID (today) | - | superseded | planner experiment E1 (beam 300, wfuel 0.7, dvmax 2.5) converted: 9 craft at 2000 kg, 223 covered; never merged |
| 2026-09-07 15:31 | E3_fleet | 9 | 239 | 239 | 25.2918 | 86.291758 | VALID (today) | - | superseded | planner experiment E3 (beam 800, wrare 0.5) converted: 9 usable craft, 239 covered; never merged |
| 2026-09-07 19:10 | v3 | 11 | 301 | 294 | 29.0490 | 35.049013 | VALID (today) | - | superseded | v1 with re-inserted/mop-up chains (11 craft, 294); intermediate of v4 |
| 2026-09-07 19:28 | v4 | 12 | 305 | 298 | 30.2595 | 32.259473 | VALID | 09-07 19:35 -> 29.8765 | keep | first 298-cover: v3 + mop-up craft (12 craft); first row on the board 09-07 19:35 (shown 29.877, N 14) |
| 2026-09-10 09:29 | swarm_prelim | 15 | 280 | 280 | 23.3156 | 43.315586 | VALID (today) | - | superseded | jointly planned small-tank "patient" swarm, preliminary conversion: 15 craft, 280 covered |
| 2026-09-10 10:38 | swarm_v1 | 17 | 297 | 290 | 25.9984 | 35.998386 | VALID (today) | - | superseded | patient swarm + tail craft: 17 craft, 290 covered |
| 2026-09-10 11:14 | swarm_v2 | 16 | 300 | 291 | 25.9824 | 34.982409 | VALID (today) | - | superseded | pure swarm file: 12 main + 4 tail craft (greedy fragments), 16 craft, 291 covered, J 34.98 |
| 2026-09-10 11:38 | mix1 | 18 | 491 | 295 | 38.8430 | 43.843002 | VALID (today) | - | superseded | greedy mix of swarm + chaser fragments (order 1): 18 craft, 295 covered; never submitted |
| 2026-09-10 11:38 | mix2 | 18 | 361 | 293 | 32.6382 | 39.638242 | VALID (today) | - | superseded | greedy mix of swarm + chaser fragments (order 2): 18 craft, 293 covered; never submitted |
| 2026-09-10 11:52 | milp1 | 12 | 301 | 298 | 30.1114 | 32.111399 | VALID | - | superseded | set-cover MILP over 81 converted fragments of both families: 12 craft, 298, J 32.111 = TEAM on 09-10 |
| 2026-09-14 22:09 | milp2 | 12 | 301 | 298 | 30.1114 | 32.111399 | VALID (today) | - | derived | 09-14 re-run of the fragment MILP with the phasing-campaign fragments: same 12-craft selection and J as milp1 (bytes differ) |
| 2026-09-14 23:25 | milp3 | 16 | 382 | 294 | 24.7886 | 30.788599 | VALID | 09-15 01:28 -> 29.3109 | keep | fragment MILP allowing the patient fleets: 16 craft, 294 covered, J 30.789; submitted 09-15 01:28 (shown 29.311, N 22) |
| 2026-09-16 03:06 | cg1 | 14 | 347 | 298 | 18.5415 | 20.541506 | VALID | - | superseded | column generation + waypoint LNS fleet, fragment MILP: 14 craft, 298, J 20.54 |
| 2026-09-16 04:30 | cg3 | 14 | 347 | 298 | 17.8892 | 19.889156 | VALID | - | superseded | cg1 fragments after tank tightening: 14 craft, J 19.889 |
| 2026-09-16 05:30 | cg4 | 14 | 347 | 298 | 17.7173 | 19.717267 | VALID | 09-16 09:10 -> 18.8640 | keep | second tightening: 14 craft, J 19.717; submitted 09-16 09:10 (shown 18.864, N 16) |
| 2026-09-16 09:27 | cg5 | 13 | 322 | 297 | 16.5511 | 19.551056 | VALID | - | superseded | new LNS routes + MILP: 13 craft, 297 covered (miss 27), J 19.551 |
| 2026-09-16 11:24 | cgX | 13 | 323 | 298 | 16.5989 | 18.598869 | VALID | - | superseded | 13 craft, 298 covered, J 18.599 |
| 2026-09-16 13:25 | cgF | 13 | 323 | 298 | 16.4917 | 18.491747 | VALID | - | superseded | final colgen fleet: 13 craft, 298, J 18.4917; start point of the globalopt campaign |
| 2026-09-16 20:29 | go1 | 13 | 323 | 298 | 14.8231 | 16.823063 | VALID | 09-16 22:13 -> 16.1277 | keep | whole-trajectory min-fuel SCP (ctoc14/globalopt.py) on the cgF flyby sequences: 13 craft, J 16.823; submitted 09-16 22:13 (shown 16.128, N 15) |
| 2026-09-16 22:34 | dd2 | 13 | 298 | 298 | 14.7200 | 16.720044 | VALID | - | superseded | go1 with repeated (waypoint) flybys removed: 13 craft, J 16.720 |
| 2026-09-17 00:49 | d11 | 12 | 298 | 298 | 13.8097 | 15.809736 | VALID | 09-17 01:18 -> 15.1635 | keep | relocate + dissolve route 11: 12 craft, J 15.810; submitted 09-17 01:18 (shown 15.164, N 14) |
| 2026-09-17 01:59 | r15 | 11 | 298 | 298 | 13.0374 | 15.037440 | VALID | - | superseded | dissolve route 08 + relocate passes: 11 craft, J 15.037 |
| 2026-09-17 03:25 | r14 | 11 | 298 | 298 | 12.9501 | 14.950134 | VALID | - | superseded | relocate passes on 11 craft: J 14.950 |
| 2026-09-17 05:13 | t10a | 10 | 298 | 298 | 12.8454 | 14.845435 | VALID | - | superseded | annealed dissolve of route 07: 10 craft, J 14.845 |
| 2026-09-17 06:26 | t10b | 10 | 298 | 298 | 12.4959 | 14.495927 | VALID | - | superseded | relocate passes: 10 craft, J 14.496 |
| 2026-09-17 07:39 | t10c | 10 | 298 | 298 | 12.3950 | 14.395002 | VALID | - | superseded | relocate passes to saturation: 10 craft, J 14.395 |
| 2026-09-17 10:14 | t10d | 10 | 298 | 298 | 12.3661 | 14.366055 | VALID | 09-17 12:04 -> 13.8019 | keep | final export of the saturated 10-craft impulsive fleet (results/newgen/ifleet_t10d -> efleet_t10d): J 14.366055; submitted 09-17 12:04 (shown 13.801851, N 12), banked until 09-23 |
| 2026-09-17 10:40 | TEAM | 10 | 298 | 298 | 12.3661 | 14.366055 | VALID | - | derived | byte-identical copy of t10d (same md5 265d3ac7...), written 09-17 10:40 as the upload name; results/validator2_TEAM.log is a stale 09-10 copy of the milp1 log; STALE - never upload |
| 2026-09-23 11:41 | s15a | 9 | 298 | 298 | 11.8422 | 13.842227 | VALID | 09-23 11:52 -> 13.5948 | keep | first closed 9-craft 298 fleet: prefix-seeded trading re-plan + closer + relocate, regridded honest twin 11.8295 (results/s15b/reloc_h2/fleet), export J 13.842227; submitted 09-23 11:52 (shown 13.594796, N 11) |
| 2026-09-23 17:02 | s16E | 9 | 298 | 298 | 11.7466 | 13.746635 | VALID | 09-23 17:05 -> 13.5116 | keep | fleet E: iterated MILP single-target moves + pair exchange + multi-start polish (twin 11.7227, results/s16/best/fleet), export J 13.746635; submitted 09-23 17:05 (shown 13.511581, N 11) = FINAL board score, rank 10 |
| 2026-09-23 17:02 | Miki_s16E (~/.Trash) | 9 | 298 | 298 | 11.7466 | 13.746635 | VALID | - | derived | ~/.Trash/CTOC14_Result_Miki_s16E.txt = s16E renamed for upload (same md5 cc438c69...); not part of results/ |
| 2026-09-23 17:32 | s16a | 9 | 298 | 298 | 11.7399 | 13.739947 | VALID | - | keep | fleet E with exact-aware tighter conversions (tools/s16_xconv.py): J 13.739947 = best valid file; NOT submitted (k is taken from the last submission, it would show ~13.513 > 13.5116) |

Selection-only tags (fragment lists that were never combined into a file): cg2 (2026-09-16 03:29, `results/CTOC14_Result_cg2.txt.selection`), cg6 (2026-09-16 12:24, `results/CTOC14_Result_cg6.txt.selection`), cg7 (2026-09-16 12:50, `results/CTOC14_Result_cg7.txt.selection`).

Byte-identical pairs (md5): TEAM = t10d (265d3ac7...); ~/.Trash/CTOC14_Result_Miki_s16E.txt = s16E (cc438c69...). milp2 has the same craft/J as milp1 but different bytes (re-export).

md5 of every file is in `results/catalog/fleets.json` (`submissions[].md5`).

## 3. Twin fleets (impulsive `route_*.npz` + `fleet.json`) and exact fleets (`craft_*.npz`)

418 `fleet.json` directories were found (by stage: n8 28, n8s6 67, n9 2, newgen 39, s10 22, s11 1, s12 2, s13 76, s14 5, s15 22, s15b 42, s16 30, s18 35, s7 45, s8 2) plus 447 directories holding `route_*.npz` without a `fleet.json` (route libraries, column caches, unsaved intermediates; their n/covered/J were computed from the npz files with the run_ialns tank formula). All are listed in `results/catalog/fleets.json`; the table below shows the ones that matter (all "keep" fleets and the milestones of each stage).

| date | directory | kind | craft | flybys | covered | sum J_i | J (twin) | status | provenance |
|---|---|---|---|---|---|---|---|---|---|
| 2026-09-16 22:30 | `results/newgen/ifleet0` | twin-fleet | 13 | 298 | 298 | 14.6624 | 16.6624 | superseded | import results/newgen/fleet_dd2 |
| 2026-09-16 23:12 | `results/newgen/ifleet_d11` | twin-fleet | 12 | 298 | 298 | 13.7869 | 15.7869 | superseded | round 0 dissolved 11 |
| 2026-09-17 02:40 | `results/newgen/ifleet_r14` | twin-fleet | 11 | 298 | 298 | 12.9404 | 14.9404 | superseded | final |
| 2026-09-17 04:28 | `results/newgen/ifleet_t10a` | twin-fleet | 10 | 298 | 298 | 12.7880 | 14.7880 | superseded | round 0 dissolved 07 |
| 2026-09-17 06:31 | `results/newgen/ifleet_t10c` | twin-fleet | 10 | 298 | 298 | 12.3779 | 14.3779 | superseded | round 1 relocate |
| 2026-09-17 09:59 | `results/newgen/ifleet_t10d` | twin-fleet | 10 | 298 | 298 | 12.3574 | 14.3574 | keep | impulsive twin of t10d (final search state, sum J_i 12.3574) |
| 2026-09-17 10:14 | `results/newgen/efleet_t10d` | exact-fleet | 10 | 298 | 298 | 12.3661 | 14.3661 | keep | exact (run_gfleet) export of ifleet_t10d = CTOC14_Result_t10d.txt (craft_*.npz, m0 per craft) |
| 2026-09-17 22:23 | `results/newgen/ifleet_skel8` | twin-fleet | 8 | 265 | 265 | 13.2000 | 48.2000 | superseded | grown |
| 2026-09-19 02:20 | `results/n8/fill2` | twin-fleet | 8 | 241 | 241 | 10.2438 | 69.2438 | superseded | filled |
| 2026-09-19 06:46 | `results/n9/grow1` | twin-fleet | 9 | 263 | 263 | 11.7084 | 48.7084 | superseded | all round 5 |
| 2026-09-20 13:37 | `results/n8s6/PART/best` | twin-fleet | 8 | 234 | 234 | 10.4250 | 76.4250 | superseded | filled |
| 2026-09-20 22:30 | `results/s7/F1p/fleet` | twin-fleet | 8 | 249 | 249 | 11.5611 | 62.5611 | superseded | leftover pass |
| 2026-09-21 00:18 | `results/s7/cover1/fleet` | twin-fleet | 10 | 298 | 298 | 12.3546 | 14.3546 | superseded | set cover over all routes |
| 2026-09-21 17:01 | `results/s10/milp/fleet` | twin-fleet | 10 | 287 | 244 | 18.1516 | 74.1516 | superseded | joint fleet MILP (stage 10) |
| 2026-09-22 03:02 | `results/s13/plan13/close_N6a/closed` | twin-fleet | 8 | 280 | 280 | 11.9325 | 31.9325 | superseded | closing round 7 |
| 2026-09-22 16:40 | `results/s14/A/best_routes` | twin-fleet | 10 | 283 | 283 | 13.2254 | 30.2254 | superseded | s14 track A: best linearised-twin re-fly per t10d route (inspection only, not an export) |
| 2026-09-22 17:24 | `results/s13/gen/g02/master/N8/fleet` | twin-fleet | 8 | 282 | 282 | 11.3185 | 29.3185 | superseded | s13_master N<=8 |
| 2026-09-22 17:28 | `results/s13/gen/g02/master/N9/fleet` | twin-fleet | 9 | 290 | 290 | 11.7045 | 21.7045 | keep | stage 13 master N<=9: 290 covered at 11.70 (pool ceiling before closing) |
| 2026-09-22 23:49 | `results/s15/select/r02/N9/fleet` | twin-fleet | 9 | 296 | 296 | 12.4375 | 16.4375 | superseded | s15_select N<=9 |
| 2026-09-23 00:42 | `results/s15/twinfleet/F296a` | twin-fleet | 9 | 296 | 296 | 11.9208 | 15.9208 | keep | stage 15 twin fleet 296/298 at 11.92 (TB5f + 110 closed) |
| 2026-09-23 02:11 | `results/s15b/reloc/B297/fleet` | twin-fleet | 9 | 297 | 297 | 11.9540 | 14.9540 | superseded | s15b_relocate of results/s15b/fleets/B297 |
| 2026-09-23 07:31 | `results/s15b/best_closed/fleet` | twin-fleet | 9 | 298 | 298 | 11.9957 | 13.9957 | keep | first closed 9-craft 298 twin (11.9957, NOT regridded: optimistic tanks, see honest/) |
| 2026-09-23 11:16 | `results/s15b/honest/fleet` | twin-fleet | 9 | 298 | 298 | 12.0508 | 14.0508 | superseded | all 9 routes regridded to 20 d bins (honest per-window thrust cap) |
| 2026-09-23 11:30 | `results/s15b/reloc_h2/fleet` | twin-fleet | 9 | 298 | 298 | 11.8327 | 13.8327 | keep | honest (regridded) + relocated twin 11.8327 = source of s15a |
| 2026-09-23 11:41 | `results/s15b/efleet2` | exact-fleet | 9 | 298 | 298 | 11.8422 | 13.8422 | keep | exact export of reloc_h2 = CTOC14_Result_s15a.txt |
| 2026-09-23 14:33 | `results/s16/fleets/C` | twin-fleet | 9 | 298 | 298 | 11.7439 | 13.7439 | superseded | select C: all s16 columns after swap run 1 + polish it6 |
| 2026-09-23 14:59 | `results/s16/fleets/E` | twin-fleet | 9 | 298 | 298 | 11.7227 | 13.7227 | superseded | select E: all s16 columns after polish D (14:59) |
| 2026-09-23 15:14 | `results/s16/submit_E` | exact-fleet | 9 | 298 | 298 | 11.7466 | 13.7466 | keep | exact export of fleet E = CTOC14_Result_s16E.txt (13.7466) |
| 2026-09-23 17:04 | `results/s16/best/fleet` | twin-fleet | 9 | 298 | 298 | 11.7227 | 13.7227 | keep | fleet E twin routes (11.7227) = source of s16E and s16a |
| 2026-09-23 17:04 | `results/s16/efleet` | exact-fleet | 9 | 298 | 298 | 11.7399 | 13.7399 | keep | exact-aware conversions of s16/best = CTOC14_Result_s16a.txt (13.7399) |
| 2026-09-24 07:40 | `results/s18/g2_plan/fleet` | twin-fleet | 8 | 224 | 224 | 10.9406 | 86.9406 | keep | stage 18 RHFA (planbeam proposer) first full 8-craft run: 224/298, 3189 kg, honest |
| 2026-09-24 09:39 | `results/s18/p3_wide/fleet` | twin-fleet | 8 | 229 | 229 | 11.2222 | 82.2222 | superseded | stage 18 RHFA run p3_wide: commit mode regrid; tanks after regrid + settle |
| 2026-09-24 11:45 | `results/s18/p10_claims/fleet` | twin-fleet | 8 | 242 | 242 | 10.5457 | 68.5457 | keep | RHFA + plan claims (ttl 2): 242/298, 2836 kg (lightest 8-craft RHFA fleet) |
| 2026-09-24 12:38 | `results/s18/p13_claims_ttl4/fleet` | twin-fleet | 8 | 246 | 246 | 10.8628 | 64.8628 | keep | RHFA + claims ttl 4 bonus 1.0: 246/298, 3123 kg = best 8-craft RHFA fleet |
| 2026-09-24 13:13 | `results/s18/p14_claims_w300/fleet` | twin-fleet | 8 | 245 | 245 | 10.9634 | 65.9634 | superseded | stage 18 RHFA run p14_claims_w300: commit mode regrid; tanks after regrid + settle |
| 2026-09-24 14:19 | `results/s18/p16_claims_ttl9/fleet` | twin-fleet | 8 | 240 | 240 | 11.1897 | 71.1897 | superseded | stage 18 RHFA run p16_claims_ttl9: commit mode regrid; tanks after regrid + settle |
| 2026-09-24 15:27 | `results/s18/p19_claims_N9/fleet` | twin-fleet | 9 | 262 | 262 | 12.2924 | 50.2924 | keep | RHFA 9-craft diagnostic: 262/298, 3575 kg |

Twin J values are impulsive-model estimates (tank = 601.5 exp(dv/ve)); only the `honest: true` stage-18 fleets and the regridded s15b/s16 fleets were priced with the per-window thrust cap, so pre-09-23 twin tanks are optimistic by up to ~0.06 J per fleet (docs: memory ctoc14-stage15b-closed, EXPORT LESSON). Exact fleets (`craft_*.npz`, `kind exact-fleet`) carry the converted m0 and match the submission files.

## 4. Viewer selection (`results/catalog/viz_selection.json`)

| key | source | kind | craft | covered | raw J | date | why |
|---|---|---|---|---|---|---|---|
| b300 | `results/CTOC14_Result_b300.txt` | submission | 10 | 277 | 49.1969 | 2026-09-06 | earliest validator2-VALID file: stage-1 sequential full-tank tours (2000 kg craft) + low-thrust conversion; 23 misses; the "chaser" era |
| v4 | `results/CTOC14_Result_v4.txt` | submission | 12 | 298 | 32.2595 | 2026-09-07 | first 298-cover and first leaderboard entry (shown 29.877 at k 0.926); full-tank chasers + mop-up craft |
| mix1 | `results/CTOC14_Result_mix1.txt` | submission | 18 | 295 | 43.8430 | 2026-09-10 | most craft among VALID files: greedy mix of swarm + chaser fragments (order 1): 18 craft, 295 covered; never submitted |
| t10d | `results/CTOC14_Result_t10d.txt` | submission | 10 | 298 | 14.3661 | 2026-09-17 | impulsive fleet search (relocate + dissolve) from the 13-craft colgen fleet; 10 craft at 793-978 kg; first banked score, rank 8 on 09-17 |
| s16E | `results/CTOC14_Result_s16E.txt` | submission | 9 | 298 | 13.7466 | 2026-09-23 | stage 16 fleet E: closed 9-craft fleet after iterated MILP moves + swaps + polish; the file on the board (rank 10) |
| s16a | `results/CTOC14_Result_s16a.txt` | submission | 9 | 298 | 13.7399 | 2026-09-23 | same target sets as s16E with exact-aware tighter conversions; best raw J of the project (-0.0067 vs s16E); withheld because k is taken from the last submission |
| p13 | `results/s18/p13_claims_ttl4/fleet` | twin | 8 | 246 | 64.8628 | 2026-09-24 | stage 18 rolling-horizon fleet auction with plan claims (ttl 4): all 8 craft launched together and advanced in 540 d slices; best 8-craft RHFA fleet; 52 misses |
| p19 | `results/s18/p19_claims_N9/fleet` | twin | 9 | 262 | 50.2924 | 2026-09-24 | stage 18 RHFA 9-craft diagnostic (same claims settings as p13): coverage scales with N at ~29 flybys per craft; 36 misses |
| g2 | `results/s18/g2_plan/fleet` | twin | 8 | 224 | 86.9406 | 2026-09-24 | stage 18 RHFA first full run with the planbeam proposer (no claims): shows the lane-coherence loss that claims fixed (+22 targets in p13) |

Priority-1 entries are the seven the spec asks for (earliest valid, most craft, t10d, a submitted 9-craft file, s16a, RHFA 8-craft p13, RHFA 9-craft p19); v4 and g2 are optional extras (12-craft first board entry; RHFA without claims).

## 5. Cross-checks against the project notes

| claim (memory / docs) | measured |
|---|---|
| t10d: raw 14.366055, 10 craft, submitted 09-17 | J 14.366055, 10 craft, 298 covered, VALID, submitted 2026-09-17 12:04, shown 13.801851 |
| s15a: 13.842227, submitted 09-23 11:52 | J 13.842227, 9 craft, 298 covered, VALID, submitted 2026-09-23 11:52, shown 13.594796 |
| s16E: 13.746635, submitted 09-23 17:05, shown 13.5116 | J 13.746635, 9 craft, 298 covered, VALID, submitted 2026-09-23 17:05, shown 13.511581 |
| s16a: 13.739947 best valid | J 13.739947, 9 craft, 298 covered, VALID, submitted None, shown None |
| TEAM.txt is a stale copy of t10d | md5 TEAM 265d3ac7f310ea87b3687de09aa28e8e = t10d 265d3ac7f310ea87b3687de09aa28e8e |
| Miki_s16E is a renamed copy of s16E | not under results/; found at ~/.Trash/CTOC14_Result_Miki_s16E.txt, md5 cc438c69d975bad84646953b70bd6f03 = s16E cc438c69d975bad84646953b70bd6f03 |
| 34 submission files | 33 under results/ + the Trash copy = 34 parsed |
| results/validator2_TEAM.log | STALE: it validated the 09-10 TEAM file (= milp1, J 32.111399); the current TEAM.txt bytes are t10d's, covered by results/validator2_t10d.log |
