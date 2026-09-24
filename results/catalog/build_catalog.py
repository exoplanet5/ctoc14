"""Build results/catalog/fleets.json, docs/results_catalog.md, results/catalog/README.md, results/catalog/viz_selection.json
from results/catalog/inventory_raw.json (measured) + validator logs + provenance notes (from docs/ and memory)."""
import json, re, os, time, pathlib, glob
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
inv = json.load(open(ROOT / 'results/catalog/inventory_raw.json'))

# ------------------------------------------------------------------ board history (decoded from leaderboard/*/ctoc14.educoder.net.html, team row "Mickey")
BOARD = [
    dict(submitted='2026-09-07 19:35', shown=29.876511, N_board=14, tag='v4'),
    dict(submitted='2026-09-15 01:28', shown=29.310936, N_board=22, tag='milp3'),
    dict(submitted='2026-09-16 09:10', shown=18.863969, N_board=16, tag='cg4'),
    dict(submitted='2026-09-16 22:13', shown=16.127697, N_board=15, tag='go1'),
    dict(submitted='2026-09-17 01:18', shown=15.163503, N_board=14, tag='d11'),
    dict(submitted='2026-09-17 12:04', shown=13.801851, N_board=12, tag='t10d'),
    dict(submitted='2026-09-23 11:52', shown=13.594796, N_board=11, tag='s15a'),
    dict(submitted='2026-09-23 17:05', shown=13.511581, N_board=11, tag='s16E'),
]
T_START = time.mktime(time.strptime('2026-08-31 12:00', '%Y-%m-%d %H:%M'))
T_END = time.mktime(time.strptime('2026-09-28 12:00', '%Y-%m-%d %H:%M'))
def k_of(ts):
    t = time.mktime(time.strptime(ts, '%Y-%m-%d %H:%M'))
    return 1.0 - 0.1 * (T_END - t) / (28 * 86400.0)
for b in BOARD:
    b['k'] = round(k_of(b['submitted']), 6)
    b['raw_from_shown'] = round(b['shown'] / b['k'], 6)
board_by_tag = {b['tag']: b for b in BOARD}

# ------------------------------------------------------------------ provenance (stage, one line, status)
PROV = {
 'b300': ('stage 1 (report.md s6)', 'first complete solution: sequential full-tank beam-search tours + low-thrust conversion (tools/convert_fleet.py); replica checker PASS', 'keep'),
 'v1':   ('stage 1 (combine_v1.log)', 'b300 + mop-up tail chains merged (11 craft); validator2 VALID', 'superseded'),
 'E1':   ('stage 1 plan_E1 (run_plan_experiments.sh)', 'planner experiment E1 (beam 300, wfuel 0.7, dvmax 2.5) converted: 9 craft at 2000 kg, 223 covered; never merged', 'superseded'),
 'E3_fleet': ('stage 1 plan_E3 (campaign_E3.log)', 'planner experiment E3 (beam 800, wrare 0.5) converted: 9 usable craft, 239 covered; never merged', 'superseded'),
 'v3':   ('stage 1 (combine_v3.log)', 'v1 with re-inserted/mop-up chains (11 craft, 294); intermediate of v4', 'superseded'),
 'v4':   ('stage 1 (report.md s9)', 'first 298-cover: v3 + mop-up craft (12 craft); first row on the board 09-07 19:35 (shown 29.877, N 14)', 'keep'),
 'swarm_prelim': ('stage j20 patient swarm (j20_research.md s6)', 'jointly planned small-tank "patient" swarm, preliminary conversion: 15 craft, 280 covered', 'superseded'),
 'swarm_v1': ('stage j20 patient swarm', 'patient swarm + tail craft: 17 craft, 290 covered', 'superseded'),
 'swarm_v2': ('stage j20 patient swarm', 'pure swarm file: 12 main + 4 tail craft (greedy fragments), 16 craft, 291 covered, J 34.98', 'superseded'),
 'mix1': ('stage j20 greedy mix', 'greedy mix of swarm + chaser fragments (order 1): 18 craft, 295 covered; never submitted', 'superseded'),
 'mix2': ('stage j20 greedy mix', 'greedy mix of swarm + chaser fragments (order 2): 18 craft, 293 covered; never submitted', 'superseded'),
 'milp1': ('stage j20 fragment MILP (tools/select_frags_milp.py)', 'set-cover MILP over 81 converted fragments of both families: 12 craft, 298, J 32.111 = TEAM on 09-10', 'superseded'),
 'milp2': ('stage phasing (phasing_campaign.md)', '09-14 re-run of the fragment MILP with the phasing-campaign fragments: same 12-craft selection and J as milp1 (bytes differ)', 'derived'),
 'milp3': ('stage phasing', 'fragment MILP allowing the patient fleets: 16 craft, 294 covered, J 30.789; submitted 09-15 01:28 (shown 29.311, N 22)', 'keep'),
 'cg1':  ('stage colgen (colgen_campaign.md)', 'column generation + waypoint LNS fleet, fragment MILP: 14 craft, 298, J 20.54', 'superseded'),
 'cg3':  ('stage colgen', 'cg1 fragments after tank tightening: 14 craft, J 19.889', 'superseded'),
 'cg4':  ('stage colgen', 'second tightening: 14 craft, J 19.717; submitted 09-16 09:10 (shown 18.864, N 16)', 'keep'),
 'cg5':  ('stage colgen', 'new LNS routes + MILP: 13 craft, 297 covered (miss 27), J 19.551', 'superseded'),
 'cgX':  ('stage colgen', '13 craft, 298 covered, J 18.599', 'superseded'),
 'cgF':  ('stage colgen', 'final colgen fleet: 13 craft, 298, J 18.4917; start point of the globalopt campaign', 'superseded'),
 'go1':  ('stage globalopt (globalopt_campaign.md)', 'whole-trajectory min-fuel SCP (ctoc14/globalopt.py) on the cgF flyby sequences: 13 craft, J 16.823; submitted 09-16 22:13 (shown 16.128, N 15)', 'keep'),
 'dd2':  ('stage globalopt', 'go1 with repeated (waypoint) flybys removed: 13 craft, J 16.720', 'superseded'),
 'd11':  ('stage globalopt / impulsive fleet search (tools/run_ialns.py)', 'relocate + dissolve route 11: 12 craft, J 15.810; submitted 09-17 01:18 (shown 15.164, N 14)', 'keep'),
 'r15':  ('stage globalopt', 'dissolve route 08 + relocate passes: 11 craft, J 15.037', 'superseded'),
 'r14':  ('stage globalopt', 'relocate passes on 11 craft: J 14.950', 'superseded'),
 't10a': ('stage globalopt', 'annealed dissolve of route 07: 10 craft, J 14.845', 'superseded'),
 't10b': ('stage globalopt', 'relocate passes: 10 craft, J 14.496', 'superseded'),
 't10c': ('stage globalopt', 'relocate passes to saturation: 10 craft, J 14.395', 'superseded'),
 't10d': ('stage globalopt', 'final export of the saturated 10-craft impulsive fleet (results/newgen/ifleet_t10d -> efleet_t10d): J 14.366055; submitted 09-17 12:04 (shown 13.801851, N 12), banked until 09-23', 'keep'),
 'TEAM': ('copy', 'byte-identical copy of t10d (same md5 265d3ac7...), written 09-17 10:40 as the upload name; results/validator2_TEAM.log is a stale 09-10 copy of the milp1 log; STALE - never upload', 'derived'),
 's15a': ('stage 15b (memory ctoc14-stage15b-closed)', 'first closed 9-craft 298 fleet: prefix-seeded trading re-plan + closer + relocate, regridded honest twin 11.8295 (results/s15b/reloc_h2/fleet), export J 13.842227; submitted 09-23 11:52 (shown 13.594796, N 11)', 'keep'),
 's16E': ('stage 16 (results/s16/PLAN.txt)', 'fleet E: iterated MILP single-target moves + pair exchange + multi-start polish (twin 11.7227, results/s16/best/fleet), export J 13.746635; submitted 09-23 17:05 (shown 13.511581, N 11) = FINAL board score, rank 10', 'keep'),
 's16a': ('stage 16', 'fleet E with exact-aware tighter conversions (tools/s16_xconv.py): J 13.739947 = best valid file; NOT submitted (k is taken from the last submission, it would show ~13.513 > 13.5116)', 'keep'),
 'Miki_s16E': ('copy (outside results/)', '~/.Trash/CTOC14_Result_Miki_s16E.txt = s16E renamed for upload (same md5 cc438c69...); not part of results/', 'derived'),
}

subs = inv['submissions']
for s in subs:
    tag = s['tag']
    st, line, status = PROV.get(tag, ('?', '', 'superseded'))
    s['stage'] = st; s['provenance'] = line; s['status'] = status
    b = board_by_tag.get(tag)
    s['submitted'] = b['submitted'] if b else None
    s['shown'] = b['shown'] if b else None
    s['N_board'] = b['N_board'] if b else None
    v = s.get('validator')
    if v:
        v['source'] = 'contest-time log' if v['log'].startswith('results/validator2_') else 'run 2026-09-24 by the CATALOG agent (no contest-time log)'
        # cross-check validator numbers vs our parse
        if 'J' in v:
            s['validator_J_diff'] = round(v['J'] - s['J_raw'], 6)
            s['validator_cov_diff'] = v['covered'] - s['covered']
    s['valid'] = bool(v and v.get('verdict') == 'VALID')
    if v and 'J' in v and abs(v['J'] - s['J_raw']) > 1e-3:
        # the log was written for an earlier file of the same name (TEAM was re-copied several times)
        s['stale_validator_log'] = dict(v, note='STALE: this log validated an earlier file of the same name '
                                        f"(J {v['J']:.6f}, {v['covered']} covered); current bytes differ")
        s['validator'] = None; s['valid'] = False
    # duplicates
same = {}
for s in subs:
    same.setdefault(s['md5'], []).append(s['tag'])
for s in subs:
    s['same_md5_as'] = [t for t in same[s['md5']] if t != s['tag']]
    if not s.get('validator') and s['same_md5_as']:
        src = next((x for x in subs if x['tag'] in s['same_md5_as'] and x.get('validator')), None)
        if src:
            s['validator'] = dict(src['validator'], source=f"inherited from {src['tag']} (identical md5)")
            s['valid'] = src['valid']

# ------------------------------------------------------------------ twin fleets: curated status
tw = inv['twin_fleets']
KEEP_TWINS = {
 'results/newgen/ifleet_t10d': 'impulsive twin of t10d (final search state, sum J_i 12.3574)',
 'results/newgen/efleet_t10d': 'exact (run_gfleet) export of ifleet_t10d = CTOC14_Result_t10d.txt (craft_*.npz, m0 per craft)',
 'results/s15b/best_closed/fleet': 'first closed 9-craft 298 twin (11.9957, NOT regridded: optimistic tanks, see honest/)',
 'results/s15b/reloc_h2/fleet': 'honest (regridded) + relocated twin 11.8327 = source of s15a',
 'results/s15b/efleet2': 'exact export of reloc_h2 = CTOC14_Result_s15a.txt',
 'results/s16/best/fleet': 'fleet E twin routes (11.7227) = source of s16E and s16a',
 'results/s16/efleet': 'exact-aware conversions of s16/best = CTOC14_Result_s16a.txt (13.7399)',
 'results/s16/submit_E': 'exact export of fleet E = CTOC14_Result_s16E.txt (13.7466)',
 'results/s18/g2_plan/fleet': 'stage 18 RHFA (planbeam proposer) first full 8-craft run: 224/298, 3189 kg, honest',
 'results/s18/p10_claims/fleet': 'RHFA + plan claims (ttl 2): 242/298, 2836 kg (lightest 8-craft RHFA fleet)',
 'results/s18/p13_claims_ttl4/fleet': 'RHFA + claims ttl 4 bonus 1.0: 246/298, 3123 kg = best 8-craft RHFA fleet',
 'results/s18/p19_claims_N9/fleet': 'RHFA 9-craft diagnostic: 262/298, 3575 kg',
 'results/s13/gen/g02/master/N9/fleet': 'stage 13 master N<=9: 290 covered at 11.70 (pool ceiling before closing)',
 'results/s15/twinfleet/F296a': 'stage 15 twin fleet 296/298 at 11.92 (TB5f + 110 closed)',
}
DERIVED_PAT = [r'/snap_\d+$', r'^results/s16/hist/', r'^results/newgen/fleet_team$', r'^results/newgen/test_dissolve$',
               r'^results/s16/fleets/Ex$', r'^results/s15b/resettle/C3b/fleet$', r'^results/s18/_', r'^results/s18/planbeam_test/']
for t in tw:
    d = t['dir']
    if d in KEEP_TWINS:
        t['status'] = 'keep'; t['provenance'] = KEEP_TWINS[d]
    elif any(re.search(p, d) for p in DERIVED_PAT):
        t['status'] = 'derived'; t['provenance'] = 'snapshot / copy / fixture of a neighbouring fleet'
    else:
        t['status'] = 'superseded'; t['provenance'] = t.get('note', '')
    t['kind'] = 'exact-fleet' if t.get('format') == 'craft' else 'twin-fleet'

# ------------------------------------------------------------------ selection for the viewer
def sub(tag):
    return next(s for s in subs if s['tag'] == tag)
def twin(d):
    return next(t for t in tw if t['dir'] == d)
valid_subs = [s for s in subs if s['valid'] and s['tag'] not in ('TEAM', 'Miki_s16E')]
earliest = min(valid_subs, key=lambda s: s['mtime'])
most = max(valid_subs, key=lambda s: (s['n_craft'], s['mtime']))
def sub_entry(key, s, title, note, priority=1):
    return dict(key=key, source=s['file'], kind='submission', n_craft=s['n_craft'], covered=s['covered'], missed=s['missed'],
                n_flybys=s['n_flybys'], J_raw=round(s['J_raw'], 6), sumJi=round(s['sumJi'], 6), date=s['mtime'][:10],
                file_mtime=s['mtime'], md5=s['md5'], submitted=s['submitted'], shown=s['shown'],
                valid=s['valid'], validator_log=s['validator']['log'] if s.get('validator') else None,
                m0=[round(c['m0'], 3) for c in s['craft']], fuel_kg=round(s['fuel_total'], 1), title=title, note=note, priority=priority)
def twin_entry(key, d, title, note, priority=1):
    t = twin(d)
    return dict(key=key, source=d, kind='twin', n_craft=t['n'], covered=t['covered'], missed=t.get('misses'),
                n_flybys=t['flybys'], J_raw=round(t['J'], 6), sumJi=round(t['sumJi'], 6), date=t['mtime'][:10], file_mtime=t['mtime'],
                md5=None, submitted=None, shown=None, valid=None, validator_log=None,
                m0=[round(x, 3) for x in t['m0']], fuel_kg=round(t['fuel'], 1) if t.get('fuel') else round(sum(t['m0']) - 600 * t['n'], 1),
                honest=t.get('honest'), route_files=sorted(p.name for p in (ROOT / d).glob('route_*.npz')), title=title, note=note, priority=priority)

sel = []
sel.append(sub_entry(earliest['tag'], earliest,
    f"{earliest['tag']}: {earliest['n_craft']} craft, {earliest['covered']} targets, raw J {earliest['J_raw']:.3f} (first complete solution, {earliest['mtime'][:10]})",
    'earliest validator2-VALID file: stage-1 sequential full-tank tours (2000 kg craft) + low-thrust conversion; 23 misses; the "chaser" era'))
sel.append(sub_entry('v4', sub('v4'),
    f"v4: 12 craft, 298 targets, raw J {sub('v4')['J_raw']:.3f} (first on the board, 09-07)",
    'first 298-cover and first leaderboard entry (shown 29.877 at k 0.926); full-tank chasers + mop-up craft', priority=2))
sel.append(sub_entry(most['tag'], most,
    f"{most['tag']}: {most['n_craft']} craft, {most['covered']} targets, raw J {most['J_raw']:.3f} (most craft, {most['mtime'][:10]})",
    'most craft among VALID files: ' + most['provenance']))
sel.append(sub_entry('t10d', sub('t10d'),
    f"t10d: 10 craft, 298 targets, raw J {sub('t10d')['J_raw']:.6f} (submitted 09-17, banked 13.8019)",
    'impulsive fleet search (relocate + dissolve) from the 13-craft colgen fleet; 10 craft at 793-978 kg; first banked score, rank 8 on 09-17'))
sel.append(sub_entry('s16E', sub('s16E'),
    f"s16E: 9 craft, 298 targets, raw J {sub('s16E')['J_raw']:.6f} (submitted 09-23 17:05, final board score 13.5116)",
    'stage 16 fleet E: closed 9-craft fleet after iterated MILP moves + swaps + polish; the file on the board (rank 10)'))
sel.append(sub_entry('s16a', sub('s16a'),
    f"s16a: 9 craft, 298 targets, raw J {sub('s16a')['J_raw']:.6f} (best valid, not submitted)",
    'same target sets as s16E with exact-aware tighter conversions; best raw J of the project (-0.0067 vs s16E); withheld because k is taken from the last submission'))
sel.append(twin_entry('p13', 'results/s18/p13_claims_ttl4/fleet',
    "p13: RHFA 8 craft, 246/298 targets, honest fuel 3123 kg (twin, 09-24)",
    'stage 18 rolling-horizon fleet auction with plan claims (ttl 4): all 8 craft launched together and advanced in 540 d slices; best 8-craft RHFA fleet; 52 misses'))
sel.append(twin_entry('p19', 'results/s18/p19_claims_N9/fleet',
    "p19: RHFA 9 craft, 262/298 targets, honest fuel 3575 kg (twin, 09-24)",
    'stage 18 RHFA 9-craft diagnostic (same claims settings as p13): coverage scales with N at ~29 flybys per craft; 36 misses'))
sel.append(twin_entry('g2', 'results/s18/g2_plan/fleet',
    "g2_plan: RHFA 8 craft, 224/298 targets, honest fuel 3189 kg (twin, 09-24)",
    'stage 18 RHFA first full run with the planbeam proposer (no claims): shows the lane-coherence loss that claims fixed (+22 targets in p13)', priority=2))

selection = dict(generated=time.strftime('%Y-%m-%d %H:%M:%S'),
                 note='priority 1 = required by docs/viz_spec.md section 4; priority 2 = optional extras (v4: 12-craft first board entry; g2: RHFA without claims). '
                      'Submission files: parse Event 0/1/3 rows (km, s) -> AU, days; twin fleets: route_*.npz in run_ialns ist format, load with tools/run_ialns.ipr.',
                 frame='heliocentric ecliptic J2000, AU, days since 2030-01-01 00:00 UTC (MJD 62502.0)',
                 fleets=sel)
json.dump(selection, open(ROOT / 'results/catalog/viz_selection.json', 'w'), indent=1, ensure_ascii=False)

# ------------------------------------------------------------------ fleets.json
out = dict(generated=time.strftime('%Y-%m-%d %H:%M:%S'), team_row='Mickey (清华大学) on ctoc14.educoder.net',
           cost_model='J_i = 1 + x + x^2, x = (m0 - 600)/1400 from the Event 0 mass; raw J = sum J_i + (300 - covered); shown = k * raw, k = 1 - 0.1 (2026-09-28 12:00 - t_submit)/28 d',
           counts=dict(submission_files_under_results=sum(1 for s in subs if s['file'].startswith('results/')),
                       submission_files_total_incl_trash_copy=len(subs), selection_only=len(inv['selection_only']),
                       twin_fleet_json=len(tw), npz_dirs_without_fleet_json=len(inv['npz_dirs_without_fleet_json'])),
           board_history=BOARD, submissions=sorted(subs, key=lambda s: s['mtime']), selection_only=inv['selection_only'],
           twin_fleets=sorted(tw, key=lambda t: t['mtime']), npz_dirs_without_fleet_json=inv['npz_dirs_without_fleet_json'],
           viz_selection=[dict(key=e['key'], source=e['source'], kind=e['kind'], n_craft=e['n_craft'], covered=e['covered'], J_raw=e['J_raw'], priority=e['priority']) for e in sel])
json.dump(out, open(ROOT / 'results/catalog/fleets.json', 'w'), indent=1, ensure_ascii=False)

# ------------------------------------------------------------------ docs/results_catalog.md
L = []
L.append('# CTOC14 results catalog (measured 2026-09-24)\n')
L.append('Every submission file `results/CTOC14_Result_*.txt` and every twin-fleet directory under `results/` was parsed by the CATALOG agent '
         '(script output: `results/catalog/inventory_raw.json`, merged catalog `results/catalog/fleets.json`, viewer selection '
         '`results/catalog/viz_selection.json`). Numbers are re-computed from the files, not copied from notes: craft = distinct SC ids, '
         'flybys = Event 3 rows, covered = distinct asteroid ids, m0 from the Event 0 row, J_i = 1 + x + x^2 with x = (m0 - 600)/1400, '
         'raw J = sum J_i + (300 - covered). Validity is the `VERDICT` line of `results/validator2_<tag>.log` (contest-time log) or, for the '
         '10 files that had none, of `results/catalog/validator2_<tag>.log` (run today with `tools/validator2.py --quiet`). '
         'The validator agreed with the parse on covered and J for every file (max |dJ| = %.1e), except that `results/validator2_TEAM.log` '
         'is stale (see section 5).\n' % max(abs(s.get('validator_J_diff', 0)) for s in subs if s.get('validator')))
L.append('Status tags: **keep** = submitted to the contest, best valid, or chosen for the viewer; **derived** = byte copy / re-export with identical content; '
         '**superseded** = intermediate step replaced by a later, better file. Nothing was deleted or moved.\n')
L.append('## 1. Board history (team row "Mickey", decoded from leaderboard/*/ctoc14.educoder.net.html; snapshot dirs are UTC)\n')
L.append('| submitted (CST) | file | shown J | board N (craft + misses) | k | raw J from shown | raw J of the file |')
L.append('|---|---|---|---|---|---|---|')
for b in BOARD:
    s = sub(b['tag'])
    L.append(f"| {b['submitted']} | {b['tag']} | {b['shown']:.6f} | {b['N_board']} | {b['k']:.4f} | {b['raw_from_shown']:.4f} | {s['J_raw']:.6f} |")
L.append('\nThe shown/k values reproduce each file\'s raw J to < 1e-3, which identifies the submitted files unambiguously. '
         'Final board state (09-24 07:24 UTC): rank 10, shown 13.511581 = s16E.\n')
L.append('## 2. Submission files (sorted by file date)\n')
L.append('| date | tag | craft | flybys | covered | sum J_i | raw J | validator2 | submitted -> shown | status | provenance |')
L.append('|---|---|---|---|---|---|---|---|---|---|---|')
for s in sorted(subs, key=lambda s: s['mtime']):
    v = s.get('validator') or {}
    vs = (v.get('verdict') or 'no log') + ('' if v.get('log', '').startswith('results/validator2_') else ' (today)' if v else '')
    subm = f"{s['submitted'][5:]} -> {s['shown']:.4f}" if s['submitted'] else '-'
    tagcell = s['tag'] if s['file'].startswith('results/') else f"{s['tag']} (~/.Trash)"
    L.append(f"| {s['mtime']} | {tagcell} | {s['n_craft']} | {s['n_flybys']} | {s['covered']} | {s['sumJi']:.4f} | {s['J_raw']:.6f} | {vs} | {subm} | {s['status']} | {s['provenance']} |")
L.append('\nSelection-only tags (fragment lists that were never combined into a file): ' +
         ', '.join(f"{x['tag']} ({x['mtime']}, `{x['file']}`)" for x in inv['selection_only']) + '.\n')
L.append('Byte-identical pairs (md5): TEAM = t10d (265d3ac7...); ~/.Trash/CTOC14_Result_Miki_s16E.txt = s16E (cc438c69...). '
         'milp2 has the same craft/J as milp1 but different bytes (re-export).\n')
L.append('md5 of every file is in `results/catalog/fleets.json` (`submissions[].md5`).\n')
L.append('## 3. Twin fleets (impulsive `route_*.npz` + `fleet.json`) and exact fleets (`craft_*.npz`)\n')
from collections import Counter
cnt = Counter(t['stage'] for t in tw)
L.append(f"{len(tw)} `fleet.json` directories were found (by stage: " + ', '.join(f"{k} {v}" for k, v in sorted(cnt.items())) +
         f") plus {len(inv['npz_dirs_without_fleet_json'])} directories holding `route_*.npz` without a `fleet.json` (route libraries, column caches, "
         'unsaved intermediates; their n/covered/J were computed from the npz files with the run_ialns tank formula). All are listed in '
         '`results/catalog/fleets.json`; the table below shows the ones that matter (all "keep" fleets and the milestones of each stage).\n')
L.append('| date | directory | kind | craft | flybys | covered | sum J_i | J (twin) | status | provenance |')
L.append('|---|---|---|---|---|---|---|---|---|---|')
MILESTONES = ['results/newgen/ifleet0', 'results/newgen/ifleet_d11', 'results/newgen/ifleet_r14', 'results/newgen/ifleet_t10a',
              'results/newgen/ifleet_t10c', 'results/newgen/ifleet_skel8', 'results/n8/fill2', 'results/n9/grow1', 'results/n8s6/PART/best',
              'results/s7/F1p/fleet', 'results/s7/cover1/fleet', 'results/s10/milp/fleet', 'results/s13/plan13/close_N6a/closed',
              'results/s13/gen/g02/master/N8/fleet', 'results/s14/A/best_routes', 'results/s15/select/r02/N9/fleet',
              'results/s15b/reloc/B297/fleet', 'results/s15b/honest/fleet', 'results/s16/fleets/C', 'results/s16/fleets/E',
              'results/s18/p3_wide/fleet', 'results/s18/p14_claims_w300/fleet', 'results/s18/p16_claims_ttl9/fleet']
shown_tw = [t for t in tw if t['status'] == 'keep' or t['dir'] in MILESTONES]
for t in sorted(shown_tw, key=lambda t: t['mtime']):
    L.append(f"| {t['mtime']} | `{t['dir']}` | {t['kind']} | {t['n']} | {t.get('flybys')} | {t['covered']} | {t['sumJi']:.4f} | {t['J']:.4f} | {t['status']} | {t['provenance']} |")
L.append('\nTwin J values are impulsive-model estimates (tank = 601.5 exp(dv/ve)); only the `honest: true` stage-18 fleets and the regridded '
         's15b/s16 fleets were priced with the per-window thrust cap, so pre-09-23 twin tanks are optimistic by up to ~0.06 J per fleet '
         '(docs: memory ctoc14-stage15b-closed, EXPORT LESSON). Exact fleets (`craft_*.npz`, `kind exact-fleet`) carry the converted m0 and match the submission files.\n')
L.append('## 4. Viewer selection (`results/catalog/viz_selection.json`)\n')
L.append('| key | source | kind | craft | covered | raw J | date | why |')
L.append('|---|---|---|---|---|---|---|---|')
for e in sel:
    L.append(f"| {e['key']} | `{e['source']}` | {e['kind']} | {e['n_craft']} | {e['covered']} | {e['J_raw']:.4f} | {e['date']} | {e['note']} |")
L.append('\nPriority-1 entries are the seven the spec asks for (earliest valid, most craft, t10d, a submitted 9-craft file, s16a, RHFA 8-craft p13, RHFA 9-craft p19); '
         'v4 and g2 are optional extras (12-craft first board entry; RHFA without claims).\n')
L.append('## 5. Cross-checks against the project notes\n')
L.append('| claim (memory / docs) | measured |')
L.append('|---|---|')
for tag, claim in [('t10d', 'raw 14.366055, 10 craft, submitted 09-17'), ('s15a', '13.842227, submitted 09-23 11:52'),
                   ('s16E', '13.746635, submitted 09-23 17:05, shown 13.5116'), ('s16a', '13.739947 best valid')]:
    s = sub(tag)
    L.append(f"| {tag}: {claim} | J {s['J_raw']:.6f}, {s['n_craft']} craft, {s['covered']} covered, {s['validator']['verdict']}, submitted {s['submitted']}, shown {s['shown']} |")
L.append(f"| TEAM.txt is a stale copy of t10d | md5 TEAM {sub('TEAM')['md5']} = t10d {sub('t10d')['md5']} |")
L.append(f"| Miki_s16E is a renamed copy of s16E | not under results/; found at ~/.Trash/CTOC14_Result_Miki_s16E.txt, md5 {sub('Miki_s16E')['md5']} = s16E {sub('s16E')['md5']} |")
L.append(f"| 34 submission files | 33 under results/ + the Trash copy = 34 parsed |")
tm = sub('TEAM')
if tm.get('stale_validator_log'):
    L.append(f"| results/validator2_TEAM.log | STALE: it validated the 09-10 TEAM file (= milp1, J {tm['stale_validator_log']['J']:.6f}); the current TEAM.txt bytes are t10d's, covered by results/validator2_t10d.log |")
open(ROOT / 'docs/results_catalog.md', 'w').write('\n'.join(L) + '\n')

# ------------------------------------------------------------------ README
R = f"""# results/catalog -- sorted inventory of every CTOC14 result (written 2026-09-24 by the CATALOG agent)

Nothing here is a contest record; the records are `results/CTOC14_Result_*.txt` and `results/validator2_*.log`, which this
directory only reads. Layout:

- `inventory_raw.json` -- raw measurements: every `results/CTOC14_Result_*.txt` (craft, flybys, covered, missed, per-craft m0 /
  J_i / launch and end day / fuel, md5, size, date, validator verdict) and every `fleet.json` / `route_*.npz` directory under
  `results/` (n, covered, flybys, tanks, sum J_i, J). Produced by `inventory.py` (here); re-run it to refresh
  (`nice -n 5 ~/.venvs/astro313/bin/python results/catalog/inventory.py`, ~1 min).
- `fleets.json` -- the catalog: the same records with stage, provenance, keep/derived/superseded status, the decoded board
  history (`board_history`), duplicate detection (`same_md5_as`), stale-log detection (`stale_validator_log`) and the viewer
  selection summary. Produced by `build_catalog.py` (here) from `inventory_raw.json`; it also writes `viz_selection.json`,
  `README.md` and `docs/results_catalog.md`.
- `viz_selection.json` -- the {len(sel)} fleets chosen for the 3D viewer (keys: {', '.join(e['key'] for e in sel)}), each with
  source path, kind (submission | twin), craft count, covered, missed, raw J, date, m0 per craft, title and note. The FLEETS
  agent converts these into `results/viz/data/fleets/<key>.json` (contract: docs/viz_spec.md section 2).
- `validator2_<tag>.log` -- `tools/validator2.py --quiet` output for the 10 submission files that had no contest-time log
  (b300, E1, E3_fleet, milp2, mix1, mix2, swarm_prelim, swarm_v1, swarm_v2, v3), run 2026-09-24.
- Human-readable table: `docs/results_catalog.md`.

Conventions: dates are file mtimes (local, CST); "submitted" times come from the leaderboard snapshots (`leaderboard/<UTC stamp>/`);
raw J = sum J_i + (300 - covered); twin J = impulsive-model estimate, exact J = submission file value.
"""
open(ROOT / 'results/catalog/README.md', 'w').write(R)
print('earliest valid:', earliest['tag'], earliest['mtime'], '| most craft:', most['tag'], most['n_craft'])
print('selection:', [(e['key'], e['n_craft'], e['covered'], e['J_raw']) for e in sel])
print('subs status:', Counter(s['status'] for s in subs), '| twins status:', Counter(t['status'] for t in tw))
