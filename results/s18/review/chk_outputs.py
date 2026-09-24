"""Reviewer checks on a finished/partial s18 run: honesty of craft npz vs state tanks, lookahead never committed,
column epochs == twin epochs, disjointness, edf/quota consistency, potential recomputed from ckpt.pkl."""
import sys, json, pathlib, pickle, glob
sys.path.insert(0, '/Users/mickey/solarsystem/ctoc14/tools'); sys.path.insert(0, '/Users/mickey/solarsystem/ctoc14')
import numpy as np
import s18_common as C
rd = C.run_dir(sys.argv[1]); S = C.load_state(rd); prm = S['params']; H = C.d2s(prm['H']); L = C.d2s(prm['L'])
W = C.load_windows(rd / 'windows.json')
print('run', rd.name, 'slice', S['slice'], 'T', C.s2d(S['T']))
# 1 honesty: recompute from craft npz
for c in S['craft']:
    if not c['launched']: continue
    st = C.load_twin(C.craft_npz(rd, c['i'])); ip = C.twin_of(st); hc = C.honest_check(ip)
    o = C.order_of(st); tg = [int(st['asts'][i]) for i in o]; ep = [float(st['tf'][i]) for i in o]
    print(f"  c{c['i']}: state tank {c['tank']:.3f} npz tank {hc['tank']:.3f} d {c['tank']-hc['tank']:+.4f} irregular {hc['irregular']} "
          f"wincap {hc['max_window_cap']} impcap {hc['max_imp_cap']} miss {hc['miss_km']} ok {hc['ok']} targets== {tg==c['targets']} epochs== {np.allclose(ep, c['epochs'])} "
          f"t_end {C.s2d(max(ep)):.1f} d <= T {C.s2d(S['T']):.0f}: {max(ep) <= S['T']}")
# 2 per slice: columns -> twin epochs; lookahead never in twin; chosen columns consistent
for sd in sorted((rd / 'slices').glob('s*')):
    cf = sd / 'columns.json'
    if not cf.exists(): continue
    sin = C.jload(sd / 'in' / 'state_in.json'); T = float(sin['T']); rem = set(sin['remaining'])
    cols = C.jload(cf); bad = 0; nlook = 0; n_chk = 0; pot_pos = 0
    for col in cols:
        st = C.load_twin(rd / col['npz']); tf = np.asarray(st['tf'], float); asts = [int(a) for a in st['asts']]
        n_chk += 1
        if tf.max() > T + H + 1e-6: nlook += 1
        o = C.order_of(st)
        new = [(int(st['asts'][i]), float(st['tf'][i])) for i in o if float(st['tf'][i]) > T]
        if [a for a, _ in new] != col['targets'] or not np.allclose([e for _, e in new], col['epochs']): bad += 1
        if not set(col['targets']) <= rem: bad += 1; print('   column with non-remaining target', col['id'])
        if len(set(asts)) != len(asts): bad += 1; print('   dup in twin', col['id'])
        if col['potential'] > 0: pot_pos += 1
    ch = C.jload(sd / 'auction.json') if (sd / 'auction.json').exists() else None
    cm = C.jload(sd / 'commit.json') if (sd / 'commit.json').exists() else None
    print(f"  {sd.name}: T {C.s2d(T):.0f} d, {n_chk} columns checked, mismatches {bad}, twins with flyby > T+H {nlook}, potential>0 {pot_pos}, "
          f"auction {ch['status'] if ch else '-'} chosen {len(ch['chosen']) if ch else '-'}, commit acc {len(cm['accepted']) if cm else '-'} rej {len(cm['rejected']) if cm else '-'}")
    if cm:
        for a in cm['accepted']:
            r = cm['results'].get(f"{a['craft']}:{a['column_id']}", {})
            print(f"     acc c{a['craft']} {a['column_id']} n_new {a['n_new']} tank {a['tank']:.2f} (in {r.get('tank_in', -1):.2f}, over_in {r.get('over_in', -1):.3f} -> {r.get('over_out', -1):.3f}) miss {a['miss']} honest {a['honest']} secs {a['secs']}")
        for rj in cm['rejected']:
            print('     rej', rj)
# 3 potential recomputed from a ckpt
for ck in sorted(rd.glob('slices/s*/*/ckpt.pkl'))[:3]:
    with open(ck, 'rb') as f: K = pickle.load(f)
    sd = ck.parent.parent; sin = C.jload(sd / 'in' / 'state_in.json'); T = float(sin['T'])
    res = C.jload(ck.parent / 'result.json')
    states = K['states']; by = {s['id']: s for s in states}; kids = {}
    for s in states:
        if s['parent'] is not None: kids.setdefault(s['parent'], []).append(s['id'])
    def look(s): return sum(1 for t in s['tfs'] if T + H < float(t) <= T + H + L)
    def pot(i):
        return max([look(by[i])] + [pot(c) for c in kids.get(i, [])])
    mism = 0
    for col in res['columns']:
        key = frozenset(col['targets'])
        # states with same new-target set
        cand = [s for s in states if s['t_end'] <= T + H + 1e-6 and frozenset(int(a) for a in s['asts'] if float(e) > T for a, e in [(a, e)] ) == key] if False else None
        ps = [pot(s['id']) for s in states if s['t_end'] <= T + H + 1e-6 and frozenset(int(a) for a, e in zip(s['asts'], s['tfs']) if float(e) > T) == key]
        if not ps or max(ps) != col['potential']: mism += 1
    print(f"  ckpt {ck.parent.relative_to(rd)}: {len(states)} states, {len(res['columns'])} columns, potential mismatches {mism}, max look of any state {max(look(s) for s in states)}")
