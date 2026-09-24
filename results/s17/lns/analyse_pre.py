"""Scan the prefix-job checkpoints/fronts: states keeping ALL own suffix targets, with r9 targets added; price vs the
original route tank (pre-registered metric, THRESHOLDS.txt)."""
import sys, json, pickle, pathlib
sys.path.insert(0, '/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0, '/Users/mickey/solarsystem/ctoc14/tools')
import numpy as np
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
an = json.load(open(ROOT / 'results/s17/lns/anatomy.json'))
r9 = set(an['r9']['asts'])
rep = {}
for tag in sys.argv[1:]:
    J = json.load(open(ROOT / f'results/s17/lns/job_{tag}.json'))
    host = J['route'].split('route_')[1][:-4]; own = set(J['own']); k = J['k']
    tank0 = an[host]['tank']; prefix = set(an[host]['asts'][:k])
    C = pickle.load(open(ROOT / J['out'] / 's1' / 'ckpt.pkl', 'rb'))
    states = []
    for n, f in C['front'].items():
        for kind, (v, ts) in f.items():
            states.append((kind, ts))
    for ts in C['beam']:
        states.append(('beam', ts))
    rows = []; seen = set()
    for kind, ts in states:
        a = set(int(x) for x in ts['asts']); key = (frozenset(a), round(ts['tank'], 1))
        if key in seen: continue
        seen.add(key)
        nown = len(a & own); add = sorted(a & r9)
        rows.append(dict(kind=kind, n=len(a), own=nown, own_missing=sorted(own - a), add=add, tank=round(ts['tank'], 1),
                         dkg=round(ts['tank'] - tank0, 1), t_end=round(ts['t_end'] / 86400)))
    full = [r for r in rows if r['own'] == len(own)]
    best_own_only = min([r['tank'] for r in full if not r['add']], default=None)
    withS = sorted([r for r in full if r['add']], key=lambda r: (-len(r['add']), r['tank']))
    trade = sorted([r for r in rows if r['add'] and r['own'] < len(own)], key=lambda r: (-(len(r['add']) - len(r['own_missing'])), r['tank']))
    print(f'== {tag}: host {host} orig tank {tank0:.1f}, k {k}, own suffix {len(own)}; level {C["level"]} end {C["end"]}; '
          f'settles {C["settles"]} ok {C["settle_ok"]}; wall {C.get("wall",0):.0f} s')
    print(f'   best all-own, no-r9 state: {best_own_only} kg (re-plan penalty {None if best_own_only is None else round(best_own_only - tank0, 1)} kg)')
    for r in withS[:8]:
        print(f'   ALL-OWN + r9 {r["add"]}: tank {r["tank"]} ({r["dkg"]:+.1f} kg vs orig, {r["dkg"] / len(r["add"]):.1f} kg/target) t_end {r["t_end"]}')
    for r in trade[:8]:
        print(f'   TRADE +{r["add"]} -{r["own_missing"]}: n {r["n"]} tank {r["tank"]} ({r["dkg"]:+.1f} kg) t_end {r["t_end"]}')
    deep = max(rows, key=lambda r: (r['n'], -r['tank']))
    print(f'   deepest: n {deep["n"]} (orig {len(an[host]["asts"])}) own {deep["own"]}/{len(own)} add {deep["add"]} tank {deep["tank"]}')
    rep[tag] = dict(host=host, tank0=tank0, best_own_only=best_own_only, all_own_with_r9=withS[:10], trades=trade[:10],
                    hist=C['hist'], end=C['end'])
json.dump(rep, open(ROOT / 'results/s17/lns/suffix_test.json', 'w'), indent=1, default=str)
