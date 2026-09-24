"""Compile the before/after table of results/review2/<tag>_scN_info.json. Usage: review2_table.py tag1 tag2 ..."""
import sys, json, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
D = ROOT / 'results/review2'
tags = sys.argv[1:] or ['base', 'all']
scs = (3, 8, 7)
print('| run | ' + ' | '.join(f'SC{s} flybys (fuel kg, s)' for s in scs) + ' | total flybys | valid | events |')
print('|---|' + '---|' * (len(scs) + 3))
for tag in tags:
    cells = []; tot = 0; ok = True; ev = []
    for s in scs:
        f = D / f'{tag}_sc{s}_info.json'
        if not f.exists():
            cells.append('-'); ok = False; continue
        r = json.load(open(f))
        tot += r['n_flybys']; ok = ok and r['ok']
        cells.append(f"{r['n_flybys']}/{r['planned']} ({r['fuel']:.0f}, {r['runtime']:.0f})")
        kinds = {}
        for e in r.get('events', []):
            k = e['kind'] + (':' + e['how'] if e['kind'] == 'replan' else '')
            kinds[k] = kinds.get(k, 0) + 1
        ev.append(f"SC{s}: " + ', '.join(f'{k}x{v}' for k, v in kinds.items()) if kinds else f'SC{s}: none')
    print(f'| {tag} | ' + ' | '.join(cells) + f' | {tot} | {"PASS" if ok else "FAIL/missing"} | ' + '; '.join(ev) + ' |')
