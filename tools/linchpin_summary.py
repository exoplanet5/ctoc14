"""Stage-13 linchpin: tabulate linchpin_replan result.json files (one line per run) and, per condition, the fleet level
(best settled candidate per route across builders, ranked by pool targets flown then tank).
Usage: linchpin_summary.py dir [dir ...] [--lam 0.15]"""
import sys, json, pathlib, argparse
import numpy as np

ap = argparse.ArgumentParser(); ap.add_argument('dirs', nargs='+'); ap.add_argument('--lam', type=float, default=0.15)
ap.add_argument('--by', default='k', choices=['k', 'value'])
a = ap.parse_args()
cost = lambda t: 1 + (t - 600) / 1400 + ((t - 600) / 1400) ** 2
rows = []
for d in a.dirs:
    for f in sorted(pathlib.Path(d).glob('*/result.json')):
        r = json.load(open(f)); r['dir'] = f.parent.name; rows.append(r)
print(f'{"run":28s} {"pool":>4s} {"deep":>4s} | {"k":>3s} {"n":>3s} {"tank":>7s} {"kg/fb":>6s} {"kg/pool":>7s} | src {"fb":>3s} '
      f'{"kg/fb":>6s} | own/ext flown | settles ok/tried | t_build')
for r in rows:
    ok = [s for s in r['settled'] if s.get('ok')]
    b = (max(ok, key=lambda s: (s['k'], -s['tank'])) if a.by == 'k' else min(ok, key=lambda s: s['value'])) if ok else None
    src_kgfb = (r['source_tank'] - 600) / r['source_fb']
    if b:
        print(f'{r["dir"]:28s} {r["pool_size"]:4d} {r["deepest"]:4d} | {b["k"]:3d} {b["n"]:3d} {b["tank"]:7.1f} {b["kg_per_fb"]:6.2f} '
              f'{(b["tank"] - 600) / b["k"]:7.2f} | {r["source_fb"]:6d} {src_kgfb:6.2f} | {b["own_flown"]:3d}/{b["extra_flown"]:<3d}      '
              f'| {len(ok)}/{len(r["settled"])} | {r["t_build_s"]:.0f} s')
    else:
        print(f'{r["dir"]:28s} {r["pool_size"]:4d} {r["deepest"]:4d} | nothing settled ({len(r["settled"])} tried)        '
              f'| src {r["source_fb"]} {src_kgfb:.2f}')
