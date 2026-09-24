"""Rebuild a skeleton dir from a partition.json (stage 6).

`skel_partition.py` writes its best partition as {route: [hard target ids]}; the epochs are not in that file because
a move never changes them -- every hard target keeps the t10d crossing it was found at.  This script looks each
epoch up in a reference skeleton and writes a skeleton IFleet, so a search can be resumed from its own best result
(or that result can be re-filled at production settings).

Usage: skel_frompart.py partition.json out_dir [--ref results/n8/skel8b]
"""
import sys, json, pathlib, argparse
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from ctoc14.constants import DAY

ap = argparse.ArgumentParser(); ap.add_argument('part'); ap.add_argument('out')
ap.add_argument('--ref', default='results/n8/skel8b')
a = ap.parse_args()

ref = RI.IFleet(a.ref)
epoch = {}
st0 = {}
for n, r in ref.routes.items():
    st0[n] = r['st']
    for x, t in zip(r['st']['asts'], r['st']['tf']): epoch[int(x)] = float(t)
part = json.load(open(a.part))
miss = [t for v in part.values() for t in v if t not in epoch]
if miss:
    print(f'no epoch in {a.ref} for {miss}'); sys.exit(1)

out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
fl = RI.IFleet()
for n, ts in part.items():
    wps = sorted((epoch[t], int(t)) for t in ts)
    st = dict(st0[n])
    st['asts'] = np.array([x for _, x in wps], float)
    st['tf'] = np.array([t for t, _ in wps], float)
    fl.routes[n] = dict(st=st, tank=ref.routes[n].get('tank', 700.0))
    print(f'  {n}: {len(wps)} waypoints, {wps[0][0]/DAY:.0f}..{wps[-1][0]/DAY:.0f} d')
fl.save(out, note='from partition.json')
print(f'-> {out}')
