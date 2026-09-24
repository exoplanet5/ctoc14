"""Promote a filled fleet's LEFTOVER targets to mandatory waypoints, then let the fill re-plan around them.

Stage 5 inserted leftovers into a FIXED tour and paid 100-280 kg each (the route was already time-full).  Stage 5 §1
also showed the opposite: a target a route is PLANNED AROUND costs the same as any other (0.51 vs 0.48 km/s).  So the
cheap way to place a leftover is not insertion but re-planning: give it to a route as a mandatory time window and
re-run the segmented fill, which rebuilds the whole easy sequence around the enlarged waypoint set.

Epochs come from `run_ialns.w_cands`: the local minima of the distance between the leftover and the route's current
trajectory, i.e. the times that route is already near it.  Each leftover goes to the route whose closest approach is
smallest, subject to --cap extra waypoints per route; leftovers with no approach within --dmax are reported and left
out.  The output is a SKELETON fleet (waypoint ids + epochs only, no settling) for `skel_fill --mode segmented`.

Usage: skel_addwp.py filled_dir skel_dir out_dir [--dmax 0.5] [--cap 8] [--nproc 8]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, pathlib, argparse, collections, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from ctoc14.search import UNREACHABLE
from ctoc14.constants import DAY

ap = argparse.ArgumentParser(); ap.add_argument('filled'); ap.add_argument('skel'); ap.add_argument('out')
ap.add_argument('--dmax', type=float, default=0.5, help='AU: max closest approach to accept an epoch')
ap.add_argument('--cap', type=int, default=8, help='max extra waypoints per route')
ap.add_argument('--nproc', type=int, default=8)
ap.add_argument('--min-gap', type=float, default=40.0, help='days: keep new waypoints this far from existing ones')
a = ap.parse_args()
ALL = [t for t in range(1, 301) if t not in UNREACHABLE]

filled = RI.IFleet(a.filled); skel = RI.IFleet(a.skel)
cov = set(filled.coverage()); left = [t for t in ALL if t not in cov]
print(f'filled fleet: {filled.summary()}')
print(f'{len(left)} leftovers: {left}')
if not left:
    print('nothing to do'); sys.exit(0)

with mp.get_context('fork').Pool(a.nproc) as pool:
    jobs = [(n, r['st'], left, a.dmax) for n, r in filled.routes.items()]
    res = pool.map(RI.w_cands, jobs)
cands = collections.defaultdict(list)
for lst in res:
    for c in lst: cands[c['ast']].append(c)
for t in cands: cands[t].sort(key=lambda c: c['dist'])

# existing waypoint epochs per route (hard skeleton), used to keep new windows from colliding
wp = {n: [(float(t), int(x)) for x, t in zip(r['st']['asts'], r['st']['tf'])] for n, r in skel.routes.items()}
extra = {n: [] for n in skel.routes}
noc = [t for t in left if not cands[t]]
order = sorted((t for t in left if cands[t]), key=lambda t: (len(cands[t]), cands[t][0]['dist']))
placed = {}
for t in order:
    for c in cands[t]:
        n = c['host']
        if n not in extra or len(extra[n]) >= a.cap: continue
        ts = [x[0] for x in wp[n]] + [x[0] for x in extra[n]]
        if any(abs(c['t'] - u) < a.min_gap * DAY for u in ts): continue
        extra[n].append((c['t'], t)); placed[t] = (n, c['dist']); break
un = [t for t in left if t not in placed]
print(f'\nplaced {len(placed)} of {len(left)}: ' + ' '.join(f'{n}:{len(extra[n])}' for n in sorted(extra)))
print(f'no approach within {a.dmax} AU ({len(noc)}): {noc}')
print(f'unplaced (cap/gap) ({len(un)}): {un}')
d = [v[1] for v in placed.values()]
if d: print(f'closest approach of the placed: min {min(d):.3f} med {np.median(d):.3f} max {max(d):.3f} AU')

out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
new = RI.IFleet()
for n, r in skel.routes.items():
    wps = sorted(wp[n] + extra[n])
    st = dict(r['st'])
    st['asts'] = np.array([x for _, x in wps], float)
    st['tf'] = np.array([t for t, _ in wps], float)
    new.routes[n] = dict(st=st, tank=r.get('tank', 700.0))
    print(f'  {n}: {len(wp[n])} hard + {len(extra[n])} promoted = {len(wps)} waypoints'
          + (f'  (+{sorted(x for _, x in extra[n])})' if extra[n] else ''))
new.save(out, note='hard + promoted leftovers')
json.dump(dict(placed={str(k): v[0] for k, v in placed.items()}, unplaced=un, no_candidate=noc),
          open(out / 'promote.json', 'w'), indent=1)
print(f'-> {out}')
