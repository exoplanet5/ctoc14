"""Scarcity prizes and fill order from a column library (stage 6 C).

Sequential filling starves the last routes.  fill2 filled by waypoint count, so skeleton 13 -- which reaches 154 of
the 231 easy targets -- went third while skeleton 10, which reaches only 82, went seventh and got 23 flybys.  Two
fixes come straight out of the library's reach matrix:

  order   : fill the NARROWEST skeletons first (fewest reachable easy targets), so a route with few options takes
            them before a broad route spends them.
  prizes  : prize(t) = min(cap, (nskel / reach(t)) ** gamma) -- a target only one skeleton can fly is worth more
            than one any of eight can.  Whoever can take it does so early, and broad routes drift to common targets.
            (Stage 4's "hardness prizes" boosted GEOMETRICALLY hard targets and lost depth; this boosts targets that
            are scarce in the MEASURED reach matrix, which is a different and much better conditioned signal.)

Usage: scarcity_prizes.py 'results/n8s6/B/cols_*.jsonl' prizes.json [--gamma 0.6] [--cap 2.0]
"""
import sys, json, glob, argparse, pathlib, collections
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
from ctoc14.search import UNREACHABLE

ap = argparse.ArgumentParser(); ap.add_argument('cols', nargs='+'); ap.add_argument('out')
ap.add_argument('--gamma', type=float, default=0.6); ap.add_argument('--cap', type=float, default=2.0)
ap.add_argument('--skel-dir', default='results/n8/skel8b')
ap.add_argument('--unreached', type=float, default=2.0, help='prize for easy targets no column reached at all')
a = ap.parse_args()

cols = []
for pat in a.cols:
    for f in sorted(glob.glob(pat)): cols += [json.loads(l) for l in open(f)]
import run_ialns as RI
skel = RI.IFleet(a.skel_dir)
hard = {int(t): n for n, r in skel.routes.items() for t in r['st']['asts']}
ALL = [t for t in range(1, 301) if t not in UNREACHABLE]
easy = [t for t in ALL if t not in hard]
skels = sorted(skel.routes)
R = {t: set() for t in easy}; per = {s: set() for s in skels}
for c in cols:
    for t in c['targets']:
        if t not in hard:
            R[t].add(c['skel']); per[c['skel']].add(t)

PZ = np.ones(300)
for t in easy:
    n = len(R[t])
    PZ[t - 1] = a.unreached if n == 0 else min(a.cap, (len(skels) / n) ** a.gamma)
hist = collections.Counter(len(R[t]) for t in easy)
print('reach histogram (skeletons per easy target): ' + ' '.join(f'{k}:{hist[k]}' for k in sorted(hist)))
print('prize by reach count: ' + ' '.join(
    f'{k}->{min(a.cap,(len(skels)/max(k,1))**a.gamma) if k else a.unreached:.2f}' for k in sorted(hist)))
order = sorted(skels, key=lambda s: len(per[s]))
print('reach breadth per skeleton: ' + ' '.join(f'{s}:{len(per[s])}' for s in order))
print('suggested fill order (narrowest first): ' + ','.join(order))
json.dump([float(x) for x in PZ], open(a.out, 'w'))
print(f'-> {a.out}  (easy prizes {PZ[[t-1 for t in easy]].min():.2f}..{PZ[[t-1 for t in easy]].max():.2f})')
