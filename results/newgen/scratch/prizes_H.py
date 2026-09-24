"""R0.1 -- hardness prizes for the N=8 fleet tree (docs/stage4_n8_search_tree.md section 2).

prize_t = clip( 1 + 1.2 * (1 - f_t/median_f)+ + 0.6 * o_t , 1.0, 3.0 ),  unreachable -> 0
  f_t = share of the >=34-target columns of results/newgen/pool16.jsonl that contain target t
  o_t = 1 if t is in the greedy-cover orphan list (results/newgen/gc16, the 67 left over), else 0
Writes results/n8/prizes_H.json (300 floats, index t-1) and prints the histogram + the hard list (prize>=1.5).
"""
import sys, json, re, pathlib, argparse
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from ctoc14.search import UNREACHABLE

# cs was 1.2 in the doc, but the deep-column frequency is heavily skewed (80 reachable targets never appear in a deep
# column, p25(f)=0), so cs=1.2 flags 156 hard -- which auto-trips branch-A prune (iii) (hard_left > 12*(8-d), budget 96
# at the root, calibrated for ~87 hard).  cs=0.45 keeps |prize>=1.5| = the 67 structural orphans while still lifting the
# f=0 tail to ~1.45 (preferred by the beam, not hard-flagged).  Branch A's runtime remedy = drop cs to ~0.22.
ap = argparse.ArgumentParser()
ap.add_argument('--cs', type=float, default=0.45, help='scarcity coefficient')
ap.add_argument('--orphan', type=float, default=0.6, help='bonus for greedy-cover orphans')
a = ap.parse_args()

ALL = [t for t in range(1, 301) if t not in UNREACHABLE]
POOL = ROOT / 'results/newgen/pool16.jsonl'
GC = ROOT / 'results/newgen/gc16/log.txt'
OUT = ROOT / 'results/n8/prizes_H.json'

# --- f_t : column frequency among deep (>=34) columns
cnt = np.zeros(301, int); ncol = 0
for line in open(POOL):
    line = line.strip()
    if not line:
        continue
    r = json.loads(line)
    tg = r.get('targets') or []
    if len(tg) < 34:
        continue
    ncol += 1
    for t in tg:
        cnt[t] += 1
f = cnt / max(ncol, 1)
fr = np.array([f[t] for t in ALL])
med = float(np.median(fr))

# --- o_t : greedy-cover orphans
m = re.search(r'uncovered \(\d+\): \[([0-9,\s]+)\]', open(GC).read())
orphans = set(int(x) for x in m.group(1).split(',')) if m else set()

# --- prize
pz = np.zeros(300)
for t in ALL:
    scarce = max(0.0, 1.0 - f[t] / med) if med > 0 else 0.0
    pz[t - 1] = float(np.clip(1.0 + a.cs * scarce + a.orphan * (1.0 if t in orphans else 0.0), 1.0, 3.0))

OUT.parent.mkdir(parents=True, exist_ok=True)
json.dump([round(float(x), 4) for x in pz], open(OUT, 'w'))

vals = np.array([pz[t - 1] for t in ALL])
hard = sorted(t for t in ALL if pz[t - 1] >= 1.5)
print(f'{ncol} deep (>=34) columns in {POOL.name}; median f = {med:.3f}; {len(orphans)} orphans')
print(f'reachable prizes: min {vals.min():.2f}  median {np.median(vals):.2f}  mean {vals.mean():.2f}  max {vals.max():.2f}')
edges = [1.0, 1.01, 1.25, 1.5, 1.75, 2.0, 2.5, 3.001]
h, _ = np.histogram(vals, bins=edges)
labels = ['=1.00', '1.01-1.25', '1.25-1.50', '1.50-1.75', '1.75-2.00', '2.00-2.50', '2.50-3.00']
for lab, c in zip(labels, h):
    print(f'  {lab:>11}: {c:3d}  ' + '#' * int(c))
print(f'hard (prize>=1.5): {len(hard)} targets')
print(hard)
print(f'wrote {OUT}')
