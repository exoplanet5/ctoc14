"""Neighbourhood matrix: for each s16 route, approach minima (run_ialns.w_cands, 2-d sampled track, no 10-d
exclusion) to every other covered target; counts by owner at 0.02/0.035/0.06/0.10 AU.  -> nbr.json"""
import sys, json, pathlib, collections
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
RI.eph()
F = {}
for f in sorted((ROOT / 'results/s16/best/fleet').glob('route_*.npz')):
    z = np.load(f); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL']); F[f.stem[6:]] = st
owner = {int(a): k for k, st in F.items() for a in st['asts']}
allT = sorted(owner)
R = {}
for k, st in F.items():
    apps = RI.w_cands((k, st, allT, 0.10))
    best = {}
    for a in apps:
        if a['ast'] not in best or a['dist'] < best[a['ast']]:
            best[a['ast']] = a['dist']
    row = {}
    for th in (0.02, 0.035, 0.06, 0.10):
        c = collections.Counter(owner[x] for x, d in best.items() if d <= th)
        row[str(th)] = dict(c)
    R[k] = dict(row=row, best={int(x): round(d, 4) for x, d in best.items()})
    print(k, len(st['asts']), {th: sum(v.values()) for th, v in row.items()}, 'at 0.035:', dict(sorted(row['0.035'].items())))
# targets with no route within 0.035 / 0.06 except the owner
iso = {}
for th in (0.035, 0.06, 0.10):
    cnt = sum(1 for x in allT if not any(x in R[k]['best'] and R[k]['best'][x] <= th for k in F if k != owner[x]))
    by = collections.Counter(owner[x] for x in allT if not any(x in R[k]['best'] and R[k]['best'][x] <= th for k in F if k != owner[x]))
    iso[str(th)] = dict(n=cnt, by_owner=dict(sorted(by.items())))
print('targets with NO non-owner route within th:', iso)
json.dump(dict(routes=R, isolated=iso), open(ROOT / 'results/s17/redteam/nbr.json', 'w'))
