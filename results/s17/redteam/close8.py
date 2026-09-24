"""Leftover pricing of the pool-derived 8-craft skeleton: r1-r4 (+17 lin-screened absorptions, not applied to the
twins) + 4 pool routes maximising residual coverage (milp4_resid.txt).  For each leftover: approaches (<=0.3 AU) to all
8 routes and lin_price (as s15b_close._price_host)."""
import sys, json, csv, pathlib, collections
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
sys.path.insert(0, str(ROOT / 'results/s17/redteam'))
import numpy as np
import run_ialns as RI
import dissolve_price as D
RI.eph()
T = {int(r['ast']): r for r in csv.DictReader(open(ROOT / 'results/s17/redteam/targets.csv'))}
tail = sorted(k for k, v in T.items() if v['cls'] == 'tail')
absorb = [4, 300, 265, 142, 176, 14, 182, 130, 39, 184, 213, 18, 230, 56, 242, 26, 6]
F = D.load()
H = {k: F[k] for k in ('r1', 'r2', 'r3', 'r4', 'r8')}
for key in ('b33a7d00d0639bf6', '007b9f038060badd', '1f61b93d5dbe0ab6'):
    z = np.load(ROOT / f'results/s16/honest/route_{key}.npz'); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL']); H['p' + key[:4]] = st
cov = set(int(a) for st in H.values() for a in st['asts']) | set(absorb)
left = sorted(set(T) - cov - {131, 144})
left = [x for x in left if T[x]['cls'] != 'miss']
print('covered', len(cov & set(k for k, v in T.items() if v['cls'] != 'miss')), 'leftovers', len(left), left)
rows = []
for h, st in H.items():
    r = D.job(('LEFT', h, st, left, 0.3)); rows += r['rows']
by = collections.defaultdict(list)
for r in rows: by[r['ast']].append(r)
tot = 0
for x in left:
    L = sorted(by.get(x, []), key=lambda r: r['lin_dJ'])
    nd = min((r['dist'] for r in L), default=9)
    tot += L[0]['lin_dJ'] if L else 1.0
    print(f'  {x:4d} owner {T[x]["owner"]} nearest {nd:.3f} AU; cheapest: ' + ' | '.join(f"{r['host']} d{r['dist']:.3f} dJ{r['lin_dJ']:.3f}" for r in L[:3]))
print('sum of cheapest lin dJ (1.0 if none):', round(tot, 3))
json.dump(dict(left=left, rows=rows), open(ROOT / 'results/s17/redteam/close8.json', 'w'))
