import os
for _v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, pathlib
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT/'tools')): sys.path.insert(0, _p)
os.chdir(ROOT)
import numpy as np, run_ialns as RI
F = RI.IFleet('results/s16/best/fleet')
r9 = set(int(a) for a in F.routes['r9']['st']['asts'])
cen = json.load(open('results/s17/lns/census_r9.json'))
best = {}
for c in cen['cands']:
    if c['ast'] not in best or c['lin_kg'] < best[c['ast']]['lin_kg']: best[c['ast']] = c
rows = []
for tag, host in (('L1r6k17S','r6'), ('L1r8k13S','r8')):
    own = set(int(a) for a in F.routes[host]['st']['asts']); tank0 = F.routes[host]['tank']
    for f in sorted((ROOT/'results/s17/lns/pre'/tag).glob('route_*.npz')):
        z = np.load(f); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL'])
        a = set(int(x) for x in st['asts']); tk = float(RI.ipr(st).tank())
        add = sorted(a & r9); drop = sorted(own - a)
        rows.append(dict(file=str(f.relative_to(ROOT)), host=host, n=len(a), tank=round(tk,1), dkg=round(tk-tank0,1), add=add,
                         add_census_kg={x: round(best[x]['lin_kg'],1) if x in best else None for x in add}, drop=drop))
for r in rows:
    if r['add']: print(r)
print(sum(1 for r in rows), 'columns;', sum(1 for r in rows if r['add']), 'with r9')
json.dump(rows, open('results/s17/redteam/r2/trades.json','w'), indent=1)
