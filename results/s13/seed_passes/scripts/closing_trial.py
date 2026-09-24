"""Read-only probe: REAL twin insertion prices of a pass's uncovered targets into its beam-built routes."""
import os, sys, json, time, pathlib
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
SP = pathlib.Path(__file__).resolve().parent
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
os.chdir(ROOT)
import numpy as np, multiprocessing as mp
import run_ialns as RI
from insert_model import predict_dJ, predict_ok
from route_cover import ALL
V = sys.argv[1]; ntgt = int(sys.argv[2]); per = int(sys.argv[3])
fl = RI.IFleet(SP / 'passes' / V)
cov = fl.coverage(); U = sorted(set(ALL) - set(cov))
with mp.get_context('fork').Pool(6) as pl:
    res = pl.map(RI.w_cands, [(n, r['st'], U, 0.3) for n, r in fl.routes.items()])
cands = []
for n, lst in zip(fl.routes, res):
    r = fl.routes[n]; depth = len(r['st']['asts']); tank = r['tank']; tf = np.sort(np.asarray(r['st']['tf'], float))
    for c in lst:
        marg = float(np.min(np.abs(tf - c['t']))) / 86400.0
        p = predict_ok(c['dist'], depth, tank, marg); dj = predict_dJ(c['dist'], depth, tank, marg)
        cands.append(dict(host=n, ast=c['ast'], t=c['t'], d=c['dist'], depth=depth, tank=tank, p=p, pdj=dj, e=p * dj + (1 - p)))
by = {}
for c in sorted(cands, key=lambda c: c['e']):
    by.setdefault(c['ast'], [])
    if len(by[c['ast']]) < per: by[c['ast']].append(c)
tg = sorted(by, key=lambda a: by[a][0]['e'])[:ntgt]
jobs = [(c['host'], fl.routes[c['host']]['st'], c['ast'], c['t']) for a in tg for c in by[a]]
meta = {(c['host'], c['ast'], c['t']): c for a in tg for c in by[a]}
print(f'{V}: {len(U)} uncovered; trying {len(jobs)} insertions for {len(tg)} targets', flush=True)
tic = time.time()
with mp.get_context('fork').Pool(6, maxtasksperchild=4) as pl:
    out = pl.map(RI.w_insert, jobs)
rows = []
for r in out:
    m = meta[(r['name'], r['ast'], r['t'])]
    dj = (RI.cost(r['tank']) - RI.cost(fl.routes[r['name']]['tank'])) if r.get('ok') else None
    rows.append(dict(ast=r['ast'], host=r['name'], d=round(m['d'], 3), depth=m['depth'], ok=bool(r.get('ok')),
                     kg=None if not r.get('ok') else round(r['tank'] - fl.routes[r['name']]['tank'], 1),
                     dJ=None if dj is None else round(dj, 4), pred_dJ=round(m['pdj'], 4), pred_p=round(m['p'], 2)))
best = {}
for x in rows:
    if x['ok'] and (x['ast'] not in best or x['dJ'] < best[x['ast']]['dJ']): best[x['ast']] = x
print(f'{len(out)} trials in {time.time()-tic:.0f} s; {sum(x["ok"] for x in rows)} settled; targets placed {len(best)}/{len(tg)}')
b = sorted(best.values(), key=lambda x: x['dJ'])
print('best real dJ per placed target:', [x['dJ'] for x in b])
print('their predicted dJ           :', [x['pred_dJ'] for x in b])
print('real kg                      :', [x['kg'] for x in b])
print('distances                    :', [x['d'] for x in b])
json.dump(rows, open(SP / 'passes' / f'closing_{V}.json', 'w'), indent=1)
