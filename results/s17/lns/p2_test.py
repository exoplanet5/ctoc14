"""P2 test: first steps of a lambda-ramped dissolution of r9 by point moves.
Truncations of r9 (remove + settle) and independent twin insertions of the removed targets into their census-best
host (s15b_close._insert: w_insert or w_insert_long, SIGALRM timeout).  Independent insertions into the same host do
not compound, so the insertion sum is OPTIMISTIC.  Writes results/s17/lns/p2_test.json."""
import os
for _v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, multiprocessing as mp
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT/'tools')):
    sys.path.insert(0, _p)
os.chdir(ROOT)
import numpy as np
import run_ialns as RI
import s15b_close as CL
F = RI.IFleet('results/s16/best/fleet')
cheap = [265, 176, 142, 208, 130, 6, 184]; mid = [243, 258, 180, 75]
cen = json.load(open('results/s17/lns/census_r9.json'))['cands']
best = {}
for c in cen:
    if c['ast'] in cheap + mid and (c['ast'] not in best or c['lin_kg'] < best[c['ast']]['lin_kg']):
        best[c['ast']] = c
r9 = F.routes['r9']
jobs_rem = [('r9', r9['st'], cheap), ('r9', r9['st'], cheap + mid)]
jobs_ins = [(best[x]['host'], F.routes[best[x]['host']]['st'], x, best[x]['t'], best[x]['dist'], 300.0,
             F.routes[best[x]['host']]['tank']) for x in cheap + mid]
tic = time.time()
with mp.get_context('fork').Pool(2, maxtasksperchild=2) as pool:
    ins = pool.map(CL._insert, jobs_ins, chunksize=1)
    rem = pool.map(RI.w_remove, jobs_rem, chunksize=1)
out = dict(wall=time.time() - tic, r9_tank=r9['tank'], rem=[], ins=[])
for j, r in zip(jobs_rem, rem):
    out['rem'].append(dict(removed=j[2], ok=r['ok'], tank=r.get('tank'), dJ=(RI.cost(r['tank']) - RI.cost(r9['tank'])) if r['ok'] else None))
for x, r in zip(cheap + mid, ins):
    out['ins'].append(dict(ast=x, host=r['host'], dist=r['dist'], lin_kg=best[x]['lin_kg'], ok=r['ok'], dkg=r.get('dkg'), dJ=r.get('dJ'), sec=r['sec']))
json.dump(out, open('results/s17/lns/p2_test.json', 'w'), indent=1, default=float)
for r in out['rem']: print('remove', len(r['removed']), 'ok', r['ok'], 'tank', r['tank'], 'dJ', r['dJ'])
for r in out['ins']: print('insert', r['ast'], '->', r['host'], f"dist {r['dist']:.3f} lin {r['lin_kg']:.0f} kg", 'ok', r['ok'], 'dkg', r['dkg'], 'dJ', r['dJ'], 'sec', r['sec'])
print('wall', out['wall'])
