"""Two-carrier test: insert a SECOND carrier's whole base into a settled carrier route as one block."""
import os
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'): os.environ.setdefault(v,'1')
import sys, json, itertools, multiprocessing as mp; sys.path.insert(0,'.'); sys.path.insert(0,'tools')
import numpy as np, run_ialns as RI
from ctoc14.globalopt import insert_block
from ctoc14.constants import VE, DAY

fl = RI.IFleet('results/s8/J9base')
D = json.load(open('results/s8/J9b/best.json'))
bases = {f'{k:02d}': r['base'] for k, r in enumerate(D['routes'])}

def job(p):
    host, don = p
    st = fl.routes[host]['st']; tank0 = fl.routes[host]['tank']
    have = set(int(x) for x in st['asts'])
    items = [(int(t), float(tt)) for t, tt, _ in bases[don] if int(t) not in have]
    tf = np.asarray(st['tf'])
    items = [(t, tt) for t, tt in items if np.abs(tf - tt).min() > 10 * DAY]
    if len(items) < 8: return None
    ip = RI.ipr(st)
    try:
        d = insert_block(ip, items, log=RI.QUIET)
        if d.max() > 1e4: return (host, don, len(items), None, None)
        miss = RI.settle(ip, 100)
        if miss > 150: return (host, don, len(items), None, None)
    except Exception:
        return (host, don, len(items), None, None)
    tk = float(ip.tank()); n = len(ip.asts)
    return (host, don, len(items), n, tk)

pairs = [(h, d) for h in fl.routes for d in bases if h != d]
np.random.default_rng(0).shuffle(pairs)
with mp.get_context('fork').Pool(8) as pl:
    res = [r for r in pl.imap_unordered(job, pairs[:24]) if r]
ok = [r for r in res if r[3]]
print(f'{len(ok)} of {len(res)} block insertions settled\n')
for h, d, ni, n, tk in sorted(ok, key=lambda r: RI.cost(r[4]) / r[3]):
    t0 = fl.routes[h]['tank']; n0 = len(fl.routes[h]['st']['asts'])
    print(f'  {h} ({n0} fb @ {t0:.0f}) + carrier {d} block of {ni} -> {n} fb @ {tk:.0f} kg  '
          f'({(tk-t0)/max(1,n-n0):.0f} kg per new target, {RI.cost(tk)/n:.4f} J/fb)')
print('\n  reference: single-carrier base 13 kg/target; greedy insertion 45-90 kg; t10d fleet 0.0415 J/fb')
