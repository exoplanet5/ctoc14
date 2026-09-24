"""Settle the base routes of a tiling (tools/carrier_tile.py select) into a fleet + per-route pools.
Usage: carrier_build.py tile.json out_dir"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, pathlib, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI, carrier_dp as CD
from ctoc14.globalopt import insert_block


def w(job):
    k, c = job
    base = [b for b in c['base'] if b[0] in set(c['pool'])]
    if len(base) < 3: return k, None
    full = base
    for drop in [None] + list(range(len(full))):               # leave-one-out retries when the full base does not settle
        base = [b for i, b in enumerate(full) if i != drop]
        t = np.array([b[1] for b in base]); lag = np.array([b[2] for b in base])
        try:
            ip = CD.seed_ip(RI.eph(), c['tL'], np.array(c['vinf']), t, lag)
            d = insert_block(ip, [(b[0], b[1]) for b in base], log=RI.QUIET)
            if d.max() > 1e4: continue
            if RI.settle(ip, 100) <= 150: return k, (RI.ist(ip), float(ip.tank()))
        except Exception:
            continue
    return k, None


tile = json.load(open(sys.argv[1])); out = pathlib.Path(sys.argv[2])
with mp.get_context('fork').Pool(8) as pl: res = dict(pl.map(w, list(enumerate(tile['carriers']))))
fl = RI.IFleet(); pools = {}
for k, c in enumerate(tile['carriers']):
    n = f'{k:02d}'
    if res[k] is None: print(f'{n}: base did not settle'); continue
    fl.routes[n] = dict(st=res[k][0], tank=res[k][1]); pools[n] = c['pool']
    print(f'{n}: base {len(res[k][0]["asts"])} @ {res[k][1]:.0f} kg, pool {len(c["pool"])}')
fl.save(out, note='tiled carrier bases'); json.dump(pools, open(out / 'pools.json', 'w'))
print(f'covered {len(fl.coverage())}, pools union {len(set(t for p in pools.values() for t in p))}')
