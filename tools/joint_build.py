"""Settle the routes of a joint design (tools/carrier_joint.py best.json) into an IFleet. Usage: joint_build.py best.json out"""
import os
for _v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'): os.environ.setdefault(_v,'1')
import sys, json, pathlib, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np, run_ialns as RI, carrier_dp as CD
from ctoc14.globalopt import insert_block
CD.K_SLOPE = CD.K_PHYS * 1.5


def w(job):
    k, c = job
    full = c['base']
    for drop in [None] + list(range(len(full))):
        base = [b for i, b in enumerate(full) if i != drop]
        if len(base) < 4: break
        try:
            ip = CD.seed_ip(RI.eph(), c['tL'], np.array(c['vinf']),
                            np.array([b[1] for b in base]), np.array([b[2] for b in base]))
            if insert_block(ip, [(b[0], b[1]) for b in base], log=RI.QUIET).max() > 1e4: continue
            if RI.settle(ip, 100) <= 150: return k, (RI.ist(ip), float(ip.tank()))
        except Exception:
            continue
    return k, None


d = json.load(open(sys.argv[1])); out = pathlib.Path(sys.argv[2])
with mp.get_context('fork').Pool(8) as pl: res = dict(pl.map(w, list(enumerate(d['routes']))))
fl = RI.IFleet()
for k, r in sorted(res.items()):
    if r is None: print(f'{k:02d}: did not settle'); continue
    fl.routes[f'{k:02d}'] = dict(st=r[0], tank=r[1])
    print(f'{k:02d}: {len(r[0]["asts"])} fb @ {r[1]:.0f} kg')
fl.save(out, note='joint design bases'); print('covered', len(fl.coverage()))
