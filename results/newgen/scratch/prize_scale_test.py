"""How does the pricing beam's depth / reduced cost depend on the prize scale? Uses the run's current duals."""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'): os.environ.setdefault(_v, '1')
import sys, json, glob, pathlib, time, numpy as np
sys.path.insert(0, '.')
from ctoc14.kepler import Ephemeris
from ctoc14.colgen import CostModel, ColumnStore, PricingSpec, price, TARGETS
d = json.load(open(sys.argv[1])); pi = np.zeros(300)
for k, v in d['pi'].items(): pi[int(k) - 1] = v
mu = float(sys.argv[2]); scales = [float(x) for x in sys.argv[3].split(',')]
eph = Ephemeris(); store = ColumnStore(CostModel(s=0.70, reserve=2.0))
for sc in scales:
    spec = PricingSpec(tag=f's{sc}', launch=(0.0, 320.0, 20.0), beam=300, nproc=4, m0=1000.0, seed=1)
    n0 = len(store)
    added, st = price(eph, store, pi * sc, pi, mu, (), spec, log=print)
    n = np.array([store.n[j] for j in range(n0, len(store))]); rc = np.array([store.src[j][1].get('rc', 0) for j in range(n0, len(store))])
    if len(n):
        o = np.argsort(rc)[:5]
        print(f'  scale {sc}: added {len(n)}, depth p50 {np.median(n):.0f} max {n.max()}, best five rc/depth: ' + ' '.join(f'{rc[i]:.2f}/{n[i]}' for i in o), flush=True)
