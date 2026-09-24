"""How dense are the pool's columns? (targets per column, and the best coverage-per-cost columns)"""
import sys, glob, pathlib, numpy as np
sys.path.insert(0, '.')
from ctoc14.colgen import CostModel, ColumnStore
cm = CostModel(s=0.70, reserve=2.0); store = ColumnStore(cm)
for pat in ['results/phasing/camp*/pool.jsonl', 'results/phasing/pool/*.jsonl']:
    for f in sorted(glob.glob(pat)): store.load_jsonl(f, tag=pathlib.Path(f).parent.name)
for f in sorted(glob.glob('results/colgen/run*/columns.jsonl')) + sorted(glob.glob('results/colgen/*/columns.jsonl')):
    store.load_columns_file(f)
n = np.array([len(t) for t in store.tidx]); c = np.array(store.cost, float)
print(f'{len(n)} columns: targets per column p50 {np.percentile(n,50):.0f} p90 {np.percentile(n,90):.0f} '
      f'p99 {np.percentile(n,99):.0f} max {n.max()}')
ok = c < 9
print(f'columns with finite cost: {ok.sum()}; among them targets p90 {np.percentile(n[ok],90):.0f} max {n[ok].max()}')
for lo in (20, 25, 30, 35):
    m = ok & (n >= lo)
    print(f'  >= {lo} targets: {m.sum():6d} columns, best cost {c[m].min() if m.any() else float("nan"):.3f}, '
          f'best targets/cost {np.max(n[m]/c[m]) if m.any() else float("nan"):.1f}')
print(f'the flown 10-craft fleet: 24-37 targets per route at cost 1.156-1.341 (targets/cost 20.8-29.0)')
