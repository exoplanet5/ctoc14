import sys, glob, time, pathlib, numpy as np
sys.path.insert(0, '.')
from ctoc14.colgen import CostModel, ColumnStore, RMP, TARGETS
s, R = float(sys.argv[1]), float(sys.argv[2])
cm = CostModel(s=s, reserve=R); store = ColumnStore(cm); tic = time.time()
for pat in ['results/phasing/camp*/pool.jsonl', 'results/phasing/pool/*.jsonl']:
    for f in sorted(glob.glob(pat)): store.load_jsonl(f, tag=pathlib.Path(f).parent.name)
for f in sorted(glob.glob('results/colgen/run*/columns.jsonl')) + sorted(glob.glob('results/colgen/*/columns.jsonl')):
    store.load_columns_file(f)
print(f's={s} R={R}: store {len(store)} ({time.time()-tic:.0f} s)', flush=True)
rows = TARGETS - 1
for N in [13, 12, 11, 10, 9]:
    rmp = RMP(store); r = rmp.solve(rows, N)
    print(f'  N={N}: LP {r.value:.3f}, sum y {r.y.sum():.2f}, mu {r.mu:.3f}', flush=True)
