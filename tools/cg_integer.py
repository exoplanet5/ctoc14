"""Integer fleets from a column-generation run WITHOUT new pricing (fast, for comparison with the priced dive):
  1. LP at N over the whole store (exact), 2. MIP restricted to the LP support (+ the K most negative-rc neighbours),
  3. LP-guided dive without re-pricing (fix best look-ahead route, remove its targets, re-solve).
Usage: cg_integer.py cgdir --pools ... [--N 10] [--support-extra 300] [--mip-time 120] [--tag int]
Writes cgdir/<tag>_supmip_N{N}/ and cgdir/<tag>_lpdive_N{N}/ (tour_sc*.json + summary.json)."""
import sys, json, time, pathlib, argparse, glob
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from ctoc14.colgen import CostModel, ColumnStore, RMP, fleet_summary, write_fleet, TARGETS

ap = argparse.ArgumentParser(); ap.add_argument('cgdir'); ap.add_argument('--pools', nargs='*', default=[])
ap.add_argument('--N', type=int, default=10); ap.add_argument('--support-extra', type=int, default=300)
ap.add_argument('--mip-time', type=float, default=120.0); ap.add_argument('--tag', default='int')
ap.add_argument('--cost-s', type=float, default=1.15); ap.add_argument('--reserve', type=float, default=20.0)
a = ap.parse_args(); d = pathlib.Path(a.cgdir); tic = time.time()
store = ColumnStore(CostModel(s=a.cost_s, reserve=a.reserve))
for pat in a.pools:
    for f in sorted(glob.glob(pat)): store.load_jsonl(f, tag=pathlib.Path(f).parent.name + '/' + pathlib.Path(f).stem)
if (d / 'columns.jsonl').exists(): store.load_columns_file(d / 'columns.jsonl')
print(f'store {len(store)} columns ({time.time()-tic:.0f} s)', flush=True)
rmp = RMP(store); rows = TARGETS - 1
res = rmp.solve(rows, a.N)
print(f'LP N={a.N}: {res.value:.4f}, sum y {res.y.sum():.2f}, support {(res.y > 1e-6).sum()}', flush=True)
# 2. MIP on the LP support + nearest-rc columns
C = store.csr(); rc = store.costs() - C @ res.pi + res.mu
sup = res.cols[res.y > 1e-6]
extra = np.argsort(rc)[:a.support_extra + len(sup)]
cols = np.unique(np.concatenate([sup, extra]))
J, sel, missed = rmp.mip(rows, a.N, time_limit=a.mip_time, cols=cols)
s1 = write_fleet(store, sel, d / f'{a.tag}_supmip_N{a.N}', extra=dict(lp=res.value, ncols=int(len(cols))))
print(f'support MIP over {len(cols)} columns: J {J:.3f}, {s1["n"]} craft, covered {s1["covered"]}, sum J_i {s1["sumJ"]:.3f}, missed {s1["missed"]} ({time.time()-tic:.0f} s)', flush=True)
# 3. LP-guided dive without pricing
rows_d = rows.copy(); fixed = []; fc = 0.0
for level in range(a.N):
    r = rmp.solve(rows_d, a.N - level, fixed_cost=fc)
    order = np.argsort(-r.y)[:4]; best = None
    for i in order:
        if r.y[i] <= 1e-6: continue
        j = int(r.cols[i]); tj = np.array(store.targets(j)) - 1; rj = np.setdiff1d(rows_d, tj)
        v = rmp.lp_value_fixed(rj, a.N - level - 1, fc + store.cost[j])
        if best is None or v < best[0]: best = (v, j, rj)
    if best is None: break
    v, j, rows_d = best; fixed.append(j); fc += store.cost[j]
    print(f'  dive level {level}: fix {j} [{store.tag[j]}] {store.n[j]} targets, look-ahead LP {v:.3f}, rows left {len(rows_d)}', flush=True)
    if len(rows_d) == 0: break
s2 = write_fleet(store, fixed, d / f'{a.tag}_lpdive_N{a.N}')
print(f'LP dive: J {s2["J"]:.3f}, {s2["n"]} craft, covered {s2["covered"]}, sum J_i {s2["sumJ"]:.3f}, missed {s2["missed"]} ({time.time()-tic:.0f} s)')
