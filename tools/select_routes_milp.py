"""Set-cover MILP over PLANNED routes (tour JSONs and/or route-pool jsonl files) -> a fleet of tours to convert.
   min sum_r c_r y_r + sum_a (1 - z_a),  z_a <= sum_{r covers a} y_r,  optional N <= --max-craft.
   c_r = cost_sc(600 + fuel_est*(1+--fuel-margin) + --reserve) (right-sized tank, capped at 2000).
Memory-lean: routes are held as packed bitmasks (38 bytes each); only the selected routes are re-read from disk.
Column filter: --per-target-k K keeps, for every target, the K routes with the lowest cost per target that contain it
(union over targets), so the MILP sees <= 298*K columns instead of the whole pool.
Usage: select_routes_milp.py outdir src [src ...] [--fuel-margin 0.3] [--reserve 20] [--max-craft 0] [--exclude ids.json]
       [--min-targets 2] [--per-target-k 0] [--time-limit 600]"""
import sys, json, glob, os, pathlib, argparse, time, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import csr_matrix
from ctoc14.constants import cost_sc, M_DRY, M0_MAX, DAY
ap = argparse.ArgumentParser(); ap.add_argument('outdir'); ap.add_argument('src', nargs='+'); ap.add_argument('--fuel-margin', type=float, default=0.3)
ap.add_argument('--reserve', type=float, default=20.0); ap.add_argument('--max-craft', type=int, default=0); ap.add_argument('--exclude', default=None)
ap.add_argument('--min-targets', type=int, default=2); ap.add_argument('--time-limit', type=float, default=600); ap.add_argument('--per-target-k', type=int, default=0)
ap.add_argument('--mip-gap', type=float, default=1e-3); ap.add_argument('--colgen', type=int, default=0, help='column-generation rounds: LP duals price all routes, add the most negative reduced-cost columns'); ap.add_argument('--colgen-add', type=int, default=1500)
a = ap.parse_args(); tic = time.time()
def tank(fuel): return min(M0_MAX, M_DRY + fuel * (1 + a.fuel_margin) + a.reserve)
masks, costs, nts, locs = [], [], [], []          # locs: (file, byte offset or -1 for whole file)
def add_route(t, loc):
    ids = [l['ast'] for l in t['legs']]
    if len(set(ids)) < a.min_targets: return
    m = np.zeros(300, np.uint8); m[[i - 1 for i in ids]] = 1
    masks.append(np.packbits(m)); costs.append(float(cost_sc(tank(t['fuel_est'])))); nts.append(len(set(ids))); locs.append(loc)
for s in a.src:
    if os.path.isdir(s):
        for f in sorted(glob.glob(f'{s}/tour_*.json')): add_route(json.load(open(f)), (f, -1))
    elif s.endswith('.jsonl'):
        with open(s, 'rb') as fh:
            off = 0
            for line in fh:
                if line.strip(): add_route(json.loads(line), (s, off))
                off += len(line)
    else: add_route(json.load(open(s)), (s, -1))
M = np.unpackbits(np.array(masks), axis=1)[:, :300].astype(bool); costs = np.array(costs); nts = np.array(nts)
print(f'{len(costs)} routes loaded ({time.time()-tic:.0f} s)', flush=True)
# exact-duplicate target sets: keep the cheapest
key = np.array([m.tobytes() for m in masks]); order = np.lexsort((costs, key)); keep = np.ones(len(costs), bool)
for i in range(1, len(order)):
    if key[order[i]] == key[order[i - 1]]: keep[order[i]] = False
excl = set(json.load(open(a.exclude))) if a.exclude else set()
targets = np.array([t for t in range(1, 301) if t not in (131, 144) and t not in excl]); tcol = targets - 1; nT = len(targets)
if a.per_target_k:
    cpt = costs / np.maximum(nts, 1); sel = np.zeros(len(costs), bool)
    for c in tcol:
        idx = np.where(M[:, c] & keep)[0]
        if len(idx): sel[idx[np.argsort(cpt[idx])[:a.per_target_k]]] = True
    keep &= sel
if a.colgen:
    from scipy.optimize import linprog
    Mall = M[:, tcol].astype(np.float64); cur = np.where(keep)[0] if a.per_target_k else np.array([], int)
    if len(cur) == 0:
        cpt = costs / np.maximum(nts, 1); cur = np.argsort(cpt)[:a.colgen_add]
    for rnd in range(a.colgen):
        n0 = len(cur); Mc = M[cur][:, tcol]
        rr, cc = np.nonzero(Mc.T); nrow = nT + (1 if a.max_craft else 0)
        rows = list(rr) + list(range(nT)) + ([nT] * n0 if a.max_craft else []); cols = list(cc) + list(range(n0, n0 + nT)) + (list(range(n0)) if a.max_craft else [])
        vals = [-1.0] * len(rr) + [1.0] * nT + ([1.0] * n0 if a.max_craft else [])
        A = csr_matrix((vals, (rows, cols)), shape=(nrow, n0 + nT)); hi = np.zeros(nrow)
        if a.max_craft: hi[nT] = a.max_craft
        lp = linprog(np.concatenate([costs[cur], -np.ones(nT)]), A_ub=A, b_ub=hi, bounds=[(0, 1)] * (n0 + nT), method='highs')
        duals = -lp.ineqlin.marginals[:nT]; mu = -lp.ineqlin.marginals[nT] if a.max_craft else 0.0
        rc = costs - Mall @ duals + mu                      # reduced cost of every route (uncovered -> keep only non-dup)
        rc[~keep] = np.inf; rc[cur] = np.inf
        cand = np.where(rc < -1e-6)[0]
        if len(cand) == 0: print(f'colgen round {rnd}: LP {lp.fun + nT + 2:.3f}, no negative reduced costs'); break
        addn = cand[np.argsort(rc[cand])[:a.colgen_add]]; cur = np.concatenate([cur, addn])
        print(f'colgen round {rnd}: LP {lp.fun + nT + 2:.3f}, {len(cand)} negative columns, added {len(addn)} -> {len(cur)}', flush=True)
    keep = np.zeros(len(costs), bool); keep[cur] = True
idx = np.where(keep)[0]; n = len(idx); Mk = M[idx][:, tcol]; ck = costs[idx]
print(f'{n} columns after dedup/filter, {nT} targets ({time.time()-tic:.0f} s)', flush=True)
# variables y (n), z (nT): min c.y - 1.z ; z_a - sum_{r∋a} y_r <= 0 ; sum y <= max_craft
rows = []; cols = []; vals = []
rr, cc = np.nonzero(Mk.T)                         # rr = target row, cc = route col
rows += list(rr); cols += list(cc); vals += [-1.0] * len(rr)
rows += list(range(nT)); cols += list(range(n, n + nT)); vals += [1.0] * nT
nrow = nT + (1 if a.max_craft else 0)
if a.max_craft: rows += [nT] * n; cols += list(range(n)); vals += [1.0] * n
A = csr_matrix((vals, (rows, cols)), shape=(nrow, n + nT))
lo = np.full(nrow, -np.inf); hi = np.zeros(nrow)
if a.max_craft: hi[nT] = a.max_craft
c = np.concatenate([ck, -np.ones(nT)])
res = milp(c, constraints=LinearConstraint(A, lo, hi), integrality=np.concatenate([np.ones(n), np.zeros(nT)]), bounds=Bounds(0, 1),
           options=dict(time_limit=a.time_limit, mip_rel_gap=a.mip_gap))
if res.x is None: print('MILP failed', res.message); sys.exit(1)
y = np.round(res.x[:n]).astype(int); sel = [j for j in range(n) if y[j]]
cov = np.zeros(300, bool)
for j in sel: cov |= M[idx[j]]
ncov = int(cov[tcol].sum()); J = float(ck[sel].sum()) + (nT - ncov) + 2 + len(excl)
out = pathlib.Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
def load(loc):
    f, off = loc
    if off < 0: return json.load(open(f))
    with open(f, 'rb') as fh: fh.seek(off); return json.loads(fh.readline())
for k, j in enumerate(sorted(sel, key=lambda j: -nts[idx[j]]), 1):
    t = load(locs[idx[j]]); t = dict(t); t['m0'] = float(tank(t['fuel_est'])); t['src'] = locs[idx[j]][0]
    json.dump(t, open(out / f'tour_sc{k}.json', 'w'), indent=1)
    print(f'  SC{k:2d}: {nts[idx[j]]:2d} targets, fuel_est {t["fuel_est"]:.0f} kg -> m0 {t["m0"]:.0f}, J_i {ck[j]:.3f}, launch {t["t_launch"]/DAY:.0f} d  [{os.path.basename(t["src"])}]')
missed = [int(t) for t in targets if not cov[t - 1]]
print(f'SELECTED {len(sel)} craft, covered {ncov}/{nT}, planned J (incl. 131/144{"+excluded" if excl else ""}) = {J:.3f}, gap {res.mip_gap if hasattr(res, "mip_gap") else "?"}, missed {missed}  ({time.time()-tic:.0f} s)')
json.dump(dict(J=J, covered=[int(t) for t in targets if cov[t - 1]], missed=missed, n=len(sel)), open(out / 'summary.json', 'w'))
