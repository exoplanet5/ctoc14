"""Partition MIP: pick N columns (routes) that cover as many targets as possible, nearly disjointly, at low cost.
Columns come from twin-recost files (`--twin-files`: ok records with tank + inline/npz state) and/or from planner
pools (`--plan-files`: tour_json records, cost from CostModel(s, reserve) -- calibrated to the twin within 2 %).
max sum_t w_t c_t - sum_j cost_j y_j  s.t.  c_t <= sum_{j in t} y_j <= 1 + s_t (s_t = 0 on hard targets), sum y = N.
Usage: partition_mip.py out_dir [--N 8] [--twin-files f,f] [--plan-files f,f] [--w-hard .6] [--time 600]
"""
import sys, json, glob, pathlib, argparse, numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix, csr_matrix, vstack, hstack, identity
from ctoc14.constants import cost_sc
from ctoc14.colgen import CostModel
from ctoc14.search import UNREACHABLE

ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--N', default='8')
ap.add_argument('--twin-files', default='results/newgen/recost35.jsonl')
ap.add_argument('--plan-files', default='')
ap.add_argument('--w-hard', type=float, default=0.6); ap.add_argument('--w-easy', type=float, default=0.03)
ap.add_argument('--overlap-pen', type=float, default=0.05); ap.add_argument('--time', type=float, default=600.0)
ap.add_argument('--max-tank', type=float, default=1050.0); ap.add_argument('--min-targets', type=int, default=0)
ap.add_argument('--min-hard', type=int, default=0); ap.add_argument('--hard-thr', type=int, default=5)
ap.add_argument('--easy-slack', type=float, default=3.0); ap.add_argument('--save-states', action='store_true')
a = ap.parse_args()
out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
ALL = [t for t in range(1, 301) if t not in UNREACHABLE]

cnt = np.zeros(301, int)
for line in open(ROOT / 'results/newgen/recost35.jsonl'):
    r = json.loads(line)
    if r.get('ok'):
        for t in r['targets']: cnt[t] += 1
hard = set(t for t in ALL if cnt[t] <= a.hard_thr); easy = set(ALL) - hard
print(f'{len(hard)} hard targets, {len(easy)} easy')

cm = CostModel(s=0.70, reserve=2.0)
cols = {}
def add(targets, cost, src, tank=None):
    k = frozenset(int(t) for t in targets if int(t) not in UNREACHABLE)
    if not k or cost >= 9.0: return
    if k not in cols or cost < cols[k]['cost']:
        cols[k] = dict(targets=sorted(k), cost=float(cost), src=src, tank=(float(tank) if tank else None))

for f in [x for x in a.twin_files.split(',') if x]:
    for g in sorted(glob.glob(str(ROOT / f))) or [f]:
        n0 = len(cols)
        for line in open(g):
            r = json.loads(line)
            if r.get('ok') is False or not r.get('tank') or r['tank'] > a.max_tank: continue
            st = r.get('st')
            if st is None:
                q = pathlib.Path(g).parent / f'twin_states/col_{r["col"]}.npz'
                st = ('npz', str(q)) if q.exists() else None
            else:
                st = ('inline', st)
            if st is None: continue
            add(r['targets'], cost_sc(r['tank']), st, r['tank'])
        print(f'  twin {pathlib.Path(g).name}: {len(cols)-n0} new')
for f in [x for x in a.plan_files.split(',') if x]:
    for g in sorted(glob.glob(str(ROOT / f))) or [f]:
        n0 = len(cols)
        for line in open(g):
            r = json.loads(line)
            tg = r.get('targets') or [l['ast'] for l in r['legs']]
            add(tg, float(cm.cost(r['fuel_est'], r['m0'])), ('plan', r), None)
        print(f'  plan {pathlib.Path(g).name}: {len(cols)-n0} new')

C = [c for c in cols.values() if len(c['targets']) >= a.min_targets and len(set(c['targets']) & hard) >= a.min_hard]
n = len(C); print(f'{n} distinct columns (of {len(cols)}) after filters')
T = sorted(hard | easy); tix = {t: i for i, t in enumerate(T)}; m = len(T)
A = lil_matrix((m, n))
for j, c in enumerate(C):
    for t in c['targets']: A[tix[t], j] = 1
A = A.tocsr()
w = np.array([a.w_hard if t in hard else a.w_easy for t in T])
cost = np.array([c['cost'] for c in C])
nv = n + 2 * m
obj = np.concatenate([cost, -w, np.full(m, a.overlap_pen)])
R1 = hstack([-A, identity(m, format='csr'), csr_matrix((m, m))])
R2 = hstack([A, csr_matrix((m, m)), -identity(m, format='csr')])
R3 = csr_matrix((np.ones(n), (np.zeros(n, int), np.arange(n))), shape=(1, nv))
for N in [int(x) for x in str(a.N).split(',')]:
    lo = np.concatenate([np.full(m, -np.inf), np.full(m, -np.inf), [N]])
    hi = np.concatenate([np.zeros(m), np.ones(m), [N]])
    ub = np.concatenate([np.ones(n), np.ones(m), np.array([0.0 if t in hard else a.easy_slack for t in T])])
    res = milp(obj, constraints=LinearConstraint(vstack([R1, R2, R3]).tocsr(), lo, hi),
               integrality=np.concatenate([np.ones(n), np.zeros(2 * m)]), bounds=Bounds(0, ub),
               options=dict(time_limit=a.time, mip_rel_gap=1e-3))
    if res.x is None:
        print(f'N={N}: no solution ({res.message})'); continue
    y = np.round(res.x[:n]).astype(int); sel = np.where(y > 0)[0]
    cov = {}
    for j in sel:
        for t in C[j]['targets']: cov.setdefault(t, []).append(int(j))
    hc = sum(1 for t in hard if t in cov); ec = sum(1 for t in easy if t in cov)
    ov = sum(1 for t, v in cov.items() if len(v) > 1)
    print(f'N={N} status {res.status}: {len(sel)} routes, covered {len(cov)}/298 (hard {hc}/{len(hard)}, easy {ec}), '
          f'overlaps {ov}, sum J_i {cost[sel].sum():.3f}, J if all placed {cost[sel].sum()+2:.3f}')
    for k, j in enumerate(sel, 1):
        c = C[j]
        print(f'   route {k}: {len(c["targets"]):2d} targets ({len(set(c["targets"])&hard):2d} hard) '
              f'cost {c["cost"]:.3f} tank {c["tank"] or float("nan"):.0f} src {c["src"][0]}')
    d = out / f'N{N}'; d.mkdir(parents=True, exist_ok=True)
    json.dump(dict(N=N, hard=sorted(hard), covered=sorted(cov), sum_cost=float(cost[sel].sum()),
                   uncovered_hard=sorted(t for t in hard if t not in cov),
                   uncovered_easy=sorted(t for t in easy if t not in cov),
                   routes=[dict(targets=C[j]['targets'], cost=C[j]['cost'], tank=C[j]['tank'],
                                src=C[j]['src'][0]) for j in sel],
                   overlap={int(t): v for t, v in cov.items() if len(v) > 1}), open(d / 'partition.json', 'w'), indent=1)
    for k, j in enumerate(sel, 1):
        kind, st = C[j]['src']
        if kind == 'plan':
            json.dump(st, open(d / f'tour_{k:02d}.json', 'w'))
        else:
            if kind == 'npz':
                z = np.load(st); st = {kk: z[kk] for kk in z.files}
            st = {kk: np.asarray(v) for kk, v in st.items()}
            np.savez(d / f'route_{k:02d}.npz', **st)
    print(f'   -> {d}')
