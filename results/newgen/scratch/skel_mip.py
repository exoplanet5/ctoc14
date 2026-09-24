"""Skeleton selection for the 8-craft fleet: choose N twin-costed columns (with saved twin states) that cover the hard
targets, nearly disjointly, at low cost. MIP: max sum_t w_t c_t - sum_j cost_j y_j, c_t <= sum_{j in t} y_j, c_t <= 1,
sum_{j in t} y_j <= 1 + slack for easy targets (overlap allowed only on easy targets, penalised), sum y = N.
Usage: skel_mip.py out_dir [--N 8] [--w-hard 0.3] [--w-easy 0.02] [--overlap-pen 0.05] [--time 120]"""
import sys, json, glob, pathlib, argparse, numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix, csr_matrix, vstack, hstack, identity
from ctoc14.constants import cost_sc

ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--N', type=int, default=8)
ap.add_argument('--w-hard', type=float, default=0.3); ap.add_argument('--w-easy', type=float, default=0.02)
ap.add_argument('--overlap-pen', type=float, default=0.05); ap.add_argument('--time', type=float, default=120.0)
ap.add_argument('--max-tank', type=float, default=1000.0); ap.add_argument('--min-hard', type=int, default=0)
a = ap.parse_args()
out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)

# hardness
cnt = np.zeros(301, int)
for line in open(ROOT / 'results/newgen/recost35.jsonl'):
    r = json.loads(line)
    if r['ok']:
        for t in r['targets']: cnt[t] += 1
hard = set(t for t in range(1, 301) if cnt[t] <= 5) - {131, 144}
easy = set(range(1, 301)) - hard - {131, 144}
print(f'{len(hard)} hard targets, {len(easy)} easy')

# columns with states
cols = {}   # key: frozenset targets -> dict
def add(targets, tank, state_src):
    k = frozenset(int(t) for t in targets)
    if tank > a.max_tank: return
    if k not in cols or tank < cols[k]['tank']:
        cols[k] = dict(targets=sorted(k), tank=float(tank), src=state_src)
for line in open(ROOT / 'results/newgen/recost35.jsonl'):
    r = json.loads(line)
    if r['ok']: add(r['targets'], r['tank'], ('inline', r['st']))
for f in glob.glob(str(ROOT / 'results/colgen/n8a/twin_run1.jsonl')) + glob.glob(str(ROOT / 'results/colgen/n8a/twin.jsonl')):
    for line in open(f):
        r = json.loads(line)
        if r.get('tank') and (ROOT / f'results/colgen/n8a/twin_states/col_{r["col"]}.npz').exists():
            add(r['targets'], r['tank'], ('npz', str(ROOT / f'results/colgen/n8a/twin_states/col_{r["col"]}.npz')))
C = list(cols.values())
C = [c for c in C if len(set(c['targets']) & hard) >= a.min_hard]
n = len(C); print(f'{n} distinct twin-costed columns with states')
T = sorted(hard | easy); tix = {t: i for i, t in enumerate(T)}; m = len(T)
A = lil_matrix((m, n))
for j, c in enumerate(C):
    for t in c['targets']:
        if t in tix: A[tix[t], j] = 1
A = A.tocsr()
w = np.array([a.w_hard if t in hard else a.w_easy for t in T])
cost = np.array([cost_sc(c['tank']) for c in C])
# variables: y (n) binary, c (m) in [0,1], s (m) >= 0 overlap slack (easy only)
nv = n + 2 * m
obj = np.concatenate([cost, -w, np.full(m, a.overlap_pen)])
# c_t - sum_j A_tj y_j <= 0
R1 = hstack([-A, identity(m, format='csr'), csr_matrix((m, m))])
# sum_j A_tj y_j - s_t <= 1  (hard: s fixed 0 -> disjoint)
R2 = hstack([A, csr_matrix((m, m)), -identity(m, format='csr')])
R3 = csr_matrix((np.ones(n), (np.zeros(n, int), np.arange(n))), shape=(1, nv))
Acon = vstack([R1, R2, R3]).tocsr()
lo = np.concatenate([np.full(m, -np.inf), np.full(m, -np.inf), [a.N]]); hi = np.concatenate([np.zeros(m), np.ones(m), [a.N]])
ub = np.concatenate([np.ones(n), np.ones(m), np.array([0.0 if t in hard else 3.0 for t in T])])
res = milp(obj, constraints=LinearConstraint(Acon, lo, hi), integrality=np.concatenate([np.ones(n), np.zeros(2 * m)]),
           bounds=Bounds(0, ub), options=dict(time_limit=a.time, mip_rel_gap=1e-3))
y = np.round(res.x[:n]).astype(int); sel = np.where(y > 0)[0]
cov = {}
for j in sel:
    for t in C[j]['targets']: cov.setdefault(t, []).append(int(j))
hc = sum(1 for t in hard if t in cov); ec = sum(1 for t in easy if t in cov); ov = sum(1 for t, v in cov.items() if len(v) > 1)
print(f'MIP status {res.status}: {len(sel)} skeletons, cover hard {hc}/{len(hard)}, easy {ec}/{len(easy)}, overlapping targets {ov}, sum J_i {cost[sel].sum():.3f}')
for k, j in enumerate(sel, 1):
    c = C[j]; print(f'  skel {k}: {len(c["targets"])} targets ({len(set(c["targets"]) & hard)} hard), tank {c["tank"]:.0f}, J_i {cost[j]:.3f}, src {c["src"][0]}')
    st = c['src'][1]
    if c['src'][0] == 'npz':
        z = np.load(st); st = {kk: z[kk] for kk in z.files}
    else:
        st = {kk: np.asarray(v) for kk, v in st.items()}
    np.savez(out / f'route_{k:02d}.npz', **st)
json.dump(dict(hard=sorted(hard), covered=sorted(cov), uncovered_hard=sorted(t for t in hard if t not in cov),
               uncovered_easy=sorted(t for t in easy if t not in cov), overlap={int(t): v for t, v in cov.items() if len(v) > 1}),
          open(out / 'skeletons.json', 'w'), indent=1)
print('uncovered hard:', sorted(t for t in hard if t not in cov))
