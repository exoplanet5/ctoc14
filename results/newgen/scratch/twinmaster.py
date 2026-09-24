"""Master problem over twin-costed columns (gate F follow-up): set cover with a craft cap, columns = recosted pool
columns (recost_pool.py output) + the flown impulsive fleet's routes. Cost of a column = J_i(tank) from the twin.
Usage: twinmaster.py recost.jsonl [more.jsonl ...] [--fleet results/newgen/ifleet_t10d] [--N 8,9,10] [--mip 120]"""
import sys, json, glob, pathlib, argparse, numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
from ctoc14.colgen import CostModel, ColumnStore, RMP, TARGETS
from ctoc14.constants import cost_sc

ap = argparse.ArgumentParser(); ap.add_argument('files', nargs='+'); ap.add_argument('--fleet', default='results/newgen/ifleet_t10d')
ap.add_argument('--N', default='8,9,10'); ap.add_argument('--mip', type=float, default=120.0); ap.add_argument('--max-tank', type=float, default=1900.0)
a = ap.parse_args()

store = ColumnStore(CostModel())
def add(targets, tank, src, tag):
    j, new = store.add(targets, 0.0, 1000.0, 0.0, ('inline', src), tag)
    if j >= 0:
        c = float(cost_sc(tank))
        if new or c < store.cost[j]:
            store.cost[j] = c; store.fuel[j] = float(tank)      # fuel slot holds the tank
    return j

n_ok = 0; n_fail = 0
for f in a.files:
    for line in open(f):
        r = json.loads(line)
        if not r.get('ok') or r['tank'] > a.max_tank:
            n_fail += 1; continue
        add(r['targets'], r['tank'], dict(col=r['col'], tank=r['tank'], dv=r['dv']), r.get('tag', 'pool')); n_ok += 1
fleet_ids = []
if a.fleet:
    for f in sorted(glob.glob(str(ROOT / a.fleet / 'route_*.npz'))):
        z = np.load(f); tank = float(json.load(open(ROOT / a.fleet / 'fleet.json'))['routes'][pathlib.Path(f).stem[6:]]['tank'])
        fleet_ids.append(add([int(x) for x in z['asts']], tank, dict(route=pathlib.Path(f).stem), 'flown'))
n = np.array(store.n); c = np.array(store.cost)
print(f'{n_ok} recosted columns ok ({n_fail} failed/skipped) + {len(fleet_ids)} flown routes -> {len(store)} distinct columns')
cov = np.zeros(300, bool)
for t in store.tidx: cov[t] = True
rows = TARGETS - 1
print(f'union coverage {cov[rows].sum()} / 298; missing: {[int(i) + 1 for i in rows if not cov[i]]}')
for lo in (30, 33, 35, 37):
    m = n >= lo
    if m.any():
        dvpf = np.array([store.src[j][1].get('dv', np.nan) / n[j] for j in np.where(m)[0]])
        print(f'  >= {lo} targets: {m.sum():5d} cols, J_i min {c[m].min():.3f} p25 {np.percentile(c[m], 25):.3f} median {np.median(c[m]):.3f}; '
              f'dv/flyby min {np.nanmin(dvpf):.3f} median {np.nanmedian(dvpf):.3f}')
# how often is each target covered by a cheap column?
C = store.csr(); cnt = np.asarray(C.sum(axis=0)).ravel()
print(f'targets covered by <= 3 columns: {int((cnt[rows] <= 3).sum())}; by 0: {int((cnt[rows] == 0).sum())}')
rmp = RMP(store)
for N in [int(x) for x in a.N.split(',')]:
    rmp.active = np.arange(len(store))
    res = rmp.solve(rows, N, verbose=None)
    frac = int(np.sum((res.y > 1e-6) & (res.y < 1 - 1e-6)))
    print(f'N={N}: LP {res.value:.4f} (sum y {res.y.sum():.2f}, misses {res.w.sum():.2f}, mu {res.mu:.3f}, fractional {frac}, '
          f'#pi>0.99 {(res.pi > 0.99).sum()}, pi median {np.median(res.pi[rows]):.3f})')
    Jm, sel, missed = rmp.mip(rows, N, time_limit=a.mip)
    print(f'      MIP {Jm:.4f}: {len(sel)} craft, missed {len(missed)}; columns: ' +
          ' '.join(f'{store.tag[j][:5]}:{n[j]}@{store.fuel[j]:.0f}' for j in sorted(sel, key=lambda j: -n[j])))
    json.dump(dict(N=N, lp=res.value, mip=Jm, sel=[dict(tag=store.tag[j], n=int(n[j]), tank=store.fuel[j], src=store.src[j][1], targets=store.targets(j)) for j in sel],
                   missed=missed, pi={int(i) + 1: float(res.pi[i]) for i in rows}), open(ROOT / f'results/newgen/scratch/twinmaster_N{N}.json', 'w'))
