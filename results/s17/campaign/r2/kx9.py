"""R2 CAMPAIGN (same check, part 2): exact k-exchange from s16a at N=9, full 298 cover.
For each set D of k s16a routes: U = targets of D not covered by the kept 9-k; solve min sum cost_j over <= k columns
(parents + stitched children, dominance-pruned) covering U; improvement = cost(D) - optimum.  k = 1, 2, 3.
Child tanks are unregridded estimates (R1 X2: 5/5 settled at -15..+1 kg).  Output r2/kx9.json."""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, itertools, pathlib, multiprocessing as mp
sys.path.insert(0, '/Users/mickey/solarsystem/ctoc14/results/s17/campaign')
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import csr_matrix, coo_matrix
import xover_scan as X
import run_ialns as RI
an = json.load(open('/Users/mickey/solarsystem/ctoc14/results/s17/lns/anatomy.json'))
BASE = {f'r{i}': (frozenset(a for a in an[f'r{i}']['asts'] if a in X.TI), RI.cost(an[f'r{i}']['tank'])) for i in range(1, 10)}
COLS = X.load_cols(True)
SETS = [c['set'] & set(X.TI) for c in COLS]; COST = np.array([c['cost'] for c in COLS])


def solve(D):
    keep = [r for r in BASE if r not in D]
    covk = set().union(*[BASE[r][0] for r in keep])
    U = sorted(set().union(*[BASE[r][0] for r in D]) - covk)
    ui = {t: i for i, t in enumerate(U)}
    J = [j for j, s in enumerate(SETS) if s & set(U)]
    R = []; C = []
    for q, j in enumerate(J):
        for t in SETS[j]:
            if t in ui:
                R.append(ui[t]); C.append(q)
    A = coo_matrix((np.ones(len(R)), (R, C)), shape=(len(U), len(J))).tocsr()
    cons = [LinearConstraint(A, 1, np.inf), LinearConstraint(csr_matrix(np.ones((1, len(J)))), -np.inf, len(D))]
    r = milp(COST[J], constraints=cons, integrality=np.ones(len(J)), bounds=Bounds(0, 1),
             options=dict(time_limit=90, mip_rel_gap=1e-6))
    base = sum(BASE[x][1] for x in D)
    if r.x is None:
        return dict(D=list(D), U=len(U), status=str(r.message)[:40])
    pick = [J[q] for q in range(len(J)) if r.x[q] > 0.5]
    return dict(D=list(D), U=len(U), base=round(base, 4), opt=round(float(COST[pick].sum()), 4),
                gain=round(base - float(COST[pick].sum()), 4), kinds=[COLS[j]['kind'] for j in pick],
                n=[len(SETS[j]) for j in pick], tanks=[round(COLS[j]['tank'], 1) for j in pick], msg=str(r.message)[:30])


if __name__ == '__main__':
    jobs = [D for k in (1, 2, 3) for D in itertools.combinations(sorted(BASE), k)]
    print(len(COLS), 'cols;', len(jobs), 'exchanges', flush=True)
    out = []; t0 = time.time()
    with mp.get_context('fork').Pool(2) as P:
        for res in P.imap_unordered(solve, jobs):
            out.append(res)
            if res.get('gain', 0) > 0.005:
                print('IMPROVE', json.dumps(res), flush=True)
    out.sort(key=lambda r: -r.get('gain', -9))
    json.dump(out, open(X.OUT / 'r2' / 'kx9.json', 'w'), indent=1)
    print(f'done {time.time() - t0:.0f} s; best', json.dumps(out[:3]), flush=True)
