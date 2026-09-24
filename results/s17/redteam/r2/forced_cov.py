"""R2 check (REDTEAM): is the 8-craft hard core SPECIFIC or does it ROTATE?
For each target X of the 285-optimum's miss set, solve the P1 max-coverage MILP (N<=8, lin arcs <= thr, cap 150 kg)
with y_X forced to 1.  Rotation = coverage stays ~285 but a DIFFERENT target drops out.  Reuses p1_milp.load (read-only).
usage: forced_cov.py OUT --thr 0.08 --tlim 240 --nproc 2 [--targets 17,36,...]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, argparse, pathlib, multiprocessing as mp
sys.path.insert(0, '/Users/mickey/solarsystem/ctoc14/results/s17/redteam/p1')
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix, csr_matrix
import p1_milp as PM
R2 = pathlib.Path('/Users/mickey/solarsystem/ctoc14/results/s17/redteam/r2')
MISS = [17, 36, 75, 108, 118, 180, 216, 219, 220, 243, 247, 258, 288]
COLS = None


def solve(X, N, tlim):
    cols = COLS; TI = PM.TI; NT = PM.NT
    nC = len(cols); arcs = [(j, a) for j, c in enumerate(cols) for a in c['arcs']]; nA = len(arcs); nV = nC + nA + NT
    obj = np.concatenate([1e-3 * np.array([c['c'] for c in cols] + [a['lin_dJ'] for _, a in arcs]), -np.ones(NT)])
    rows, colx, val, lo, hi = [], [], [], [], []
    r = 0
    for t in range(NT):
        rows.append(r); colx.append(nC + nA + t); val.append(1.0); r += 1
    for j, c in enumerate(cols):
        for x in c['S']:
            rows.append(TI[x]); colx.append(j); val.append(-1.0)
    for k, (j, a) in enumerate(arcs):
        rows.append(TI[a['ast']]); colx.append(nC + k); val.append(-1.0)
    lo += [-np.inf] * NT; hi += [0.0] * NT
    for k, (j, a) in enumerate(arcs):
        rows += [r, r]; colx += [nC + k, j]; val += [1.0, -1.0]; lo.append(-np.inf); hi.append(0.0); r += 1
    byj = {}
    for k, (jj, _) in enumerate(arcs):
        byj.setdefault(jj, []).append(k)
    for j in range(nC):
        ks = byj.get(j, [])
        if not ks:
            continue
        for k in ks:
            rows.append(r); colx.append(nC + k); val.append(arcs[k][1]['lin_kg'])
        rows.append(r); colx.append(j); val.append(-150.0); lo.append(-np.inf); hi.append(0.0); r += 1
    for j in range(nC):
        rows.append(r); colx.append(j); val.append(1.0)
    lo.append(-np.inf); hi.append(N); r += 1
    A = csr_matrix(coo_matrix((val, (rows, colx)), shape=(r, nV)))
    lb = np.zeros(nV); ub = np.ones(nV)
    if X is not None:
        lb[nC + nA + TI[X]] = 1.0
    tic = time.time()
    res = milp(obj, constraints=LinearConstraint(A, np.array(lo), np.array(hi)), integrality=np.ones(nV),
               bounds=Bounds(lb, ub), options=dict(time_limit=tlim, mip_rel_gap=1e-4, disp=False))
    out = dict(X=X, N=N, status=int(res.status), sec=round(time.time() - tic),
               dual_cov=None if getattr(res, 'mip_dual_bound', None) is None else round(-float(res.mip_dual_bound), 2))
    if res.x is None:
        return out
    x = res.x
    pick = [j for j in range(nC) if x[j] > 0.5]; used = [arcs[k] for k in range(nA) if x[nC + k] > 0.5]
    cov = set().union(*[cols[j]['S'] for j in pick]) | set(a['ast'] for _, a in used)
    add = {}
    for j, a in used:
        add[j] = add.get(j, 0.0) + a['lin_kg']
    out.update(covered=len(cov), misses=sorted(set(PM.REACH) - cov),
               sumJi_add=round(sum(cols[j]['c'] for j in pick) + sum(a['lin_dJ'] for _, a in used), 4),
               cols=[(cols[j]['key'], cols[j]['n'], round(cols[j]['tank'], 1), round(add.get(j, 0.0), 1)) for j in pick],
               X_by=[cols[j]['key'] for j in pick if X in cols[j]['S']] + [f"arc:{cols[j]['key']}:{a['lin_kg']:.0f}kg" for j, a in used if a['ast'] == X])
    return out


def _job(a):
    return solve(*a)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--thr', type=float, default=0.08)
    ap.add_argument('--tlim', type=float, default=240.0); ap.add_argument('--nproc', type=int, default=2)
    ap.add_argument('--targets', default=''); ap.add_argument('--N', type=int, default=8); ap.add_argument('--nmin', type=int, default=0); ap.add_argument('--first', type=int, default=0)
    a = ap.parse_args()
    COLS = PM.load(a.thr, a.nmin)
    if a.first:
        keep = set(json.loads(l)['key'] for i, l in zip(range(a.first), open(PM.P1 / 'cols.jsonl')))
        COLS = [c for c in COLS if c['key'] in keep]
    print(f'{len(COLS)} columns after pruning, {sum(len(c["arcs"]) for c in COLS)} arcs', flush=True)
    T = [int(x) for x in a.targets.split(',')] if a.targets else [None] + MISS
    R = []
    with mp.get_context('fork').Pool(a.nproc) as pool:
        for r in pool.imap_unordered(_job, [(X, a.N, a.tlim) for X in T]):
            R.append(r); print(json.dumps({k: v for k, v in r.items() if k != 'cols'}), flush=True)
            json.dump(R, open(R2 / a.out, 'w'), indent=1)
    print('all done', flush=True)
