"""P1 MILP: closure-aware fleet selection = facility location with lin-priced insertion arcs.

    min  sum_j c_j x_j + sum_a dJ_a z_a
    s.t. y_t <= sum_{j: t in S_j} x_j + sum_{a: ast(a)=t} z_a        (coverage)
         z_a <= x_{host(a)}                                           (arc needs its host)
         sum_{a on j} kg_a z_a <= CAP x_j                              (<= CAP kg of insertions per host)
         sum_j x_j <= N,  sum_t y_t >= K,   x, z, y binary
c_j = cost(honest twin tank), arcs from p1/cols.jsonl with lin dJ <= THR.  ADDITIVE arc prices on the base tank =
OPTIMISTIC (no interaction, convexity ignored); each solution is also re-priced with the convex host cost
cost(tank + sum kg) - cost(tank) ('convex').  Safe dominance pruning: column A is dropped when some B has
S_A U arcs(A) subset of S_B and c_B <= c_A.
usage: p1_milp.py OUT_JSON --N 8,9 --K 0,290,294,296,298 --thr 0.08 [--cap 150] [--tlim 900] [--nproc 2]   (K = 0: max coverage)"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, argparse, pathlib, multiprocessing as mp
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix, csr_matrix, vstack

P1 = pathlib.Path('/Users/mickey/solarsystem/ctoc14/results/s17/redteam/p1')
REACH = [x for x in range(1, 301) if x not in (131, 144)]
TI = {x: i for i, x in enumerate(REACH)}; NT = len(REACH)


def cost(t):
    x = (t - 600.0) / 1400.0
    return 1.0 + x + x * x


def load(thr, nmin=0):
    cols = []
    for l in open(P1 / 'cols.jsonl'):
        r = json.loads(l)
        if not r.get('ok', True) or r['n'] < nmin:
            continue
        arcs = [a for a in r['arcs'] if a['lin_dJ'] <= thr and a['ast'] in TI and a['ast'] not in r['asts']]
        cols.append(dict(key=r['key'], file=r['file'], n=r['n'], tank=r['tank'], c=cost(r['tank']), S=set(r['asts']) & set(TI),
                         arcs=arcs))
    # safe dominance pruning
    cols.sort(key=lambda c: c['c'])
    keep = []
    for A in cols:
        reach = A['S'] | set(a['ast'] for a in A['arcs'])
        if any(B['c'] <= A['c'] and reach <= B['S'] for B in keep):
            continue
        keep.append(A)
    return keep


def solve(cols, N, K, cap, tlim, maxcov=False):
    nC = len(cols)
    arcs = [(j, a) for j, c in enumerate(cols) for a in c['arcs']]
    nA = len(arcs)
    nV = nC + nA + NT
    obj = np.concatenate([[c['c'] for c in cols], [a['lin_dJ'] for _, a in arcs], np.zeros(NT)])
    if maxcov:   # K = 0: maximise coverage first, then min cost (weight 1e-3 < 1 target)
        obj = np.concatenate([1e-3 * obj[:nC + nA], -np.ones(NT)])
    rows, colx, val = [], [], []
    r = 0
    # coverage rows: y_t - sum x_j - sum z_a <= 0
    for t in range(NT):
        rows.append(r); colx.append(nC + nA + t); val.append(1.0); r += 1
    for j, c in enumerate(cols):
        for x in c['S']:
            rows.append(TI[x]); colx.append(j); val.append(-1.0)
    for k, (j, a) in enumerate(arcs):
        rows.append(TI[a['ast']]); colx.append(nC + k); val.append(-1.0)
    lo = [-np.inf] * NT; hi = [0.0] * NT
    # z_a - x_j <= 0
    for k, (j, a) in enumerate(arcs):
        rows += [r, r]; colx += [nC + k, j]; val += [1.0, -1.0]; lo.append(-np.inf); hi.append(0.0); r += 1
    # capacity: sum kg z - cap x_j <= 0
    byj = {}
    for k, (jj, _) in enumerate(arcs):
        byj.setdefault(jj, []).append(k)
    for j, c in enumerate(cols):
        ks = byj.get(j, [])
        if not ks:
            continue
        for k in ks:
            rows.append(r); colx.append(nC + k); val.append(arcs[k][1]['lin_kg'])
        rows.append(r); colx.append(j); val.append(-cap); lo.append(-np.inf); hi.append(0.0); r += 1
    # sum x <= N
    for j in range(nC):
        rows.append(r); colx.append(j); val.append(1.0)
    lo.append(-np.inf); hi.append(N); r += 1
    # sum y >= K
    for t in range(NT):
        rows.append(r); colx.append(nC + nA + t); val.append(1.0)
    lo.append(K); hi.append(np.inf); r += 1
    A = csr_matrix(coo_matrix((val, (rows, colx)), shape=(r, nV)))
    tic = time.time()
    res = milp(obj, constraints=LinearConstraint(A, np.array(lo), np.array(hi)), integrality=np.ones(nV),
               bounds=Bounds(0, 1), options=dict(time_limit=tlim, mip_rel_gap=1e-4, disp=False))
    out = dict(maxcov=maxcov, N=N, K=K, status=int(res.status), message=str(res.message)[:80], sec=round(time.time() - tic),
               dual_bound=None if getattr(res, 'mip_dual_bound', None) is None else round(float(res.mip_dual_bound), 4),
               nC=nC, nA=nA)
    if res.x is None:
        return out
    x = res.x
    pick = [j for j in range(nC) if x[j] > 0.5]
    used = [arcs[k] for k in range(nA) if x[nC + k] > 0.5]
    cov = set().union(*[cols[j]['S'] for j in pick]) | set(a['ast'] for _, a in used)
    add = {}
    for j, a in used:
        add[j] = add.get(j, 0.0) + a['lin_kg']
    base = sum(cols[j]['c'] for j in pick)
    convex = sum(cost(cols[j]['tank'] + add.get(j, 0.0)) for j in pick)
    out.update(obj=round(float(sum(cols[j]['c'] for j in pick) + sum(a['lin_dJ'] for _, a in used)), 4), covered=len(cov), base_sumJi=round(base, 4), arcs_dJ=round(sum(a['lin_dJ'] for _, a in used), 4),
               convex_sumJi=round(convex, 4), misses=sorted(set(REACH) - cov),
               columns=[dict(key=cols[j]['key'], file=cols[j]['file'], n=cols[j]['n'], tank=round(cols[j]['tank'], 1),
                             J_i=round(cols[j]['c'], 4), add_kg=round(add.get(j, 0.0), 1)) for j in pick],
               arcs=[dict(host=cols[j]['key'], ast=a['ast'], t_day=a['t_day'], dist=a['dist'], lin_kg=a['lin_kg'], lin_dJ=a['lin_dJ'],
                          res_km=a['res_km'], untrusted=a['lin_dJ'] > 0.06) for j, a in sorted(used, key=lambda u: -u[1]['lin_dJ'])],
               overlap=sum(len(cols[j]['S']) for j in pick) - len(set().union(*[cols[j]['S'] for j in pick])))
    return out


def _job(args):
    thr, N, K, cap, tlim, nmin = args
    cols = load(thr, nmin)
    r = solve(cols, N, K, cap, tlim, maxcov=(K == 0)); r['thr'] = thr; r['cap'] = cap
    return r


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--N', default='8,9'); ap.add_argument('--K', default='290,294,296,298')
    ap.add_argument('--thr', default='0.08'); ap.add_argument('--cap', default='150'); ap.add_argument('--tlim', type=float, default=900.0)
    ap.add_argument('--nproc', type=int, default=2); ap.add_argument('--nmin', type=int, default=0)
    a = ap.parse_args()
    jobs = [(float(th), int(N), int(K), float(cp), a.tlim, a.nmin) for th in a.thr.split(',') for cp in a.cap.split(',')
            for N in a.N.split(',') for K in a.K.split(',')]
    with mp.get_context('fork').Pool(a.nproc) as pool:
        R = []
        for r in pool.imap_unordered(_job, jobs):
            R.append(r)
            print(json.dumps({k: v for k, v in r.items() if k not in ('columns', 'arcs')}), flush=True)
            json.dump(sorted(R, key=lambda r: (r['thr'], r['cap'], r['N'], r['K'])), open(P1 / a.out, 'w'), indent=1)
    print('all done', flush=True)
