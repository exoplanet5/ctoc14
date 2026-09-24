"""Stage 17 CAMPAIGN X1 solve: max coverage N<=8 over parents + crossover children, done properly.
1. LP bound (HiGHS); 2. reduced-cost fixing against the parents-only incumbent (286); 3. MILP on the survivors;
4. independent 1-swap / 2-swap local search from the parents-only pick over ALL columns.
Output results/s17/campaign/solve.json"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, itertools
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix, csr_matrix, hstack, identity, vstack
import xover_scan as X

OUT = X.OUT; N = int(sys.argv[1]) if len(sys.argv) > 1 else 8; MODE = sys.argv[2] if len(sys.argv) > 2 else 'all'


def mats(cols):
    R = []; C = []
    for j, c in enumerate(cols):
        for t in c['set']:
            if t in X.TI:
                R.append(X.TI[t]); C.append(j)
    return coo_matrix((np.ones(len(R)), (R, C)), shape=(X.NT, len(cols))).tocsr()


def cov(cols, pick):
    return len(set().union(*[cols[j]['set'] for j in pick]) & set(X.ALL)) if pick else 0


def local(cols, pick, say, tlim=900):
    t0 = time.time(); sets = [np.array([X.TI[t] for t in c['set'] if t in X.TI], int) for c in cols]
    M = np.zeros((len(cols), X.NT), bool)
    for j, s in enumerate(sets):
        M[j, s] = True
    pick = list(pick); best = cov(cols, pick); improved = True
    while improved and time.time() - t0 < tlim:
        improved = False
        # 1-swap and 2-swap: drop set D, add greedily |D| columns maximising new coverage
        for k in (1, 2):
            for D in itertools.combinations(range(len(pick)), k):
                keep = [pick[i] for i in range(len(pick)) if i not in D]
                c = M[keep].any(0) if keep else np.zeros(X.NT, bool)
                add = []
                for _ in range(k):
                    gain = (M & ~c).sum(1); j = int(gain.argmax()); add.append(j); c = c | M[j]
                v = int(c.sum())
                if v > best:
                    pick = keep + add; best = v; improved = True
                    say(f'local: {k}-swap -> {best}'); break
            if improved:
                break
    return pick, best


def main():
    say = lambda s: X.say('[solve] ' + s)
    cols = X.load_cols(True); nr = len(cols)
    par = X.load_cols(False)
    Mi = json.load(open(OUT / 'milp.json'))[f'parents_N{N}']
    pset = {frozenset(p['set']) for p in Mi['picks']}
    inc = [j for j, c in enumerate(cols) if frozenset(sorted(c['set'])) in pset]
    # a parent pick may have been dominated away by a child: then take the dominating column
    for s in pset:
        if not any(frozenset(cols[j]['set']) == s for j in inc):
            inc.append(max(range(nr), key=lambda j: len(cols[j]['set'] & s) - 1e-3 * len(cols[j]['set'])))
    L0 = cov(cols, inc); say(f'N<={N}: {nr} columns, incumbent (parents) {L0} with {len(inc)} picks')
    res = dict(N=N, n_cols=nr, incumbent=L0)
    A = mats(cols); cost = np.array([c['cost'] for c in cols]); eps = 1e-4
    if MODE in ('all', 'lp'):
        t0 = time.time()
        Aub = vstack([hstack([-A, identity(X.NT, format='csr')]), csr_matrix(np.concatenate([np.ones(nr), np.zeros(X.NT)])[None, :])]).tocsr()
        bub = np.concatenate([np.zeros(X.NT), [N]])
        r = linprog(np.concatenate([eps * cost, -np.ones(X.NT)]), A_ub=Aub, b_ub=bub, bounds=(0, 1), method='highs')
        U = -r.fun; rc = r.lower.marginals[:nr]
        res['lp_bound'] = float(U); say(f'LP bound {U:.3f} ({time.time()-t0:.0f} s)')
        keep = sorted(set(np.nonzero(rc <= U - (L0 + 1) + 1e-6)[0].tolist()) | set(inc))
        res['rc_kept'] = len(keep); say(f'reduced-cost fixing keeps {len(keep)} columns (rc <= {U - L0 - 1:.3f})')
        sub = [cols[j] for j in keep]; As = mats(sub); nsub = len(sub)
        Mc = hstack([As, -identity(X.NT, format='csr')]).tocsr()
        cons = [LinearConstraint(Mc, 0, np.inf), LinearConstraint(csr_matrix(np.concatenate([np.ones(nsub), np.zeros(X.NT)])[None, :]), -np.inf, N),
                LinearConstraint(csr_matrix(np.concatenate([np.zeros(nsub), np.ones(X.NT)])[None, :]), L0, np.inf)]
        t0 = time.time()
        m = milp(np.concatenate([eps * np.array([c['cost'] for c in sub]), -np.ones(X.NT)]), constraints=cons,
                 integrality=np.concatenate([np.ones(nsub), np.zeros(X.NT)]), bounds=Bounds(0, 1),
                 options=dict(time_limit=900, mip_rel_gap=1e-7))
        if m.x is not None:
            pick = [j for j in range(nsub) if m.x[j] > 0.5]
            res['milp'] = dict(covered=cov(sub, pick), bound=float(-m.mip_dual_bound) if m.mip_dual_bound is not None else None,
                               msg=m.message[:60], sec=round(time.time() - t0), sumJi=float(sum(sub[j]['cost'] for j in pick)),
                               picks=[dict(kind=sub[j]['kind'], depth=len(sub[j]['set']), tank=round(sub[j]['tank'], 1),
                                           **({k: sub[j][k] for k in ('a', 'b', 'k', 'L', 'dvj', 'npre', 'nsuf')} if sub[j]['kind'] == 'child' else dict(f=sub[j]['f'])),
                                           set=sorted(sub[j]['set'])) for j in pick],
                               misses=sorted(set(X.ALL) - set().union(*[sub[j]['set'] for j in pick])))
            say(f"MILP restricted: covered {res['milp']['covered']} bound {res['milp']['bound']} ({res['milp']['msg']}, {res['milp']['sec']} s)")
        else:
            res['milp'] = dict(fail=m.message[:80]); say('MILP fail ' + m.message[:80])
        json.dump(res, open(OUT / f'solve_N{N}.json', 'w'), indent=1)
    if MODE in ('all', 'local'):
        pk, b = local(cols, inc, say)
        res['local'] = dict(covered=b, kinds=[cols[j]['kind'] for j in pk], depths=[len(cols[j]['set']) for j in pk])
        say(f'local search from parents incumbent: {b}')
        json.dump(res, open(OUT / f'solve_N{N}_{MODE}.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
