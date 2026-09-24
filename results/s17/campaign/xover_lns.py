"""Stage 17 CAMPAIGN X1: exact-refill LNS for max coverage N<=8 over parents + children.
Drop k of the 8 picks (all combinations), refill k columns EXACTLY (MILP over the columns' traces on the uncovered set,
deduped), accept improvements; k = 2, 3, 4.  Also reports the parents-only LP bound.  Output lns_N8.json"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, itertools
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix, csr_matrix, hstack, identity, vstack
import xover_scan as X

say = lambda s: X.say('[lns] ' + s)
N = 8


def lp_bound(cols):
    R = []; C = []
    for j, c in enumerate(cols):
        for t in c['set']:
            if t in X.TI:
                R.append(X.TI[t]); C.append(j)
    nr = len(cols); A = coo_matrix((np.ones(len(R)), (R, C)), shape=(X.NT, nr)).tocsr()
    Aub = vstack([hstack([-A, identity(X.NT, format='csr')]), csr_matrix(np.concatenate([np.ones(nr), np.zeros(X.NT)])[None, :])]).tocsr()
    r = linprog(np.concatenate([np.zeros(nr), -np.ones(X.NT)]), A_ub=Aub, b_ub=np.concatenate([np.zeros(X.NT), [N]]), bounds=(0, 1), method='highs')
    return -r.fun


def refill(M, U, k, tlim=60):
    """best k columns maximising coverage of the boolean target mask U; returns (value, column indices)."""
    T = M[:, U]; nz = T.any(1)
    idx = np.nonzero(nz)[0]
    if len(idx) == 0:
        return 0, []
    keys = {}
    for j in idx:
        key = T[j].tobytes()
        if key not in keys:
            keys[key] = j
    cand = np.array(list(keys.values())); Tc = T[cand]
    # dominance on traces (vectorised): drop a trace that is a subset of a kept one
    order = np.argsort(-Tc.sum(1)); Tc = Tc[order]; cand = cand[order]
    keep = []; K = np.zeros((0, Tc.shape[1]), bool)
    for i in range(len(cand)):
        if len(keep) and ((Tc[i][None, :] & ~K).sum(1) == 0).any():
            continue
        keep.append(i); K = np.vstack([K, Tc[i][None, :]])
    Tc = Tc[keep]; cand = cand[keep]; nr = len(cand); nu = Tc.shape[1]
    A = csr_matrix(Tc.T.astype(float))
    cons = [LinearConstraint(hstack([A, -identity(nu, format='csr')]).tocsr(), 0, np.inf),
            LinearConstraint(csr_matrix(np.concatenate([np.ones(nr), np.zeros(nu)])[None, :]), -np.inf, k)]
    m = milp(np.concatenate([np.zeros(nr), -np.ones(nu)]), constraints=cons, integrality=np.concatenate([np.ones(nr), np.zeros(nu)]),
             bounds=Bounds(0, 1), options=dict(time_limit=tlim))
    if m.x is None:
        return 0, []
    pick = [int(cand[j]) for j in range(nr) if m.x[j] > 0.5]
    return int(round(-m.fun)), pick


def main():
    t00 = time.time()
    par = X.load_cols(False); res = dict(lp_parents=lp_bound(par)); say(f"LP bound parents only {res['lp_parents']:.3f}")
    cols = X.load_cols(True)
    M = np.zeros((len(cols), X.NT), bool)
    for j, c in enumerate(cols):
        for t in c['set']:
            if t in X.TI:
                M[j, X.TI[t]] = True
    pset = {frozenset(p['set']) for p in json.load(open(X.OUT / 'milp.json'))['parents_N8']['picks']}
    pick = [j for j, c in enumerate(cols) if frozenset(c['set']) in pset]
    for s in pset:
        if not any(frozenset(cols[j]['set']) == s for j in pick):
            pick.append(int(np.argmax([len(c['set'] & s) for c in cols])))
    best = int(M[pick].any(0).sum()); say(f'start {best} ({len(pick)} picks)'); hist = [best]
    for k in (2, 3, 4):
        improved = True
        while improved and time.time() - t00 < 1500:
            improved = False
            for D in itertools.combinations(range(len(pick)), k):
                keep = [pick[i] for i in range(len(pick)) if i not in D]
                c = M[keep].any(0); U = ~c
                v, add = refill(M, U, k)
                if int(c.sum()) + v > best:
                    pick = keep + add; best = int(M[pick].any(0).sum()); hist.append(best); improved = True
                    say(f'k={k}: -> {best}  kinds {[cols[j]["kind"][0] for j in pick]}'); break
        say(f'k={k} done: {best} ({time.time()-t00:.0f} s)')
    res.update(best=best, hist=hist, picks=[dict(kind=cols[j]['kind'], depth=len(cols[j]['set']), tank=round(cols[j]['tank'], 1),
                                                   **({q: cols[j][q] for q in ('a', 'b', 'k', 'L', 'dvj', 'npre', 'nsuf')} if cols[j]['kind'] == 'child' else dict(f=cols[j]['f'])),
                                                   set=sorted(cols[j]['set'])) for j in pick],
               misses=sorted(set(X.ALL) - set().union(*[cols[j]['set'] for j in pick])), sec=round(time.time() - t00))
    json.dump(res, open(X.OUT / 'lns_N8.json', 'w'), indent=1)
    say(f'final {best}, misses {res["misses"]}')


if __name__ == '__main__':
    main()
