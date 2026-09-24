"""Stage 15 [M']: exact N<=9 selection that finishes.  s13_master's plain HiGHS MILP stopped at its 300 s limit on
the round-1 pool (1245 columns) with the OLD incumbent (J 21.70) against an LP bound of 13.98, so this selector adds:

  1. load  = s13_master.load (miss-checked columns, dedupe by target set, lightest first), then DOMINANCE pruning
             (drop c2 when some c1 covers a superset of c2's targets at <= c2's cost);
  2. LP    = the master LP (sum y <= N, miss variables z at cost 1);
  3. DIVE  = LP-guided rounding: fix the column with the largest fractional y (ties: more new targets), re-solve the
             LP over what is left, repeat until N columns are fixed or everything is covered; several dives with
             different tie-breaks give an incumbent UB;
  4. LOCAL = 1-swap improvement of the incumbent (drop one column, add the best column for the uncovered set);
  5. FIX   = reduced-cost fixing: a column whose LP reduced cost exceeds UB - LB can be in no better solution;
  6. MILP  = HiGHS on the reduced column set, time limit --tlim, with the incumbent known.
J = sum cost(tank) + #uncovered (of the 298) + 2.   Outputs OUT/N<N>/fleet (RI.IFleet), OUT/select.json, OUT/log.txt.

usage: s15_select.py OUT PAT [PAT ...] [--N 9] [--tlim 1800] [--nproc 4] [--dives 6]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, glob, time, fnmatch, pathlib, argparse, collections, multiprocessing as mp
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
_CWD0 = os.getcwd()
os.chdir(ROOT)
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix, csr_matrix, vstack
import run_ialns as RI
import s13_master as SM

ALL = SM.ALL
TI = {t: i for i, t in enumerate(ALL)}
NT = len(ALL)


def load_columns(pats, nproc, say, exclude=('*/relocated/*',), max_tank=1400.0):
    files = sorted(set(f for p in pats for f in glob.glob(p, recursive=True)
                       if not any(fnmatch.fnmatch(f, e) for e in exclude)))
    with mp.get_context('fork').Pool(nproc) as pool:
        res = pool.map(SM.load, [(f, max_tank) for f in files], chunksize=8)
    why = collections.Counter(r for _, r in res)
    best = {}
    for c, _ in res:
        if c is None:
            continue
        k = frozenset(c['asts'])
        if k not in best or c['tank'] < best[k]['tank']:
            best[k] = c
    cols = sorted(best.values(), key=lambda x: (-len(x['asts']), x['tank']))
    for c in cols:
        c['cost'] = RI.cost(c['tank']); c['set'] = frozenset(c['asts'])
    # dominance: c2 dropped if some c1 has set >= set2 and cost <= cost2 (c1 != c2)
    keep = []
    by_t = collections.defaultdict(list)
    for i, c in enumerate(cols):
        for t in c['asts']:
            by_t[t].append(i)
    dom = 0
    for j, c2 in enumerate(cols):
        t0 = min(c2['asts'], key=lambda t: len(by_t[t]))
        if any(i != j and cols[i]['cost'] <= c2['cost'] + 1e-12 and len(cols[i]['set']) >= len(c2['set'])
               and c2['set'] <= cols[i]['set'] and (cols[i]['cost'] < c2['cost'] - 1e-12 or len(cols[i]['set']) > len(c2['set']) or i < j)
               for i in by_t[t0]):
            dom += 1; continue
        keep.append(c2)
    say(f'load: {len(files)} files ({dict(why)}), {len(cols)} distinct, {dom} dominated -> {len(keep)} columns')
    return keep


def matrices(cols):
    R = []; C = []
    for j, c in enumerate(cols):
        for t in c['asts']:
            R.append(TI[t]); C.append(j)
    A = coo_matrix((np.ones(len(R)), (R, C)), shape=(NT, len(cols))).tocsr()
    return A, np.array([c['cost'] for c in cols])


def lp(cols, N, fixed=(), banned=(), covered_fixed=None):
    """LP over y (columns) and z (misses): min c y + 1 z, A y + z >= 1, sum y <= N, fixed y = 1, banned y = 0.
    Returns (value incl. +2, y, reduced costs of y) or None."""
    A, c = matrices(cols); nr = len(cols)
    Aub = vstack([-_hcat(A), csr_matrix(np.concatenate([np.ones(nr), np.zeros(NT)])[None, :])]).tocsr()
    bub = np.concatenate([-np.ones(NT), [N]])
    lo = np.zeros(nr + NT); hi = np.ones(nr + NT)
    for j in fixed:
        lo[j] = 1.0
    for j in banned:
        hi[j] = 0.0
    obj = np.concatenate([c, np.ones(NT)])
    r = linprog(obj, A_ub=Aub, b_ub=bub, bounds=list(zip(lo, hi)), method='highs')
    if r.status != 0:
        return None
    y = r.x[:nr]
    dual = -r.ineqlin.marginals[:NT]                 # >= 0 cover duals
    mu = -r.ineqlin.marginals[NT]                    # >= 0 card dual
    rc = c - A.T @ dual + mu
    return float(r.fun) + 2.0, y, rc


def _hcat(A):
    from scipy.sparse import hstack, identity
    return hstack([A, identity(NT, format='csr')]).tocsr()


def evaluate(cols, pick):
    cov = set()
    for j in pick:
        cov |= cols[j]['set']
    sj = sum(cols[j]['cost'] for j in pick)
    return sj + (NT - len(cov)) + 2.0, len(cov), sj


def dive(cols, N, say, tie=0, max_lp=40):
    fixed = []; cov = set(); nlp = 0
    while len(fixed) < N and len(cov) < NT and nlp < max_lp:
        r = lp(cols, N, fixed=fixed); nlp += 1
        if r is None:
            break
        val, y, rc = r
        cand = [j for j in range(len(cols)) if j not in fixed and y[j] > 1e-6]
        if not cand:
            break
        def key(j):
            new = len(cols[j]['set'] - cov)
            return (y[j], new) if tie == 0 else ((new, y[j]) if tie == 1 else (y[j] * new, y[j]))
        jb = max(cand, key=key)
        fixed.append(jb); cov |= cols[jb]['set']
    # fill remaining slots greedily if coverage incomplete
    while len(fixed) < N and len(cov) < NT:
        jb = max((j for j in range(len(cols)) if j not in fixed),
                 key=lambda j: (len(cols[j]['set'] - cov) - cols[j]['cost']))
        if len(cols[jb]['set'] - cov) - cols[jb]['cost'] <= 0:
            break
        fixed.append(jb); cov |= cols[jb]['set']
    return fixed


def local_swap(cols, pick, N, say, rounds=4):
    best = list(pick); bJ = evaluate(cols, best)[0]
    for it in range(rounds):
        improved = False
        for i in range(len(best)):
            rest = best[:i] + best[i + 1:]
            covr = set()
            for j in rest:
                covr |= cols[j]['set']
            # best replacement: min cost + uncovered
            bj = None; bv = None
            for j in range(len(cols)):
                if j in rest:
                    continue
                v = cols[j]['cost'] + (NT - len(covr | cols[j]['set']))
                if bv is None or v < bv:
                    bv = v; bj = j
            cand = rest + [bj]
            J = evaluate(cols, cand)[0]
            if J < bJ - 1e-9:
                best, bJ, improved = cand, J, True
        # also try adding a column if fewer than N
        if len(best) < N:
            covb = set()
            for j in best:
                covb |= cols[j]['set']
            bj = min((j for j in range(len(cols)) if j not in best),
                     key=lambda j: cols[j]['cost'] - len(cols[j]['set'] - covb))
            if cols[bj]['cost'] - len(cols[bj]['set'] - covb) < 0:
                best = best + [bj]; bJ = evaluate(cols, best)[0]; improved = True
        if not improved:
            break
    return best, bJ


def exact(cols, N, tlim, say, gap=1e-4):
    A, c = matrices(cols); nr = len(cols)
    M = _hcat(A)
    cons = [LinearConstraint(M, np.ones(NT), np.full(NT, np.inf)),
            LinearConstraint(csr_matrix(np.concatenate([np.ones(nr), np.zeros(NT)])[None, :]), -np.inf, N)]
    r = milp(np.concatenate([c, np.ones(NT)]), constraints=cons, integrality=np.concatenate([np.ones(nr), np.zeros(NT)]),
             bounds=Bounds(0, 1), options=dict(time_limit=tlim, mip_rel_gap=gap, disp=False))
    if r.x is None:
        say(f'  MILP: {r.message}'); return None, r.message
    pick = [j for j in range(nr) if r.x[j] > 0.5]
    return pick, r.message


def frontier(cols, N, Ks, tlim, say):
    """Min sum J_i with sum y <= N and coverage >= K (c_t <= sum_{r: t in r} y_r, sum c_t >= K), per K."""
    A, c = matrices(cols); nr = len(cols)
    from scipy.sparse import hstack, identity
    Mcov = hstack([A, -identity(NT, format='csr')]).tocsr()          # sum y - c_t >= 0
    out = {}
    for K in Ks:
        cons = [LinearConstraint(Mcov, np.zeros(NT), np.full(NT, np.inf)),
                LinearConstraint(csr_matrix(np.concatenate([np.ones(nr), np.zeros(NT)])[None, :]), -np.inf, N),
                LinearConstraint(csr_matrix(np.concatenate([np.zeros(nr), np.ones(NT)])[None, :]), K, np.inf)]
        r = milp(np.concatenate([c, np.zeros(NT)]), constraints=cons, integrality=np.ones(nr + NT), bounds=Bounds(0, 1),
                 options=dict(time_limit=tlim, mip_rel_gap=1e-4))
        if r.x is None:
            say(f'  frontier K>={K}: {r.message[:60]}'); out[K] = None; continue
        pick = [j for j in range(nr) if r.x[j] > 0.5]
        J, cv, sj = evaluate(cols, pick)
        say(f'  frontier K>={K}: {len(pick)} craft {cv} covered sum J_i {sj:.4f} ({r.message[:30]}); '
            f'misses {sorted(set(ALL) - set().union(*[cols[j]["set"] for j in pick]))}; '
            f'depths {sorted([len(cols[j]["asts"]) for j in pick], reverse=True)}')
        out[K] = dict(covered=cv, sumJi=round(sj, 4), picks=[cols[j]['f'] for j in pick],
                      misses=sorted(set(ALL) - set().union(*[cols[j]['set'] for j in pick])))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('pats', nargs='+')
    ap.add_argument('--N', type=int, default=9); ap.add_argument('--tlim', type=float, default=1800.0)
    ap.add_argument('--nproc', type=int, default=4); ap.add_argument('--dives', type=int, default=3)
    ap.add_argument('--frontier', default='', help='also min sum J_i at coverage >= K for these K, e.g. 291,292,294,296,298')
    ap.add_argument('--frontier-only', action='store_true'); ap.add_argument('--ftlim', type=float, default=600.0)
    a = ap.parse_args()
    out = pathlib.Path(a.out) if os.path.isabs(a.out) else pathlib.Path(_CWD0) / a.out
    out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt'); tic = time.time()
    RI.eph()
    cols = load_columns(a.pats, a.nproc, say)
    if a.frontier:
        fr = frontier(cols, a.N, [int(x) for x in a.frontier.split(',')], a.ftlim, say)
        json.dump({str(k): v for k, v in fr.items()}, open(out / 'frontier.json', 'w'), indent=1)
        if a.frontier_only:
            return
    r = lp(cols, a.N)
    LB, y, rc = r
    say(f'LP (sum y <= {a.N}): {LB:.4f}; {int((y > 1e-6).sum())} fractional columns')
    best = None
    for tie in range(a.dives):
        p = dive(cols, a.N, say, tie=tie)
        J, cv, sj = evaluate(cols, p)
        p2, J2 = local_swap(cols, p, a.N, say)
        J2, cv2, sj2 = evaluate(cols, p2)
        say(f'  dive {tie}: {len(p)} craft {cv} covered sum J_i {sj:.4f} J {J:.4f} -> swap: {cv2} covered {sj2:.4f} J {J2:.4f}')
        if best is None or J2 < best[0]:
            best = (J2, p2)
    UB, inc = best
    gap = UB - LB
    keep = [j for j in range(len(cols)) if rc[j] <= gap + 1e-9 or j in inc]
    say(f'incumbent J {UB:.4f}; LB {LB:.4f}; gap {gap:.4f}: reduced-cost fixing keeps {len(keep)} of {len(cols)} columns')
    sub = [cols[j] for j in keep]
    pick, msg = exact(sub, a.N, a.tlim, say)
    if pick is not None:
        Jm, cvm, sjm = evaluate(sub, pick)
        say(f'MILP ({msg[:50]}): {len(pick)} craft, {cvm} covered, sum J_i {sjm:.4f}, J {Jm:.4f}')
        if Jm < UB - 1e-9:
            UB = Jm; final = [sub[j] for j in pick]
        else:
            final = [cols[j] for j in inc]
    else:
        final = [cols[j] for j in inc]
    Jf, cvf, sjf = evaluate(final, list(range(len(final))))
    fl = RI.IFleet()
    for i, c in enumerate(sorted(final, key=lambda c: (-len(c['asts']), c['tank']))):
        st = SM._state(c['f']); fl.routes[f'{i + 1:02d}'] = dict(st=st, tank=float(c['tank']))
    fl.save(out / f'N{a.N}' / 'fleet', note=f's15_select N<={a.N}')
    cov = set(fl.coverage())
    rep = dict(N=a.N, LB=round(LB, 4), J=round(Jf, 4), covered=cvf, sumJi=round(sjf, 4), misses=sorted(set(ALL) - cov),
               depths=[len(c['asts']) for c in sorted(final, key=lambda c: -len(c['asts']))],
               picks=[dict(f=c['f'], n=len(c['asts']), tank=round(c['tank'], 2)) for c in final],
               columns=len(cols), kept=len(keep), milp=msg if pick is not None else None, wall_s=round(time.time() - tic),
               fleet=str(out / f'N{a.N}' / 'fleet'))
    json.dump(rep, open(out / 'select.json', 'w'), indent=1)
    say(f'FINAL N<={a.N}: {cvf} covered, sum J_i {sjf:.4f}, J {Jf:.4f} (LB {LB:.4f}); misses {rep["misses"]}; '
        f'depths {rep["depths"]} ({rep["wall_s"]} s)')


if __name__ == '__main__':
    main()
