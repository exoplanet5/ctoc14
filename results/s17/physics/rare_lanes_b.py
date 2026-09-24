"""T1b (exploratory follow-up, not pre-registered): lane column generation for the rare set + exact max-coverage MILP (N=8, 9)."""
import os
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'): os.environ.setdefault(v,'1')
import sys, json, time, pathlib, multiprocessing as mp
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0, str(ROOT/'results/s17/physics'))
import numpy as np
import rare_lanes as RL
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix, csr_matrix
if __name__ == '__main__':
    t0 = time.time(); rng = np.random.default_rng(17); n = 40000
    tLs = rng.uniform(0, 730*RL.DAY, n)
    d = rng.normal(size=(n,3)); d /= np.linalg.norm(d, axis=1)[:,None]
    vinf = d*(RL.VINF_MAX*rng.uniform(0,1,n)**(1/3))[:,None]
    RV = np.array([np.concatenate(RL.E.earth_state(x)) for x in tLs]); r = RV[:,:3]; v = RV[:,3:]+vinf
    pa, na, Ea = RL.frame(RL.elem[:, :5]); cnt = np.zeros(n, int)
    for c0 in range(0, n, 10000):
        s = slice(c0, c0+10000); pc, nc, Ec, ac = RL.rv2frame(r[s], v[s]); g, _ = RL.gaps(pc, nc, Ec, pa, na, Ea)
        cnt[s] = (np.abs(g).min(2) < RL.EPS).sum(1)
    top = np.argsort(-cnt)[:3000]
    cols = []
    with mp.get_context('fork').Pool(2) as pl:
        w = np.ones(301)
        for it in range(6):
            L = pl.map(RL.lane, [(float(tLs[c]), r[c], v[c], w, False) for c in top], chunksize=20)
            for c, x in zip(top, L):
                if x: cols.append((frozenset(x['tg']), x['dv'], int(c)))
            # upweight rare targets that few columns carry (column generation by coverage deficit)
            cc = np.zeros(301)
            for s_, _, _ in cols:
                for X in s_: cc[X] += 1
            w = np.ones(301); low = [X for X in RL.RARE if cc[X] < np.percentile(cc[RL.RARE], 25)]
            w[low] = 3.0 + rng.uniform(0, 2, len(low))
            print(f'gen {it}: columns {len(cols)}; rare never in a lane: {[X for X in RL.RARE if cc[X]==0]}  ({time.time()-t0:.0f}s)', flush=True)
    uniq = {}
    for s_, dv, c in cols:
        if s_ not in uniq or dv < uniq[s_][0]: uniq[s_] = (dv, c)
    cols = [(s_, dv, c) for s_, (dv, c) in uniq.items()]
    RI = {X: i for i, X in enumerate(RL.RARE)}; nR = len(RL.RARE); m = len(cols)
    A = lil_matrix((nR, m))
    for j, (s_, _, _) in enumerate(cols):
        for X in s_: A[RI[X], j] = 1
    A = csr_matrix(A); out = {}
    for N in (8, 9):
        nv = m + nR; cobj = np.zeros(nv); cobj[m:] = -1.0; cobj[:m] = 1e-4*np.array([dv for _, dv, _ in cols])
        M = lil_matrix((nR, nv)); M[:, :m] = -A
        for i in range(nR): M[i, m+i] = 1
        S = np.zeros((1, nv)); S[0, :m] = 1
        res = milp(cobj, constraints=[LinearConstraint(csr_matrix(M), -np.inf, 0), LinearConstraint(S, -np.inf, N)],
                   integrality=np.ones(nv), bounds=Bounds(0, 1), options=dict(time_limit=300))
        y = np.where(np.round(res.x[:m]))[0]
        cov = set().union(*[cols[j][0] for j in y])
        out[N] = dict(covered=len(cov), sizes=[len(cols[j][0]) for j in y], dv=[round(cols[j][1],1) for j in y], missing=sorted(set(RL.RARE)-cov), msg=res.message[:40])
        print(f'N={N}: max rare coverage {len(cov)}/{nR} with lanes of {out[N]["sizes"]} rare, model dv {out[N]["dv"]}; missing {out[N]["missing"]}  ({time.time()-t0:.0f}s)', flush=True)
    json.dump(out, open(RL.OUT/'rare_lanes_b.json','w'), indent=1)
