import sys, glob, time, pathlib, json, numpy as np
sys.path.insert(0, '.')
from scipy.optimize import milp, LinearConstraint, Bounds, linprog
from scipy.sparse import hstack, vstack, identity, csr_matrix
from ctoc14.colgen import CostModel, ColumnStore, RMP, TARGETS
cm = CostModel(s=0.70, reserve=2.0); store = ColumnStore(cm)
for pat in ['results/phasing/camp*/pool.jsonl', 'results/phasing/pool/*.jsonl']:
    for f in sorted(glob.glob(pat)): store.load_jsonl(f, tag=pathlib.Path(f).parent.name)
for f in sorted(glob.glob('results/colgen/run*/columns.jsonl')) + sorted(glob.glob('results/colgen/*/columns.jsonl')):
    store.load_columns_file(f)
rows = TARGETS - 1; C = store.csr(); cost = store.costs(); nt = np.asarray(C.sum(axis=1)).ravel()
out = {}
for pen in [float(x) for x in sys.argv[1].split(',')]:
    for N in [int(x) for x in sys.argv[2].split(',')]:
        # LP with miss penalty pen over the whole store (column generation on the store)
        act = np.argsort(cost / np.maximum(nt, 1))[:3000]
        for it in range(30):
            Cs = C[act][:, rows]; n, nR = Cs.shape
            A = vstack([hstack([-Cs.T.tocsr(), -identity(nR, format='csr')]), csr_matrix((np.ones(n), (np.zeros(n, int), np.arange(n))), shape=(1, n + nR))]).tocsr()
            r = linprog(np.concatenate([cost[act], np.full(nR, pen)]), A_ub=A, b_ub=np.concatenate([-np.ones(nR), [N]]), bounds=(0, 1), method='highs')
            m = -r.ineqlin.marginals; pi = np.zeros(300); pi[rows] = m[:nR]; mu = m[nR]
            rc = cost - C @ pi + mu; rc[act] = np.inf; neg = np.where(rc < -1e-6)[0]
            if len(neg) == 0: break
            act = np.concatenate([act, neg[np.argsort(rc[neg])[:3000]]])
        lpv = r.fun + 2
        # MIP over LP support + best reduced-cost columns
        rc_all = cost - C @ pi + mu
        cand = np.unique(np.concatenate([act[r.x[:len(act)] > 1e-6], np.argsort(rc_all)[:4000]]))
        Cs = C[cand][:, rows]; n, nR = Cs.shape
        A = vstack([hstack([-Cs.T.tocsr(), -identity(nR, format='csr')]), csr_matrix((np.ones(n), (np.zeros(n, int), np.arange(n))), shape=(1, n + nR))]).tocsr()
        res = milp(np.concatenate([cost[cand], np.full(nR, pen)]), constraints=LinearConstraint(A, -np.inf, np.concatenate([-np.ones(nR), [N]])),
                   integrality=np.concatenate([np.ones(n), np.zeros(nR)]), bounds=Bounds(0, 1), options=dict(time_limit=240, mip_rel_gap=1e-3))
        sel = cand[np.round(res.x[:n]) > 0]
        cov = set()
        for j in sel: cov |= set((store.tidx[j] + 1).tolist())
        miss = sorted(set(TARGETS.tolist()) - cov)
        J_pen = float(cost[sel].sum()) + pen * len(miss) + 2
        print(f'pen {pen} N {N}: LP {lpv:.3f}; MIP {len(sel)} routes, sum c {cost[sel].sum():.3f}, covered {len(cov)}, misses {len(miss)} -> J(pen) {J_pen:.3f}, J(all missed=1) {cost[sel].sum()+len(miss)+2:.3f}', flush=True)
        out[f'{pen}_{N}'] = dict(sel=[int(j) for j in sel], tags=[store.tag[j] for j in sel], src=[store.src[j] if store.src[j][0] != 'inline' else ('inline', None) for j in sel], miss=miss, sumc=float(cost[sel].sum()))
json.dump(out, open('results/newgen/mip_pen.json', 'w'), indent=1, default=str)
