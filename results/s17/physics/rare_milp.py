"""T0: can 8 deep pool routes carry pairwise-DISJOINT rare sets that cover the 95 rare targets, and what union do they give?"""
import json, numpy as np, pathlib, time
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix, csr_matrix
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); OUT = ROOT/'results/s17/physics'
ev = json.load(open(OUT/'events.json'))
nev = np.zeros(301, int)
for e in ev:
    if e['d'] < 0.05: nev[e['ast']] += 1
TG = [X for X in range(1,301) if X not in (131,144)]; TI = {X: i for i, X in enumerate(TG)}
rare = set(X for X in TG if nev[X] <= 4)
cols = []; seen = set()
for line in open(ROOT/'results/s16/honest/index.jsonl'):
    d = json.loads(line)
    if d.get('status') == 'rejected': continue
    fn = ROOT/'results/s16/honest'/f"route_{d['key']}.npz"
    if not fn.exists(): continue
    a = tuple(sorted(int(x) for x in np.load(fn)['asts']))
    if len(a) < 36 or a in seen: continue
    seen.add(a); tank = d.get('tank_out') or d.get('tank_in')
    cols.append(dict(key=d['key'], asts=a, tank=tank))
print('deep distinct columns', len(cols))
R = len(cols); NT = len(TG)
A = lil_matrix((NT, R))
for j, c in enumerate(cols):
    for X in c['asts']:
        if X in TI: A[TI[X], j] = 1
A = csr_matrix(A)
rare_idx = [TI[X] for X in TG if X in rare]; freq_idx = [TI[X] for X in TG if X not in rare]
def solve(N, w_freq, disjoint_rare=True, tl=300):
    # vars: y (R), c (NT) coverage indicators
    nv = R + NT
    cobj = np.zeros(nv)
    cobj[R + np.array(rare_idx)] = -1.0; cobj[R + np.array(freq_idx)] = -w_freq
    cons = []
    # c_t <= sum_r A_tr y_r
    M = lil_matrix((NT, nv)); M[:, :R] = -A; 
    for i in range(NT): M[i, R+i] = 1
    cons.append(LinearConstraint(csr_matrix(M), -np.inf, 0))
    if disjoint_rare:
        D = lil_matrix((len(rare_idx), nv)); D[:, :R] = A[rare_idx]
        cons.append(LinearConstraint(csr_matrix(D), -np.inf, 1))
    S = np.zeros((1, nv)); S[0, :R] = 1; cons.append(LinearConstraint(S, -np.inf, N))
    t0 = time.time()
    r = milp(cobj, constraints=cons, integrality=np.ones(nv), bounds=Bounds(0, 1), options=dict(time_limit=tl, disp=False))
    y = np.round(r.x[:R]).astype(int) if r.x is not None else np.zeros(R, int)
    pick = [cols[j] for j in np.where(y)[0]]
    cov = set(X for c in pick for X in c['asts'])
    rc = len(cov & rare); uc = len(cov)
    slots = sum(len(c['asts']) for c in pick)
    ov_rare = slots_r = sum(1 for c in pick for X in c['asts'] if X in rare)
    fuel = sum(c['tank'] - 600 for c in pick)
    return dict(N=N, w=w_freq, disjoint=disjoint_rare, status=r.message[:40], secs=round(time.time()-t0,1), rare_cov=rc, union=uc, slots=slots,
                rare_slots=slots_r, fuel=round(fuel), depths=[len(c['asts']) for c in pick], keys=[c['key'] for c in pick])
res = []
for N in (8, 9):
    for w in (0.0, 0.3):
        s = solve(N, w); res.append(s)
        print({k: v for k, v in s.items() if k != 'keys'}, flush=True)
s = solve(8, 0.3, disjoint_rare=False); res.append(s); print('no-disjoint', {k: v for k, v in s.items() if k != 'keys'}, flush=True)
json.dump(res, open(OUT/'rare_milp.json', 'w'), indent=1)
