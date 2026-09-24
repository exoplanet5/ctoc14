"""Stage 17 CAMPAIGN test X1: state-proximal recombination (2-segment crossover children) of the honest pool.

Child(A, t, B, T) = A's flybys before t  +  Lambert junction from A's twin state at t to B's twin state at t+T  +
B's flybys after t+T (B's impulses unchanged).  dv_child = dvcum_A(t) + dv_J + (dv_B - dvcum_B(t+T)).
Parents: results/s16/colcache.pkl (read-only), deduped by target set; junction search over Jaccard<0.8 representatives.

usage: xover_scan.py stage   (stage = states | pairs | milp | all)
Outputs (results/s17/campaign/): states.npz, children.pkl, milp.json, log.txt
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pickle, pathlib, collections
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)
import numpy as np
import warnings; warnings.filterwarnings('ignore')
from scipy.spatial import cKDTree
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix, csr_matrix, hstack, identity, vstack
from ctoc14.kepler import propagate_batch
from ctoc14.lambert import lambert
from ctoc14.constants import DAY, AU, VE, T_MISSION
import run_ialns as RI

OUT = ROOT / 'results/s17/campaign'
UNREACH = {131, 144}
ALL = [t for t in range(1, 301) if t not in UNREACH]; TI = {t: i for i, t in enumerate(ALL)}; NT = len(ALL)
GSTEP = 10 * DAY
GRID = np.arange(0, T_MISSION, GSTEP)
LEADS = (3, 5, 8, 11, 15)            # x 10 d
DVJ_MAX = 1.0; TANK_MAX = 1250.0; JAC = 0.8


def say(*a):
    s = time.strftime('%H:%M:%S ') + ' '.join(str(x) for x in a)
    print(s, flush=True)
    with open(OUT / 'log.txt', 'a') as f:
        f.write(s + '\n')


def parents():
    c = pickle.load(open(ROOT / 'results/s16/colcache.pkl', 'rb'))
    best = {}
    for v in c.values():
        col = v[0]
        if col is None:
            continue
        k = frozenset(col['asts'])
        if k not in best or col['tank'] < best[k]['tank']:
            best[k] = col
    cols = sorted(best.values(), key=lambda x: (-len(x['asts']), x['tank']))
    reps = []; rsets = []
    for cc in cols:
        s = frozenset(cc['asts'])
        if all(len(s & r) / len(s | r) < JAC for r in rsets):
            reps.append(cc); rsets.append(s)
    return cols, reps


def stage_states():
    cols, reps = parents()
    say(f'parents {len(cols)} distinct sets, {len(reps)} reps (Jaccard<{JAC})')
    R = len(reps); K = len(GRID)
    S = np.full((R, K, 6), np.nan); DVC = np.full((R, K), np.nan); TK = np.zeros(R); DVT = np.zeros(R)
    fl_t = []; fl_a = []; imp_t = []; imp_n = []; tL = np.zeros(R)
    r0s = []; v0s = []; dts = []; idx = []
    for i, cc in enumerate(reps):
        z = np.load(cc['f']); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL'])
        ip = RI.ipr(st); Yf, Ys = ip.integrate()
        o = np.argsort(ip.tf); tf = np.asarray(ip.tf)[o]; asts = np.asarray(ip.asts)[o]
        fl_t.append(tf); fl_a.append(asts)
        ts = np.asarray(ip.ts); Ts = np.asarray(ip.Ts); nrm = np.linalg.norm(Ts, axis=1)
        imp_t.append(ts); imp_n.append(nrm); tL[i] = ip.tL
        TK[i] = ip.tank(); DVT[i] = nrm.sum()
        post_r = np.vstack([ip.rE[None], Ys[:, :3]]); post_v = np.vstack([(ip.vE + ip.vinf)[None], Ys[:, 3:6] + Ts])
        node_t = np.concatenate([[ip.tL], ts]); cum = np.concatenate([[0.0], np.cumsum(nrm)])
        ks = np.nonzero(GRID >= ip.tL)[0]
        n = np.searchsorted(node_t, GRID[ks], side='right') - 1
        r0s.append(post_r[n]); v0s.append(post_v[n]); dts.append(GRID[ks] - node_t[n]); idx.append((i, ks))
        DVC[i, ks] = cum[n]
    r0 = np.vstack(r0s); v0 = np.vstack(v0s); dt = np.concatenate(dts)
    say(f'propagating {len(dt)} grid states')
    rr = np.zeros_like(r0); vv = np.zeros_like(v0)
    for a in range(0, len(dt), 200000):
        rr[a:a + 200000], vv[a:a + 200000] = propagate_batch(r0[a:a + 200000], v0[a:a + 200000], dt[a:a + 200000])
    p = 0
    for i, ks in idx:
        S[i, ks, :3] = rr[p:p + len(ks)]; S[i, ks, 3:] = vv[p:p + len(ks)]; p += len(ks)
    pickle.dump(dict(reps=[dict(f=c['f'], tank=c['tank'], asts=list(c['asts'])) for c in reps], fl_t=fl_t, fl_a=fl_a,
                     imp_t=imp_t, imp_n=imp_n, tL=tL, TK=TK, DVT=DVT), open(OUT / 'reps.pkl', 'wb'))
    np.savez(OUT / 'states.npz', S=S, DVC=DVC)
    say('states saved')


def stage_pairs():
    D = pickle.load(open(OUT / 'reps.pkl', 'rb')); Z = np.load(OUT / 'states.npz'); S = Z['S']; DVC = Z['DVC']
    R, K, _ = S.shape
    fl_t, fl_a, imp_t, imp_n, DVT = D['fl_t'], D['fl_a'], D['imp_t'], D['imp_n'], D['DVT']
    children = {}; stats = collections.Counter(); jn = []
    for L in LEADS:
        T = L * GSTEP; rho = min(0.06, DVJ_MAX * T / AU)
        valid = ~np.isnan(S[:, :, 0])
        # ballistic propagation of every (route, k) by T
        ii, kk = np.nonzero(valid[:, :K - L])
        P = np.zeros((len(ii), 3))
        for a in range(0, len(ii), 200000):
            P[a:a + 200000], _ = propagate_batch(S[ii[a:a + 200000], kk[a:a + 200000], :3],
                                                  S[ii[a:a + 200000], kk[a:a + 200000], 3:], np.full(min(200000, len(ii) - a), T))
        A_l = []; B_l = []; k_l = []
        byk = collections.defaultdict(list)
        for q, k in enumerate(kk):
            byk[k].append(q)
        for k, qs in byk.items():
            bs = np.nonzero(valid[:, k + L])[0]
            if len(bs) == 0:
                continue
            tree = cKDTree(S[bs, k + L, :3] / AU)
            qs = np.array(qs)
            hits = tree.query_ball_point(P[qs] / AU, r=rho)
            for q, h in zip(qs, hits):
                a = ii[q]
                for hb in h:
                    b = bs[hb]
                    if b != a:
                        A_l.append(a); B_l.append(b); k_l.append(k)
        A_l = np.array(A_l, int); B_l = np.array(B_l, int); k_l = np.array(k_l, int)
        stats[f'cand_L{L}'] = len(A_l)
        if len(A_l) == 0:
            continue
        # Lambert pricing, in chunks
        dvj = np.full(len(A_l), np.inf)
        for a in range(0, len(A_l), 100000):
            s = slice(a, a + 100000)
            r1 = S[A_l[s], k_l[s], :3]; r2 = S[B_l[s], k_l[s] + L, :3]
            try:
                v1, v2 = lambert(r1, r2, np.full(len(r1), T))
                d = np.linalg.norm(v1 - S[A_l[s], k_l[s], 3:], axis=1) + np.linalg.norm(S[B_l[s], k_l[s] + L, 3:] - v2, axis=1)
                dvj[s] = np.where(np.isfinite(d), d, np.inf)
            except Exception as e:
                say('lambert fail', e)
        ok = dvj <= DVJ_MAX
        stats[f'junc_L{L}'] = int(ok.sum())
        for a, b, k, d in zip(A_l[ok], B_l[ok], k_l[ok], dvj[ok]):
            t1 = GRID[k]; t2 = GRID[k + L]
            pre = fl_a[a][fl_t[a] < t1]; suf = fl_a[b][fl_t[b] > t2]
            dv = DVC[a, k] + d + (DVT[b] - DVC[b, k + L])
            tank = 601.5 * np.exp(dv / VE)
            jn.append((int(a), int(b), int(k), L, float(d), len(pre), len(suf)))
            if tank > TANK_MAX:
                stats['tank'] += 1; continue
            tcap = 0.43 * 20 * DAY / (tank * 1e3)
            sufimp = imp_n[b][imp_t[b] > t2]
            if len(sufimp) and sufimp.max() > 1.05 * tcap:
                stats['tcap'] += 1; continue
            s = frozenset(int(x) for x in np.concatenate([pre, suf]) if int(x) not in UNREACH)
            if len(s) < 5:
                continue
            if s not in children or tank < children[s]['tank']:
                children[s] = dict(tank=float(tank), a=int(a), b=int(b), k=int(k), L=L, dvj=float(d),
                                   npre=len(pre), nsuf=len(suf), dup=len(pre) + len(suf) - len(s))
        say(f'L={L} ({L*10} d, rho {rho:.3f} AU): cand {stats[f"cand_L{L}"]}, junctions <= {DVJ_MAX} km/s '
            f'{stats[f"junc_L{L}"]}, children so far {len(children)}')
    jn = np.array(jn, float) if jn else np.zeros((0, 7))
    np.save(OUT / 'junctions.npy', jn)
    pickle.dump(children, open(OUT / 'children.pkl', 'wb'))
    say('pairs done', dict(stats), 'children', len(children))


def load_cols(with_children):
    cols, _ = parents()
    out = [dict(f=c['f'], tank=c['tank'], set=frozenset(c['asts']), kind='parent') for c in cols]
    if with_children:
        ch = pickle.load(open(OUT / 'children.pkl', 'rb'))
        pset = {c['set']: c['tank'] for c in out}
        for s, c in ch.items():
            if s in pset and pset[s] <= c['tank']:
                continue
            out.append(dict(c, set=s, kind='child'))
    for c in out:
        c['cost'] = RI.cost(c['tank'])
    # dominance pruning: drop c2 if some c1 has superset targets at <= cost
    out.sort(key=lambda c: (-len(c['set']), c['cost']))
    by_t = collections.defaultdict(list); keep = []
    for c in out:
        t0 = min(c['set'], key=lambda t: len(by_t[t])) if c['set'] else None
        if t0 is not None and any(keep[i]['cost'] <= c['cost'] + 1e-12 and c['set'] <= keep[i]['set'] for i in by_t[t0]):
            continue
        j = len(keep); keep.append(c)
        for t in c['set']:
            by_t[t].append(j)
    return keep


def maxcov(cols, N, tlim, eps=1e-3):
    R = []; C = []
    for j, c in enumerate(cols):
        for t in c['set']:
            if t in TI:
                R.append(TI[t]); C.append(j)
    nr = len(cols)
    A = coo_matrix((np.ones(len(R)), (R, C)), shape=(NT, nr)).tocsr()
    cost = np.array([c['cost'] for c in cols])
    M = hstack([A, -identity(NT, format='csr')]).tocsr()
    cons = [LinearConstraint(M, 0, np.inf),
            LinearConstraint(csr_matrix(np.concatenate([np.ones(nr), np.zeros(NT)])[None, :]), -np.inf, N)]
    obj = np.concatenate([eps * cost, -np.ones(NT)])
    r = milp(obj, constraints=cons, integrality=np.concatenate([np.ones(nr), np.zeros(NT)]), bounds=Bounds(0, 1),
             options=dict(time_limit=tlim, mip_rel_gap=1e-6))
    if r.x is None:
        return None, r.message, None
    pick = [j for j in range(nr) if r.x[j] > 0.5]
    cov = set().union(*[cols[j]['set'] for j in pick]) & set(ALL)
    return pick, r.message, dict(covered=len(cov), sumJi=float(sum(cols[j]['cost'] for j in pick)),
                                  bound=float(-r.mip_dual_bound) if getattr(r, 'mip_dual_bound', None) is not None else None)


def stage_milp():
    res = {}
    for tag, wc in (('parents', False), ('with_children', True)):
        cols = load_cols(wc)
        nchild = sum(c['kind'] == 'child' for c in cols)
        say(f'{tag}: {len(cols)} columns after dominance ({nchild} children)')
        for N in (8, 9):
            t0 = time.time()
            pick, msg, info = maxcov(cols, N, tlim=600 if N == 8 else 300)
            if pick is None:
                say(tag, N, 'fail', msg); continue
            info['msg'] = msg[:60]; info['sec'] = round(time.time() - t0)
            info['picks'] = [dict(kind=cols[j]['kind'], depth=len(cols[j]['set']), tank=round(cols[j]['tank'], 1),
                                  **({k: cols[j][k] for k in ('a', 'b', 'k', 'L', 'dvj', 'npre', 'nsuf', 'dup')}
                                     if cols[j]['kind'] == 'child' else dict(f=cols[j]['f'])),
                                  set=sorted(cols[j]['set'])) for j in pick]
            info['misses'] = sorted(set(ALL) - set().union(*[cols[j]['set'] for j in pick]))
            res[f'{tag}_N{N}'] = info
            say(f'{tag} N<={N}: covered {info["covered"]} sumJi {info["sumJi"]:.4f} ({info["msg"]}, {info["sec"]} s) '
                f'depths {[p["depth"] for p in info["picks"]]} kinds {[p["kind"][0] for p in info["picks"]]}')
            json.dump(res, open(OUT / 'milp.json', 'w'), indent=1)


if __name__ == '__main__':
    st = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if st in ('states', 'all'):
        stage_states()
    if st in ('pairs', 'all'):
        stage_pairs()
    if st in ('milp', 'all'):
        stage_milp()
