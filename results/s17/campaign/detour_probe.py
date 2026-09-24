"""Stage 17 CAMPAIGN X3: 3-segment 'detour' children on the 286 max-coverage 8-route selection.

Detour(A, B) = A's flybys before t1 + junction A(t1) -> B(t1+T1) + B's flybys in (t1+T1, t2) + junction B(t2) -> A(t2+T2)
+ A's flybys after t2+T2.  Net coverage change of the fleet = (B-window targets that are fleet misses or A-window uniques)
- (A-window targets covered by no other pick).  Extra dv = dv1 + dv2 + dv_B(window) - dv_A(window).
Output detour.json.  Reads milp.json (parents_N8), reps.pkl, states.npz (read-only products of xover_scan)."""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pickle, pathlib, collections
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import numpy as np
from scipy.spatial import cKDTree
import xover_scan as X
from ctoc14.kepler import propagate_batch
from ctoc14.lambert import lambert
from ctoc14.constants import DAY, AU, VE
import run_ialns as RI

say = lambda s: X.say('[detour] ' + s)
GRID = X.GRID; K = len(GRID); LEADS = X.LEADS; GS = X.GSTEP


def route_states(files):
    S = np.full((len(files), K, 6), np.nan); DVC = np.full((len(files), K), np.nan); FL = []; DVT = []; TK = []
    for i, f in enumerate(files):
        z = np.load(f); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL'])
        ip = RI.ipr(st); Yf, Ys = ip.integrate(); nrm = np.linalg.norm(ip.Ts, axis=1)
        pr = np.vstack([ip.rE[None], Ys[:, :3]]); pv = np.vstack([(ip.vE + ip.vinf)[None], Ys[:, 3:6] + ip.Ts])
        nt = np.concatenate([[ip.tL], ip.ts]); cum = np.concatenate([[0.0], np.cumsum(nrm)])
        ks = np.nonzero(GRID >= ip.tL)[0]; m = np.searchsorted(nt, GRID[ks], side='right') - 1
        r, v = propagate_batch(pr[m], pv[m], GRID[ks] - nt[m]); S[i, ks, :3] = r; S[i, ks, 3:] = v; DVC[i, ks] = cum[m]
        o = np.argsort(ip.tf); FL.append((np.asarray(ip.tf)[o], np.asarray(ip.asts)[o])); DVT.append(nrm.sum()); TK.append(ip.tank())
    return S, DVC, FL, np.array(DVT), np.array(TK)


def junctions(S1, S2, dvmax=1.0):
    """(i1, i2, k, L, dv) for Lambert junctions S1[i1] at k -> S2[i2] at k+L."""
    out = []
    for L in LEADS:
        T = L * GS; rho = min(0.06, dvmax * T / AU)
        v1ok = ~np.isnan(S1[:, :K - L, 0]); ii, kk = np.nonzero(v1ok)
        P, _ = propagate_batch(S1[ii, kk, :3], S1[ii, kk, 3:], np.full(len(ii), T))
        A = []; B = []; Kk = []
        for k in np.unique(kk):
            q = np.nonzero(kk == k)[0]; bs = np.nonzero(~np.isnan(S2[:, k + L, 0]))[0]
            if len(bs) == 0:
                continue
            tree = cKDTree(S2[bs, k + L, :3] / AU)
            for qq, h in zip(q, tree.query_ball_point(P[q] / AU, r=rho)):
                for hb in h:
                    A.append(ii[qq]); B.append(bs[hb]); Kk.append(k)
        if not A:
            continue
        A = np.array(A); B = np.array(B); Kk = np.array(Kk)
        v1, v2 = lambert(S1[A, Kk, :3], S2[B, Kk + L, :3], np.full(len(A), T))
        d = np.linalg.norm(v1 - S1[A, Kk, 3:], axis=1) + np.linalg.norm(S2[B, Kk + L, 3:] - v2, axis=1)
        ok = np.isfinite(d) & (d <= dvmax)
        out.append(np.column_stack([A[ok], B[ok], Kk[ok], np.full(ok.sum(), L), d[ok]]))
    return np.vstack(out) if out else np.zeros((0, 5))


def main():
    t0 = time.time()
    M = json.load(open(X.OUT / 'milp.json'))['parents_N8']
    files = [p['f'] for p in M['picks']]; misses = set(M['misses'])
    SA, DA, FA, DVA, TKA = route_states(files)
    D = pickle.load(open(X.OUT / 'reps.pkl', 'rb')); SB = np.load(X.OUT / 'states.npz')['S']; DB = np.load(X.OUT / 'states.npz')['DVC']
    reps = D['reps']
    Bidx = [j for j, r in enumerate(reps) if misses & set(r['asts'])]
    say(f'8 picks, {len(misses)} misses {sorted(misses)}; {len(Bidx)} rep families contain a miss')
    SBm = SB[Bidx]
    Jab = junctions(SA, SBm); Jba = junctions(SBm, SA)
    say(f'junctions A->B {len(Jab)}, B->A {len(Jba)} ({time.time()-t0:.0f} s)')
    cover_count = collections.Counter(t for f in FA for t in f[1])
    best = []
    for a in range(len(files)):
        tfA, asA = FA[a]
        uniqA = np.array([cover_count[int(t)] == 1 for t in asA])
        for jb, b in enumerate(Bidx):
            d1 = Jab[(Jab[:, 0] == a) & (Jab[:, 1] == jb)]; d2 = Jba[(Jba[:, 0] == jb) & (Jba[:, 1] == a)]
            if len(d1) == 0 or len(d2) == 0:
                continue
            tfB = D['fl_t'][b]; asB = D['fl_a'][b]
            isM = np.array([int(t) in misses for t in asB])
            # prefix counts on the grid: B misses with tf < GRID[k]; A uniques with tf < GRID[k] (+ <= for window ends)
            nMB = np.searchsorted(np.sort(tfB[isM]), GRID, side='left') if isM.any() else np.zeros(K, int)
            nUA = np.searchsorted(np.sort(tfA[uniqA]), GRID, side='right') if uniqA.any() else np.zeros(K, int)
            nUA0 = np.searchsorted(np.sort(tfA[uniqA]), GRID, side='left') if uniqA.any() else np.zeros(K, int)
            k1 = d1[:, 2].astype(int); L1 = d1[:, 3].astype(int); k2 = d2[:, 2].astype(int); L2 = d2[:, 3].astype(int)
            arr = k1 + L1; ret = k2 + L2
            gain = nMB[k2][None, :] - nMB[np.minimum(arr, K - 1)][:, None]          # misses of B in (t_arr, t2)
            loss = nUA[np.minimum(ret, K - 1)][None, :] - nUA0[k1][:, None]        # A uniques in [t_dep, t_ret]
            valid = k2[None, :] > arr[:, None]
            net = np.where(valid, gain - loss, -99)
            if net.max() <= 0:
                continue
            for i1, i2 in zip(*np.nonzero(net == net.max())):
                q1 = d1[i1]; q2 = d2[i2]
                kk1, LL1, dv1 = int(q1[2]), int(q1[3]), q1[4]; kk2, LL2, dv2 = int(q2[2]), int(q2[3]), q2[4]
                ta = GRID[kk1 + LL1]; td = GRID[kk1]; t2 = GRID[kk2]; tr = GRID[kk2 + LL2]
                winA = (tfA >= td) & (tfA <= tr)
                lostset = set(int(t) for t in asA[winA & uniqA])
                midB = set(int(t) for t in asB[(tfB > ta) & (tfB < t2)])
                ddv = dv1 + dv2 + (DB[b, kk2] - DB[b, kk1 + LL1]) - (DA[a, kk2 + LL2] - DA[a, kk1])
                tank = 601.5 * np.exp((DVA[a] + ddv) / VE)
                best.append(dict(a=a, b=int(b), k1=kk1, L1=LL1, k2=kk2, L2=LL2, dvj=round(float(dv1 + dv2), 3),
                                 ddv=round(float(ddv), 3), tank=round(float(tank), 1), tank_A=round(float(TKA[a]), 1),
                                 net=int(len(midB & (misses | lostset)) - len(lostset)), gained=sorted(midB & misses),
                                 lost=sorted(lostset - midB), nB_mid=len(midB)))
    say(f'{len(best)} positive detours ({time.time()-t0:.0f} s)')
    # per miss: cheapest detour that recovers it with net >= 1 and tank <= 1250
    okd = [d for d in best if d['tank'] <= 1250]
    per_miss = {}
    for m in sorted(misses):
        c = [d for d in okd if m in d['gained']]
        per_miss[m] = dict(n=len(c), best=min(c, key=lambda d: (d['tank'] - d['tank_A'])) if c else None)
    # greedy fleet: one detour per pick, maximise total net, recompute conflicts
    chosen = []; usedA = set(); got = set()
    for d in sorted(okd, key=lambda d: (-d['net'], d['tank'] - d['tank_A'])):
        if d['a'] in usedA or set(d['gained']) & got:
            continue
        chosen.append(d); usedA.add(d['a']); got |= set(d['gained'])
    res = dict(misses=sorted(misses), n_positive=len(best), n_positive_tank_ok=len(okd),
               per_miss={str(k): v for k, v in per_miss.items()},
               greedy=dict(n=len(chosen), net=sum(d['net'] for d in chosen), recovered=sorted(got),
                           extra_kg=round(sum(d['tank'] - d['tank_A'] for d in chosen), 1), moves=chosen),
               sec=round(time.time() - t0))
    json.dump(res, open(X.OUT / 'detour.json', 'w'), indent=1)
    say(f"per-miss detours: {{{', '.join(f'{m}:{v['n']}' for m, v in per_miss.items())}}}; greedy one-per-pick: net +{res['greedy']['net']} "
        f"(recovers {res['greedy']['recovered']}) for +{res['greedy']['extra_kg']} kg")


if __name__ == '__main__':
    main()
