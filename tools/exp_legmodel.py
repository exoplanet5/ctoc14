"""Leg-model diagnostic: Lambert-junction (planner rule) vs linearised continuous-thrust min-fuel legs (L2 / IRLS-L1),
for the states of a given tour and its target pool. Usage: exp_legmodel.py tour.json [exclude.json] [tofmax_days]"""
import sys, json, pathlib, time, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ctoc14.kepler import Ephemeris, propagate_twobody
from ctoc14.lambert import lambert
from ctoc14.constants import DAY, TMAX, VE
eph = Ephemeris(); tour = json.load(open(sys.argv[1]))
excl = set(json.load(open(sys.argv[2]))) if len(sys.argv) > 2 and sys.argv[2] != '-' else set()
TOFMAX = float(sys.argv[3]) if len(sys.argv) > 3 else 1000
legs = tour['legs']; n = len(legs); ETA = 0.6; UB = 0.8; EPS = 1e-4; DTAU = 10 * DAY
TOFS = np.arange(20, TOFMAX + 1e-9, 10) * DAY
taus = (np.arange(int(TOFMAX / 10)) + 0.5) * DTAU
# states after each flyby (Lambert arrival velocity, linear mass model)
states = []; m = tour['m0']; fuel_per = tour['fuel_est'] / n
rE, vE = eph.earth_state(tour['t_launch']); r_prev = rE; visited = set(excl) | {131, 144}
for l in legs:
    R = eph.ast_state(l['ast'] - 1, l['t_flyby'])[0]
    v1, v2 = lambert(r_prev[None], R[None], np.array([l['tof']]))
    m -= fuel_per; visited.add(l['ast']); states.append((l['t_flyby'], R, v2[0], m, set(visited))); r_prev = R

def linear_legs(t, r, v, m, mask):
    """Best linearised cost per target over TOFS. Returns dict tof-> (cost[300], feas[300]) summarised as best per target."""
    a = TMAX / m * 1e-3
    rk, vk = propagate_twobody(r, v, taus)
    rc, vc = propagate_twobody(r, v, TOFS)
    # Phi_rv[j, k, :, i]
    nJ, nK = len(TOFS), len(taus)
    Phi = np.zeros((nJ, nK, 3, 3))
    for k in range(nK):
        jj = TOFS > taus[k]
        if not jj.any(): break
        for i in range(3):
            dvp = np.zeros(3); dvp[i] = EPS
            rp, _ = propagate_twobody(rk[k], vk[k] + dvp, TOFS[jj] - taus[k])
            Phi[jj, k, :, i] = (rp - rc[jj]) / EPS * DTAU
    best_l1 = np.full(300, np.inf); best_l2 = np.full(300, np.inf); best_tof = np.full(300, np.nan)
    for j, tof in enumerate(TOFS):
        if t + tof > 473364000: break
        K = int(np.sum(taus < tof))
        if K < 2: continue
        Mk = Phi[j, :K]                                  # (K,3,3) blocks: dr = sum_k Mk[k] @ u_k
        Rall = eph.all_ast_states(t + tof)[0]; dR = Rall - rc[j]   # (300,3)
        MMk = np.einsum('kij,kjl->kil', Mk, Mk.transpose(0, 2, 1))  # (K,3,3) = M_k M_k^T
        w = np.ones((300, K))
        for it in range(5):
            G = np.einsum('tk,kij->tij', 1 / w, MMk)
            y = np.linalg.solve(G, dR[..., None])[..., 0]   # (300,3)
            U = np.einsum('tk,kji,tj->tki', 1 / w, Mk, y)   # u_tk = (1/w) M_k^T y   (300,K,3)
            un = np.linalg.norm(U, axis=-1)              # (300,K)
            if it == 0:
                l2_cost = un.sum(1) * DTAU; l2_max = un.max(1)
            w = np.maximum(un, 1e-9)
        l1_cost = un.sum(1) * DTAU; l1_max = un.max(1)
        cost = np.where(l1_max <= UB * a, l1_cost, np.where(l2_max <= UB * a, l2_cost, np.inf))
        cost[~mask] = np.inf
        upd = cost < best_l1; best_l1[upd] = cost[upd]; best_tof[upd] = tof / DAY
    return best_l1, best_tof

def lambert_legs(t, r, v, m, mask):
    a = TMAX / m * 1e-3
    best = np.full(300, np.inf); best_tof = np.full(300, np.nan)
    for tof in TOFS:
        if t + tof > 473364000: break
        Rall = eph.all_ast_states(t + tof)[0]
        v1, _ = lambert(np.broadcast_to(r, (300, 3)), Rall, np.full(300, tof))
        dv = np.linalg.norm(v1 - v, axis=-1)
        feas = mask & np.isfinite(dv) & (dv <= ETA * a * tof) & (dv <= 2.5)
        kappa = 1 + dv / (2 * a * tof); cost = np.where(feas, dv * kappa, np.inf)
        upd = cost < best; best[upd] = cost[upd]; best_tof[upd] = tof / DAY
    return best, best_tof

rows = []; tic = time.time()
for si, (t, r, v, m, vis) in enumerate(states[:-1]):
    if t > 473364000 - 400 * DAY: break
    mask = np.ones(300, bool); mask[[i - 1 for i in vis]] = False
    bl, tl = lambert_legs(t, r, v, m, mask); bq, tq = linear_legs(t, r, v, m, mask)
    both = np.isfinite(bl) & np.isfinite(bq)
    rows.append(dict(state=si, t_yr=t / DAY / 365.25, pool=int(mask.sum()),
                     lam=[int((bl <= c).sum()) for c in (0.5, 1.0, 1.5, 2.5)], lin=[int((bq <= c).sum()) for c in (0.5, 1.0, 1.5, 2.5)],
                     ratio=float(np.median(bl[both] / bq[both])) if both.any() else np.nan,
                     lin_only_le1=int(((bq <= 1.0) & ~(bl <= 1.0)).sum()), tof_lin_med=float(np.nanmedian(tq[bq <= 1.0])) if (bq <= 1).any() else np.nan,
                     tof_lam_med=float(np.nanmedian(tl[bl <= 1.0])) if (bl <= 1).any() else np.nan))
    print(f"state {si:2d} t={rows[-1]['t_yr']:.1f}yr pool={rows[-1]['pool']:3d} | Lambert n(<=0.5/1/1.5/2.5 km/s)={rows[-1]['lam']} | linear={rows[-1]['lin']} | "
          f"lin-only<=1: {rows[-1]['lin_only_le1']}  med ratio lam/lin={rows[-1]['ratio']:.2f}  medTOF(<=1) lam {rows[-1]['tof_lam_med']:.0f} lin {rows[-1]['tof_lin_med']:.0f} d  ({time.time()-tic:.0f}s)", flush=True)
L = np.array([r['lam'] for r in rows]); Q = np.array([r['lin'] for r in rows])
print('MEAN per state  Lambert', L.mean(0).round(1), ' linear', Q.mean(0).round(1), ' lin-only<=1 per state', np.mean([r['lin_only_le1'] for r in rows]).round(1))
json.dump(rows, open(pathlib.Path(sys.argv[1]).stem + '_legmodel.json', 'w'))
