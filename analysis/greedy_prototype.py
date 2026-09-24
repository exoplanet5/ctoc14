"""Quick greedy multi-flyby tour prototype (Lambert-junction model) to estimate flybys per full-tank spacecraft."""
import sys, pathlib, numpy as np, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ctoc14.kepler import Ephemeris
from ctoc14.lambert import lambert
from ctoc14.constants import VE, TMAX, T_MISSION, DAY, M_DRY

eph = Ephemeris()
EXCL = {131, 144}
TOFS = np.arange(15, 401, 5) * DAY          # candidate leg durations
ETA = float(sys.argv[2]) if len(sys.argv) > 2 else 0.6

def greedy_tour(t_launch, m0=2000.0, eta=ETA, w_time=1.0, verbose=False):
    visited = set(EXCL); m = m0; t = t_launch
    rE, vE = eph.earth_state(t)
    # first leg: pick target/TOF minimising excess over v_inf 4, prefer short TOF
    best = None
    for tof in TOFS:
        R, _ = eph.all_ast_states(t + tof)
        v1, v2 = lambert(rE, R, tof)
        vinf = np.linalg.norm(v1 - vE, axis=1)
        for k in np.argsort(vinf)[:5]:
            if k + 1 in visited or not np.isfinite(vinf[k]): continue
            excess = max(0.0, vinf[k] - 4.0)
            score = tof / DAY + 200 * excess
            if excess < 0.05 and (best is None or score < best[0]):
                best = (score, k, tof, v1[k], v2[k], R[k])
    if best is None: return []
    _, k, tof, v1, v2, R = best
    tour = [(k + 1, t + tof, 0.0)]; visited.add(k + 1); t += tof; r = R; v = v2
    while True:
        a = TMAX / m * 1e-3   # km/s^2
        R_all = np.stack([eph.all_ast_states(t + tof)[0] for tof in TOFS])   # (nT, 300, 3)
        v1, v2 = lambert(r[None, None, :], R_all, np.broadcast_to(TOFS[:, None], R_all.shape[:2]))
        dv = np.linalg.norm(v1 - v, axis=-1)                                  # (nT, 300)
        feas = dv <= eta * a * TOFS[:, None]
        for j in visited: feas[:, j - 1] = False
        feas &= (t + TOFS[:, None]) <= T_MISSION
        if not feas.any(): break
        # cost proxy: time + fuel-equivalent time (dv / (a) is the burn time)  -> minimise tof + dv/a * w
        cost = np.where(feas, TOFS[:, None] / DAY + w_time * (dv / a) / DAY, np.inf)
        idx = np.unravel_index(np.argmin(cost), cost.shape)
        tof = TOFS[idx[0]]; k = idx[1]; ddv = dv[idx]
        kappa = 1 + ddv / (2 * a * tof)
        dm = m * (1 - np.exp(-kappa * ddv / VE))
        if m - dm < M_DRY: break
        m -= dm; t += tof; r = R_all[idx[0], k]; v = v2[idx[0], k]
        tour.append((k + 1, t, ddv)); visited.add(k + 1)
        if verbose: print(f'  flyby {len(tour):3d}: ast {k+1:3d} tof {tof/DAY:5.0f} d  dv {ddv:.3f} km/s  m {m:7.1f}  t {t/DAY/365.25:5.2f} yr')
    return tour, m

if __name__ == '__main__':
    tl = float(sys.argv[1]) * DAY if len(sys.argv) > 1 else 30 * DAY
    tic = time.time()
    tour, m = greedy_tour(tl, verbose=True)
    print(f'launch day {tl/DAY:.0f}: {len(tour)} flybys, fuel used {2000 - m:.0f} kg, dv_sum {sum(x[2] for x in tour):.1f} km/s, end at {tour[-1][1]/DAY/365.25:.2f} yr, eta={ETA}, runtime {time.time()-tic:.0f}s')
