"""M1 check: along tour_sc1, count targets reachable within short TOFs under (a) the planner's Lambert/eta rule and
(b) a continuous-thrust linearised reach test (least-norm acceleration profile, |u_k| <= a).  Output: results/exp_j20/reach.txt"""
import sys, json, pathlib, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ctoc14.kepler import Ephemeris, propagate_twobody
from ctoc14.lambert import lambert
from ctoc14.constants import DAY, TMAX
eph = Ephemeris(); tour = json.load(open('results/fleet_b300/tour_sc1.json'))
legs = tour['legs']; n = len(legs); ETA = 0.6; K = 12; EPS = 1e-4
# reconstruct states after each flyby
states = []; m = tour['m0']; fuel_per = tour['fuel_est'] / n
rE, vE = eph.earth_state(tour['t_launch']); r_prev, t_prev = rE, tour['t_launch']; visited = set()
for l in legs:
    R = eph.ast_state(l['ast'] - 1, l['t_flyby'])[0]
    v1, v2 = lambert(r_prev[None], R[None], np.array([l['tof']]))
    m -= fuel_per; visited.add(l['ast']); states.append((l['t_flyby'], R, v2[0], m, set(visited)))
    r_prev = R
TOFS = np.array([30, 45, 60, 90, 120]) * DAY
out = []; agg = {tof: [0, 0, 0, [], []] for tof in TOFS}
for (t, r, v, m, vis) in states[:-1]:
    a = TMAX / m * 1e-3
    mask = np.ones(300, bool); mask[[i - 1 for i in vis]] = False; mask[[130, 143]] = False
    for tof in TOFS:
        Rall = eph.all_ast_states(t + tof)[0]
        # (a) Lambert / eta rule
        v1, _ = lambert(np.broadcast_to(r, (300, 3)), Rall, np.full(300, tof))
        dv = np.linalg.norm(v1 - v, axis=-1)
        fa = mask & np.isfinite(dv) & (dv <= ETA * a * tof) & (dv <= 2.5)
        # (b) linearised continuous-thrust reach
        rc, vc = propagate_twobody(r, v, np.array([tof])); rc = rc[0]
        taus = (np.arange(K) + 0.5) * tof / K; dtau = tof / K
        rk, vk = propagate_twobody(r, v, taus)
        M = np.zeros((3, 3 * K))
        for k in range(K):
            for i in range(3):
                dvp = np.zeros(3); dvp[i] = EPS
                rp, _ = propagate_twobody(rk[k], vk[k] + dvp, np.array([tof - taus[k]]))
                M[:, 3 * k + i] = (rp[0] - rc) / EPS * dtau
        dR = Rall - rc
        G = np.linalg.inv(M @ M.T)
        U = (dR @ G) @ M                      # (300, 3K) least-norm accelerations km/s^2
        umax = np.linalg.norm(U.reshape(300, K, 3), axis=-1).max(axis=1)
        dvlt = np.linalg.norm(U.reshape(300, K, 3), axis=-1).sum(axis=1) * dtau
        fb = mask & (umax <= a); fb6 = mask & (umax <= ETA * a)
        agg[tof][0] += fa.sum(); agg[tof][1] += fb.sum(); agg[tof][2] += fb6.sum()
        agg[tof][3] += list(dv[fa]); agg[tof][4] += list(dvlt[fb])
lines = ['TOF[d]  eta-rule  LT(|u|<=a)  LT(|u|<=0.6a)   mean dv Lambert  mean dv LT   (sums over %d states of tour_sc1)' % (n - 1)]
for tof in TOFS:
    A = agg[tof]
    lines.append(f'{tof/DAY:5.0f}   {A[0]:6d}    {A[1]:6d}      {A[2]:6d}        {np.mean(A[3]) if A[3] else 0:.2f}            {np.mean(A[4]) if A[4] else 0:.2f}')
print('\n'.join(lines)); open('results/exp_j20/reach.txt', 'w').write('\n'.join(lines) + '\n')
