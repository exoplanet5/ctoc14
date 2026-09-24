"""Is the fleet's out-of-plane Delta-v near the minimum for its flyby sequences?

Each flyby forces the craft's orbit plane to contain the asteroid: i_vec . u(theta) = z / r, a LINE in the
inclination-vector plane. The cheapest plane history visiting those lines in time order is the shortest path touching
them in order (convex). Compare its length (x 29.8 km/s per rad) with the route's actual normal Delta-v."""
import sys, numpy as np
from scipy.optimize import minimize
sys.path.insert(0, '.'); sys.path.insert(0, 'tools')
from run_ialns import IFleet, ipr, eph
from ctoc14.phasemodel import elements, K_PLANE, K_TOTAL
from ctoc14.constants import AU
E = eph(); F = IFleet('results/newgen/ifleet10')
print('route | actual plane dv | shortest touring path | saving | max |i| deg')
tot_a = tot_b = 0.0
for n, r in sorted(F.routes.items()):
    ip = ipr(r['st']); Yf, Ys = ip.integrate(); o = np.argsort(ip.tf)
    ra, _ = E.ast_states_at(np.array(ip.asts)[o] - 1, ip.tf[o])
    th = np.arctan2(ra[:, 1], ra[:, 0]); rr = np.linalg.norm(ra, axis=1) / AU; z = ra[:, 2] / AU
    u = np.stack([np.sin(th), -np.cos(th)], axis=1)            # i_vec . u = z / r
    c = z / rr
    p0 = u * c[:, None]                                        # closest point of each line to the origin
    d = np.stack([-u[:, 1], u[:, 0]], axis=1)                  # line direction
    def length(sv):
        x = p0 + d * sv[:, None]
        return np.sqrt(((np.diff(x, axis=0)) ** 2).sum(1) + 1e-12).sum()
    best = None
    for s0 in (np.zeros(len(c)), np.full(len(c), 0.01), np.full(len(c), -0.01)):
        res = minimize(length, s0, method='L-BFGS-B', options=dict(maxiter=4000))
        if best is None or res.fun < best.fun:
            best = res
    x = p0 + d * best.x[:, None]
    el = elements(Ys[:, :3], Ys[:, 3:6])
    act = K_TOTAL * K_PLANE * np.linalg.norm(np.diff(el['i_vec'], axis=0), axis=1).sum()
    opt = K_TOTAL * K_PLANE * best.fun
    tot_a += act; tot_b += opt
    print(f'{n} | {act:6.2f} | {opt:6.2f} | {act - opt:5.2f} | {np.degrees(np.linalg.norm(x, axis=1)).max():.2f}')
print(f'fleet plane dv: actual {tot_a:.1f}, lower bound for the SAME sequences {tot_b:.1f} km/s '
      f'(potential saving {tot_a - tot_b:.1f} km/s = {100 * (tot_a - tot_b) / tot_a:.0f} %)')
