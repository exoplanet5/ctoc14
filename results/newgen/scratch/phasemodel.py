"""Is a route's Delta-v predicted by the slope variation of its phase curve?

Model: every flyby happens at r ~ 1 AU. Write the craft's mean longitude lambda(t). Between flybys the craft coasts on a
fixed semi-major axis, so lambda(t) is a straight line with slope s = n(a) - n_E [deg/yr] relative to Earth's frame; a
manoeuvre changes the slope. Changing the drift by ds costs dv = k |ds| with k = 0.028 km/s per deg/yr
(dv = mu da / (2 v a^2), ds/n = -1.5 da/a). Out of plane the craft's z at the flyby must match the asteroid's.
Predicted cost: k * sum |slope change| over the route (launch slope free: v_inf) + out-of-plane term."""
import sys, numpy as np
sys.path.insert(0, '.'); sys.path.insert(0, 'tools')
from run_ialns import IFleet, ipr, eph
from ctoc14.constants import AU, DAY, MU
YR = 365.25 * DAY
E = eph(); F = IFleet('results/newgen/ifleet10')
def wrap(x): return (x + np.pi) % (2 * np.pi) - np.pi
rE0, vE0 = E.earth_state(0.0); nE = np.sqrt(MU / (1.00091750 * AU) ** 3)
K = 0.028
print('route | flybys | actual dv | sum|d slope| (deg/yr) | k*sum | out-of-plane sum |z| (AU) | pred total')
rows = []
for n, r in sorted(F.routes.items()):
    ip = ipr(r['st']); o = np.argsort(ip.tf); tf = ip.tf[o]; asts = np.array(ip.asts)[o]
    ra, va = E.ast_states_at(asts - 1, tf)
    lam = np.unwrap(np.arctan2(ra[:, 1], ra[:, 0]))
    # Earth-frame phase: subtract Earth's mean motion
    lamE = nE * tf
    ph = np.degrees(np.unwrap(lam - lamE))          # craft phase relative to Earth at each flyby [deg]
    t = tf / YR
    slopes = np.diff(ph) / np.diff(t)               # deg/yr between consecutive flybys
    dsl = np.abs(np.diff(slopes)).sum()
    z = np.abs(ra[:, 2]) / AU
    # out-of-plane: cost of matching z sequence with one inclination is ~0 if the z's fit a sinusoid; take the
    # residual around the best-fit sinusoid in lambda as the part that must be paid
    Aq = np.stack([np.sin(lam), np.cos(lam), np.ones_like(lam)], axis=1)
    coef, *_ = np.linalg.lstsq(Aq, ra[:, 2] / AU, rcond=None)
    zres = np.abs(ra[:, 2] / AU - Aq @ coef)
    dv_z = 29.8 * zres.sum() * 2                     # crude: each residual paid twice (out and back)
    pred = K * dsl + dv_z
    rows.append((ip.dv(), K * dsl, dv_z, pred))
    print(f'{n} | {len(tf):2d} | {ip.dv():6.2f} | {dsl:7.1f} | {K * dsl:6.2f} | {zres.sum():.3f} -> {dv_z:5.2f} | {pred:6.2f}')
rows = np.array(rows)
c = np.corrcoef(rows[:, 0], rows[:, 3])[0, 1]
print(f'correlation actual vs predicted: {c:.3f}; mean ratio pred/actual {np.mean(rows[:, 3] / rows[:, 0]):.2f}')
print(f'phase-only: corr {np.corrcoef(rows[:, 0], rows[:, 1])[0, 1]:.3f}, ratio {np.mean(rows[:, 1] / rows[:, 0]):.2f}')
