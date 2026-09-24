"""Reduced (element-space) cost model of a route, calibrated on the optimised fleet.

Every flyby of this problem happens near r = 1 AU, so a craft is a slowly moving point in orbital-element space and its
propellant is the path length of that motion in the Delta-v metric (Gauss' equations at r ~ a ~ 1 AU, v = 29.8 km/s):

    dv = (v / 2) |d a| / a  +  (v / 2) |d e_vec|  +  v |d i_vec|,
    e_vec = (e cos w, e sin w),  i_vec = (i cos W, i sin W)   (inclination vector, rad)

Measured on the 10-craft fleet: the tangential part equals 0.028 km/s per deg/yr of drift-rate change (= (v/2)|da|/a
written in drift units, corr 0.994) and the normal part equals 29.8 km/s per rad of orbit-normal change (corr 0.997).

Phase: the craft's mean longitude L = W + w + M moves at n(a); relative to Earth the drift is n(a) - n_E, so a piecewise
linear phase curve in (t, L - L_E) whose slope changes cost K_DRIFT per deg/yr is the in-plane part of a route.
"""
import numpy as np
from .constants import MU, AU, DAY

YR = 365.25 * DAY
V_C = 29.784                      # km/s, circular speed at 1 AU
N_E = np.degrees(np.sqrt(MU / (1.00091750 * AU) ** 3)) * YR      # deg/yr, Earth mean motion
K_DRIFT = 0.028                   # km/s per (deg/yr) of drift-rate change  [= (V_C/2) * (2/3) / N_E]
K_PLANE = V_C                     # km/s per rad of inclination-vector change
K_ECC = V_C / 2                   # km/s per unit eccentricity-vector change
K_TOTAL = 0.82                    # measured |sum T| / (sum|T_t| + sum|T_n|): components overlap in each impulse


def elements(r, v):
    """Osculating (a, e_vec, i_vec, mean longitude) for state arrays r, v (..., 3). Angles in rad, a in km."""
    r = np.atleast_2d(r); v = np.atleast_2d(v)
    rn = np.linalg.norm(r, axis=1); v2 = (v * v).sum(1)
    a = 1.0 / (2.0 / rn - v2 / MU)
    h = np.cross(r, v); hn = np.linalg.norm(h, axis=1)
    ev = (np.cross(v, h) / MU) - r / rn[:, None]
    e = np.linalg.norm(ev, axis=1)
    hh = h / hn[:, None]
    inc = np.arccos(np.clip(hh[:, 2], -1, 1))
    # node direction: h x z_hat gives the line of nodes
    nod = np.stack([-hh[:, 1], hh[:, 0], np.zeros_like(hh[:, 0])], axis=1)
    nn = np.linalg.norm(nod, axis=1)
    Om = np.where(nn > 1e-12, np.arctan2(hh[:, 0], -hh[:, 1]), 0.0)
    i_vec = np.stack([inc * np.cos(Om), inc * np.sin(Om)], axis=1)
    # eccentricity vector in the ecliptic frame (small inclination): (e cos w~, e sin w~) with w~ = W + w
    ec = np.stack([ev[:, 0], ev[:, 1]], axis=1)
    # mean longitude via the eccentric anomaly
    cosE = np.clip((1 - rn / a) / np.maximum(e, 1e-12), -1, 1)
    sinE = np.clip((r * v).sum(1) / np.sqrt(MU * a) / np.maximum(e, 1e-12), -1, 1)
    EA = np.arctan2(sinE, cosE); M = EA - e * np.sin(EA)
    lon_p = np.arctan2(ev[:, 1], ev[:, 0])
    L = M + lon_p
    return dict(a=a, e_vec=ec, i_vec=i_vec, L=L, e=e, inc=inc)


def drift(a):
    """Phase drift rate relative to Earth [deg/yr] for semi-major axis a [km]."""
    return np.degrees(np.sqrt(MU / a ** 3)) * YR - N_E


def path_cost(el, k_total=K_TOTAL):
    """Predicted Delta-v [km/s] of a sampled element path (dict from `elements`).

    Only two independent terms: a tangential impulse moves a and e TOGETHER (d(a)/a = 2 a v dv / mu and de = 2 dv / v
    are the same number at r = a), so charging both double counts; the normal impulse moves i_vec. Measured on the
    10-craft fleet: ratio predicted/actual 1.00 +- 0.03.
    """
    d_a = K_DRIFT * np.abs(np.diff(drift(el['a']))).sum()
    d_i = K_PLANE * np.linalg.norm(np.diff(el['i_vec'], axis=0), axis=1).sum()
    return k_total * (d_a + d_i), dict(drift=d_a, plane=d_i)


def phase_of(r, t, eph):
    """Craft/asteroid phase relative to Earth [deg]: true longitude minus Earth's true longitude."""
    rE, _ = eph.earth_state(t)
    th = np.arctan2(r[..., 1], r[..., 0]); thE = np.arctan2(rE[..., 1], rE[..., 0])
    return np.degrees((th - thE + np.pi) % (2 * np.pi) - np.pi)


def bump_cost(dphi_deg, window_yr):
    """Delta-v of a phase bump: gain dphi over the free window and return to the original drift.
    Triangular slope profile: total variation 8 |dphi| / window."""
    return K_DRIFT * 8.0 * abs(dphi_deg) / max(window_yr, 1e-3)


K_Z = 8.8            # km/s per AU of |z| at a flyby, CALIBRATED on the fleet (62.6 km/s of plane Delta-v over
                     # 298 flybys of median |z| 0.024 AU). A naive round trip (2 K_PLANE |dz| / r = 59.6 per AU) is
                     # ~7x too pessimistic: the inclination vector follows one smooth path serving many flybys
                     # (each flyby only constrains it to a LINE i_vec . u(theta) = z / r), not an excursion per flyby.


def plane_cost(dz_au, r_au=1.0):
    """Delta-v attributable to flying by at |z| = dz (calibrated marginal rate, not a round trip)."""
    return K_Z * abs(dz_au) / max(r_au, 1e-6)
