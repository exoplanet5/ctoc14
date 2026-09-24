"""
CTOC14 Keplerian ephemeris (contest definition), vectorised numpy.

Units: km, s, rad internally.  Times are given as *days after t0* (t0 = MJD 62502.0,
2030-01-01 00:00 UTC) unless stated otherwise.

Contest definition (problem statement eqs. 4-5):
    M(t) = M0 + n (t - t_eph),   n = sqrt(mu_sun / (a*AU)^3)  [rad/s],  t in seconds
Earth elements at MJD 60676.0, asteroid elements (MEA.txt) at MJD 61200.0.
"""
import os
import numpy as np

MU_SUN = 1.32712440018e11        # km^3/s^2
AU = 149597870.7                 # km
DAY = 86400.0                    # s
MJD_T0 = 62502.0                 # mission start 2030-01-01
T_MAX_DAYS = 5478.75             # 15 yr = 4.73364e8 s
MJD_EPOCH_EARTH = 60676.0
MJD_EPOCH_AST = 61200.0
ISP = 4000.0
G0 = 9.80665e-3                  # km/s^2
TMAX = 0.5e-3                    # kN  (0.5 N)
M_DRY, M_MAX = 600.0, 2000.0

# Earth elements (deg for angles)
EARTH_ELEMENTS = np.array([1.0009175020, 0.017566762041, 0.002976847126,
                           189.953211282428, 273.196254000254, 357.4135031077])

_HERE = os.path.dirname(os.path.abspath(__file__))
MEA_PATH = os.path.normpath(os.path.join(_HERE, '..', '..', 'MEA.txt'))


def load_mea(path=MEA_PATH):
    """Return (ids, elements) with elements = (N,6) array [a AU, e, i deg, Om deg, w deg, M0 deg]."""
    d = np.loadtxt(path)
    return d[:, 0].astype(int), d[:, 1:7].copy()


def kepler_E(M, e, tol=1e-13, itmax=40):
    """Solve Kepler's equation E - e sin E = M (elliptic), vectorised Newton.
    M, e broadcastable arrays (rad).  Starting value pi for high e guarantees convergence."""
    M = np.mod(M, 2 * np.pi)
    e = np.broadcast_to(e, M.shape).astype(float)
    E = np.where(e < 0.8, M + e * np.sin(M), np.pi * np.ones_like(M))
    for _ in range(itmax):
        f = E - e * np.sin(E) - M
        fp = 1.0 - e * np.cos(E)
        dE = f / fp
        E = E - dE
        if np.max(np.abs(dE)) < tol:
            break
    return E


def _rot_pqw(elements_deg):
    """Perifocal->inertial rotation matrices (N,3,3) from (i, Om, w) in deg."""
    inc = np.radians(elements_deg[..., 2])
    Om = np.radians(elements_deg[..., 3])
    w = np.radians(elements_deg[..., 4])
    cO, sO = np.cos(Om), np.sin(Om)
    ci, si = np.cos(inc), np.sin(inc)
    cw, sw = np.cos(w), np.sin(w)
    P = np.stack([cO * cw - sO * sw * ci, sO * cw + cO * sw * ci, sw * si], axis=-1)
    Q = np.stack([-cO * sw - sO * cw * ci, -sO * sw + cO * cw * ci, cw * si], axis=-1)
    return P, Q


def state_from_elements(elements_deg, t_days, epoch_mjd):
    """Heliocentric state (km, km/s) of body/bodies on fixed Keplerian orbits.

    elements_deg : (6,) or (N,6) [a AU, e, i, Om, w, M0] (deg)
    t_days       : scalar or array of days after mission t0 (MJD 62502.0); broadcast against N.
    Returns r, v with shape broadcast(N, t) x 3.
    """
    el = np.atleast_2d(np.asarray(elements_deg, float))
    t_days = np.asarray(t_days, float)
    a = el[:, 0] * AU
    e = el[:, 1]
    M0 = np.radians(el[:, 5])
    n = np.sqrt(MU_SUN / a ** 3)
    # broadcast: elements along leading axes, times along trailing
    shp = np.broadcast_shapes(a.shape + (1,) * t_days.ndim, (1,) + t_days.shape) if t_days.ndim else a.shape
    dt = (t_days + (MJD_T0 - epoch_mjd)) * DAY
    if t_days.ndim:
        a_ = a.reshape(a.shape + (1,) * t_days.ndim)
        e_ = e.reshape(a_.shape)
        M = M0.reshape(a_.shape) + n.reshape(a_.shape) * dt
    else:
        a_, e_ = a, e
        M = M0 + n * dt
    E = kepler_E(M, e_)
    cE, sE = np.cos(E), np.sin(E)
    b = a_ * np.sqrt(1 - e_ ** 2)
    r_p = a_ * (cE - e_)
    r_q = b * sE
    rmag = a_ * (1 - e_ * cE)
    fac = np.sqrt(MU_SUN * a_) / rmag
    v_p = -fac * sE
    v_q = fac * np.sqrt(1 - e_ ** 2) * cE
    P, Q = _rot_pqw(el)
    if t_days.ndim:
        P = P.reshape(P.shape[:1] + (1,) * t_days.ndim + (3,))
        Q = Q.reshape(P.shape)
    r = r_p[..., None] * P + r_q[..., None] * Q
    v = v_p[..., None] * P + v_q[..., None] * Q
    if el.shape[0] == 1 and np.ndim(elements_deg) == 1:
        r, v = r[0], v[0]
    return r, v


def earth_state(t_days):
    return state_from_elements(EARTH_ELEMENTS, t_days, MJD_EPOCH_EARTH)


def asteroid_state(elements_deg, t_days):
    return state_from_elements(elements_deg, t_days, MJD_EPOCH_AST)


# ---------------------------------------------------------------------------
# Universal-variable two-body propagator (for validating Lambert solutions)
# ---------------------------------------------------------------------------
def _stumpff(z):
    z = np.asarray(z, float)
    C = np.empty_like(z)
    S = np.empty_like(z)
    pos = z > 1e-8
    neg = z < -1e-8
    mid = ~(pos | neg)
    sz = np.sqrt(np.abs(z))
    C[pos] = (1 - np.cos(sz[pos])) / z[pos]
    S[pos] = (sz[pos] - np.sin(sz[pos])) / sz[pos] ** 3
    C[neg] = (np.cosh(sz[neg]) - 1) / (-z[neg])
    S[neg] = (np.sinh(sz[neg]) - sz[neg]) / sz[neg] ** 3
    C[mid] = 0.5 - z[mid] / 24
    S[mid] = 1 / 6 - z[mid] / 120
    return C, S


def propagate_universal(r0, v0, dt, mu=MU_SUN, itmax=60, tol=1e-11):
    """Vectorised universal-variable Kepler propagation (Vallado alg. 8). r0,v0 (N,3), dt (N,) s."""
    r0 = np.atleast_2d(np.asarray(r0, float))
    v0 = np.atleast_2d(np.asarray(v0, float))
    dt = np.broadcast_to(np.asarray(dt, float), r0.shape[:1]).astype(float)
    R0 = np.linalg.norm(r0, axis=1)
    V0sq = np.sum(v0 * v0, axis=1)
    rdotv = np.sum(r0 * v0, axis=1)
    alpha = 2.0 / R0 - V0sq / mu          # 1/a
    smu = np.sqrt(mu)
    chi = np.where(alpha > 1e-12, smu * dt * alpha,
                   np.where(alpha < -1e-12,
                            np.sign(dt) * np.sqrt(-1 / np.where(alpha < -1e-12, alpha, -1)) *
                            np.log(np.abs(-2 * mu * alpha * dt /
                                          (rdotv + np.sign(dt) * np.sqrt(-mu / np.where(alpha < -1e-12, alpha, -1)) *
                                           (1 - R0 * alpha)) + 1e-300)),
                            smu * dt / R0))
    for _ in range(itmax):
        z = alpha * chi ** 2
        C, S = _stumpff(z)
        chi2 = chi ** 2
        r = chi2 * C + rdotv / smu * chi * (1 - z * S) + R0 * (1 - z * C)
        F = rdotv / smu * chi2 * C + (1 - alpha * R0) * chi ** 3 * S + R0 * chi - smu * dt
        dchi = F / r
        chi = chi - dchi
        if np.max(np.abs(dchi)) < tol:
            break
    z = alpha * chi ** 2
    C, S = _stumpff(z)
    chi2 = chi ** 2
    r = chi2 * C + rdotv / smu * chi * (1 - z * S) + R0 * (1 - z * C)
    f = 1 - chi2 / R0 * C
    g = dt - chi ** 3 / smu * S
    fdot = smu / (r * R0) * chi * (z * S - 1)
    gdot = 1 - chi2 / r * C
    r1 = f[:, None] * r0 + g[:, None] * v0
    v1 = fdot[:, None] * r0 + gdot[:, None] * v0
    return r1, v1


def elements_from_state(r, v, mu=MU_SUN):
    """Return a [km], e, i [rad], q, Q [km] of heliocentric orbit(s); (N,3) inputs."""
    r = np.atleast_2d(r); v = np.atleast_2d(v)
    R = np.linalg.norm(r, axis=1)
    V2 = np.sum(v * v, axis=1)
    h = np.cross(r, v)
    hm = np.linalg.norm(h, axis=1)
    inc = np.arccos(np.clip(h[:, 2] / hm, -1, 1))
    evec = (np.cross(v, h) / mu) - r / R[:, None]
    e = np.linalg.norm(evec, axis=1)
    energy = V2 / 2 - mu / R
    with np.errstate(divide='ignore'):
        a = -mu / (2 * energy)
    p = hm ** 2 / mu
    q = p / (1 + e)
    with np.errstate(divide='ignore', invalid='ignore'):
        Q = np.where(e < 1, p / (1 - e), np.inf)
    return a, e, inc, q, Q


def kepler_H(N, e, tol=1e-13, itmax=60):
    """Solve hyperbolic Kepler equation e sinh H - H = N (vectorised Newton)."""
    N = np.asarray(N, float)
    e = np.broadcast_to(e, N.shape).astype(float)
    H = np.where(np.abs(N) > 6 * e, np.sign(N) * np.log(2 * np.abs(N) / e + 1.8), N / (e - 1 + 1e-300))
    H = np.clip(H, -50, 50)
    for _ in range(itmax):
        f = e * np.sinh(H) - H - N
        fp = e * np.cosh(H) - 1.0
        dH = f / fp
        dH = np.clip(dH, -2.0, 2.0)
        H = H - dH
        if np.max(np.abs(dH)) < tol:
            break
    return H


def propagate_kepler(r0, v0, dt, mu=MU_SUN):
    """Robust two-body propagation via classical elements (elliptic and hyperbolic), vectorised.
    r0, v0 (N,3) km, km/s; dt (N,) s.  Exact for all conics except the parabola (|e-1| < 1e-9 -> NaN)."""
    r0 = np.atleast_2d(np.asarray(r0, float)); v0 = np.atleast_2d(np.asarray(v0, float))
    dt = np.broadcast_to(np.asarray(dt, float), r0.shape[:1]).astype(float)
    R0 = np.linalg.norm(r0, axis=1)
    h = np.cross(r0, v0); hm = np.linalg.norm(h, axis=1)
    evec = np.cross(v0, h) / mu - r0 / R0[:, None]
    e = np.linalg.norm(evec, axis=1)
    p = hm ** 2 / mu
    a = p / (1 - e ** 2)
    # perifocal basis
    P = evec / e[:, None]
    W = h / hm[:, None]
    Q = np.cross(W, P)
    cnu0 = np.clip(np.sum(r0 * P, axis=1) / R0, -1, 1)
    snu0 = np.sum(r0 * Q, axis=1) / R0
    nu0 = np.arctan2(snu0, cnu0)
    ell = e < 1 - 1e-9
    hyp = e > 1 + 1e-9
    r1 = np.full_like(r0, np.nan); v1 = np.full_like(v0, np.nan)
    # elliptic
    if ell.any():
        ee, aa = e[ell], a[ell]
        E0 = 2 * np.arctan2(np.sqrt(1 - ee) * np.sin(nu0[ell] / 2), np.sqrt(1 + ee) * np.cos(nu0[ell] / 2))
        M0 = E0 - ee * np.sin(E0)
        n = np.sqrt(mu / aa ** 3)
        E = kepler_E(M0 + n * dt[ell], ee)
        cE, sE = np.cos(E), np.sin(E)
        rp = aa * (cE - ee); rq = aa * np.sqrt(1 - ee ** 2) * sE
        rm = aa * (1 - ee * cE)
        fac = np.sqrt(mu * aa) / rm
        vp = -fac * sE; vq = fac * np.sqrt(1 - ee ** 2) * cE
        r1[ell] = rp[:, None] * P[ell] + rq[:, None] * Q[ell]
        v1[ell] = vp[:, None] * P[ell] + vq[:, None] * Q[ell]
    if hyp.any():
        eh, ah = e[hyp], -a[hyp]           # ah > 0
        F0 = 2 * np.arctanh(np.clip(np.sqrt((eh - 1) / (eh + 1)) * np.tan(nu0[hyp] / 2), -1 + 1e-16, 1 - 1e-16))
        N0 = eh * np.sinh(F0) - F0
        n = np.sqrt(mu / ah ** 3)
        F = kepler_H(N0 + n * dt[hyp], eh)
        cF, sF = np.cosh(F), np.sinh(F)
        rp = ah * (eh - cF); rq = ah * np.sqrt(eh ** 2 - 1) * sF
        rm = ah * (eh * cF - 1)
        fac = np.sqrt(mu * ah) / rm
        vp = -fac * sF; vq = fac * np.sqrt(eh ** 2 - 1) * cF
        r1[hyp] = rp[:, None] * P[hyp] + rq[:, None] * Q[hyp]
        v1[hyp] = vp[:, None] * P[hyp] + vq[:, None] * Q[hyp]
    return r1, v1
