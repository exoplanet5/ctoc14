"""
CTOC14 Keplerian ephemeris module (numpy only).

Implements the contest's exact two-body propagation:
    n = sqrt(mu_sun / (a*AU)^3)   [rad/s]
    M(t) = M0 + n * (t - t_eph)   [t in seconds; MJD difference * 86400]
Frame: Heliocentric Ecliptic Inertial J2000 (HEI). Units km, km/s, s.

Epochs:
    asteroids : MJD 61200.0 (2026-06-09 00:00 UTC), elements from MEA.txt
    Earth     : MJD 60676.0 (2025-01-01 00:00 UTC), elements from Table 3 of the problem PDF
Mission window: t0 = MJD 62502.0 (2030-01-01), duration 15 yr = 4.73364e8 s = 5478.75 d.
"""
import numpy as np

MU_SUN = 1.32712440018e11      # km^3/s^2
AU = 149597870.7               # km
DAY = 86400.0                  # s
MJD_T0 = 62502.0               # mission start
T_MAX_S = 4.73364e8            # 15 yr in s
T_MAX_D = T_MAX_S / DAY        # 5478.75 d
MJD_END = MJD_T0 + T_MAX_D     # 67980.75
MJD_EPH_AST = 61200.0
MJD_EPH_EARTH = 60676.0
YEAR_D = 365.25

# Earth heliocentric ecliptic elements (a[AU], e, i, Omega, omega, M0 [deg]) at MJD 60676.0
EARTH_ELEMENTS = np.array([1.0009175020, 0.017566762041, 0.002976847126,
                           189.953211282428, 273.196254000254, 357.4135031077])


def load_mea(path):
    """Return (ids[int array], elements[N,6] = a[AU], e, i, Omega, omega, M0 [deg])."""
    rows = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            rows.append([float(x) for x in s.split()])
    arr = np.array(rows)
    return arr[:, 0].astype(int), arr[:, 1:7]


def mean_motion(a_au):
    """Mean motion in rad/s from a in AU (contest formula)."""
    return np.sqrt(MU_SUN / (np.asarray(a_au, dtype=float) * AU) ** 3)


def period_days(a_au):
    return 2.0 * np.pi / mean_motion(a_au) / DAY


def solve_kepler(M, e, tol=1e-14, maxit=60):
    """Vectorised Newton solver for E - e sin E = M (elliptic). M, e broadcastable."""
    M = np.mod(np.asarray(M, dtype=float), 2.0 * np.pi)
    e = np.asarray(e, dtype=float)
    # Danby starter, robust for high e
    E = M + 0.85 * e * np.sign(np.sin(M))
    E = np.where(np.sin(M) == 0, M, E)
    E, e, M = np.broadcast_arrays(E, e, M)
    E = E.copy()
    for _ in range(maxit):
        f = E - e * np.sin(E) - M
        fp = 1.0 - e * np.cos(E)
        dE = -f / fp
        E += dE
        if np.max(np.abs(dE)) < tol:
            break
    return E


def mean_anomaly_at(elements, t_eph_mjd, t_mjd):
    """M(t) [rad] for elements (...,6) at times t_mjd (broadcast)."""
    a = elements[..., 0]
    M0 = np.radians(elements[..., 5])
    n = mean_motion(a)
    return M0 + n * (np.asarray(t_mjd) - t_eph_mjd) * DAY


def rotation_matrix(i_deg, Om_deg, w_deg):
    """Perifocal -> HEI rotation matrices, shape (...,3,3)."""
    i = np.radians(i_deg); Om = np.radians(Om_deg); w = np.radians(w_deg)
    cO, sO = np.cos(Om), np.sin(Om)
    ci, si = np.cos(i), np.sin(i)
    cw, sw = np.cos(w), np.sin(w)
    R = np.empty(np.broadcast(cO, ci, cw).shape + (3, 3))
    R[..., 0, 0] = cO * cw - sO * sw * ci
    R[..., 0, 1] = -cO * sw - sO * cw * ci
    R[..., 0, 2] = sO * si
    R[..., 1, 0] = sO * cw + cO * sw * ci
    R[..., 1, 1] = -sO * sw + cO * cw * ci
    R[..., 1, 2] = -cO * si
    R[..., 2, 0] = sw * si
    R[..., 2, 1] = cw * si
    R[..., 2, 2] = ci
    return R


def state_from_E(elements, E):
    """Position [km] and velocity [km/s] in HEI given eccentric anomaly E (broadcast with elements[...,6])."""
    a = elements[..., 0] * AU
    e = elements[..., 1]
    cE, sE = np.cos(E), np.sin(E)
    r = a * (1.0 - e * cE)
    b = a * np.sqrt(1.0 - e ** 2)
    xp = a * (cE - e)
    yp = b * sE
    sqrt_mua = np.sqrt(MU_SUN * a)
    vxp = -sqrt_mua * sE / r
    vyp = sqrt_mua * np.sqrt(1.0 - e ** 2) * cE / r
    R = rotation_matrix(elements[..., 2], elements[..., 3], elements[..., 4])
    pos = R[..., :, 0] * xp[..., None] + R[..., :, 1] * yp[..., None]
    vel = R[..., :, 0] * vxp[..., None] + R[..., :, 1] * vyp[..., None]
    return pos, vel


def propagate(elements, t_eph_mjd, t_mjd):
    """Propagate elements (...,6) to times t_mjd (broadcastable with elements[...,0]).
    Returns pos (...,3) km, vel (...,3) km/s.
    Typical use: elements[:,None,:] with t_mjd[None,:] -> (N,T,3)."""
    elements = np.asarray(elements, dtype=float)
    M = mean_anomaly_at(elements, t_eph_mjd, t_mjd)
    E = solve_kepler(M, elements[..., 1])
    return state_from_E(elements, E)


def propagate_t_s(elements, t_eph_mjd, t_s):
    """Same as propagate but time given in seconds since mission start t0 (MJD 62502.0)."""
    return propagate(elements, t_eph_mjd, MJD_T0 + np.asarray(t_s) / DAY)


def earth_state(t_mjd):
    return propagate(EARTH_ELEMENTS, MJD_EPH_EARTH, t_mjd)


def true_to_mean_anomaly(nu, e):
    """nu [rad] -> M [rad] in [0, 2pi)."""
    E = 2.0 * np.arctan2(np.sqrt(1.0 - e) * np.sin(nu / 2.0), np.sqrt(1.0 + e) * np.cos(nu / 2.0))
    M = E - e * np.sin(E)
    return np.mod(M, 2.0 * np.pi)


def mjd_to_date(mjd):
    """MJD -> 'YYYY-MM-DD' (UTC, proleptic Gregorian)."""
    jd = float(mjd) + 2400000.5 + 0.5
    Z = int(jd); F = jd - Z
    if Z < 2299161:
        A = Z
    else:
        alpha = int((Z - 1867216.25) / 36524.25)
        A = Z + 1 + alpha - alpha // 4
    B = A + 1524
    C = int((B - 122.1) / 365.25)
    D = int(365.25 * C)
    Ee = int((B - D) / 30.6001)
    day = B - D - int(30.6001 * Ee) + F
    month = Ee - 1 if Ee < 14 else Ee - 13
    year = C - 4716 if month > 2 else C - 4715
    return f"{year:04d}-{month:02d}-{int(day):02d}"


def mjd_to_decimal_year(mjd):
    """Approximate decimal year (for plotting)."""
    return 2030.0 + (np.asarray(mjd) - MJD_T0) / YEAR_D
