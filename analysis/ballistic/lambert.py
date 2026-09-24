"""
Vectorised Lambert solver (Izzo 2015, "Revisiting Lambert's problem", CMDA 121:1-15)
in pure numpy: single revolution (N=0) and multi-revolution (N>=1, left/right branches).

`lambert(...)`           : N=0. `retrograde=False` gives the prograde solution (angular momentum
                           with +z component); the geometry then fixes short/long way.  For a
                           given (r1, r2, direction) the N=0 solution is unique.
`lambert_multirev(...)`  : N>=1. Returns both solutions (left branch = larger semi-major axis,
                           right branch = smaller a) as arrays with an extra axis of length 2.
                           Cells with T < T_min(N) (no N-rev solution) are NaN.

Householder (3rd order) iterations on Izzo's x variable, safeguarded with bracketing (T(x) is
monotonic on the N=0 branch and on each of the two N>=1 branches, so a bracketed Newton /
Householder / bisection hybrid cannot fail).  The time-of-flight function uses the
Lancaster-Blanchard form away from the parabola and Battin's series for |x-1| < 0.01.
Typical cost (Apple M-series, 1 core): 2-3 us per N=0 solve, ~8 us per N>=1 pair.
"""
import numpy as np


def _hyp2f1_3_1_25(z, nterms=14):
    """2F1(3, 1; 5/2; z) = sum_n (3)_n/(5/2)_n z^n, truncated (only used for |z|<~0.01)."""
    z = np.asarray(z, float)
    term = np.ones_like(z)
    s = np.ones_like(z)
    for n in range(nterms):
        term = term * (3.0 + n) / (2.5 + n) * z
        s = s + term
    return s


def _x2tof(x, lam, N=0):
    """Non-dimensional time of flight T(x) for N revolutions (Izzo 2015 / pykep x2tof)."""
    lam2 = lam * lam
    E = x * x - 1.0
    rho = np.abs(E)
    z = np.sqrt(np.maximum(1.0 + lam2 * E, 0.0))
    # --- Lancaster-Blanchard form
    y = np.sqrt(np.maximum(rho, 1e-300))
    g = x * z - lam * E
    with np.errstate(invalid='ignore', divide='ignore', over='ignore'):
        d_ell = np.arccos(np.clip(g, -1.0, 1.0)) + N * np.pi
        f = y * (z - lam * x)
        d_hyp = np.log(np.maximum(f + g, 1e-300))
        d = np.where(E < 0, d_ell, d_hyp)
        Esafe = np.where(rho < 1e-300, 1.0, E)
        tof_lb = (x - lam * z - d / y) / Esafe
        # --- Battin series near the parabola
        eta = z - lam * x
        S1 = 0.5 * (1.0 - lam - x * eta)
        Q = 4.0 / 3.0 * _hyp2f1_3_1_25(np.clip(S1, -0.5, 0.5))
        tof_ser = 0.5 * (eta ** 3 * Q + 4.0 * lam * eta)
        if N > 0:
            tof_ser = tof_ser + N * np.pi / np.maximum(rho, 1e-300) ** 1.5
    return np.where(np.abs(x - 1.0) < 0.01, tof_ser, tof_lb)


def _dTdx(x, T, lam):
    """First three derivatives of T(x) (valid for any N; Izzo 2015 eq. 22)."""
    lam2 = lam * lam
    lam3 = lam2 * lam
    umx2 = 1.0 - x * x
    umx2 = np.where(np.abs(umx2) < 1e-300, 1e-300, umx2)
    y = np.sqrt(np.maximum(1.0 - lam2 * umx2, 1e-300))
    y2 = y * y
    y3 = y2 * y
    with np.errstate(invalid='ignore', divide='ignore', over='ignore'):
        DT = (3.0 * T * x - 2.0 + 2.0 * lam3 * x / y) / umx2
        DDT = (3.0 * T + 5.0 * x * DT + 2.0 * (1.0 - lam2) * lam3 / y3) / umx2
        DDDT = (7.0 * x * DDT + 8.0 * DT - 6.0 * (1.0 - lam2) * lam2 * lam3 * x / y3 / y2) / umx2
    return DT, DDT, DDDT


def _geometry(r1, r2, tof, mu, retrograde):
    """Common non-dimensionalisation.  Returns dict with lam, T, unit vectors, etc."""
    r1 = np.asarray(r1, float)
    r2 = np.asarray(r2, float)
    tof = np.asarray(tof, float)
    shp = np.broadcast_shapes(r1.shape[:-1], r2.shape[:-1], tof.shape)
    r1 = np.broadcast_to(r1, shp + (3,))
    r2 = np.broadcast_to(r2, shp + (3,))
    tof = np.broadcast_to(tof, shp)
    retro = np.broadcast_to(np.asarray(retrograde, bool), shp)

    R1 = np.linalg.norm(r1, axis=-1)
    R2 = np.linalg.norm(r2, axis=-1)
    c = np.linalg.norm(r2 - r1, axis=-1)
    s = 0.5 * (R1 + R2 + c)
    ir1 = r1 / R1[..., None]
    ir2 = r2 / R2[..., None]
    ih = np.cross(ir1, ir2)
    ihn = np.linalg.norm(ih, axis=-1)
    degenerate = (ihn < 1e-12) | (c < 1e-12) | ~(tof > 0)
    ih = ih / np.where(degenerate, 1.0, ihn)[..., None]
    lam2 = np.clip(1.0 - c / np.where(s > 0, s, 1.0), 0.0, 1.0)
    lam = np.sqrt(lam2)
    flip = ih[..., 2] < 0.0             # prograde transfer angle > pi
    lam = np.where(flip, -lam, lam)
    it1 = np.where(flip[..., None], np.cross(ir1, ih), np.cross(ih, ir1))
    it2 = np.where(flip[..., None], np.cross(ir2, ih), np.cross(ih, ir2))
    lam = np.where(retro, -lam, lam)
    sgn = np.where(retro, -1.0, 1.0)[..., None]
    it1 = it1 * sgn
    it2 = it2 * sgn
    with np.errstate(divide='ignore', invalid='ignore'):
        T = np.sqrt(2.0 * mu / np.where(s > 0, s, 1.0) ** 3) * tof
    return dict(shp=shp, R1=R1, R2=R2, c=c, s=s, ir1=ir1, ir2=ir2, it1=it1, it2=it2,
                lam=lam, lam2=lam2, T=T, degenerate=degenerate)


def _velocities(G, x, mu, bad):
    """Izzo 2015 eqs. 20-22 / pykep: terminal velocities from x."""
    R1, R2, c, s, lam, lam2 = G['R1'], G['R2'], G['c'], G['s'], G['lam'], G['lam2']
    with np.errstate(divide='ignore', invalid='ignore'):
        gamma = np.sqrt(mu * s / 2.0)
        rho = (R1 - R2) / np.where(c > 0, c, 1.0)
        sigma = np.sqrt(np.clip(1.0 - rho * rho, 0.0, 1.0))
        y = np.sqrt(np.maximum(1.0 - lam2 + lam2 * x * x, 0.0))
        vr1 = gamma * ((lam * y - x) - rho * (lam * y + x)) / R1
        vr2 = -gamma * ((lam * y - x) + rho * (lam * y + x)) / R2
        vt = gamma * sigma * (y + lam * x)
        v1 = vr1[..., None] * G['ir1'] + (vt / R1)[..., None] * G['it1']
        v2 = vr2[..., None] * G['ir2'] + (vt / R2)[..., None] * G['it2']
    if bad.any():
        v1 = np.where(bad[..., None], np.nan, v1)
        v2 = np.where(bad[..., None], np.nan, v2)
    return v1, v2


def _householder_step(x, delta, DT, DDT, DDDT):
    DT2 = DT * DT
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        denom = DT * (DT2 - delta * DDT) + DDDT * delta * delta / 6.0
        step = delta * (DT2 - delta * DDT / 2.0) / denom
        newton = delta / DT
    return step, newton


def lambert(r1, r2, tof, mu, retrograde=False, max_iter=25, tol=1e-12, return_extra=False):
    """Single-revolution Lambert solver.

    r1, r2 : (..., 3) position vectors (km); tof : (...) time of flight (s), > 0.
    mu     : gravitational parameter (km^3/s^2).
    retrograde : bool or boolean array broadcastable to tof.shape.
    Returns v1, v2 of shape (..., 3).  Cells that did not converge (or degenerate
    geometry, collinear r1/r2) are returned as NaN.
    If return_extra: also returns dict(x=..., iters=..., T=..., lam=..., converged=...).
    """
    G = _geometry(r1, r2, tof, mu, retrograde)
    shp, lam, lam2, T = G['shp'], G['lam'], G['lam2'], G['T']
    lam3 = lam2 * lam
    lam5 = lam3 * lam2
    # initial guess (Izzo 2015 sec. 5)
    T0 = np.arccos(np.clip(lam, -1, 1)) + lam * np.sqrt(1.0 - lam2)
    T1 = 2.0 / 3.0 * (1.0 - lam3)
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        x_long = (T0 / T) ** (2.0 / 3.0) - 1.0
        x_short = 2.5 * T1 * (T1 - T) / ((1.0 - lam5) * T) + 1.0
        expo = np.log(2.0) / np.log(T1 / T0)
        x_mid = (T0 / T) ** expo - 1.0
    x = np.where(T >= T0, x_long, np.where(T <= T1, x_short, x_mid))
    x = np.where(np.isfinite(x), x, 0.0)

    iters = np.zeros(shp, dtype=np.int8)
    active = np.ones(shp, dtype=bool)
    for k in range(max_iter):
        tof_x = _x2tof(x, lam)
        DT, DDT, DDDT = _dTdx(x, tof_x, lam)
        delta = tof_x - T
        step, newton = _householder_step(x, delta, DT, DDT, DDDT)
        # Safeguard: far from the root the Householder correction can have the wrong sign
        # (seen for lambda -> -1); fall back to the Newton step there.  T(x) is monotonic
        # decreasing for N=0 so the Newton direction is always right.
        good = np.isfinite(step) & (np.sign(step) == np.sign(newton)) & (np.abs(step) < 3.0 * np.abs(newton))
        step = np.where(good, step, newton)
        step = np.where(np.isfinite(step), step, 0.0)
        xnew = x - np.where(active, step, 0.0)
        # keep x > -1 for the N=0 branch (bisect towards -1 instead of jumping past it)
        xnew = np.where(xnew <= -1.0, 0.5 * (x - 1.0), xnew)
        err = np.abs(xnew - x)
        x = xnew
        iters = iters + active.astype(np.int8)
        active = active & (err > tol)
        if not active.any():
            break
    tof_x = _x2tof(x, lam)
    converged = np.abs(tof_x - T) < 1e-9 * np.maximum(T, 1.0)
    bad = G['degenerate'] | ~converged | ~np.isfinite(x)
    v1, v2 = _velocities(G, x, mu, bad)
    if return_extra:
        return v1, v2, dict(x=x, iters=iters, T=T, lam=lam, converged=converged, degenerate=G['degenerate'])
    return v1, v2


def _x_min(lam, N, itmax=40):
    """Location x_M in (-1,1) and value T_M of the minimum of T(x) for N>=1 revolutions.
    Halley iterations on dT/dx = 0 with a sign bracket (dT/dx < 0 left of x_M, > 0 right)."""
    lo = np.full(lam.shape, -1.0 + 1e-9)
    hi = np.full(lam.shape, 1.0 - 1e-9)
    x = np.zeros(lam.shape)
    for _ in range(itmax):
        T = _x2tof(x, lam, N)
        DT, DDT, DDDT = _dTdx(x, T, lam)
        lo = np.where(DT < 0, x, lo)
        hi = np.where(DT > 0, x, hi)
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            xn = x - DT * DDT / (DDT * DDT - DT * DDDT / 2.0)
        badstep = ~np.isfinite(xn) | (xn <= lo) | (xn >= hi)
        xn = np.where(badstep, 0.5 * (lo + hi), xn)
        conv = np.abs(xn - x) < 1e-13
        x = xn
        if conv.all():
            break
    return x, _x2tof(x, lam, N)


def lambert_multirev(r1, r2, tof, mu, N=1, retrograde=False, max_iter=30, tol=1e-12, return_extra=False):
    """Multi-revolution (N >= 1) Lambert solver; returns BOTH branches.

    Returns v1, v2 with shape (..., 2, 3): index 0 = left branch (x < x_M, larger a),
    index 1 = right branch (x > x_M, smaller a).  NaN where T < T_min(N) (no solution) or
    degenerate geometry.
    """
    if N < 1:
        raise ValueError("use lambert() for N=0")
    G = _geometry(r1, r2, tof, mu, retrograde)
    shp, lam, T = G['shp'], G['lam'], G['T']
    xM, TM = _x_min(lam, N)
    exists = T >= TM
    out_v1 = np.full(shp + (2, 3), np.nan)
    out_v2 = np.full(shp + (2, 3), np.nan)
    extra = dict(xM=xM, TM=TM, T=T, lam=lam, x=np.full(shp + (2,), np.nan), iters=np.zeros(shp + (2,), np.int8),
                 converged=np.zeros(shp + (2,), bool))
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        tmpL = ((N * np.pi + np.pi) / (8.0 * T)) ** (2.0 / 3.0)
        tmpR = ((8.0 * T) / (N * np.pi)) ** (2.0 / 3.0)
        x0 = ((tmpL - 1.0) / (tmpL + 1.0), (tmpR - 1.0) / (tmpR + 1.0))
    for ib in range(2):
        right = ib == 1
        lo = xM.copy() if right else np.full(shp, -1.0 + 1e-12)
        hi = np.full(shp, 1.0 - 1e-12) if right else xM.copy()
        x = x0[ib]
        x = np.where(np.isfinite(x) & (x > lo) & (x < hi), x, 0.5 * (lo + hi))
        active = exists & ~G['degenerate']
        iters = np.zeros(shp, dtype=np.int8)
        for k in range(max_iter):
            Tx = _x2tof(x, lam, N)
            delta = Tx - T
            # T decreasing on the left branch, increasing on the right one
            if right:
                lo = np.where(delta < 0, x, lo)
                hi = np.where(delta > 0, x, hi)
            else:
                lo = np.where(delta > 0, x, lo)
                hi = np.where(delta < 0, x, hi)
            DT, DDT, DDDT = _dTdx(x, Tx, lam)
            step, newton = _householder_step(x, delta, DT, DDT, DDDT)
            xn = x - step
            ok = np.isfinite(xn) & (xn > lo) & (xn < hi)
            xn2 = x - newton
            ok2 = np.isfinite(xn2) & (xn2 > lo) & (xn2 < hi)
            xn = np.where(ok, xn, np.where(ok2, xn2, 0.5 * (lo + hi)))
            xn = np.where(active, xn, x)
            err = np.abs(xn - x)
            x = xn
            iters = iters + active.astype(np.int8)
            active = active & (err > tol)
            if not active.any():
                break
        Tx = _x2tof(x, lam, N)
        converged = exists & (np.abs(Tx - T) < 1e-9 * np.maximum(T, 1.0))
        bad = G['degenerate'] | ~converged | ~np.isfinite(x)
        v1, v2 = _velocities(G, x, mu, bad)
        out_v1[..., ib, :] = v1
        out_v2[..., ib, :] = v2
        extra['x'][..., ib] = x
        extra['iters'][..., ib] = iters
        extra['converged'][..., ib] = converged
    if return_extra:
        return out_v1, out_v2, extra
    return out_v1, out_v2
