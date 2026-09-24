"""Vectorised single-revolution Lambert solver (universal variables, Curtis Alg. 5.2 with a bracketed Newton/bisection
iteration in z). Inputs broadcast: r1, r2 (...,3) km, tof (...) s. Returns v1, v2 (...,3) km/s. NaN where no solution.
prograde=True selects the transfer with the ecliptic-z angular momentum component >= 0 (short/long way chosen automatically)."""
import numpy as np
from .constants import MU

def _stumpff(z):
    C = np.empty_like(z); S = np.empty_like(z)
    pos = z > 1e-6; neg = z < -1e-6; mid = ~(pos | neg)
    sz = np.sqrt(np.abs(z))
    C[pos] = (1 - np.cos(sz[pos])) / z[pos]; S[pos] = (sz[pos] - np.sin(sz[pos])) / sz[pos] ** 3
    C[neg] = (np.cosh(sz[neg]) - 1) / (-z[neg]); S[neg] = (np.sinh(sz[neg]) - sz[neg]) / sz[neg] ** 3
    C[mid] = 0.5 - z[mid] / 24 + z[mid] ** 2 / 720; S[mid] = 1 / 6 - z[mid] / 120 + z[mid] ** 2 / 5040
    return C, S

def lambert(r1, r2, tof, mu=MU, prograde=True, long_way=None, maxit=60, tol=1e-10):
    r1 = np.asarray(r1, float); r2 = np.asarray(r2, float); tof = np.asarray(tof, float)
    r1, r2, tof = np.broadcast_arrays(r1, r2, tof[..., None])
    tof = tof[..., 0]
    shape = tof.shape
    r1n = np.linalg.norm(r1, axis=-1); r2n = np.linalg.norm(r2, axis=-1)
    cr = np.cross(r1, r2)
    cosdnu = np.clip(np.einsum('...i,...i', r1, r2) / (r1n * r2n), -1, 1)
    dnu = np.arccos(cosdnu)
    if long_way is None:
        # prograde: if z-component of r1 x r2 < 0 the short way is retrograde -> take the long way
        lw = (cr[..., 2] < 0) if prograde else (cr[..., 2] >= 0)
    else:
        lw = np.broadcast_to(np.asarray(long_way, bool), shape)
    dnu = np.where(lw, 2 * np.pi - dnu, dnu)
    sindnu = np.sin(dnu)
    A = sindnu * np.sqrt(r1n * r2n / (1 - cosdnu + 1e-300))
    bad = (np.abs(A) < 1e-12) | (tof <= 0)
    A = np.where(bad, 1.0, A)
    def y_of(z):
        C, S = _stumpff(z)
        return r1n + r2n + A * (z * S - 1) / np.sqrt(C)
    def F_of(z):
        """F(z), y(z) and the Stumpff values C, S (evaluated once; same arithmetic as y_of)."""
        C, S = _stumpff(z)
        y = r1n + r2n + A * (z * S - 1) / np.sqrt(C)
        yy = np.where(y > 0, y, np.nan)
        return (yy / C) ** 1.5 * S + A * np.sqrt(yy) - np.sqrt(mu) * tof, y, C, S
    # bracket: z in (z_lo, z_hi); F is increasing in z. Lower bound: y(z) > 0 constraint (for A<0 / long way).
    z_lo = np.full(shape, -4 * np.pi ** 2); z_hi = np.full(shape, 4 * np.pi ** 2)  # single rev: z < (2pi)^2
    # ensure y>0 at z_lo by raising it
    for _ in range(80):
        y = y_of(z_lo); need = y <= 0
        if not np.any(need): break
        z_lo = np.where(need, z_lo + 0.5 * (z_hi - z_lo), z_lo)
    z = 0.5 * (z_lo + z_hi)
    for _ in range(maxit):
        F, y, C, S = F_of(z)
        # derivative dF/dz (Curtis 5.43)
        with np.errstate(all='ignore'):
            dF = np.where(np.abs(z) < 1e-6,
                          np.sqrt(2) / 40 * y ** 1.5 + A / 8 * (np.sqrt(y) + A * np.sqrt(1 / (2 * y))),
                          (y / C) ** 1.5 * (1 / (2 * z) * (C - 1.5 * S / C) + 0.75 * S ** 2 / C) + A / 8 * (3 * S / C * np.sqrt(y) + A * np.sqrt(C / y)))
        # update bracket
        pos = F > 0
        z_hi = np.where(pos, z, z_hi); z_lo = np.where(~pos & np.isfinite(F), z, z_lo)
        with np.errstate(all='ignore'):
            zn = z - F / dF
        inside = (zn > z_lo) & (zn < z_hi) & np.isfinite(zn)
        znew = np.where(inside, zn, 0.5 * (z_lo + z_hi))
        conv = np.abs(znew - z) < tol * np.maximum(1, np.abs(z))
        z = znew
        if np.all(conv | bad): break
    F, y, C, S = F_of(z)
    f = 1 - y / r1n; g = A * np.sqrt(y / mu); gdot = 1 - y / r2n
    v1 = (r2 - f[..., None] * r1) / g[..., None]
    v2 = (gdot[..., None] * r2 - r1) / g[..., None]
    ok = ~bad & np.isfinite(F) & (np.abs(F) < 1e-6 * np.sqrt(mu) * np.maximum(tof, 1.0)) & (y > 0)
    v1 = np.where(ok[..., None], v1, np.nan); v2 = np.where(ok[..., None], v2, np.nan)
    return v1, v2


# ================================================================================================ multi-revolution
# Same universal-variable formulation (z = chi^2 / a = (Delta E)^2 for ellipses) extended to z > (2 pi)^2.
# For nrev >= 1 the variable lies in ((2 pi nrev)^2, (2 pi (nrev+1))^2), where the time-of-flight curve T(z) blows up at
# both ends (C(z) -> 0) and has a single interior minimum T_min: two solutions ('left', z < z_min, and 'right', z > z_min)
# when tof > T_min, none otherwise. The left solution covers the smaller change of eccentric anomaly in the same time,
# i.e. it has the longer period / larger semi-major axis: 'left' = high-energy, 'right' = low-energy solution.
# For nrev = 0 the same machinery reproduces lambert(): T(z) is monotone on (z_min, (2 pi)^2) with z_min the lower end of the
# domain (y > 0), so only the 'right' branch exists.
BRANCH_LEFT, BRANCH_RIGHT = 0, 1
_FOUR_PI2 = 4.0 * np.pi ** 2
_Z_EDGE = 1e-6          # distance kept from the poles of T(z) at z = (2 pi k)^2


def _stumpff_mr(z):
    """Stumpff C, S; for z > 0 uses C = 2 sin^2(sqrt z / 2) / z, which stays accurate next to the poles z = (2 pi k)^2."""
    z = np.asarray(z, float)
    C = np.empty_like(z); S = np.empty_like(z)
    pos = z > 1e-6; neg = z < -1e-6; mid = ~(pos | neg)
    sz = np.sqrt(np.abs(z))
    sp = sz[pos]; sn = sz[neg]
    C[pos] = 2.0 * np.sin(0.5 * sp) ** 2 / z[pos]; S[pos] = (sp - np.sin(sp)) / sp ** 3
    C[neg] = (np.cosh(sn) - 1.0) / (-z[neg]); S[neg] = (np.sinh(sn) - sn) / sn ** 3
    C[mid] = 0.5 - z[mid] / 24 + z[mid] ** 2 / 720; S[mid] = 1 / 6 - z[mid] / 120 + z[mid] ** 2 / 5040
    return C, S


def _psi(z):
    """(z S - 1) / sqrt(C) in closed form: -sqrt2 cosh(sqrt(-z)/2) for z < 0, -sqrt2 cos(sqrt z/2) sgn(sin(sqrt z/2)) for z > 0.
    y(z) = r1 + r2 + A psi(z); the closed form avoids the 0/0 of the Stumpff expression at the poles."""
    z = np.asarray(z, float)
    s = np.sqrt(np.abs(z))
    out = np.empty_like(z)
    neg = z < 0
    out[neg] = -np.sqrt(2.0) * np.cosh(0.5 * s[neg])
    h = 0.5 * s[~neg]
    sg = np.where(np.sin(h) < 0, -1.0, 1.0)
    out[~neg] = -np.sqrt(2.0) * np.cos(h) * sg
    return out


class _MRGeom:
    """Broadcast geometry of a batch of Lambert problems (shared by the branch/nrev evaluations)."""
    __slots__ = ('r1', 'r2', 'tof', 'r1n', 'r2n', 'A', 'bad', 'shape', 'sqmu_tof', 'mu')

    def __init__(self, r1, r2, tof, mu, prograde, long_way):
        r1 = np.asarray(r1, float); r2 = np.asarray(r2, float); tof = np.asarray(tof, float)
        r1, r2, tof = np.broadcast_arrays(r1, r2, tof[..., None])
        tof = tof[..., 0]
        self.r1 = r1; self.r2 = r2; self.tof = tof; self.shape = tof.shape; self.mu = mu
        r1n = np.linalg.norm(r1, axis=-1); r2n = np.linalg.norm(r2, axis=-1)
        cr = np.cross(r1, r2)
        cosdnu = np.clip(np.einsum('...i,...i', r1, r2) / (r1n * r2n), -1, 1)
        dnu = np.arccos(cosdnu)
        if long_way is None:
            lw = (cr[..., 2] < 0) if prograde else (cr[..., 2] >= 0)
        else:
            lw = np.broadcast_to(np.asarray(long_way, bool), self.shape)
        dnu = np.where(lw, 2 * np.pi - dnu, dnu)
        A = np.sin(dnu) * np.sqrt(r1n * r2n / (1 - cosdnu + 1e-300))
        self.bad = (np.abs(A) < 1e-12) | (tof <= 0) | ~np.isfinite(tof)
        self.A = np.where(self.bad, 1.0, A)
        self.r1n = r1n; self.r2n = r2n
        self.sqmu_tof = np.sqrt(mu) * tof

    def y_of(self, z):
        return self.r1n + self.r2n + self.A * _psi(z)

    def T_of(self, z):
        """sqrt(mu) * tof(z) and its derivative dT/dz (Curtis 5.43), plus y, C, S. NaN where y <= 0."""
        C, S = _stumpff_mr(z)
        y = self.y_of(z)
        with np.errstate(all='ignore'):
            yy = np.where(y > 0, y, np.nan)
            T = (yy / C) ** 1.5 * S + self.A * np.sqrt(yy)
            A = self.A
            dT = np.where(np.abs(z) < 1e-6,
                          np.sqrt(2) / 40 * yy ** 1.5 + A / 8 * (np.sqrt(yy) + A * np.sqrt(1 / (2 * yy))),
                          (yy / C) ** 1.5 * (1 / (2 * z) * (C - 1.5 * S / C) + 0.75 * S ** 2 / C)
                          + A / 8 * (3 * S / C * np.sqrt(yy) + A * np.sqrt(C / yy)))
        return T, dT, y, C, S

    def velocities(self, z, ok):
        C, S = _stumpff_mr(z)
        y = self.y_of(z)
        with np.errstate(all='ignore'):
            f = 1 - y / self.r1n; g = self.A * np.sqrt(y / self.mu); gdot = 1 - y / self.r2n
            v1 = (self.r2 - f[..., None] * self.r1) / g[..., None]
            v2 = (gdot[..., None] * self.r2 - self.r1) / g[..., None]
        ok = ok & np.isfinite(v1).all(axis=-1) & np.isfinite(v2).all(axis=-1)
        v1 = np.where(ok[..., None], v1, np.nan); v2 = np.where(ok[..., None], v2, np.nan)
        return v1, v2, ok


def _mr_bounds(G, nrev):
    """Domain (z_lo, z_hi) of the nrev-revolution branch (nrev broadcast to G.shape). For nrev = 0 the lower end is the
    point y = 0 (A > 0: z_y = -(2 arccosh((r1+r2)/(sqrt2 A)))^2, where T = 0) or, for A < 0 (y > 0 everywhere), a value
    with T(z_lo) below the target time (extended downwards by factors of 4, at most to -4 pi^2 * 4^5)."""
    nrev = np.broadcast_to(np.asarray(nrev, int), G.shape)
    z_hi = _FOUR_PI2 * (nrev + 1) ** 2 - _Z_EDGE
    z_lo = np.where(nrev >= 1, _FOUR_PI2 * nrev ** 2 + _Z_EDGE, -_FOUR_PI2)
    if np.any(nrev == 0):
        with np.errstate(all='ignore'):
            rho = (G.r1n + G.r2n) / (np.sqrt(2.0) * G.A)              # >= 1 for A > 0
            z_y = -(2.0 * np.arccosh(np.maximum(rho, 1.0))) ** 2
        posA = (nrev == 0) & (G.A > 0)
        z_lo = np.where(posA, z_y + 1e-9 * np.maximum(1.0, np.abs(z_y)), z_lo)
        for _ in range(5):
            T, _, y, _, _ = G.T_of(z_lo)
            need = (nrev == 0) & ~posA & (T > G.sqmu_tof)
            if not np.any(need): break
            z_lo = np.where(need, 4.0 * z_lo, z_lo)
    return z_lo, z_hi, nrev


def _mr_minimum(G, z_lo, z_hi, nrev, maxit=60):
    """z_min, T_min of the unimodal T(z) on (z_lo, z_hi) for nrev >= 1 (bisection on the sign of dT/dz);
    nrev = 0 rows get z_min = z_lo, T_min = -inf (monotone branch, existence decided by the root residual)."""
    lo = z_lo.copy(); hi = z_hi.copy()
    multi = nrev >= 1
    for _ in range(maxit):
        z = 0.5 * (lo + hi)
        _, dT, _, _, _ = G.T_of(z)
        up = dT > 0
        hi = np.where(multi & up, z, hi); lo = np.where(multi & ~up, z, lo)
        if np.all(((hi - lo) <= 1e-13 * np.maximum(1.0, np.abs(hi))) | ~multi): break
    z_min = np.where(multi, 0.5 * (lo + hi), z_lo)
    T_min, _, _, _, _ = G.T_of(z_min)
    T_min = np.where(multi, T_min, -np.inf)       # nrev = 0: existence decided by the residual check of the root
    return z_min, T_min


def _mr_root(G, z_a, z_b, increasing, maxit=60, tol=1e-12):
    """Bracketed Newton for T(z) = sqrt(mu) tof on (z_a, z_b), T monotone there (increasing flag per row)."""
    lo = z_a.copy(); hi = z_b.copy()
    z = 0.5 * (lo + hi)
    F = np.full(G.shape, np.nan)
    for _ in range(maxit):
        T, dT, y, _, _ = G.T_of(z)
        F = T - G.sqmu_tof
        fin = np.isfinite(F)
        above = F > 0
        # root is to the left of z where (F > 0 and increasing) or (F < 0 and decreasing); y <= 0 (F NaN) only occurs
        # below the y = 0 point of the nrev = 0 domain, i.e. left of the root
        left = (above == increasing) & fin
        hi = np.where(left, z, hi); lo = np.where(~left, z, lo)
        with np.errstate(all='ignore'):
            zn = z - F / dT
        inside = (zn > lo) & (zn < hi) & np.isfinite(zn)
        znew = np.where(inside, zn, 0.5 * (lo + hi))
        conv = np.abs(znew - z) < tol * np.maximum(1.0, np.abs(z))
        z = znew
        if np.all(conv | G.bad): break
    T, _, y, _, _ = G.T_of(z)
    F = T - G.sqmu_tof
    ok = ~G.bad & np.isfinite(F) & (np.abs(F) < 1e-6 * np.maximum(G.sqmu_tof, np.sqrt(G.mu))) & (y > 0)
    return z, ok


def _branch_flags(branch, shape):
    if isinstance(branch, str):
        b = {'left': BRANCH_LEFT, 'l': BRANCH_LEFT, 'low': BRANCH_LEFT, 'right': BRANCH_RIGHT, 'r': BRANCH_RIGHT,
             'high': BRANCH_RIGHT}[branch.lower()]
        return np.broadcast_to(np.asarray(b == BRANCH_RIGHT), shape)
    b = np.asarray(branch)
    if b.dtype.kind in 'US':
        b = np.vectorize(lambda s: s.lower().startswith('r'))(b)
    return np.broadcast_to(b.astype(int) == BRANCH_RIGHT, shape)


def _mr_solve_branches(G, nrev, z_lo, z_hi, z_min, T_min, want_right):
    """Solve the requested branch per row given the minimum. Returns z, ok."""
    exists = ~G.bad & (G.sqmu_tof > T_min)          # NaN T_min compares False
    multi = nrev >= 1
    # right branch: (z_min, z_hi), T increasing; left branch (nrev >= 1 only): (z_lo, z_min), T decreasing
    za = np.where(want_right, z_min, z_lo); zb = np.where(want_right, z_hi, z_min)
    z, ok = _mr_root(G, za, zb, want_right)
    ok = ok & exists & (want_right | multi)
    return z, ok


def lambert_multirev(r1, r2, tof, nrev=1, branch='right', mu=MU, prograde=True, long_way=None, maxit=60, tol=1e-12,
                     return_info=False):
    """Multi-revolution Lambert solver (universal variables, vectorised). Inputs broadcast: r1, r2 (...,3) km, tof (...) s,
    nrev (int, broadcastable), branch 'left' / 'right' (or an array of BRANCH_LEFT / BRANCH_RIGHT, broadcastable).
    Returns v1, v2 (...,3) km/s, NaN where the requested solution does not exist (tof below the minimum of the nrev branch,
    'left' with nrev = 0, r1 parallel to r2, tof <= 0). nrev = 0 gives the single-revolution solution ('right' only),
    identical to lambert() to round-off. With return_info=True also returns dict(z, z_min, tof_min[s], ok).
    'left' (z < z_min): high-energy solution (larger semi-major axis); 'right' (z > z_min): low-energy solution."""
    G = _MRGeom(r1, r2, tof, mu, prograde, long_way)
    z_lo, z_hi, nrev = _mr_bounds(G, nrev)
    want_right = _branch_flags(branch, G.shape)
    z_min, T_min = _mr_minimum(G, z_lo, z_hi, nrev, maxit=maxit)
    z, ok = _mr_solve_branches(G, nrev, z_lo, z_hi, z_min, T_min, want_right)
    v1, v2, ok = G.velocities(z, ok)
    if return_info:
        return v1, v2, dict(z=np.where(ok, z, np.nan), z_min=z_min, tof_min=T_min / np.sqrt(mu), ok=ok)
    return v1, v2


def lambert_tof_min(r1, r2, nrev=1, mu=MU, prograde=True, long_way=None):
    """Minimum time of flight [s] of the nrev-revolution (nrev >= 1) transfer r1 -> r2 (below it no solution exists)."""
    G = _MRGeom(r1, r2, np.ones(np.broadcast_shapes(np.shape(r1)[:-1], np.shape(r2)[:-1])), mu, prograde, long_way)
    z_lo, z_hi, nrev = _mr_bounds(G, nrev)
    _, T_min = _mr_minimum(G, z_lo, z_hi, nrev)
    return T_min / np.sqrt(mu)


def lambert_all(r1, r2, tof, nrev_max=2, mu=MU, prograde=True, long_way=None):
    """All solutions of the batch: v1, v2 of shape (K, ..., 3) with K = 1 + 2 nrev_max, ordered
    (0,right), (1,left), (1,right), (2,left), (2,right), ...; NaN where absent. Also returns (nrevs[K], branches[K])."""
    G = _MRGeom(r1, r2, tof, mu, prograde, long_way)
    V1 = []; V2 = []; nrevs = []; branches = []
    for n in range(nrev_max + 1):
        z_lo, z_hi, nr = _mr_bounds(G, n)
        z_min, T_min = _mr_minimum(G, z_lo, z_hi, nr)
        for b in ((BRANCH_RIGHT,) if n == 0 else (BRANCH_LEFT, BRANCH_RIGHT)):
            want_right = np.broadcast_to(np.asarray(b == BRANCH_RIGHT), G.shape)
            z, ok = _mr_solve_branches(G, nr, z_lo, z_hi, z_min, T_min, want_right)
            v1, v2, ok = G.velocities(z, ok)
            V1.append(v1); V2.append(v2); nrevs.append(n); branches.append(b)
    return np.array(V1), np.array(V2), np.array(nrevs), np.array(branches)


def lambert_best(r1, r2, tof, v_prev, nrev_max=2, mu=MU, prograde=True, long_way=None):
    """Per row, the Lambert solution (over nrev 0..nrev_max and both branches) minimising the junction dv |v1 - v_prev|.
    Returns v1, v2 (...,3), nrev (...), branch (..., BRANCH_LEFT / BRANCH_RIGHT). Rows without any solution: NaN, nrev 0."""
    V1, V2, nrevs, branches = lambert_all(r1, r2, tof, nrev_max=nrev_max, mu=mu, prograde=prograde, long_way=long_way)
    v_prev = np.broadcast_to(np.asarray(v_prev, float), V1.shape[1:])
    dv = np.linalg.norm(V1 - v_prev[None], axis=-1)
    dv = np.where(np.isfinite(dv), dv, np.inf)
    k = np.argmin(dv, axis=0)
    v1 = np.take_along_axis(V1, k[None, ..., None], axis=0)[0]
    v2 = np.take_along_axis(V2, k[None, ..., None], axis=0)[0]
    return v1, v2, nrevs[k], branches[k]
