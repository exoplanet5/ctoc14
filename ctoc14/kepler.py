"""Keplerian ephemerides (Earth + 300 asteroids) and two-body propagation, exactly as defined in the CTOC14 rules.

Time convention: `t` is seconds since t0 = MJD 62502.0 (2030-01-01 00:00 UTC) unless stated otherwise.
All vectors are heliocentric ecliptic J2000 (HEI), km and km/s.
"""
import numpy as np
from pathlib import Path
from .constants import MU, AU, EARTH_ELEMENTS, T_OFF_AST, T_OFF_EARTH

DEG = np.pi / 180.0

def load_mea(path=None):
    """Return (ids[int], elements[N,6]) with columns a[AU] e i Omega omega M0 [deg] (epoch MJD 61200)."""
    if path is None:
        path = Path(__file__).resolve().parents[1] / 'MEA.txt'
    dat = np.loadtxt(path)
    return dat[:, 0].astype(int), dat[:, 1:7]

def solve_kepler(M, e, tol=1e-13, maxit=60):
    """Solve E - e sin E = M for elliptical orbits (vectorised). M in rad (any range), e < 1."""
    M = np.asarray(M, dtype=float); e = np.asarray(e, dtype=float)
    M = np.mod(M, 2 * np.pi)
    # initial guess (Danby / high-e safe)
    E = np.where(e < 0.8, M + e * np.sin(M), np.pi * np.ones_like(M))
    for _ in range(maxit):
        f = E - e * np.sin(E) - M
        fp = 1.0 - e * np.cos(E)
        dE = -f / fp
        # Halley-type safeguard for large steps
        dE = np.clip(dE, -1.0, 1.0)
        E = E + dE
        if np.all(np.abs(dE) < tol):
            break
    return E

class Body:
    """A Keplerian body defined by classical elements at epoch; precomputes rotation and mean motion."""
    __slots__ = ('a', 'e', 'i', 'Om', 'w', 'M0', 'n', 'p', 'R', 't_off', 'sqrt_mu_p')
    def __init__(self, elem, t_off):
        a, e, i, Om, w, M0 = elem
        self.a = a * AU; self.e = e; self.i = i * DEG; self.Om = Om * DEG; self.w = w * DEG; self.M0 = M0 * DEG
        self.n = np.sqrt(MU / self.a ** 3)          # rad/s, exactly eq. (5)
        self.p = self.a * (1 - e * e)
        self.sqrt_mu_p = np.sqrt(MU / self.p)
        self.t_off = t_off                          # seconds from element epoch to t0
        cO, sO = np.cos(self.Om), np.sin(self.Om); ci, si = np.cos(self.i), np.sin(self.i); cw, sw = np.cos(self.w), np.sin(self.w)
        # perifocal -> inertial rotation matrix R3(-Om) R1(-i) R3(-w)
        self.R = np.array([[cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si],
                           [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si],
                           [sw * si, cw * si, ci]])
    def state(self, t):
        """Position (km) and velocity (km/s) at mission time t [s] (scalar or array). Returns (r[...,3], v[...,3])."""
        t = np.asarray(t, dtype=float)
        M = self.M0 + self.n * (t + self.t_off)
        E = solve_kepler(M, self.e)
        cE, sE = np.cos(E), np.sin(E)
        r = self.a * (1 - self.e * cE)
        x = self.a * (cE - self.e); y = self.a * np.sqrt(1 - self.e ** 2) * sE
        # velocity in perifocal frame
        fac = np.sqrt(MU * self.a) / r
        vx = -fac * sE; vy = fac * np.sqrt(1 - self.e ** 2) * cE
        rp = np.stack([x, y, np.zeros_like(x)], axis=-1); vp = np.stack([vx, vy, np.zeros_like(x)], axis=-1)
        return rp @ self.R.T, vp @ self.R.T

class Ephemeris:
    """Earth + asteroid ephemerides. `ast[k]` is asteroid with 1-based ID ids[k] (k = ID-1 for MEA.txt)."""
    def __init__(self, mea_path=None):
        self.ids, self.elem = load_mea(mea_path)
        self.earth = Body(EARTH_ELEMENTS, T_OFF_EARTH)
        self.ast = [Body(el, T_OFF_AST) for el in self.elem]
        self.n_ast = len(self.ast)
        # vectorised arrays over asteroids for batch evaluation
        self._a = np.array([b.a for b in self.ast]); self._e = np.array([b.e for b in self.ast])
        self._n = np.array([b.n for b in self.ast]); self._M0 = np.array([b.M0 for b in self.ast])
        self._R = np.array([b.R for b in self.ast])  # (N,3,3)
    def earth_state(self, t):
        return self.earth.state(t)
    def ast_state(self, k, t):
        """State of asteroid index k (0-based) at time(s) t."""
        return self.ast[k].state(t)
    def ast_states_at(self, idx, t):
        """States of asteroids idx (0-based int array) at times t (array, broadcast against idx) in one vectorised pass.
        Same formulas as Body.state; agrees with ast_state to rounding. Returns (r[...,3], v[...,3])."""
        idx, t = np.broadcast_arrays(np.asarray(idx, dtype=int), np.asarray(t, dtype=float))
        a = self._a[idx]; e = self._e[idx]
        M = self._M0[idx] + self._n[idx] * (t + T_OFF_AST)
        E = solve_kepler(M, e)
        cE, sE = np.cos(E), np.sin(E)
        r = a * (1 - e * cE)
        se = np.sqrt(1 - e ** 2)
        x = a * (cE - e); y = a * se * sE
        fac = np.sqrt(MU * a) / r
        vx = -fac * sE; vy = fac * se * cE
        rp = np.stack([x, y, np.zeros_like(x)], axis=-1); vp = np.stack([vx, vy, np.zeros_like(x)], axis=-1)
        R = self._R[idx]
        return np.einsum('...ij,...j->...i', R, rp), np.einsum('...ij,...j->...i', R, vp)
    def all_ast_states(self, t):
        """States of all asteroids at scalar time t. Returns (r[N,3], v[N,3])."""
        M = self._M0 + self._n * (t + T_OFF_AST)
        E = solve_kepler(M, self._e)
        cE, sE = np.cos(E), np.sin(E)
        r = self._a * (1 - self._e * cE)
        x = self._a * (cE - self._e); y = self._a * np.sqrt(1 - self._e ** 2) * sE
        fac = np.sqrt(MU * self._a) / r
        vx = -fac * sE; vy = fac * np.sqrt(1 - self._e ** 2) * cE
        rp = np.stack([x, y, np.zeros_like(x)], axis=-1); vp = np.stack([vx, vy, np.zeros_like(x)], axis=-1)
        return np.einsum('nij,nj->ni', self._R, rp), np.einsum('nij,nj->ni', self._R, vp)

# ---------------------------------------------------------------- two-body propagation (universal variables)
def _stumpff(z):
    z = np.asarray(z, dtype=float)
    C = np.empty_like(z); S = np.empty_like(z)
    pos = z > 1e-8; neg = z < -1e-8; mid = ~(pos | neg)
    sz = np.sqrt(np.abs(z))
    C[pos] = (1 - np.cos(sz[pos])) / z[pos]
    S[pos] = (sz[pos] - np.sin(sz[pos])) / sz[pos] ** 3
    C[neg] = (np.cosh(sz[neg]) - 1) / (-z[neg])
    S[neg] = (np.sinh(sz[neg]) - sz[neg]) / sz[neg] ** 3
    C[mid] = 0.5 - z[mid] / 24.0; S[mid] = 1.0 / 6.0 - z[mid] / 120.0
    return C, S

def propagate_twobody(r0, v0, dt, mu=MU, tol=1e-12, maxit=100):
    """Propagate state (r0[3], v0[3]) by dt seconds (scalar or array) under two-body gravity. Vectorised over dt.
    Returns (r[...,3], v[...,3]). Universal-variable formulation (Vallado alg. 8), robust for ellipse/parabola/hyperbola."""
    r0 = np.asarray(r0, dtype=float); v0 = np.asarray(v0, dtype=float); dt = np.asarray(dt, dtype=float)
    scalar = dt.ndim == 0; dt = np.atleast_1d(dt)
    r0n = np.linalg.norm(r0); v0n2 = v0 @ v0; rdotv = r0 @ v0
    alpha = 2.0 / r0n - v0n2 / mu                     # 1/a
    sqmu = np.sqrt(mu)
    # initial guess for chi
    chi = np.where(np.abs(alpha) > 1e-12, sqmu * dt * alpha, 0.0)
    par = np.abs(alpha) <= 1e-12
    if np.any(par):
        h = np.cross(r0, v0); p = (h @ h) / mu
        s = 0.5 * np.arctan(1.0 / (3.0 * np.sqrt(mu / p ** 3) * dt[par])); w = np.arctan(np.tan(s) ** (1 / 3))
        chi[par] = np.sqrt(p) * 2.0 / np.tan(2 * w)
    hyp = alpha < -1e-12
    if np.any(hyp):
        a = 1.0 / alpha
        sgn = np.sign(dt[hyp])
        chi[hyp] = sgn * np.sqrt(-a) * np.log(-2 * mu * alpha * dt[hyp] / (rdotv + sgn * np.sqrt(-mu * a) * (1 - r0n * alpha)))
    for _ in range(maxit):
        z = chi * chi * alpha
        C, S = _stumpff(z)
        r = chi * chi * C + rdotv / sqmu * chi * (1 - z * S) + r0n * (1 - z * C)
        F = rdotv / sqmu * chi * chi * C + (1 - alpha * r0n) * chi ** 3 * S + r0n * chi - sqmu * dt
        dchi = -F / r
        chi = chi + dchi
        if np.all(np.abs(dchi) < tol * np.maximum(1.0, np.abs(chi))):
            break
    z = chi * chi * alpha
    C, S = _stumpff(z)
    r = chi * chi * C + rdotv / sqmu * chi * (1 - z * S) + r0n * (1 - z * C)
    f = 1 - chi * chi / r0n * C; g = dt - chi ** 3 / sqmu * S
    gdot = 1 - chi * chi / r * C; fdot = sqmu / (r * r0n) * chi * (z * S - 1)
    rv = f[:, None] * r0 + g[:, None] * v0; vv = fdot[:, None] * r0 + gdot[:, None] * v0
    if scalar:
        return rv[0], vv[0]
    return rv, vv

def propagate_batch(r0, v0, dt, mu=MU, tol=1e-12, maxit=100):
    """Batched two-body propagation: row q propagates (r0[q], v0[q]) by dt[q]. r0, v0: (N,3); dt: (N,).
    Same universal-variable arithmetic as propagate_twobody (per row); agrees with it to rounding. Returns (r[N,3], v[N,3])."""
    r0 = np.atleast_2d(np.asarray(r0, dtype=float)); v0 = np.atleast_2d(np.asarray(v0, dtype=float)); dt = np.atleast_1d(np.asarray(dt, dtype=float))
    r0, v0 = np.broadcast_arrays(r0, v0); N = len(dt)
    if N == 0:
        return np.zeros((0, 3)), np.zeros((0, 3))
    r0n = np.linalg.norm(r0, axis=1); v0n2 = np.einsum('ij,ij->i', v0, v0); rdotv = np.einsum('ij,ij->i', r0, v0)
    alpha = 2.0 / r0n - v0n2 / mu
    sqmu = np.sqrt(mu)
    chi = np.where(np.abs(alpha) > 1e-12, sqmu * dt * alpha, 0.0)
    par = np.abs(alpha) <= 1e-12
    if np.any(par):
        h = np.cross(r0[par], v0[par]); p = np.einsum('ij,ij->i', h, h) / mu
        s = 0.5 * np.arctan(1.0 / (3.0 * np.sqrt(mu / p ** 3) * dt[par])); w = np.arctan(np.tan(s) ** (1 / 3))
        chi[par] = np.sqrt(p) * 2.0 / np.tan(2 * w)
    hyp = alpha < -1e-12
    if np.any(hyp):
        a = 1.0 / alpha[hyp]
        sgn = np.sign(dt[hyp])
        chi[hyp] = sgn * np.sqrt(-a) * np.log(-2 * mu * alpha[hyp] * dt[hyp] / (rdotv[hyp] + sgn * np.sqrt(-mu * a) * (1 - r0n[hyp] * alpha[hyp])))
    for _ in range(maxit):
        z = chi * chi * alpha
        C, S = _stumpff(z)
        r = chi * chi * C + rdotv / sqmu * chi * (1 - z * S) + r0n * (1 - z * C)
        F = rdotv / sqmu * chi * chi * C + (1 - alpha * r0n) * chi ** 3 * S + r0n * chi - sqmu * dt
        dchi = -F / r
        chi = chi + dchi
        if np.all(np.abs(dchi) < tol * np.maximum(1.0, np.abs(chi))):
            break
    z = chi * chi * alpha
    C, S = _stumpff(z)
    r = chi * chi * C + rdotv / sqmu * chi * (1 - z * S) + r0n * (1 - z * C)
    f = 1 - chi * chi / r0n * C; g = dt - chi ** 3 / sqmu * S
    gdot = 1 - chi * chi / r * C; fdot = sqmu / (r * r0n) * chi * (z * S - 1)
    return f[:, None] * r0 + g[:, None] * v0, fdot[:, None] * r0 + gdot[:, None] * v0
