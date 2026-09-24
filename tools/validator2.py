#!/usr/bin/env python
"""validator2.py -- independent checker for CTOC14 Problem A submission files.

Written from CTOC14_problem.pdf (Sections 3, 4, 6, 7, 8) only; it shares no code with the
ctoc14 solver package (constants were re-typed from Tables 2/3 of the PDF).

Usage:
    python tools/validator2.py SUBMISSION.txt [--mea MEA.txt] [--json OUT.json] [--rtol 1e-12]
                               [--crosscheck N] [--quiet]

Checks (PDF Section 7):
  format     blank / '#' lines ignored; 15 columns; SC_ID consecutive from 1 and contiguous;
             first row Event=0, last row Event=4; Time strictly increasing; all Time in
             [0, 473364000]; Asteroid_ID 1..300 on Event=3 rows and 0 elsewhere;
             adjacent Event=1 rows of one thrust arc >= 8640 s apart (Event=3 rows ignored)
  launch     Event=0 position within 1 km of the Keplerian Earth; v_inf <= 4 km/s (+0.01 tol)
  thrust     |T| <= 0.5 N on every Event=1 row; interpolated |T(t)| <= 0.5 N at every instant
             (exact polynomial maximum on every segment + every RHS evaluation of the integrator);
             thrust columns exactly (0,0,0) on Event=0/2/3/4 rows
  dynamics   from EACH row's state integrate eq. (1) with the interpolated thrust to the next
             row (DOP853, rtol 1e-12) and compare: |dr| <= 1 km, |dv| <= 1 m/s, |dm| <= 0.01 kg
  mass       600 <= m0 <= 2000 on the Event=0 row; m >= 600 on every row and along the integration
  flyby      every Event=3 row: distance to the named asteroid at that Time <= 1000 km
  cost       N_covered (first detection per asteroid, t <= T_max), N_miss = 300 - N_covered,
             J = sum(1 + x + x^2) + N_miss with x = (m0 - 600)/1400 from the Event=0 mass column

Exit status: 0 VALID, 1 INVALID, 2 could not parse / run.
"""
import argparse
import json
import math
import os
import sys
import time

import numpy as np
from numpy.polynomial import polynomial as P
from scipy.integrate import solve_ivp

# ----------------------------------------------------------------------------------------------
# Constants (PDF Table 2 / Table 3 / eqs. (6), (14))
# ----------------------------------------------------------------------------------------------
MU_SUN = 1.32712440018e11          # km^3/s^2
AU_KM = 149597870.7                # km
G0 = 9.80665                       # m/s^2
ISP = 4000.0                       # s
T_MAX_N = 0.5                      # N
M_DRY = 600.0                      # kg
M0_MAX = 2000.0                    # kg
VINF_MAX = 4.0                     # km/s
D_FLYBY_MAX = 1000.0               # km
DT_SAMPLE_MIN = 8640.0             # s
T_MISSION = 473364000.0            # s  (15 * 365.25 d)
T0_MJD = 62502.0                   # 2030-01-01 00:00 UTC
EARTH_EPOCH_MJD = 60676.0          # 2025-01-01
AST_EPOCH_MJD = 61200.0            # 2026-06-09
N_ASTEROIDS = 300
FUEL_NORM = 1400.0                 # kg, eq. (15)
VE_MS = ISP * G0                   # exhaust velocity, m/s = 39226.6

EARTH_ELEMENTS = (1.0009175020, 0.017566762041, 0.002976847126,
                  189.953211282428, 273.196254000254, 357.4135031077)   # a[AU] e i Om om M0 [deg]

TOL_POS_KM = 1.0
TOL_VEL_KMS = 1e-3
TOL_MASS_KG = 0.01
TOL_VINF_KMS = 0.01
TOL_LAUNCH_POS_KM = 1.0
EPS_T = 1e-12                      # grace on |T| <= 0.5 N for decimal-rounding noise
EPS_DT = 1e-6                      # grace on the 8640 s spacing (float noise in the Time column)

VALID_EVENTS = (0, 1, 2, 3, 4)


# ----------------------------------------------------------------------------------------------
# Kepler ephemerides (PDF Section 3.3, eqs. (4)-(5))
# ----------------------------------------------------------------------------------------------
def solve_kepler(M, e):
    """Eccentric anomaly from mean anomaly (rad), 0 <= e < 1. Newton with Halley-like safeguard."""
    M = math.fmod(M, 2.0 * math.pi)
    if M < 0.0:
        M += 2.0 * math.pi
    E = M + e * math.sin(M) if e < 0.8 else math.pi
    for _ in range(100):
        f = E - e * math.sin(E) - M
        fp = 1.0 - e * math.cos(E)
        dE = -f / fp
        E += dE
        if abs(dE) < 1e-15:
            break
    # one final polish
    f = E - e * math.sin(E) - M
    E -= f / (1.0 - e * math.cos(E))
    return E


def rotation_313(Om, inc, om):
    """R = Rz(Om) Rx(inc) Rz(om) mapping perifocal -> inertial."""
    cO, sO = math.cos(Om), math.sin(Om)
    ci, si = math.cos(inc), math.sin(inc)
    co, so = math.cos(om), math.sin(om)
    return np.array([
        [cO * co - sO * so * ci, -cO * so - sO * co * ci, sO * si],
        [sO * co + cO * so * ci, -sO * so + cO * co * ci, -cO * si],
        [so * si, co * si, ci],
    ])


class KeplerBody:
    """A body on a fixed heliocentric Keplerian orbit; elements (a AU, e, i, Om, om, M0 deg) at epoch_mjd."""

    def __init__(self, elements, epoch_mjd):
        a_au, e, i_deg, Om_deg, om_deg, M0_deg = elements
        self.a = a_au * AU_KM
        self.e = float(e)
        self.n = math.sqrt(MU_SUN / self.a ** 3)          # rad/s, eq. (5)
        self.M0 = math.radians(M0_deg)
        self.R = rotation_313(math.radians(Om_deg), math.radians(i_deg), math.radians(om_deg))
        self.t_offset = (T0_MJD - epoch_mjd) * 86400.0    # seconds from element epoch to t0
        self.b = self.a * math.sqrt(1.0 - self.e ** 2)

    def state(self, t):
        """Position (km) and velocity (km/s) at mission time t (s since t0)."""
        M = self.M0 + self.n * (t + self.t_offset)         # eq. (4), t - t_eph in seconds
        E = solve_kepler(M, self.e)
        cE, sE = math.cos(E), math.sin(E)
        r_pf = np.array([self.a * (cE - self.e), self.b * sE, 0.0])
        Edot = self.n / (1.0 - self.e * cE)
        v_pf = np.array([-self.a * sE * Edot, self.b * cE * Edot, 0.0])
        return self.R @ r_pf, self.R @ v_pf


def load_mea(path):
    bodies = {}
    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            s = raw.strip()
            if not s or s.startswith("#"):
                continue
            p = s.split()
            if len(p) < 7:
                raise ValueError("MEA.txt: bad line %r" % raw)
            aid = int(float(p[0]))
            bodies[aid] = KeplerBody(tuple(float(x) for x in p[1:7]), AST_EPOCH_MJD)
    if sorted(bodies) != list(range(1, N_ASTEROIDS + 1)):
        raise ValueError("MEA.txt: expected IDs 1..300, got %d bodies" % len(bodies))
    return bodies


# ----------------------------------------------------------------------------------------------
# Independent two-body propagation (universal variables) -- used only as a cross-check on coasts
# ----------------------------------------------------------------------------------------------
def _stumpff(psi):
    if psi > 1e-6:
        s = math.sqrt(psi)
        return (1.0 - math.cos(s)) / psi, (s - math.sin(s)) / (s ** 3)
    if psi < -1e-6:
        s = math.sqrt(-psi)
        return (1.0 - math.cosh(s)) / psi, (math.sinh(s) - s) / (s ** 3)
    # series
    c2 = 0.5 - psi / 24.0 + psi ** 2 / 720.0
    c3 = 1.0 / 6.0 - psi / 120.0 + psi ** 2 / 5040.0
    return c2, c3


def kepler_universal(r0, v0, dt, mu=MU_SUN):
    r0 = np.asarray(r0, float)
    v0 = np.asarray(v0, float)
    r0n = float(np.linalg.norm(r0))
    v0n2 = float(v0 @ v0)
    rv = float(r0 @ v0)
    sqmu = math.sqrt(mu)
    alpha = 2.0 / r0n - v0n2 / mu
    if abs(dt) < 1e-12:
        return r0.copy(), v0.copy()
    if alpha > 1e-12 / r0n:
        chi = sqmu * dt * alpha
    elif alpha < -1e-12 / r0n:
        a = 1.0 / alpha
        sgn = 1.0 if dt > 0 else -1.0
        num = -2.0 * mu * alpha * dt
        den = rv + sgn * math.sqrt(-mu * a) * (1.0 - r0n * alpha)
        chi = sgn * math.sqrt(-a) * math.log(num / den)
    else:
        h = np.cross(r0, v0)
        p = float(h @ h) / mu
        s = 0.5 * math.atan(1.0 / (3.0 * math.sqrt(mu / p ** 3) * dt))
        w = math.atan(math.tan(s) ** (1.0 / 3.0))
        chi = math.sqrt(p) * 2.0 / math.tan(2.0 * w)
    for _ in range(200):
        psi = chi * chi * alpha
        c2, c3 = _stumpff(psi)
        r = chi * chi * c2 + rv / sqmu * chi * (1.0 - psi * c3) + r0n * (1.0 - psi * c2)
        dchi = (sqmu * dt - chi ** 3 * c3 - rv / sqmu * chi * chi * c2 - r0n * chi * (1.0 - psi * c3)) / r
        chi += dchi
        if abs(dchi) < 1e-14 * max(1.0, abs(chi)):
            break
    psi = chi * chi * alpha
    c2, c3 = _stumpff(psi)
    f = 1.0 - chi * chi * c2 / r0n
    g = dt - chi ** 3 * c3 / sqmu
    r = f * r0 + g * v0
    rn = float(np.linalg.norm(r))
    gdot = 1.0 - chi * chi * c2 / rn
    fdot = sqmu / (rn * r0n) * chi * (psi * c3 - 1.0)
    v = fdot * r0 + gdot * v0
    return r, v


# ----------------------------------------------------------------------------------------------
# Submission parsing
# ----------------------------------------------------------------------------------------------
class Row:
    __slots__ = ("line", "sc", "event", "t", "r", "v", "m", "T", "aid", "src_line", "idx")

    def __init__(self, line, sc, event, t, r, v, m, T, aid, src_line):
        self.line, self.sc, self.event, self.t = line, sc, event, t
        self.r, self.v, self.m, self.T, self.aid = r, v, m, T, aid
        self.src_line = src_line
        self.idx = None


def _to_int(s, what, src_line, errors):
    try:
        v = float(s)
    except ValueError:
        errors.append("[format] file line %d: %s is not a number (%r)" % (src_line, what, s))
        return None
    if not v.is_integer():
        errors.append("[format] file line %d: %s must be an integer (%r)" % (src_line, what, s))
        return None
    return int(v)


def parse_submission(path, errors):
    rows = []
    with open(path, encoding="utf-8") as fh:
        for k, raw in enumerate(fh, start=1):
            s = raw.strip()
            if not s or s.startswith("#"):
                continue
            p = s.split()
            if len(p) != 15:
                errors.append("[format] file line %d: expected 15 columns, found %d" % (k, len(p)))
                continue
            line = _to_int(p[0], "Line", k, errors)
            sc = _to_int(p[1], "SC_ID", k, errors)
            ev = _to_int(p[2], "Event", k, errors)
            aid = _to_int(p[14], "Asteroid_ID", k, errors)
            try:
                vals = [float(x) for x in p[3:14]]
            except ValueError:
                errors.append("[format] file line %d: non-numeric state/thrust field" % k)
                continue
            if None in (line, sc, ev, aid):
                continue
            if not all(math.isfinite(x) for x in vals):
                errors.append("[format] file line %d: non-finite value" % k)
                continue
            rows.append(Row(line, sc, ev, vals[0], np.array(vals[1:4]), np.array(vals[4:7]),
                            vals[7], np.array(vals[8:11]), aid, k))
    return rows


# ----------------------------------------------------------------------------------------------
# Thrust arcs and the sliding-window Lagrange interpolant (PDF Section 6.1 rules 2-4)
# ----------------------------------------------------------------------------------------------
class Arc:
    def __init__(self):
        self.rows = []          # Event=1 sample rows (in order)

    @property
    def n(self):
        return len(self.rows)

    def finalize(self):
        self.t = np.array([r.t for r in self.rows])
        self.T = np.array([r.T for r in self.rows])       # (n,3) N
        self.first_idx = self.rows[0].idx
        self.last_idx = self.rows[-1].idx

    def window(self, t):
        """0-based sample indices whose Lagrange polynomial defines T at time t (t inside [t0, tn-1])."""
        n = self.n
        if n < 4:
            return list(range(n))
        j = int(np.searchsorted(self.t, t, side="right")) - 1
        j = min(max(j, 0), n - 2)
        l = min(max(j - 1, 0), n - 4)
        return [l, l + 1, l + 2, l + 3]

    def thrust_direct(self, t):
        """Direct Lagrange evaluation of the thrust vector (N) at mission time t."""
        if t < self.t[0] or t > self.t[-1]:
            return np.zeros(3)
        idx = self.window(t)
        out = np.zeros(3)
        for i in idx:
            Li = 1.0
            for k in idx:
                if k != i:
                    Li *= (t - self.t[k]) / (self.t[i] - self.t[k])
            out += Li * self.T[i]
        return out


def build_arcs(rows):
    """Maximal runs of Event=1 rows; Event=3 rows do not break a run; Event=0/2/4 end it."""
    arcs = []
    cur = None
    for r in rows:
        if r.event == 1:
            if cur is None:
                cur = Arc()
                arcs.append(cur)
            cur.rows.append(r)
        elif r.event == 3:
            continue
        else:
            cur = None
    for a in arcs:
        a.finalize()
    return arcs


def segment_polynomial(arc, tA, tB):
    """Coefficients (3 x deg+1) of the thrust vector on [tA, tB] as polynomials in u = (t-tA)/(tB-tA).

    Valid only when the whole segment lies in one interpolation interval of the arc (always true
    because every sample is itself a row)."""
    idx = arc.window(tA)
    dt = tB - tA
    coef = np.zeros((3, len(idx)))
    for i in idx:
        Li = np.array([1.0])
        for k in idx:
            if k != i:
                # (t - t_k)/(t_i - t_k) with t = tA + u*dt  ->  linear polynomial in u
                Li = P.polymul(Li, np.array([(tA - arc.t[k]), dt]) / (arc.t[i] - arc.t[k]))
        for c in range(3):
            coef[c, :len(Li)] += Li * arc.T[i][c]
    return coef


def poly_max_norm(coef):
    """Exact maximum of |T(u)| for u in [0, 1]; returns (max_norm, u_at_max)."""
    S = np.zeros(1)
    for c in range(3):
        S = P.polyadd(S, P.polymul(coef[c], coef[c]))
    dS = P.polyder(S)
    cands = [0.0, 1.0]
    if len(dS) > 1 and np.any(dS[1:] != 0.0):
        roots = P.polyroots(dS)
        for z in roots:
            if abs(z.imag) < 1e-9 and -1e-12 < z.real < 1.0 + 1e-12:
                cands.append(min(max(z.real, 0.0), 1.0))
    best, ubest = -1.0, 0.0
    for u in cands:
        val = math.sqrt(max(float(P.polyval(u, S)), 0.0))
        if val > best:
            best, ubest = val, u
    return best, ubest


# ----------------------------------------------------------------------------------------------
# Dynamics (PDF eqs. (1)-(3))
# ----------------------------------------------------------------------------------------------
class Segment:
    __slots__ = ("A", "B", "arc", "coef", "Tmax_exact", "u_at_max", "rhs_Tmax")

    def __init__(self, A, B, arc, coef):
        self.A, self.B, self.arc, self.coef = A, B, arc, coef
        self.Tmax_exact = 0.0
        self.u_at_max = 0.0
        self.rhs_Tmax = 0.0


def make_rhs(seg):
    tA = seg.A.t
    dt = seg.B.t - tA
    coef = seg.coef
    if coef is None:
        def rhs(t, y):
            r = y[0:3]
            rn = math.sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2])
            k = -MU_SUN / (rn * rn * rn)
            return [y[3], y[4], y[5], k * r[0], k * r[1], k * r[2], 0.0]
        return rhs
    cx, cy, cz = coef[0].tolist(), coef[1].tolist(), coef[2].tolist()
    deg = len(cx) - 1

    def rhs(t, y):
        u = (t - tA) / dt
        Tx = cx[deg]
        Ty = cy[deg]
        Tz = cz[deg]
        for k in range(deg - 1, -1, -1):
            Tx = Tx * u + cx[k]
            Ty = Ty * u + cy[k]
            Tz = Tz * u + cz[k]
        Tn = math.sqrt(Tx * Tx + Ty * Ty + Tz * Tz)
        if Tn > seg.rhs_Tmax:
            seg.rhs_Tmax = Tn
        r = y[0:3]
        m = y[6]
        rn = math.sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2])
        kg = -MU_SUN / (rn * rn * rn)
        inv = 1.0e-3 / m                       # N/kg = m/s^2 -> km/s^2
        return [y[3], y[4], y[5],
                kg * r[0] + Tx * inv, kg * r[1] + Ty * inv, kg * r[2] + Tz * inv,
                -Tn / VE_MS]
    return rhs


def integrate_segment(seg, rtol, atol):
    y0 = np.concatenate([seg.A.r, seg.A.v, [seg.A.m]])
    rhs = make_rhs(seg)
    sol = solve_ivp(rhs, (seg.A.t, seg.B.t), y0, method="DOP853", rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError("DOP853 failed on segment %d->%d: %s" % (seg.A.line, seg.B.line, sol.message))
    y = sol.y[:, -1]
    return y, sol.t.size


def rk4_segment(seg, nsteps):
    """Fixed-step classical RK4 re-integration (cross-check only)."""
    rhs = make_rhs(seg)
    y = np.concatenate([seg.A.r, seg.A.v, [seg.A.m]])
    t = seg.A.t
    h = (seg.B.t - seg.A.t) / nsteps
    for _ in range(nsteps):
        k1 = np.array(rhs(t, y))
        k2 = np.array(rhs(t + 0.5 * h, y + 0.5 * h * k1))
        k3 = np.array(rhs(t + 0.5 * h, y + 0.5 * h * k2))
        k4 = np.array(rhs(t + h, y + h * k3))
        y = y + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        t += h
    return y


# ----------------------------------------------------------------------------------------------
# Main validation
# ----------------------------------------------------------------------------------------------
def validate(path, mea_path, rtol=1e-12, crosscheck=20, quiet=False, log=print):
    t_start = time.time()
    errors, warnings = [], []
    report = {"file": path, "errors": errors, "warnings": warnings, "spacecraft": [], "flybys": []}

    rows = parse_submission(path, errors)
    if not rows:
        errors.append("[format] no data rows")
        report["verdict"] = "INVALID"
        return report
    for r in rows:
        if r.line <= 0:
            errors.append("[format] file line %d: Line must be a positive integer (%d)" % (r.src_line, r.line))
        if r.event not in VALID_EVENTS:
            errors.append("[format] Line %d: Event %d not in 0..4" % (r.line, r.event))
    lines = [r.line for r in rows]
    if lines != sorted(lines) or len(set(lines)) != len(lines):
        warnings.append("[format] Line column is not strictly increasing / unique")

    earth = KeplerBody(EARTH_ELEMENTS, EARTH_EPOCH_MJD)
    asteroids = load_mea(mea_path)

    # ---- group by spacecraft (contiguous blocks, ids consecutive from 1) --------------------------
    blocks = []
    for r in rows:
        if blocks and blocks[-1][0] == r.sc:
            blocks[-1][1].append(r)
        else:
            blocks.append((r.sc, [r]))
    ids = [b[0] for b in blocks]
    if ids != list(range(1, len(ids) + 1)):
        errors.append("[format] SC_ID sequence %s is not consecutive from 1 in contiguous blocks" % ids[:20])

    covered = {}                 # asteroid id -> (sc, line, t) of first counted detection
    sum_Ji = 0.0
    gmax = {"pos": (0.0, None), "vel": (0.0, None), "mass": (0.0, None)}
    gTs = (0.0, None)            # max sample |T|
    gTi = (0.0, None)            # max interpolated |T|
    n_seg_total = 0
    coast_xcheck = 0.0
    lagrange_xcheck = 0.0
    all_segments = []

    for sc, srows in blocks:
        for k, r in enumerate(srows):
            r.idx = k
        info = {"sc": sc, "rows": len(srows), "errors_before": len(errors)}
        tag = "SC %d" % sc

        # ---- structure ------------------------------------------------------------------------
        if srows[0].event != 0:
            errors.append("[format] %s: first row (Line %d) must be Event=0, found Event=%d" % (tag, srows[0].line, srows[0].event))
        if srows[-1].event != 4:
            errors.append("[format] %s: last row (Line %d) must be Event=4, found Event=%d" % (tag, srows[-1].line, srows[-1].event))
        for r in srows[1:]:
            if r.event == 0:
                errors.append("[format] %s: Line %d: extra Event=0 row (launch must be the first row only)" % (tag, r.line))
        for r in srows[:-1]:
            if r.event == 4:
                errors.append("[format] %s: Line %d: Event=4 row before the last row" % (tag, r.line))
        for a, b in zip(srows[:-1], srows[1:]):
            if not (b.t > a.t):
                errors.append("[format] %s: Time not strictly increasing at Line %d (t=%.6f) -> Line %d (t=%.6f)" % (tag, a.line, a.t, b.line, b.t))
        for r in srows:
            if r.t < 0.0 or r.t > T_MISSION:
                errors.append("[format] %s: Line %d: Time %.3f s outside [0, %.0f]" % (tag, r.line, r.t, T_MISSION))
            if r.event == 3:
                if not (1 <= r.aid <= N_ASTEROIDS):
                    errors.append("[flyby] %s: Line %d: Asteroid_ID %d not in 1..300" % (tag, r.line, r.aid))
            elif r.aid != 0:
                errors.append("[format] %s: Line %d: Asteroid_ID must be 0 on Event=%d rows (found %d)" % (tag, r.line, r.event, r.aid))
            Tn = float(np.linalg.norm(r.T))
            if r.event == 1:
                if Tn > T_MAX_N + EPS_T:
                    errors.append("[thrust] %s: Line %d: sample |T| = %.12f N > 0.5 N" % (tag, r.line, Tn))
                if Tn > gTs[0]:
                    gTs = (Tn, (sc, r.line))
            elif np.any(r.T != 0.0):
                errors.append("[thrust] %s: Line %d: Event=%d row must have thrust (0,0,0), found (%g, %g, %g)" % (tag, r.line, r.event, r.T[0], r.T[1], r.T[2]))
            if r.m < M_DRY:
                errors.append("[mass] %s: Line %d: m = %.6f kg < 600 kg" % (tag, r.line, r.m))

        # ---- launch -----------------------------------------------------------------------------
        L = srows[0]
        m0 = L.m
        if L.event == 0:
            rE, vE = earth.state(L.t)
            dpos = float(np.linalg.norm(L.r - rE))
            vinf = float(np.linalg.norm(L.v - vE))
            info.update({"t_launch": L.t, "launch_pos_err_km": dpos, "v_inf_kms": vinf})
            if dpos > TOL_LAUNCH_POS_KM:
                errors.append("[launch] %s: Line %d: position differs from Earth by %.6f km > 1 km" % (tag, L.line, dpos))
            if vinf > VINF_MAX + TOL_VINF_KMS:
                errors.append("[launch] %s: Line %d: v_inf = %.6f km/s > 4 km/s (+0.01 tol)" % (tag, L.line, vinf))
            if m0 > M0_MAX + 1e-9:
                errors.append("[mass] %s: Line %d: m0 = %.6f kg > 2000 kg" % (tag, L.line, m0))
            if m0 < M_DRY:
                errors.append("[mass] %s: Line %d: m0 = %.6f kg < 600 kg (negative fuel)" % (tag, L.line, m0))
        x = (m0 - M_DRY) / FUEL_NORM
        Ji = 1.0 + x + x * x
        sum_Ji += Ji
        info.update({"m0": m0, "J_i": Ji, "m_final": srows[-1].m, "min_row_mass": min(r.m for r in srows)})

        # ---- arcs & spacing ---------------------------------------------------------------------
        arcs = build_arcs(srows)
        for a in arcs:
            for p, q in zip(a.rows[:-1], a.rows[1:]):
                if q.t - p.t < DT_SAMPLE_MIN - EPS_DT:
                    errors.append("[format] %s: thrust samples Line %d -> Line %d only %.3f s apart (< 8640 s)" % (tag, p.line, q.line, q.t - p.t))
        info["arcs"] = len(arcs)
        info["samples"] = sum(a.n for a in arcs)
        arc_of_row = [None] * len(srows)
        for a in arcs:
            for k in range(a.first_idx, a.last_idx + 1):
                arc_of_row[k] = a

        # ---- segments -------------------------------------------------------------------------
        segs = []
        for k in range(len(srows) - 1):
            A, B = srows[k], srows[k + 1]
            arc = arc_of_row[k]
            coef = None
            if arc is not None and arc.first_idx <= k and k + 1 <= arc.last_idx and B.t > A.t:
                coef = segment_polynomial(arc, A.t, B.t)
                # self-check: polynomial equals direct Lagrange evaluation at 3 interior points
                for u in (0.2, 0.5, 0.8):
                    tt = A.t + u * (B.t - A.t)
                    dT = np.linalg.norm(P.polyval(u, coef.T) - arc.thrust_direct(tt))
                    lagrange_xcheck = max(lagrange_xcheck, float(dT))
            segs.append(Segment(A, B, arc, coef))
        all_segments.extend(segs)

        # exact interpolated-thrust maximum on every thrust segment
        for s in segs:
            if s.coef is not None:
                s.Tmax_exact, s.u_at_max = poly_max_norm(s.coef)
                if s.Tmax_exact > gTi[0]:
                    gTi = (s.Tmax_exact, (sc, s.A.line, s.A.t + s.u_at_max * (s.B.t - s.A.t)))
                if s.Tmax_exact > T_MAX_N + EPS_T:
                    errors.append("[thrust] %s: interpolated |T(t)| reaches %.9f N > 0.5 N at t = %.3f s (between Line %d and Line %d)"
                                  % (tag, s.Tmax_exact, s.A.t + s.u_at_max * (s.B.t - s.A.t), s.A.line, s.B.line))

        # dynamics: integrate every segment from the submitted state
        smax = {"pos": (0.0, None), "vel": (0.0, None), "mass": (0.0, None)}
        min_int_mass = m0
        atol = [1e-9] * 7
        for i, s in enumerate(segs):
            if s.B.t <= s.A.t:
                continue                                   # already reported as a format error
            y, nsteps = integrate_segment(s, rtol, atol)
            dr = float(np.linalg.norm(y[0:3] - s.B.r))
            dv = float(np.linalg.norm(y[3:6] - s.B.v))
            dm = abs(float(y[6]) - s.B.m)
            min_int_mass = min(min_int_mass, float(y[6]))
            if s.coef is None:
                rk, vk = kepler_universal(s.A.r, s.A.v, s.B.t - s.A.t)
                coast_xcheck = max(coast_xcheck, float(np.linalg.norm(rk - y[0:3])))
            for key, val in (("pos", dr), ("vel", dv), ("mass", dm)):
                if val > smax[key][0]:
                    smax[key] = (val, (s.A.line, s.B.line))
                if val > gmax[key][0]:
                    gmax[key] = (val, (sc, s.A.line, s.B.line))
            if dr > TOL_POS_KM:
                errors.append("[dynamics] %s: Line %d -> Line %d: position error %.6f km > 1 km" % (tag, s.A.line, s.B.line, dr))
            if dv > TOL_VEL_KMS:
                errors.append("[dynamics] %s: Line %d -> Line %d: velocity error %.6f m/s > 1 m/s" % (tag, s.A.line, s.B.line, dv * 1e3))
            if dm > TOL_MASS_KG:
                errors.append("[mass] %s: Line %d -> Line %d: mass error %.6f kg > 0.01 kg" % (tag, s.A.line, s.B.line, dm))
            if s.rhs_Tmax > T_MAX_N + EPS_T and s.Tmax_exact <= T_MAX_N + EPS_T:
                errors.append("[thrust] %s: Line %d -> Line %d: integrator saw |T| = %.9f N > 0.5 N" % (tag, s.A.line, s.B.line, s.rhs_Tmax))
            n_seg_total += 1
            if not quiet and n_seg_total % 1000 == 0:
                sys.stderr.write("  ... %d segments integrated (%.1f s)\n" % (n_seg_total, time.time() - t_start))
        if min_int_mass < M_DRY - TOL_MASS_KG:
            errors.append("[mass] %s: integrated mass drops to %.6f kg < 600 kg" % (tag, min_int_mass))
        info["min_integrated_mass"] = min_int_mass
        info["max_pos_err_km"] = smax["pos"][0]
        info["max_vel_err_kms"] = smax["vel"][0]
        info["max_mass_err_kg"] = smax["mass"][0]
        info["max_T_sample"] = max([float(np.linalg.norm(r.T)) for r in srows if r.event == 1], default=0.0)
        info["max_T_interp"] = max([s.Tmax_exact for s in segs], default=0.0)

        # ---- flybys -----------------------------------------------------------------------------
        for r in srows:
            if r.event != 3:
                continue
            if not (1 <= r.aid <= N_ASTEROIDS):
                continue
            ra, va = asteroids[r.aid].state(r.t)
            d = float(np.linalg.norm(r.r - ra))
            vrel = float(np.linalg.norm(r.v - va))
            status = "COUNTED"
            if d > D_FLYBY_MAX:
                errors.append("[flyby] %s: Line %d: distance to asteroid %d is %.3f km > 1000 km" % (tag, r.line, r.aid, d))
                status = "FAIL_DISTANCE"
            elif r.t < 0.0 or r.t > T_MISSION:
                status = "OUT_OF_WINDOW"
            elif r.aid in covered:
                status = "REPEAT"
            else:
                covered[r.aid] = (sc, r.line, r.t)
            report["flybys"].append({"sc": sc, "line": r.line, "t": r.t, "asteroid": r.aid, "d_km": d,
                                     "vrel_kms": vrel, "status": status})
        info["flyby_rows"] = sum(1 for r in srows if r.event == 3)
        info["t_end"] = srows[-1].t
        report["spacecraft"].append(info)

    # ---- cross-check the worst thrust segments with fixed-step RK4 ---------------------------------
    xc = []
    if crosscheck > 0:
        thrust_segs = [s for s in all_segments if s.coef is not None and s.B.t > s.A.t]
        thrust_segs.sort(key=lambda s: -s.Tmax_exact)
        pick = thrust_segs[:crosscheck]
        for s in pick:
            y_d, _ = integrate_segment(s, rtol, [1e-9] * 7)
            y_r = rk4_segment(s, 2000)
            xc.append(float(np.linalg.norm(y_d[0:3] - y_r[0:3])))
    report["crosscheck_rk4_max_km"] = max(xc) if xc else None

    N_cov = len(covered)
    N_miss = N_ASTEROIDS - N_cov
    J = sum_Ji + N_miss
    report.update({
        "n_rows": len(rows), "n_segments": n_seg_total, "runtime_s": time.time() - t_start,
        "N_covered": N_cov, "N_miss": N_miss, "sum_Ji": sum_Ji, "J": J,
        "max_pos_err_km": gmax["pos"][0], "max_pos_err_at": gmax["pos"][1],
        "max_vel_err_kms": gmax["vel"][0], "max_vel_err_at": gmax["vel"][1],
        "max_mass_err_kg": gmax["mass"][0], "max_mass_err_at": gmax["mass"][1],
        "max_T_sample": gTs[0], "max_T_sample_at": gTs[1],
        "max_T_interp": gTi[0], "max_T_interp_at": gTi[1],
        "coast_kepler_xcheck_km": coast_xcheck, "lagrange_poly_xcheck_N": lagrange_xcheck,
        "verdict": "VALID" if not errors else "INVALID",
    })
    return report


def print_report(rep, log=print):
    log("=" * 78)
    log("validator2  --  %s" % rep["file"])
    log("=" * 78)
    if "n_rows" not in rep:
        for e in rep["errors"]:
            log("  " + e)
        log("VERDICT: %s" % rep["verdict"])
        return
    log("rows: %d   segments integrated: %d   runtime: %.1f s" % (rep["n_rows"], rep["n_segments"], rep["runtime_s"]))
    for s in rep["spacecraft"]:
        log("SC %d: rows=%d arcs=%d samples=%d flyby_rows=%d  m0=%.6f kg  m_final=%.6f kg  J_i=%.6f"
            % (s["sc"], s["rows"], s.get("arcs", 0), s.get("samples", 0), s.get("flyby_rows", 0), s["m0"], s["m_final"], s["J_i"]))
        if "v_inf_kms" in s:
            log("      launch t=%.3f s  pos err=%.6f km  v_inf=%.6f km/s   end t=%.3f s   min row mass=%.6f kg"
                % (s["t_launch"], s["launch_pos_err_km"], s["v_inf_kms"], s["t_end"], s["min_row_mass"]))
        log("      max seg err: pos=%.6e km  vel=%.3e km/s  mass=%.3e kg   max|T| sample=%.9f interp=%.9f N"
            % (s["max_pos_err_km"], s["max_vel_err_kms"], s["max_mass_err_kg"], s["max_T_sample"], s["max_T_interp"]))
    log("global max errors: pos %.6e km at %s | vel %.3e km/s (%.3e m/s) at %s | mass %.3e kg at %s"
        % (rep["max_pos_err_km"], rep["max_pos_err_at"], rep["max_vel_err_kms"], rep["max_vel_err_kms"] * 1e3,
           rep["max_vel_err_at"], rep["max_mass_err_kg"], rep["max_mass_err_at"]))
    at = rep["max_T_interp_at"]
    at_s = "(sc %d, from Line %d, t=%.3f s)" % at if at else "-"
    log("max |T| at samples %.9f N at %s; max interpolated |T(t)| %.9f N %s"
        % (rep["max_T_sample"], rep["max_T_sample_at"], rep["max_T_interp"], at_s))
    log("cross-checks: DOP853 vs universal-variable Kepler on coasts %.3e km; polynomial vs direct Lagrange %.3e N; DOP853 vs RK4(2000 steps) on worst thrust segments %s"
        % (rep["coast_kepler_xcheck_km"], rep["lagrange_poly_xcheck_N"],
           ("%.3e km" % rep["crosscheck_rk4_max_km"]) if rep["crosscheck_rk4_max_km"] is not None else "n/a"))
    log("-" * 78)
    if rep["flybys"]:
        log("flybys (sc Line t[s] asteroid d[km] vrel[km/s] status):")
        for f in rep["flybys"]:
            log("  %2d %6d %16.3f %4d %12.3f %8.3f %s" % (f["sc"], f["line"], f["t"], f["asteroid"], f["d_km"], f["vrel_kms"], f["status"]))
        log("-" * 78)
    log("N_covered = %d   N_miss = %d   sum J_i = %.6f   J = %.6f" % (rep["N_covered"], rep["N_miss"], rep["sum_Ji"], rep["J"]))
    if rep["warnings"]:
        log("WARNINGS (%d):" % len(rep["warnings"]))
        for w in rep["warnings"]:
            log("  " + w)
    if rep["errors"]:
        log("ERRORS (%d):" % len(rep["errors"]))
        for e in rep["errors"][:200]:
            log("  " + e)
        if len(rep["errors"]) > 200:
            log("  ... %d more" % (len(rep["errors"]) - 200))
    log("VERDICT: %s" % rep["verdict"])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("submission")
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--mea", default=os.path.join(here, "..", "MEA.txt"))
    ap.add_argument("--json", default=None, help="write machine-readable report here")
    ap.add_argument("--rtol", type=float, default=1e-12)
    ap.add_argument("--crosscheck", type=int, default=20, help="re-integrate the N worst thrust segments with RK4")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)
    try:
        rep = validate(args.submission, args.mea, rtol=args.rtol, crosscheck=args.crosscheck, quiet=args.quiet)
    except Exception as exc:  # noqa: BLE001
        print("validator2: fatal: %s" % exc)
        return 2
    print_report(rep)
    if args.json:
        def _default(o):
            if isinstance(o, (np.floating,)):
                return float(o)
            if isinstance(o, (np.integer,)):
                return int(o)
            if isinstance(o, np.ndarray):
                return o.tolist()
            return str(o)
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump(rep, fh, indent=1, default=_default)
    return 0 if rep["verdict"] == "VALID" else 1


if __name__ == "__main__":
    sys.exit(main())
