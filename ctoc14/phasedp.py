"""Route design as a dynamic program on the phase plane (calibration: ctoc14/phasemodel.py).

A craft is a curve phi(t) = mean longitude relative to Earth [deg] whose slope is the drift rate s = n(a) - n_E
[deg/yr]; changing the slope costs K_DRIFT |ds| km/s (measured 0.028, corr 0.994 on the optimised 10-craft fleet).

EVENTS. Each pass of an asteroid through the craft-reachable ring (r in [0.80, 1.25] AU) gives ONE event, taken at the
sample of the pass with the smallest |z| (the cheapest plane geometry). An event carries (t, theta_rel, r, z).

REACHABILITY of event (t, theta, r, z) from craft state (phi, s) with orbit plane i_vec:
    a = a(s), e_r = (a - r) / a, e_t = (theta - phi) / 2 [rad]      (matching r and theta fixes the craft's e-vector)
    in plane:  |(e_r, e_t)| <= E_MAX
    out of plane: the craft's plane must contain the point: i_vec . u(theta) = z / r with u(theta) = (sin, -cos);
                  a mismatch dz is paid as a plane change, so events are admitted within Z_TOL and charged.
So around the craft's phase line there is a reachable band of half-width 2 sqrt(E_MAX^2 - e_r^2) [rad].

DP over (step, phi bin, slope bin), max-plus:
    1. slope change: L1 max-plus convolution (two sweeps), cost lam * K_TOTAL * K_DRIFT * |ds|
    2. prizes of events at this step reachable from (phi, s)
    3. drift: phi += s dt
This is the pricing problem of a set-covering column generation over routes (prize = LP dual of the asteroid).
"""
import numpy as np
from .constants import MU, AU, DAY, T_MISSION
from .phasemodel import K_DRIFT, K_PLANE, N_E, YR, K_TOTAL

E_MAX = 0.16              # max craft eccentricity (fleet: 0.003-0.161)
R_LO, R_HI = 0.80, 1.25   # AU, craft-reachable ring


def a_of_drift(s):
    """Semi-major axis [km] whose drift rate relative to Earth is s [deg/yr]."""
    return (MU / (np.radians(s + N_E) / YR) ** 2) ** (1.0 / 3.0)


def u_of(theta_deg):
    """Unit vector u with z = r * (i_vec . u) for a point at ecliptic longitude theta (small inclination)."""
    th = np.radians(theta_deg)
    return np.stack([np.sin(th), -np.cos(th)], axis=-1)


def events(eph, ids, dt=1.0 * DAY, r_lo=R_LO, r_hi=R_HI, z_max=0.12):
    """One event per ring pass: the sample with the smallest |z|. Returns a dict of arrays."""
    tt = np.arange(0.0, T_MISSION, dt)
    rE, _ = eph.earth_state(tt)
    thE = np.degrees(np.arctan2(rE[:, 1], rE[:, 0]))
    T, TH, R, Z, A, THabs = [], [], [], [], [], []
    for j, k in enumerate(ids):
        ra, _ = eph.ast_states_at(np.full(len(tt), k - 1), tt)
        r = np.linalg.norm(ra, axis=1) / AU; z = ra[:, 2] / AU
        m = (r >= r_lo) & (r <= r_hi)
        if not m.any():
            continue
        idx = np.nonzero(m)[0]
        brk = np.nonzero(np.diff(idx) > 1)[0]
        for run in np.split(idx, brk + 1):
            if len(run) == 0:
                continue
            b = run[np.argmin(np.abs(z[run]))]
            if abs(z[b]) > z_max:
                continue
            tha = np.degrees(np.arctan2(ra[b, 1], ra[b, 0]))
            T.append(tt[b]); TH.append((tha - thE[b] + 180) % 360 - 180); R.append(r[b]); Z.append(z[b])
            A.append(j); THabs.append(tha)
    return dict(t=np.array(T), th=np.array(TH), r=np.array(R), z=np.array(Z), ast=np.array(A),
                th_abs=np.array(THabs), ids=np.asarray(ids))


class PhaseDP:
    def __init__(self, ev, dphi=1.0, ds=1.5, s_max=90.0, dt_step=10.0 * DAY, i_vec=(0.0, 0.0), z_tol=0.010,
                 half_max=None):
        # half_max [deg] caps the phase slack credited to eccentricity. The instantaneous reachable set is
        # |(e_r, e_t)| <= E_MAX (half-width 2 sqrt(E_MAX^2 - e_r^2) rad), but the craft carries ONE eccentricity
        # vector, so that slack is not free at every event; half_max is the calibrated compromise.
        self.half_max = half_max; self.r_tol = None
        self.ev = ev; self.dphi = dphi; self.nphi = int(round(360.0 / dphi))
        self.s = np.arange(-s_max, s_max + 1e-9, ds); self.ns = len(self.s); self.ds = ds
        self.dt = dt_step; self.dt_yr = dt_step / YR
        self.nt = int(np.ceil(T_MISSION / dt_step))
        self.roll = np.rint(self.s * self.dt_yr / dphi).astype(int)
        self.a = a_of_drift(self.s) / AU
        self.set_plane(i_vec, z_tol)

    def set_plane(self, i_vec, z_tol=0.010):
        """Select the events this craft's orbit plane can reach: |z - r (i_vec . u)| <= z_tol."""
        ev = self.ev
        dz = ev['z'] - ev['r'] * (u_of(ev['th_abs']) @ np.asarray(i_vec, float))
        self.ok = np.abs(dz) <= z_tol
        self.dz = dz
        self.step_of = np.clip((ev['t'] / self.dt).astype(int), 0, self.nt - 1)
        self.cells = [[] for _ in range(self.nt)]
        for j in np.nonzero(self.ok)[0]:
            self.cells[self.step_of[j]].append(j)
        return int(self.ok.sum())

    def _bands(self, si, js):
        """phi-bin ranges reachable at slope bin si for the events js (empty when the radius is out of reach)."""
        ev = self.ev; a = self.a[si]
        e_r = (a - ev['r'][js]) / a
        good = np.abs(e_r) < (E_MAX if self.r_tol is None else self.r_tol)
        if not good.any():
            return np.zeros(0, int), np.zeros(0, int), np.zeros(0, int)
        js = np.asarray(js)[good]
        half = np.degrees(2.0 * np.sqrt(E_MAX ** 2 - e_r[good] ** 2))
        if self.half_max is not None:
            half = np.minimum(half, self.half_max)
        lo = np.ceil((ev['th'][js] - half) / self.dphi).astype(int)
        hi = np.floor((ev['th'][js] + half) / self.dphi).astype(int)
        return lo, hi, js

    def solve(self, prize, lam, dv_cap=None):
        """max sum(prize of events collected) - lam * dv over phase curves. Launch: phi = 0 at t = 0, free drift."""
        NEG = -1e18
        V = np.full((self.nphi, self.ns), NEG); V[0, :] = 0.0
        back = np.zeros((self.nt, self.nphi, self.ns), np.int16)
        kern = lam * K_TOTAL * K_DRIFT * self.ds
        for ti in range(self.nt):
            b = np.tile(np.arange(self.ns, dtype=np.int16), (self.nphi, 1))
            for i in range(1, self.ns):
                m = V[:, i - 1] - kern > V[:, i]
                V[:, i] = np.where(m, V[:, i - 1] - kern, V[:, i]); b[:, i] = np.where(m, b[:, i - 1], b[:, i])
            for i in range(self.ns - 2, -1, -1):
                m = V[:, i + 1] - kern > V[:, i]
                V[:, i] = np.where(m, V[:, i + 1] - kern, V[:, i]); b[:, i] = np.where(m, b[:, i + 1], b[:, i])
            back[ti] = b
            js = self.cells[ti]
            if js:
                for si in range(self.ns):
                    lo, hi, jj = self._bands(si, js)
                    for l_, h_, j_ in zip(lo, hi, jj):
                        p = prize[self.ev['ast'][j_]]
                        if p > 0:
                            V[np.arange(l_, h_ + 1) % self.nphi, si] += p
            for si in range(self.ns):
                if self.roll[si]:
                    V[:, si] = np.roll(V[:, si], self.roll[si])
        fi = np.unravel_index(np.argmax(V), V.shape)
        return float(V[fi]), self._trace(back, fi)

    def _trace(self, back, fi):
        phi, si = fi; path = []
        for ti in range(self.nt - 1, -1, -1):
            if self.roll[si]:
                phi = (phi - self.roll[si]) % self.nphi
            path.append((ti, phi, si))
            si = int(back[ti][phi, si])
        return path[::-1]

    def collected(self, path):
        """Events collected by a path: list of event indices (one per (asteroid, pass))."""
        out = []
        for ti, phi, si in path:
            js = self.cells[ti]
            if not js:
                continue
            lo, hi, jj = self._bands(si, js)
            for l_, h_, j_ in zip(lo, hi, jj):
                a, b = l_ % self.nphi, h_ % self.nphi
                if (a <= phi <= b) if a <= b else (phi >= a or phi <= b):
                    out.append(int(j_))
        return out

    def dv_of(self, path):
        s = np.array([self.s[si] for _, _, si in path])
        return K_TOTAL * K_DRIFT * np.abs(np.diff(s)).sum()
