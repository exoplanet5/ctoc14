"""Replica of the official CTOC14 submission checker (Section 7 of the problem statement).

Usage: python -m ctoc14.validator CTOC14_Result_TeamID.txt
Checks format/structure, launch, thrust, mass, dynamics consistency (integration of each row to the next with the
sliding-window Lagrange thrust), flyby distances, and computes N_covered, N_miss and J.
"""
import sys, numpy as np
from dataclasses import dataclass, field
from scipy.integrate import solve_ivp
from .constants import (MU, VE, TMAX, M_DRY, M0_MAX, VINF_MAX, D_FLYBY, DT_MIN_THRUST, T_MISSION,
                        TOL_POS, TOL_VEL, TOL_MASS, TOL_VINF, cost_sc)
from .kepler import Ephemeris, propagate_twobody
from .thrust import thrust_interp, arc_max_thrust

@dataclass
class Row:
    line: int; sc: int; event: int; t: float; r: np.ndarray; v: np.ndarray; m: float; T: np.ndarray; ast: int

@dataclass
class Report:
    ok: bool = True
    errors: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    n_sc: int = 0
    m0: list = field(default_factory=list)
    flybys: dict = field(default_factory=dict)      # asteroid id -> (sc, time, distance)
    max_pos_err: float = 0.0; max_vel_err: float = 0.0; max_mass_err: float = 0.0; max_thrust: float = 0.0
    J: float = float('nan'); n_covered: int = 0; n_miss: int = 300
    def err(self, msg): self.errors.append(msg); self.ok = False
    def warn(self, msg): self.warnings.append(msg)

def parse(path):
    rows = []
    with open(path, encoding='utf-8') as f:
        for ln, line in enumerate(f, 1):
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            p = s.replace('\t', ' ').split()
            if len(p) != 15:
                raise ValueError(f'line {ln}: expected 15 columns, got {len(p)}')
            rows.append(Row(int(p[0]), int(p[1]), int(p[2]), float(p[3]), np.array(p[4:7], float), np.array(p[7:10], float),
                            float(p[10]), np.array(p[11:14], float), int(p[14])))
    return rows

def _rhs_factory(ts, Ts):
    def rhs(t, y):
        r = y[0:3]; v = y[3:6]; m = y[6]
        T = thrust_interp(ts, Ts, t)
        Tn = np.linalg.norm(T)
        rn = np.linalg.norm(r)
        acc = -MU * r / rn ** 3 + (T / m) * 1e-3   # N/kg = m/s^2 -> km/s^2
        return np.concatenate([v, acc, [-Tn / (VE * 1e3)]])
    return rhs

def integrate_segment(r0, v0, m0, t0, t1, ts, Ts):
    """Integrate from (r0,v0,m0) at t0 to t1 with arc thrust samples (ts,Ts); ts=None means coast (analytic)."""
    if ts is None or len(ts) == 0:
        r1, v1 = propagate_twobody(r0, v0, t1 - t0)
        return r1, v1, m0
    sol = solve_ivp(_rhs_factory(ts, Ts), (t0, t1), np.concatenate([r0, v0, [m0]]), method='DOP853',
                    rtol=1e-12, atol=np.array([1e-6] * 3 + [1e-9] * 3 + [1e-7]), max_step=DT_MIN_THRUST / 4)
    y = sol.y[:, -1]
    return y[0:3], y[3:6], y[6]

def validate(path, eph=None, verbose=True):
    eph = eph or Ephemeris()
    rows = parse(path)
    rep = Report()
    by_sc = {}
    for r in rows:
        by_sc.setdefault(r.sc, []).append(r)
    ids = sorted(by_sc)
    if ids != list(range(1, len(ids) + 1)):
        rep.err(f'SC_ID must be consecutive from 1, got {ids}')
    # blocks must appear in ascending SC_ID order and be contiguous
    seq = [r.sc for r in rows]
    blocks = [seq[0]] + [b for a, b in zip(seq, seq[1:]) if b != a] if seq else []
    if blocks != ids:
        rep.err(f'spacecraft blocks must be contiguous and sorted by SC_ID; file order is {blocks}')
    rep.n_sc = len(ids)
    for sc in ids:
        R = by_sc[sc]
        if R[0].event != 0: rep.err(f'SC{sc}: first row must be Event=0')
        if R[-1].event != 4: rep.err(f'SC{sc}: last row must be Event=4')
        times = np.array([x.t for x in R])
        if np.any(np.diff(times) <= 0): rep.err(f'SC{sc}: Time not strictly increasing')
        if times.min() < 0 or times.max() > T_MISSION: rep.err(f'SC{sc}: event time outside [0, {T_MISSION}]')
        for x in R:
            if x.event != 1 and np.any(x.T != 0): rep.err(f'SC{sc} line {x.line}: non-thrust row must have zero thrust')
            if x.event == 1 and np.linalg.norm(x.T) > TMAX + 1e-12: rep.err(f'SC{sc} line {x.line}: |T|={np.linalg.norm(x.T):.4f} > 0.5 N')
            if x.event == 3 and not (1 <= x.ast <= 300): rep.err(f'SC{sc} line {x.line}: bad Asteroid_ID {x.ast}')
            if x.event != 3 and x.ast != 0: rep.err(f'SC{sc} line {x.line}: Asteroid_ID must be 0 on non-flyby rows')
            if x.m < M_DRY - TOL_MASS: rep.err(f'SC{sc} line {x.line}: mass {x.m} < 600 kg')
        # launch
        L = R[0]
        rE, vE = eph.earth_state(L.t)
        dpos = np.linalg.norm(L.r - rE); vinf = np.linalg.norm(L.v - vE)
        if dpos > TOL_POS: rep.err(f'SC{sc}: launch position error {dpos:.3f} km > 1 km')
        if vinf > VINF_MAX + TOL_VINF: rep.err(f'SC{sc}: launch v_inf {vinf:.4f} km/s > 4')
        if L.m > M0_MAX + 1e-9 or L.m < M_DRY: rep.err(f'SC{sc}: m0 = {L.m} outside [600, 2000]')
        rep.m0.append(L.m)
        # build arcs: maximal runs of Event=1 rows (Event=3 rows inside do not break the run)
        arcs = []  # list of (first_index, last_index) into R covering the run
        k = 0; n = len(R)
        while k < n:
            if R[k].event == 1:
                j = k; last1 = k
                while j + 1 < n and R[j + 1].event in (1, 3):
                    j += 1
                    if R[j].event == 1: last1 = j
                arcs.append((k, last1)); k = last1 + 1
            else:
                k += 1
        arc_of_row = {}
        for (a, b) in arcs:
            ts = np.array([R[i].t for i in range(a, b + 1) if R[i].event == 1])
            Ts = np.array([R[i].T for i in range(a, b + 1) if R[i].event == 1])
            if len(ts) > 1 and np.min(np.diff(ts)) < DT_MIN_THRUST - 1e-9:
                rep.err(f'SC{sc}: thrust samples spaced {np.min(np.diff(ts)):.1f} s < 8640 s in arc starting line {R[a].line}')
            tmax, tat = arc_max_thrust(ts, Ts)
            rep.max_thrust = max(rep.max_thrust, tmax)
            if tmax > TMAX + 1e-9:
                rep.err(f'SC{sc}: interpolated |T| = {tmax:.5f} N > 0.5 N at t={tat:.1f} in arc starting line {R[a].line}')
            for i in range(a, b + 1):
                arc_of_row[i] = (ts, Ts)
        # dynamics consistency row -> next row
        for i in range(n - 1):
            x0, x1 = R[i], R[i + 1]
            ts, Ts = arc_of_row.get(i, (None, None))
            # a row that is the last Event=1 of an arc followed by Event=3 rows still inside the arc is handled by arc_of_row;
            # after the arc's last sample thrust is zero anyway (outside [t1,tn]) so using the arc is exact.
            r1, v1, m1 = integrate_segment(x0.r, x0.v, x0.m, x0.t, x1.t, ts, Ts)
            ep = np.linalg.norm(r1 - x1.r); ev = np.linalg.norm(v1 - x1.v); em = abs(m1 - x1.m)
            rep.max_pos_err = max(rep.max_pos_err, ep); rep.max_vel_err = max(rep.max_vel_err, ev); rep.max_mass_err = max(rep.max_mass_err, em)
            if ep > TOL_POS: rep.err(f'SC{sc} line {x0.line}->{x1.line}: position error {ep:.3f} km > 1 km')
            if ev > TOL_VEL: rep.err(f'SC{sc} line {x0.line}->{x1.line}: velocity error {ev * 1e3:.3f} m/s > 1 m/s')
            if em > TOL_MASS: rep.err(f'SC{sc} line {x0.line}->{x1.line}: mass error {em:.4f} kg > 0.01 kg')
            if m1 < M_DRY - TOL_MASS: rep.err(f'SC{sc} line {x1.line}: integrated mass {m1:.3f} < 600 kg')
        # flybys
        for x in R:
            if x.event == 3 and 1 <= x.ast <= 300:
                ra, _ = eph.ast_state(x.ast - 1, x.t)
                d = np.linalg.norm(x.r - ra)
                if d > D_FLYBY:
                    rep.err(f'SC{sc} line {x.line}: flyby distance to asteroid {x.ast} is {d:.1f} km > 1000 km')
                elif x.t <= T_MISSION and x.ast not in rep.flybys:
                    rep.flybys[x.ast] = (sc, x.t, d)
    rep.n_covered = len(rep.flybys); rep.n_miss = 300 - rep.n_covered
    rep.J = float(np.sum(cost_sc(np.array(rep.m0)))) + rep.n_miss if rep.m0 else float('nan')
    if verbose:
        print(f'spacecraft: {rep.n_sc}, m0 = {rep.m0}')
        print(f'covered {rep.n_covered}, missed {rep.n_miss}, J = {rep.J:.4f}  (sum J_i = {rep.J - rep.n_miss:.4f})')
        print(f'max errors: pos {rep.max_pos_err:.4e} km, vel {rep.max_vel_err * 1e3:.4e} m/s, mass {rep.max_mass_err:.4e} kg, max |T| {rep.max_thrust:.5f} N')
        for e in rep.errors: print('ERROR:', e)
        for w in rep.warnings: print('WARN :', w)
        print('RESULT:', 'PASS' if rep.ok else 'FAIL')
    return rep

if __name__ == '__main__':
    rep = validate(sys.argv[1])
    sys.exit(0 if rep.ok else 1)
