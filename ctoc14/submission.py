"""Submission writer: builds per-spacecraft event rows by integrating the exact validator dynamics, then writes the
15-column CTOC14 result file."""
import numpy as np
from .constants import VE, DAY
from .kepler import propagate_twobody
from .validator import integrate_segment

class SCTrajectory:
    """Accumulates rows (Event, t, r, v, m, T, ast) for one spacecraft. All states are produced by integration so the
    file is self-consistent with the validator."""
    def __init__(self, t_launch, r, v, m0):
        self.rows = [(0, float(t_launch), np.array(r, float), np.array(v, float), float(m0), np.zeros(3), 0)]
        self.t = float(t_launch); self.r = np.array(r, float); self.v = np.array(v, float); self.m = float(m0)
    def coast(self, dt, node=False):
        """Coast for dt seconds (analytic two-body). If node=True, append an Event=2 row at the end."""
        if dt > 0:
            self.r, self.v = propagate_twobody(self.r, self.v, dt); self.t += dt
        if node:
            self.rows.append((2, self.t, self.r.copy(), self.v.copy(), self.m, np.zeros(3), 0))
    def flyby(self, ast_id):
        """Record a flyby (Event=3) at the current state/time (must be called at the flyby time)."""
        self.rows.append((3, self.t, self.r.copy(), self.v.copy(), self.m, np.zeros(3), int(ast_id)))
    def coast_to(self, t, node=False):
        self.coast(t - self.t, node)
    def thrust_arc(self, ts, Ts, flybys=()):
        """Thrust arc with samples at absolute times ts (first must be >= current time; coast to it if later) and vectors Ts.
        flybys: iterable of (t_flyby, ast_id) inside the arc (or at its ends); Event=3 rows are inserted in time order."""
        ts = np.asarray(ts, float); Ts = np.asarray(Ts, float)
        if ts[0] > self.t + 1e-9:
            self.coast(ts[0] - self.t)
        assert abs(ts[0] - self.t) < 1e-6, 'arc must start at current time'
        events = [(t, 1, k) for k, t in enumerate(ts)] + [(tf, 3, a) for (tf, a) in flybys]
        events.sort(key=lambda e: (e[0], e[1]))
        for (t, ev, payload) in events:
            if t > self.t + 1e-12:
                self.r, self.v, self.m = integrate_segment(self.r, self.v, self.m, self.t, t, ts, Ts); self.t = t
            if ev == 1:
                self.rows.append((1, self.t, self.r.copy(), self.v.copy(), self.m, Ts[payload].copy(), 0))
            else:
                self.rows.append((3, self.t, self.r.copy(), self.v.copy(), self.m, np.zeros(3), int(payload)))
    def end(self):
        self.rows.append((4, self.t, self.r.copy(), self.v.copy(), self.m, np.zeros(3), 0))

def write_submission(path, trajectories, header=None):
    """trajectories: list of SCTrajectory (SC_ID assigned 1..N in order)."""
    line = 1
    with open(path, 'w', encoding='utf-8') as f:
        if header:
            for h in header.splitlines(): f.write('# ' + h + '\n')
        for sc, traj in enumerate(trajectories, 1):
            for (ev, t, r, v, m, T, ast) in traj.rows:
                f.write(f'{line} {sc} {ev} {t:.17g} {r[0]:.17g} {r[1]:.17g} {r[2]:.17g} {v[0]:.17g} {v[1]:.17g} {v[2]:.17g} '
                        f'{m:.17g} {T[0]:.17g} {T[1]:.17g} {T[2]:.17g} {ast}\n')
                line += 1
