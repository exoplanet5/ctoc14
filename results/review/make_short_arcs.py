"""Synthetic submission exercising thrust arcs with n < 4 samples (validator2 model only).

SC 1: PDF-example launch state at t=0  ->  Event=2 coast node  ->  1-sample arc (thrust only at an
instant, i.e. effectively coast)  ->  Event=2  ->  2-sample arc (linear interpolant)  ->  Event=2  ->
3-sample arc (quadratic)  ->  Event=2  ->  4-sample arc (cubic, single window)  ->  Event=4.
States are produced with validator2's own integrator at rtol 1e-13; no flybys are declared.
"""
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import validator2 as V  # noqa: E402

DAY = 86400.0
rows = []   # (event, t, r, v, m, T)


def add(event, t, r, v, m, T=(0.0, 0.0, 0.0)):
    rows.append((event, t, np.array(r, float), np.array(v, float), float(m), np.array(T, float)))


def propagate(r, v, m, tA, tB, arc):
    """Integrate with the arc interpolant (or coast) from tA to tB."""
    if arc is None:
        def rhs(t, y):
            rr = y[0:3]
            k = -V.MU_SUN / np.linalg.norm(rr) ** 3
            return [y[3], y[4], y[5], k * rr[0], k * rr[1], k * rr[2], 0.0]
    else:
        def rhs(t, y):
            T = arc.thrust_direct(t)
            Tn = np.linalg.norm(T)
            rr = y[0:3]
            k = -V.MU_SUN / np.linalg.norm(rr) ** 3
            inv = 1e-3 / y[6]
            return [y[3], y[4], y[5], k * rr[0] + T[0] * inv, k * rr[1] + T[1] * inv, k * rr[2] + T[2] * inv, -Tn / V.VE_MS]
    sol = solve_ivp(rhs, (tA, tB), np.concatenate([r, v, [m]]), method="DOP853", rtol=1e-13, atol=1e-10)
    y = sol.y[:, -1]
    return y[0:3], y[3:6], y[6]


class FakeRow:
    def __init__(self, t, T):
        self.t, self.T = t, np.array(T, float)


def make_arc(samples):
    a = V.Arc()
    a.rows = [FakeRow(t, T) for t, T in samples]
    a.t = np.array([r.t for r in a.rows])
    a.T = np.array([r.T for r in a.rows])
    return a


r = np.array([-19500328.647, 145818483.47, -7637.2048832])
v = np.array([-29.389477847, -4.316691297, -3.9417153392])
m = 1500.0
t = 0.0
add(0, t, r, v, m)

# coast 10 d -> Event=2
r, v, m = propagate(r, v, m, t, t + 10 * DAY, None); t += 10 * DAY
add(2, t, r, v, m)

# 1-sample arc: thrust only at the instant -> next row 3 d later is a coast
add(1, t + 1 * DAY, *propagate(r, v, m, t, t + 1 * DAY, None), T=(0.3, 0.1, -0.2)); t += 1 * DAY
r, v, m = rows[-1][2], rows[-1][3], rows[-1][4]
r, v, m = propagate(r, v, m, t, t + 3 * DAY, None); t += 3 * DAY
add(2, t, r, v, m)

# 2-sample arc (linear), spacing 1 d, magnitudes 0.2 -> 0.4 N with a 30 deg direction change
s2 = [(t + 2 * DAY, (0.2, 0.0, 0.0)), (t + 3 * DAY, (0.4 * np.cos(np.radians(30)), 0.4 * np.sin(np.radians(30)), 0.0))]
r, v, m = propagate(r, v, m, t, t + 2 * DAY, None); t += 2 * DAY
add(1, t, r, v, m, T=s2[0][1])
arc2 = make_arc(s2)
r, v, m = propagate(r, v, m, t, t + 1 * DAY, arc2); t += 1 * DAY
add(1, t, r, v, m, T=s2[1][1])
r, v, m = propagate(r, v, m, t, t + 5 * DAY, None); t += 5 * DAY
add(2, t, r, v, m)

# 3-sample arc (quadratic), spacing 0.5 d (43200 s >= 8640 s)
s3 = [(t + 1 * DAY, (0.0, 0.3, 0.1)), (t + 1.5 * DAY, (0.1, 0.45, 0.0)), (t + 2 * DAY, (0.2, 0.3, -0.1))]
r, v, m = propagate(r, v, m, t, t + 1 * DAY, None); t += 1 * DAY
add(1, t, r, v, m, T=s3[0][1])
arc3 = make_arc(s3)
for k in (1, 2):
    r, v, m = propagate(r, v, m, t, s3[k][0], arc3); t = s3[k][0]
    add(1, t, r, v, m, T=s3[k][1])
r, v, m = propagate(r, v, m, t, t + 4 * DAY, None); t += 4 * DAY
add(2, t, r, v, m)

# 4-sample arc (cubic; one window), spacing 1 d, half-sine magnitude 0.45 N peak
s4 = [(t + (1 + k) * DAY, tuple(0.45 * np.sin(np.pi * (k + 0.5) / 4) * np.array([np.cos(0.2 * k), np.sin(0.2 * k), 0.1]) / np.linalg.norm([np.cos(0.2 * k), np.sin(0.2 * k), 0.1]))) for k in range(4)]
r, v, m = propagate(r, v, m, t, s4[0][0], None); t = s4[0][0]
add(1, t, r, v, m, T=s4[0][1])
arc4 = make_arc(s4)
for k in (1, 2, 3):
    r, v, m = propagate(r, v, m, t, s4[k][0], arc4); t = s4[k][0]
    add(1, t, r, v, m, T=s4[k][1])
r, v, m = propagate(r, v, m, t, t + 20 * DAY, None); t += 20 * DAY
add(4, t, r, v, m)

out = os.path.join(HERE, "synth_v2_short_arcs.txt")
with open(out, "w", encoding="utf-8") as fh:
    fh.write("# validator2 synthetic test: arcs with n = 1, 2, 3, 4 samples (m0 = 1500 kg, no flybys)\n")
    for i, (ev, tt, rr, vv, mm, TT) in enumerate(rows, start=1):
        fh.write("%d 1 %d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g 0\n"
                 % (i, ev, tt, rr[0], rr[1], rr[2], vv[0], vv[1], vv[2], mm, TT[0], TT[1], TT[2]))
print("wrote", out, "rows", len(rows), "final mass %.6f" % rows[-1][4])
