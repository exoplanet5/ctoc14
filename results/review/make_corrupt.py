"""Create three corrupted copies of results/fleet_b300/sub_sc2.txt for validator cross-checking.

  sub_sc2_flyby_moved.txt   : Event=3 row (Line 204, asteroid 278) position shifted by +2000 km in x
  sub_sc2_thrust_0p6.txt    : Event=1 row (Line 50, |T| = 0.45 N) thrust vector scaled to |T| = 0.6 N
  sub_sc2_mass_shift.txt    : Event=1 row (Line 2599) mass column increased by 0.02 kg
Every other line is copied verbatim.
"""
import math
import os

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fleet_b300", "sub_sc2.txt")
OUT = os.path.dirname(os.path.abspath(__file__))


def fmt(x):
    return "%.17g" % x


def rewrite(target_line, mutate, out_name):
    hit = 0
    with open(SRC, encoding="utf-8") as fh, open(os.path.join(OUT, out_name), "w", encoding="utf-8") as out:
        for line in fh:
            s = line.strip()
            if s and not s.startswith("#"):
                p = s.split()
                if int(p[0]) == target_line:
                    p = mutate(p)
                    line = " ".join(p) + "\n"
                    hit += 1
            out.write(line)
    assert hit == 1, (out_name, hit)
    print("wrote", out_name)


def move_flyby(p):
    assert p[2] == "3", p
    p[4] = fmt(float(p[4]) + 2000.0)
    return p


def scale_thrust(p):
    assert p[2] == "1", p
    T = [float(p[11]), float(p[12]), float(p[13])]
    n = math.sqrt(sum(t * t for t in T))
    assert n > 0
    p[11:14] = [fmt(t * 0.6 / n) for t in T]
    return p


def shift_mass(p):
    assert p[2] == "1", p
    p[10] = fmt(float(p[10]) + 0.02)
    return p


rewrite(204, move_flyby, "sub_sc2_flyby_moved.txt")
rewrite(50, scale_thrust, "sub_sc2_thrust_0p6.txt")
rewrite(2599, shift_mass, "sub_sc2_mass_shift.txt")
