#!/usr/bin/env python
"""Check that the Kepler propagation documented in docs/viz_kepler.md (the formulas the viewer's JavaScript must
implement) reproduces ctoc14.kepler for 10 asteroids at 5 epochs to better than 1e-6 AU.

The reference is ctoc14.kepler.Ephemeris (Body.state, i.e. the contest ephemeris). The candidate is `propagate()`
below, which is a literal transcription of viz_kepler.md: it reads bodies.json (degrees, AU, MJD epochs), takes time
in days since t0 and uses its own Newton solver. It shares no code with ctoc14.kepler; it only imports the constants.

Run from the project root:  nice -n 5 ~/.venvs/astro313/bin/python results/viz/data/check_bodies.py
Exit status 0 iff max error < 1e-6 AU (and the extra Earth / all-asteroid sweeps pass too).
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from ctoc14.constants import AU, DAY, MU, T0_MJD       # noqa: E402
from ctoc14.kepler import Ephemeris                    # noqa: E402

HERE = Path(__file__).resolve().parent
BODIES = json.loads((HERE / 'bodies.json').read_text())
TOL_AU = 1e-6

# ----------------------------------------------------------------------------------------------- viz_kepler.md, verbatim
def solve_kepler_newton(M, e, tol=1e-12, maxit=50):
    """E - e sin E = M by Newton's method (scalar, radians). M is reduced to [0, 2pi) first."""
    M = M % (2.0 * math.pi)
    E = M + e * math.sin(M) if e < 0.8 else math.pi
    for _ in range(maxit):
        dE = (E - e * math.sin(E) - M) / (1.0 - e * math.cos(E))
        E -= dE
        if abs(dE) < tol:
            break
    return E

def rotation(Om, i, om):
    """Perifocal -> ecliptic J2000 matrix R = R3(-Om) R1(-i) R3(-om), angles in radians (same as ctoc14/kepler.py)."""
    cO, sO = math.cos(Om), math.sin(Om)
    ci, si = math.cos(i), math.sin(i)
    cw, sw = math.cos(om), math.sin(om)
    return ((cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si),
            (sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si),
            (sw * si, cw * si, ci))

def propagate(el, t_day, mu=BODIES['mu_sun_km3s2'], au=BODIES['au_km'], t0_mjd=BODIES['epoch_mjd_t0']):
    """Heliocentric ecliptic position [AU] of a bodies.json element block at t_day days since t0 (steps 1-5 of the doc)."""
    d2r = math.pi / 180.0
    a_km = el['a'] * au
    e = el['e']
    n = math.sqrt(mu / a_km ** 3)                                   # rad/s        (step 1)
    dt_s = (t_day + (t0_mjd - el['epoch_mjd'])) * 86400.0           # s since element epoch
    M = el['M0'] * d2r + n * dt_s                                   # rad          (step 2)
    E = solve_kepler_newton(M, e)                                   #              (step 3)
    xp = a_km * (math.cos(E) - e)                                   # km, perifocal (step 4)
    yp = a_km * math.sqrt(1.0 - e * e) * math.sin(E)
    R = rotation(el['Om'] * d2r, el['i'] * d2r, el['om'] * d2r)     #              (step 5)
    return np.array([(R[k][0] * xp + R[k][1] * yp) / au for k in range(3)])

# ------------------------------------------------------------------------------------------------------------ checks
def main():
    eph = Ephemeris(ROOT / 'MEA.txt')
    ast = {b['id']: b for b in BODIES['asteroids']}
    ids = [1, 2, 4, 57, 100, 131, 144, 200, 250, 300]                 # incl. e=0.947 (144) and a=17.8 AU (131)
    epochs = [0.0, 1000.0, 2500.0, 4000.0, 5478.75]                    # days since t0 (t0 .. mission end)
    worst = (0.0, None, None)
    print(f'{"ast":>4} {"t[d]":>8} {"x_ref[AU]":>14} {"y_ref[AU]":>14} {"z_ref[AU]":>14} {"err[AU]":>10}')
    for k in ids:
        for t in epochs:
            r_ref = eph.ast_state(k - 1, t * DAY)[0] / AU
            r_doc = propagate(ast[k], t)
            err = float(np.linalg.norm(r_doc - r_ref))
            if err > worst[0]:
                worst = (err, k, t)
            print(f'{k:>4} {t:>8.2f} {r_ref[0]:>14.9f} {r_ref[1]:>14.9f} {r_ref[2]:>14.9f} {err:>10.2e}')
    print(f'\n10 asteroids x 5 epochs: max |r_doc - r_ref| = {worst[0]:.3e} AU (asteroid {worst[1]}, t = {worst[2]} d); '
          f'tolerance {TOL_AU:.0e} AU -> {"PASS" if worst[0] < TOL_AU else "FAIL"}')

    # Extra 1: Earth at the same epochs (bodies.json earth block vs Ephemeris.earth_state).
    e_worst = 0.0
    for t in epochs:
        r_ref = eph.earth_state(t * DAY)[0] / AU
        e_worst = max(e_worst, float(np.linalg.norm(propagate(BODIES['earth'], t) - r_ref)))
    print(f'Earth x 5 epochs: max error = {e_worst:.3e} AU -> {"PASS" if e_worst < TOL_AU else "FAIL"}')

    # Extra 2: every asteroid on a 50 d grid over the whole mission (Newton solver robustness for all e, all M).
    grid = np.arange(0.0, 5478.75 + 1e-9, 50.0)
    a_worst = (0.0, None, None)
    for t in grid:
        r_all = eph.all_ast_states(t * DAY)[0] / AU
        for k in range(1, 301):
            err = float(np.linalg.norm(propagate(ast[k], float(t)) - r_all[k - 1]))
            if err > a_worst[0]:
                a_worst = (err, k, float(t))
    print(f'All 300 asteroids x {len(grid)} epochs (50 d grid): max error = {a_worst[0]:.3e} AU '
          f'(asteroid {a_worst[1]}, t = {a_worst[2]} d) -> {"PASS" if a_worst[0] < TOL_AU else "FAIL"}')

    ok = worst[0] < TOL_AU and e_worst < TOL_AU and a_worst[0] < TOL_AU
    print('RESULT', 'PASS' if ok else 'FAIL', f'max_check_err_au={worst[0]:.6e}')
    return 0 if ok else 1

if __name__ == '__main__':
    sys.exit(main())
