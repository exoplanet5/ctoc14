#!/usr/bin/env python
"""Write results/viz/data/bodies.json (docs/viz_spec.md section 2) from MEA.txt + ctoc14/constants.py.

Run from the project root:  nice -n 5 ~/.venvs/astro313/bin/python results/viz/data/make_bodies.py
Reads only MEA.txt (via ctoc14.kepler.load_mea) and ctoc14.constants; writes only bodies.json next to itself.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]          # .../ctoc14
sys.path.insert(0, str(ROOT))
from ctoc14 import constants as C                   # noqa: E402
from ctoc14.kepler import load_mea                  # noqa: E402

OUT = Path(__file__).resolve().parent / 'bodies.json'
UNREACHABLE = [131, 144]                            # per docs/viz_spec.md section 2

def main():
    ids, elem = load_mea(ROOT / 'MEA.txt')
    assert len(ids) == 300 and list(ids) == list(range(1, 301)), 'MEA.txt must hold IDs 1..300 in order'
    a, e, i, Om, om, M0 = C.EARTH_ELEMENTS
    body = {
        'epoch_mjd_t0': C.T0_MJD,
        'units': {'length': 'AU', 'time': 'day', 'angle': 'deg'},
        'mu_sun_km3s2': C.MU,
        'au_km': C.AU,
        'earth': {'a': float(a), 'e': float(e), 'i': float(i), 'Om': float(Om), 'om': float(om), 'M0': float(M0),
                  'epoch_mjd': C.T_EPH_EARTH_MJD},
        # Illustrative geocentric mean lunar orbit, exactly as given in docs/viz_spec.md section 2.
        'moon': {'note': 'illustrative geocentric mean orbit, not part of the contest model',
                 'a_km': 384400, 'e': 0.0549, 'i': 5.145, 'Om': 125.08, 'om': 318.15, 'M0': 135.27,
                 'epoch_mjd': 51544.5, 'period_d': 27.321582,
                 'Om_rate_deg_per_day': -0.05295, 'om_rate_deg_per_day': 0.11140},
        'asteroids': [
            {'id': int(k), 'a': float(row[0]), 'e': float(row[1]), 'i': float(row[2]), 'Om': float(row[3]),
             'om': float(row[4]), 'M0': float(row[5]), 'epoch_mjd': C.T_EPH_AST_MJD}
            for k, row in zip(ids, elem)],
        'unreachable': UNREACHABLE,
    }
    OUT.write_text(json.dumps(body, indent=1))
    print(f'wrote {OUT} ({OUT.stat().st_size} bytes, {len(body["asteroids"])} asteroids)')

if __name__ == '__main__':
    main()
