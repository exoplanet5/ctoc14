"""Generate the VIEWER fixture: results/viz/data/fixture/{bodies.json, index.json, fleets/demo.json}.

Synthetic 2-craft, 200-day fleet built from Lambert arcs between Earth and MEA.txt asteroids, sampled with
ctoc14.kepler.propagate_twobody, so every flyby sample coincides with the asteroid position (docs/viz_spec.md section 2).
It is NOT a contest result. Run from the project root:
    nice -n 5 ~/.venvs/astro313/bin/python results/viz/data/fixture/make_fixture.py
"""
import json, sys, hashlib
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from ctoc14.constants import MU, AU, DAY, T0_MJD, T_EPH_AST_MJD, T_EPH_EARTH_MJD, EARTH_ELEMENTS, cost_sc
from ctoc14.kepler import Ephemeris, propagate_twobody
from ctoc14.lambert import lambert

OUT = Path(__file__).resolve().parent
eph = Ephemeris(ROOT / 'MEA.txt')
ids, elem = eph.ids, eph.elem
UNREACHABLE = [131, 144]


def sec(d):
    return d * DAY


def best_leg(r0, v0, t0d, t1d, exclude):
    """Cheapest Lambert leg from (r0, v0) at t0d to any asteroid at t1d; returns (ast_id, v1, v2, r2, v2ast, dv)."""
    rA, vA = eph.all_ast_states(sec(t1d))
    v1, v2 = lambert(np.broadcast_to(r0, rA.shape), rA, np.full(len(rA), sec(t1d - t0d)))
    dv = np.linalg.norm(v1 - v0, axis=1)
    dv[np.isnan(dv)] = np.inf
    for k in exclude:
        dv[k - 1] = np.inf
    k = int(np.argmin(dv))
    return ids[k], v1[k], v2[k], rA[k], vA[k], float(dv[k])


def sample_arc(r0, v0, ta, tb, grid):
    """Times in (ta, tb] on the 2-day grid (plus tb exactly) and positions in AU."""
    ts = [t for t in grid if ta < t < tb] + [tb]
    dt = np.array([sec(t - ta) for t in ts])
    r, _ = propagate_twobody(r0, v0, dt)
    return ts, (r / AU)


def build_craft(cid, tL, flyby_times, t_end, m0, exclude):
    rE, vE = eph.earth_state(sec(tL))
    grid = list(np.arange(tL, t_end + 1e-9, 2.0))
    t_out, xyz, flybys, log = [tL], [rE / AU], [], []
    r, v, t_prev = rE, vE, tL
    for tf in flyby_times:
        aid, v1, v2, rA, vA, dv = best_leg(r, v, t_prev, tf, exclude + [f['ast'] for f in flybys])
        log.append((aid, tf, dv))
        ts, ps = sample_arc(r, v1, t_prev, tf, grid)
        t_out += ts; xyz += list(ps)
        flybys.append({'t': float(tf), 'ast': int(aid)})
        r, v, t_prev = rA, v2, tf          # continue on the arrival velocity (no gravity assist, just the fixture)
    ts, ps = sample_arc(r, v, t_prev, t_end, grid)
    t_out += ts; xyz += list(ps)
    t_out = np.array(t_out)
    assert np.all(np.diff(t_out) > 0), 'non-increasing t'
    # illustrative thrust arcs and mass history (the real files have Event=1 runs)
    thrust = [[tL + 1.0, tL + 24.0]] + [[f['t'] + 3.0, f['t'] + 28.0] for f in flybys[:-1]]
    thrust = [[a, min(b, t_end)] for a, b in thrust]
    burn_rate = 0.5 / (4000 * 9.80665) * DAY   # kg/day at full thrust (1.10 kg/d)
    m = np.full(len(t_out), float(m0))
    used = np.zeros(len(t_out))
    for a, b in thrust:
        used += np.clip(np.minimum(t_out, b) - a, 0, None) * burn_rate
    m = m - used
    return {'id': cid, 'm0': float(m0), 'launch_d': float(tL), 'end_d': float(t_end), 'n_flybys': len(flybys),
            'fuel_kg': float(m0 - 600.0), 't': [round(float(x), 4) for x in t_out],
            'xyz': [[float(f'{c:.6g}') for c in p] for p in xyz], 'm': [round(float(x), 3) for x in m],
            'thrust': thrust, 'flybys': flybys}, log


def main():
    craft, logs = [], []
    c1, l1 = build_craft(1, 10.0, [88.0, 164.0], 200.0, 950.0, UNREACHABLE); craft.append(c1); logs += l1
    c2, l2 = build_craft(2, 26.0, [130.0], 200.0, 880.0, UNREACHABLE + [f['ast'] for f in c1['flybys']]); craft.append(c2); logs += l2
    flown = [f['ast'] for c in craft for f in c['flybys']]
    # 5 asteroids in the fixture: the 3 flown, one unreachable (131) and one extra reachable target that stays grey
    extra = next(int(i) for i in ids if int(i) not in flown and int(i) not in UNREACHABLE and abs(elem[i - 1, 0] - 1.1) < 0.15)
    sel = flown + [131, extra]
    bodies = {'epoch_mjd_t0': T0_MJD, 'units': {'length': 'AU', 'time': 'day', 'angle': 'deg'},
              'mu_sun_km3s2': MU, 'au_km': AU,
              'earth': dict(zip(['a', 'e', 'i', 'Om', 'om', 'M0'], [float(x) for x in EARTH_ELEMENTS]), epoch_mjd=T_EPH_EARTH_MJD),
              'moon': {'note': 'illustrative geocentric mean orbit, not part of the contest model',
                       'a_km': 384400, 'e': 0.0549, 'i': 5.145, 'Om': 125.08, 'om': 318.15, 'M0': 135.27, 'epoch_mjd': 51544.5,
                       'period_d': 27.321582, 'Om_rate_deg_per_day': -0.05295, 'om_rate_deg_per_day': 0.11140},
              'asteroids': [dict(id=int(i), **dict(zip(['a', 'e', 'i', 'Om', 'om', 'M0'], [float(x) for x in elem[i - 1]])), epoch_mjd=T_EPH_AST_MJD)
                            for i in sel],
              'unreachable': [131],
              'note': 'FIXTURE: 5 asteroids only (3 flown by demo, 131 unreachable, %d extra). Not the contest catalogue.' % extra}
    (OUT / 'fleets').mkdir(parents=True, exist_ok=True)
    (OUT / 'bodies.json').write_text(json.dumps(bodies, indent=1))
    sumJi = float(sum(cost_sc(c['m0']) for c in craft))
    reachable = [a['id'] for a in bodies['asteroids'] if a['id'] not in bodies['unreachable']]
    missed = [a for a in reachable if a not in flown]
    fleet = {'key': 'demo', 'title': 'demo fixture: 2 craft, %d flybys, 200 d (synthetic Lambert arcs)' % len(flown),
             'source': 'results/viz/data/fixture/make_fixture.py', 'kind': 'fixture', 'n_craft': 2, 'covered': len(flown),
             'missed': missed + bodies['unreachable'], 'J_raw': round(sumJi + len(missed) + len(bodies['unreachable']), 6),
             'sumJi': round(sumJi, 6), 'date': '2026-09-24', 'craft': craft}
    (OUT / 'fleets' / 'demo.json').write_text(json.dumps(fleet, separators=(',', ':')))
    index = {'fleets': [{'key': 'demo', 'file': 'fleets/demo.json', 'title': fleet['title'], 'n_craft': 2,
                         'covered': len(flown), 'J_raw': fleet['J_raw'], 'shown': None, 'submitted': None, 'kind': 'fixture',
                         'note': 'Synthetic fixture for the viewer (not a contest result). Real data: results/viz/data/index.json.'}]}
    (OUT / 'index.json').write_text(json.dumps(index, indent=1))
    # verification: flyby sample vs asteroid position from the same elements
    worst = 0.0
    for c in craft:
        for f in c['flybys']:
            k = c['t'].index(f['t'])
            rA, _ = eph.ast_state(f['ast'] - 1, sec(f['t']))
            worst = max(worst, float(np.linalg.norm(np.array(c['xyz'][k]) - rA / AU)))
    for aid, tf, dv in logs:
        print(f'leg -> ast {aid:3d} at t={tf:6.1f} d  dv={dv:6.3f} km/s')
    print('asteroids in fixture:', sel, ' unreachable:', bodies['unreachable'])
    print('samples per craft:', [len(c['t']) for c in craft], ' max flyby position error [AU]: %.2e' % worst)
    for p in ['bodies.json', 'index.json', 'fleets/demo.json']:
        b = (OUT / p).read_bytes(); print(f'{p}: {len(b)} bytes md5 {hashlib.md5(b).hexdigest()}')


if __name__ == '__main__':
    main()
