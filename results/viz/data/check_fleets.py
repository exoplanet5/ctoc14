#!/usr/bin/env python
"""Check every fleet file listed in results/viz/data/index.json against docs/viz_spec.md section 2.

For every craft and every flyby: the flyby epoch must be an exact sample epoch in `t`, and the sample position there must
match the asteroid's position recomputed with ctoc14.kepler (Body) from the elements in bodies.json (or, when bodies.json
does not exist yet, straight from MEA.txt via ctoc14.kepler.Ephemeris) to better than TOL_AU. Also checked: t strictly
increasing, xyz / m lengths, thrust arcs inside [launch_d, end_d] and ordered, n_flybys and covered / missed consistency,
J_raw = sumJi + 300 - covered with J_i from the craft m0, file size <= SIZE_LIMIT, index.json entries.

Run from the project root:  nice -n 5 ~/.venvs/astro313/bin/python results/viz/data/check_fleets.py [key ...]
Exit status 0 iff every fleet passes. Prints one line per fleet with the MEASURED max flyby position error in AU.
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from ctoc14.constants import AU, DAY, T0_MJD, cost_sc                 # noqa: E402
from ctoc14.kepler import Body, Ephemeris                             # noqa: E402

HERE = Path(__file__).resolve().parent
TOL_AU = 2e-4
SIZE_LIMIT = 1_500_000
T_END_D = 5478.75

def asteroid_bodies():
    """{id: Body} from bodies.json (elements in deg, AU, epoch_mjd) or, if absent, from MEA.txt."""
    bj = HERE / 'bodies.json'
    if bj.exists():
        b = json.loads(bj.read_text())
        src = 'bodies.json'
        bodies = {a['id']: Body([a['a'], a['e'], a['i'], a['Om'], a['om'], a['M0']], (T0_MJD - a['epoch_mjd']) * DAY)
                  for a in b['asteroids']}
    else:
        e = Ephemeris(ROOT / 'MEA.txt')
        src = 'MEA.txt (bodies.json missing)'
        bodies = {int(k): body for k, body in zip(e.ids, e.ast)}
    return bodies, src

def check_fleet(entry, bodies):
    path = HERE / entry['file']
    size = path.stat().st_size
    fleet = json.loads(path.read_text())
    errs = []
    def bad(msg):
        errs.append(msg)
    if size > SIZE_LIMIT:
        bad(f'size {size} > {SIZE_LIMIT}')
    if fleet['key'] != entry['key']:
        bad(f"key {fleet['key']} != index {entry['key']}")
    max_err = 0.0; max_at = None; n_fb = 0; n_samples = 0; worst_step = 0.0
    targets = []
    for c in fleet['craft']:
        t = np.asarray(c['t'], float); n = len(t); n_samples += n
        if n < 2:
            bad(f"craft {c['id']}: {n} samples"); continue
        if not np.all(np.diff(t) > 0):
            bad(f"craft {c['id']}: t not strictly increasing at {int(np.argmin(np.diff(t)))}")
        worst_step = max(worst_step, float(np.diff(t).max()))
        if len(c['xyz']) != n:
            bad(f"craft {c['id']}: xyz length {len(c['xyz'])} != {n}")
        if 'm' in c and len(c['m']) != n:
            bad(f"craft {c['id']}: m length {len(c['m'])} != {n}")
        if not (c['launch_d'] - 1e-9 <= t[0] and t[-1] <= c['end_d'] + 1e-9):
            bad(f"craft {c['id']}: samples [{t[0]}, {t[-1]}] outside [launch {c['launch_d']}, end {c['end_d']}]")
        if t[0] < 0 or t[-1] > T_END_D + 1e-6:
            bad(f"craft {c['id']}: samples outside the mission [{t[0]}, {t[-1]}]")
        xyz = np.asarray(c['xyz'], float)
        r = np.linalg.norm(xyz, axis=1)
        if r.min() < 0.05 or r.max() > 25:
            bad(f"craft {c['id']}: |r| range {r.min():.3f}..{r.max():.3f} AU implausible")
        for a, b in c['thrust']:
            if not (a < b):
                bad(f"craft {c['id']}: thrust arc [{a}, {b}] not ordered")
            if a < c['launch_d'] - 1e-6 or b > c['end_d'] + 1e-6:
                bad(f"craft {c['id']}: thrust arc [{a}, {b}] outside [launch, end]")
        if c['n_flybys'] != len(c['flybys']):
            bad(f"craft {c['id']}: n_flybys {c['n_flybys']} != {len(c['flybys'])}")
        tset = {float(x): k for k, x in enumerate(c['t'])}
        for f in c['flybys']:
            n_fb += 1; targets.append(int(f['ast']))
            k = tset.get(float(f['t']))
            if k is None:
                bad(f"craft {c['id']}: flyby epoch {f['t']} (ast {f['ast']}) is not a sample epoch"); continue
            body = bodies.get(int(f['ast']))
            if body is None:
                bad(f"craft {c['id']}: unknown asteroid {f['ast']}"); continue
            ra, _ = body.state(float(f['t']) * DAY)
            err = float(np.linalg.norm(xyz[k] - ra / AU))
            if err > max_err:
                max_err = err; max_at = (c['id'], int(f['ast']), float(f['t']))
            if err > TOL_AU:
                bad(f"craft {c['id']}: flyby ast {f['ast']} at {f['t']} d off by {err:.2e} AU")
    covered = len(set(targets)); missed = sorted(set(range(1, 301)) - set(targets))
    if covered != fleet['covered']:
        bad(f"covered {fleet['covered']} != {covered} distinct targets")
    if missed != fleet['missed']:
        bad('missed list inconsistent with the flybys')
    if fleet['n_craft'] != len(fleet['craft']):
        bad(f"n_craft {fleet['n_craft']} != {len(fleet['craft'])}")
    sumJi = float(sum(cost_sc(c['m0']) for c in fleet['craft']))
    if abs(sumJi + 300 - covered - fleet['J_raw']) > 1e-4:
        bad(f"J_raw {fleet['J_raw']} != sum J_i {sumJi:.6f} + {300 - covered}")
    for k in ('n_craft', 'covered', 'J_raw', 'kind', 'title'):
        if entry.get(k) != fleet.get(k):
            bad(f'index.json {k} {entry.get(k)!r} != file {fleet.get(k)!r}')
    if n_fb == 0:
        bad('no flybys')
    return dict(key=entry['key'], kind=fleet['kind'], n_craft=len(fleet['craft']), n_flybys=n_fb, covered=covered,
                samples=n_samples, size=size, max_step_d=worst_step, max_err_au=max_err, max_at=max_at, errors=errs)

def main(argv):
    bodies, src = asteroid_bodies()
    idx = json.loads((HERE / 'index.json').read_text())
    entries = [e for e in idx['fleets'] if not argv or e['key'] in argv]
    print(f'asteroid elements from {src}; tolerance {TOL_AU:.0e} AU; size limit {SIZE_LIMIT / 1e6:.1f} MB')
    print(f'{"key":>5} {"kind":>10} {"craft":>5} {"flybys":>6} {"cov":>4} {"samples":>7} {"size MB":>8} {"maxdt":>5} {"max flyby err AU":>17}  worst (craft, ast, t)  status')
    ok = True; overall = 0.0; overall_at = None; total = 0
    for e in entries:
        r = check_fleet(e, bodies)
        total += r['size']
        if r['max_err_au'] > overall:
            overall = r['max_err_au']; overall_at = (r['key'],) + tuple(r['max_at'] or ())
        status = 'PASS' if not r['errors'] else 'FAIL'
        ok &= not r['errors']
        print(f"{r['key']:>5} {r['kind']:>10} {r['n_craft']:>5} {r['n_flybys']:>6} {r['covered']:>4} {r['samples']:>7} "
              f"{r['size'] / 1e6:>8.3f} {r['max_step_d']:>5.1f} {r['max_err_au']:>17.3e}  {r['max_at']}  {status}")
        for m in r['errors']:
            print('      -', m)
    print(f'\n{len(entries)} fleets, {total / 1e6:.3f} MB total; max flyby position error {overall:.3e} AU '
          f'({overall * AU:.1f} km) at {overall_at}; RESULT {"PASS" if ok else "FAIL"} max_flyby_err_au={overall:.6e}')
    return 0 if ok else 1

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
