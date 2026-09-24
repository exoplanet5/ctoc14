"""Craft heliocentric longitude relative to Earth over time (twin trajectories) -> phase.json; bunching statistics."""
import sys, json, pathlib
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from ctoc14.constants import DAY, AU
E = RI.eph()
F = {}
for f in sorted((ROOT / 'results/s16/best/fleet').glob('route_*.npz')):
    z = np.load(f); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL']); F[f.stem[6:]] = st
grid = np.arange(0, 5470, 30.0)
ph = {}; rr = {}
for k, st in F.items():
    ip = RI.ipr(st); tt, r = RI.trajectory(ip, step=2 * DAY)
    td = tt / DAY
    rE, _ = E.earth_state(tt)
    lc = np.degrees(np.arctan2(r[:, 1], r[:, 0])); le = np.degrees(np.arctan2(rE[:, 1], rE[:, 0]))
    d = (lc - le + 180) % 360 - 180
    p = np.full(len(grid), np.nan); rad = np.full(len(grid), np.nan)
    m = (grid >= td[0]) & (grid <= td[-1])
    p[m] = np.interp(grid[m], td, np.unwrap(np.radians(d)) * 180 / np.pi)
    rad[m] = np.interp(grid[m], td, np.linalg.norm(r, axis=1) / AU)
    ph[k] = ((p + 180) % 360 - 180).tolist(); rr[k] = rad.tolist()
json.dump(dict(grid=grid.tolist(), phase=ph, r=rr), open(ROOT / 'results/s17/redteam/phase.json', 'w'))
P = np.array([ph[k] for k in sorted(F)])
print('Earth-relative longitude (deg) of each craft at yr 0.5..14:')
for yr in (0.5, 1, 2, 3, 4, 6, 8, 10, 12, 14):
    i = int(np.argmin(abs(grid - yr * 365.25)))
    print(f' yr {yr:4.1f}: ' + ' '.join(f'{k}:{P[j, i]:6.0f}' for j, k in enumerate(sorted(F))))
# phase drift range per route
for j, k in enumerate(sorted(F)):
    x = np.array(ph[k]); x = x[np.isfinite(x)]
    print(k, 'phase unwrapped span (deg):', round(float(np.ptp(np.unwrap(np.radians(x)) * 180 / np.pi))))
