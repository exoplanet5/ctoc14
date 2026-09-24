"""Summarise a joint-search fleet directory: per-craft flybys/fuel/time, marginal fuel, leftover pool. Usage: analyze_joint.py dir"""
import sys, json, glob, pathlib, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ctoc14.constants import DAY, cost_sc
d = pathlib.Path(sys.argv[1]); tours = [json.load(open(f)) for f in sorted(d.glob('tour_sc*.json'), key=lambda f: int(f.stem[7:]))]
cov = set(); tot_f = 0; tot_n = 0
for k, t in enumerate(tours, 1):
    S = {l['ast'] for l in t['legs']}; cov |= S; tot_f += t['fuel_est']; tot_n += t['n_flybys']
    dv = [l['dv_est'] for l in t['legs']]; tof = [l['tof'] / DAY for l in t['legs']]
    print(f"SC{k:2d}: n={t['n_flybys']:2d} fuel={t['fuel_est']:5.0f}/{t['m0']-600:.0f} kg  launch {t['t_launch']/DAY:4.0f} d  end {t['t_end']/DAY/365.25:5.1f} yr  dv_med {np.median(dv):.2f} max {max(dv):.2f}  tof_med {np.median(tof):4.0f} d  kg/flyby {t['fuel_est']/t['n_flybys']:.1f}")
J = sum(cost_sc(t['m0']) for t in tours) + 300 - len(cov); Jrs = sum(cost_sc(min(2000, 600 + t['fuel_est'] + 30)) for t in tours) + 300 - len(cov)
print(f"fleet: {len(tours)} craft, {tot_n} flybys, {len(cov)} covered, fuel {tot_f:.0f} kg ({tot_f/tot_n:.1f} kg/flyby), planned J {J:.2f} (right-sized tanks {Jrs:.2f})")
left = sorted(set(range(1, 301)) - cov - {131, 144}); print(f"leftover {len(left)}: {left}")
json.dump(sorted(cov), open(d / 'covered.json', 'w'))
