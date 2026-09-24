"""Copy a planned fleet (tour_sc*.json) to a conversion directory with conversion tanks m0 = 600 + k * fuel_est + c
(fuel_est at the search mass; k 1.3, c 20 = the rule under which the camp1 fleet converted 10/12 craft losslessly).
tools/convert_fleet.py then converts (pass 1) and right-sizes (pass 2: m0 = 600 + 1.03 actual + 10).
Usage: prep_convert.py fleetdir convdir [--k 1.3] [--c 20]"""
import sys, json, pathlib, argparse, shutil
ap = argparse.ArgumentParser(); ap.add_argument('fleetdir'); ap.add_argument('convdir')
ap.add_argument('--k', type=float, default=1.3); ap.add_argument('--c', type=float, default=20.0)
a = ap.parse_args()
src = pathlib.Path(a.fleetdir); dst = pathlib.Path(a.convdir); dst.mkdir(parents=True, exist_ok=True)
for f in sorted(src.glob('tour_sc*.json'), key=lambda f: int(f.stem[7:])):
    t = json.load(open(f)); t['m0_plan'] = t.get('m0'); t['m0'] = min(2000.0, 600.0 + a.k * float(t['fuel_est']) + a.c)
    json.dump(t, open(dst / f.name, 'w'), indent=1)
    print(f'{f.name}: {len(t["legs"])} flybys, fuel_est {t["fuel_est"]:.1f} kg -> conversion tank {t["m0"]:.1f} kg')
if (src / 'summary.json').exists(): shutil.copy(src / 'summary.json', dst / 'plan_summary.json')
