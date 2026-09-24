"""Sequential carrier fleet (stage 7): for each craft, DP-rank carriers on the still-uncovered targets, settle the best,
deepen it by insertion (tools/carrier_grow.py), add its targets to the covered set, repeat.

Usage: carrier_fleet.py out_dir [--ncraft 8] [--lam 1.0]
"""
import sys, json, shutil, pathlib, argparse, subprocess
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import run_ialns as RI
PY = sys.executable

ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--ncraft', type=int, default=8)
ap.add_argument('--lam', type=float, default=1.0); ap.add_argument('--n', type=int, default=60000)
ap.add_argument('--eps', type=float, default=0.012); ap.add_argument('--top', type=int, default=4)
ap.add_argument('--max-tank', type=float, default=1100.0); ap.add_argument('--max-step', type=float, default=45.0)
ap.add_argument('--target', type=int, default=60); ap.add_argument('--dmax', type=float, default=0.15)
ap.add_argument('--ntrial', type=int, default=32); ap.add_argument('--seed', type=int, default=0)
a = ap.parse_args()
out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt')
fleet = RI.IFleet(out / 'fleet') if (out / 'fleet' / 'fleet.json').exists() else RI.IFleet()
for k in range(len(fleet.routes), a.ncraft):
    cov = sorted(fleet.coverage()); ex = ','.join(str(x) for x in cov)
    d = out / f'dp_{k:02d}'
    subprocess.run([PY, str(ROOT / 'tools/carrier_dp.py'), '--n', str(a.n), '--eps', str(a.eps), '--pool', str(a.n),
                    '--top', str(a.top), '--lam', str(a.lam), '--kphase', '1.5', '--kgap', '10', '--seed', str(a.seed + k),
                    '--settle', '--fleet-out', str(d)] + (['--exclude', ex] if ex else []), capture_output=True, text=True)
    if not (d / 'fleet.json').exists(): say(f'craft {k}: no carrier route settled'); break
    info = json.load(open(d / 'fleet.json'))['routes']
    best = max(info, key=lambda n: info[n]['flybys'] - (info[n]['tank'] - 600.0) / 40.0)
    say(f'craft {k}: carrier routes ' + ' '.join(f'{info[n]["flybys"]}@{info[n]["tank"]:.0f}' for n in sorted(info)) + f' -> {best}')
    g = out / f'grow_{k:02d}'
    subprocess.run([PY, str(ROOT / 'tools/carrier_grow.py'), str(d), best, str(g), '--target', str(a.target),
                    '--max-tank', str(a.max_tank), '--max-step', str(a.max_step), '--dmax', str(a.dmax),
                    '--ntrial', str(a.ntrial)] + (['--exclude', ex] if ex else []), capture_output=True, text=True)
    src = RI.IFleet(g) if (g / 'fleet.json').exists() else RI.IFleet(d)
    fleet.routes[f'{k:02d}'] = src.routes[best]
    fleet.save(out / 'fleet', note='carrier fleet')
    r = fleet.routes[f'{k:02d}']
    say(f'craft {k}: {len(r["st"]["asts"])} flybys @ {r["tank"]:.0f} kg -> fleet covered {len(fleet.coverage())}, '
        f'sum J_i {sum(RI.cost(x["tank"]) for x in fleet.routes.values()):.3f}, J {fleet.J():.3f}')
