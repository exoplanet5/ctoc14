"""Fix-and-optimise over the route pool (stage 7, section 9).

Set cover over the pool is optimal but pool-limited, and set-cover LP duals are too degenerate to price with (adding a
column collapses its targets' duals to 0 and the bound does not move).  So drive the pool instead: keep the routes the
cover already likes best (lowest J per flyby), FREEZE their targets, and grow a fresh batch of carrier routes on the
REMAINDER only -- where they face no competition and can be deep.  Re-solve the cover over everything; the MILP is free
to keep or drop the frozen routes.  Repeat with a different keep-count so different subsets get frozen.

Usage: cover_lns.py out_dir [--rounds 6] [--keep 4,5,3,6] [--nroute 24]
"""
import sys, json, time, pathlib, argparse, subprocess
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from ctoc14.search import UNREACHABLE
PY = sys.executable
ALL = [t for t in range(1, 301) if t not in UNREACHABLE]
PATS = ['results/s7/farm*/route_*.npz', 'results/s7/cgF*/route_*.npz', 'results/s7/LNS/lnsF*/route_*.npz',
        'results/**/route_*.npz']

ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--rounds', type=int, default=6)
ap.add_argument('--keep', default='4,5,3,6,4,5'); ap.add_argument('--nroute', type=int, default=24)
ap.add_argument('--max-tank', type=float, default=1150.0); ap.add_argument('--tmax', type=float, default=700.0)
ap.add_argument('--lib', default='results/s7/lib0.pkl'); ap.add_argument('--top', type=int, default=250)
a = ap.parse_args()
out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt')
keeps = [int(x) for x in a.keep.split(',')]
best = 1e9
for rnd in range(1, a.rounds + 1):
    cov = out / f'cov{rnd}'
    r = subprocess.run([PY, str(ROOT / 'tools/route_cover.py')] + PATS + ['--out', str(cov)], capture_output=True, text=True)
    fl = RI.IFleet(cov / 'fleet'); J = fl.J()
    say(f'round {rnd}: cover = {len(fl.routes)} craft, {len(fl.coverage())} covered, J {J:.3f}'
        + ('  <-- NEW BEST' if J < best - 1e-9 else ''))
    best = min(best, J)
    k = keeps[(rnd - 1) % len(keeps)]
    order = sorted(fl.routes, key=lambda n: RI.cost(fl.routes[n]['tank']) / len(fl.routes[n]['st']['asts']))
    frozen = order[:k]
    fixed = set(int(x) for n in frozen for x in fl.routes[n]['st']['asts'])
    rest = [t for t in ALL if t not in fixed]
    say(f'  freeze {k} routes (' + ' '.join(f'{len(fl.routes[n]["st"]["asts"])}@{fl.routes[n]["tank"]:.0f}' for n in frozen)
        + f'), grow {a.nroute} routes on the remaining {len(rest)} targets')
    json.dump(rest, open(out / f'allow{rnd}.json', 'w'))
    subprocess.run([PY, str(ROOT / 'tools/carrier_farm.py'), a.lib, str(out / f'lnsF{rnd}'), '--nroute', str(a.nroute),
                    '--frac', '1.0', '--ntrial', '12', '--nfail', '3', '--tmax', str(a.tmax), '--top', str(a.top),
                    '--seed', str(200 + rnd), '--max-tank', str(a.max_tank),
                    '--allow-file', str(out / f'allow{rnd}.json')], capture_output=True, text=True)
    n = len(list((out / f'lnsF{rnd}').glob('route_*.npz'))) if (out / f'lnsF{rnd}').exists() else 0
    say(f'  round {rnd}: {n} new routes')
r = subprocess.run([PY, str(ROOT / 'tools/route_cover.py')] + PATS + ['--out', str(out / 'final')], capture_output=True, text=True)
fl = RI.IFleet(out / 'final' / 'fleet')
say(f'FINAL: {len(fl.routes)} craft, {len(fl.coverage())} covered, J {fl.J():.3f}')
