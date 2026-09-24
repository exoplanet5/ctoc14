"""Leftover pass over a carrier fleet (stage 7): grow every route on the targets no route covers yet, lightest tank
first, with the step limit a LEFTOVER deserves (the alternative is a miss = 1 J ~ 1000 kg, or a ninth craft ~ 70 kg per
target), not the 45 kg of the N trade-off.

Usage: carrier_pass.py fleet_dir out_dir [--max-step 90] [--max-tank 1250] [--dmax 0.25] [--rounds 2]
"""
import sys, pathlib, argparse, subprocess
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import run_ialns as RI
PY = sys.executable
ap = argparse.ArgumentParser(); ap.add_argument('src'); ap.add_argument('out')
ap.add_argument('--max-step', type=float, default=90.0); ap.add_argument('--max-tank', type=float, default=1250.0)
ap.add_argument('--dmax', type=float, default=0.25); ap.add_argument('--ntrial', type=int, default=40)
ap.add_argument('--rounds', type=int, default=2); ap.add_argument('--sep', type=float, default=200.0)
a = ap.parse_args()
out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt')
fleet = RI.IFleet(a.src); fleet.save(out / 'fleet', note='pass start')
say(f'start: covered {len(fleet.coverage())}, sum J_i {sum(RI.cost(r["tank"]) for r in fleet.routes.values()):.3f}')
for rd in range(a.rounds):
    gained = 0
    for n in sorted(fleet.routes, key=lambda n: fleet.routes[n]['tank']):
        ex = sorted(set(fleet.coverage()) - set(int(x) for x in fleet.routes[n]['st']['asts']))
        g = out / f'r{rd}_{n}'
        subprocess.run([PY, str(ROOT / 'tools/carrier_grow.py'), str(out / 'fleet'), n, str(g), '--target', '70',
                        '--max-tank', str(a.max_tank), '--max-step', str(a.max_step), '--dmax', str(a.dmax),
                        '--ntrial', str(a.ntrial), '--sep', str(a.sep), '--exclude', ','.join(map(str, ex))],
                       capture_output=True, text=True)
        if (g / 'fleet.json').exists():
            new = RI.IFleet(g).routes[n]; k = len(new['st']['asts']) - len(fleet.routes[n]['st']['asts'])
            if k > 0:
                say(f'  round {rd} route {n}: +{k} -> {len(new["st"]["asts"])} @ {new["tank"]:.0f} kg (was {fleet.routes[n]["tank"]:.0f})')
                fleet.routes[n] = new; gained += k; fleet.save(out / 'fleet', note='leftover pass')
    sj = sum(RI.cost(r['tank']) for r in fleet.routes.values())
    say(f'round {rd}: +{gained} -> covered {len(fleet.coverage())}, sum J_i {sj:.3f}, J {fleet.J():.3f}')
    if gained == 0: break
left = sorted(set(range(1, 301)) - set(fleet.coverage()))
say(f'left ({len(left)}): {left}')
