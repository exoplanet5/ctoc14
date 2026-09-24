"""Leave-one-out improvement of a joint (swarm) plan: re-plan each craft alone (single-craft beam search, patient settings)
on pool = its own targets + current leftovers, others' targets excluded; accept if fleet J (right-sized tanks) improves.
Usage: improve_joint.py plandir [--rounds 2] [--beam 600] [--wfuel 3] [--tofmax 500] [--nproc 8] [--m0 same|1300]"""
import sys, json, time, pathlib, argparse, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ctoc14.kepler import Ephemeris
from ctoc14.search import beam_search, Params, tour_json
from ctoc14.constants import DAY, cost_sc
ap = argparse.ArgumentParser(); ap.add_argument('plandir'); ap.add_argument('--rounds', type=int, default=2); ap.add_argument('--beam', type=int, default=600)
ap.add_argument('--wfuel', type=float, default=3.0); ap.add_argument('--tofmax', type=float, default=500); ap.add_argument('--nproc', type=int, default=8)
ap.add_argument('--m0', default='same'); ap.add_argument('--mmargin', type=float, default=20); ap.add_argument('--dvmax', type=float, default=2.5)
a = ap.parse_args(); d = pathlib.Path(a.plandir); eph = Ephemeris()
files = sorted(d.glob('tour_sc*.json'), key=lambda f: int(f.stem[7:])); tours = [json.load(open(f)) for f in files]
def tset(t): return {l['ast'] for l in t['legs']}
def rs_cost(t): return cost_sc(min(2000.0, 600 + t['fuel_est'] + 30))
def fleetJ(ts):
    cov = set().union(*[tset(t) for t in ts]); return sum(rs_cost(t) for t in ts) + 300 - len(cov), cov
J0, cov = fleetJ(tours); print(f'start: {len(tours)} craft, covered {len(cov)}, J(rs) = {J0:.2f}', flush=True)
P = Params(beam=a.beam, w_fuel=a.wfuel, m_margin=a.mmargin, dv_max=a.dvmax, tofs=np.arange(15, a.tofmax + 1e-9, 5) * DAY)
for rnd in range(a.rounds):
    improved = False
    for k in range(len(tours)):
        others = set().union(*[tset(t) for j, t in enumerate(tours) if j != k])
        m0 = tours[k]['m0'] if a.m0 == 'same' else float(a.m0)
        tL = tours[k]['t_launch']; grid = np.arange(max(0, tL - 200 * DAY), tL + 400 * DAY + 1, 20 * DAY)
        tic = time.time(); best, _ = beam_search(eph, excluded=others, m0=m0, t_launch_grid=grid, P=P, n_proc=a.nproc, verbose=False)
        if best is None: continue
        trial = list(tours); trial[k] = tour_json(best)
        J1, cov1 = fleetJ(trial)
        tag = 'accept' if J1 < J0 - 1e-6 else 'reject'
        print(f'round {rnd+1} SC{k+1}: {tours[k]["n_flybys"]} -> {best.n()} flybys, fuel {tours[k]["fuel_est"]:.0f} -> {best.fuel:.0f} kg, covered {len(cov)} -> {len(cov1)}, J {J0:.2f} -> {J1:.2f} [{tag}] ({time.time()-tic:.0f}s)', flush=True)
        if J1 < J0 - 1e-6:
            tours = trial; J0, cov = J1, cov1; improved = True
            json.dump(tours[k], open(files[k], 'w'), indent=1)
    if not improved: break
print(f'final: covered {len(cov)}, J(rs) = {J0:.2f}; leftover {sorted(set(range(1,301)) - cov - {131,144})}')
