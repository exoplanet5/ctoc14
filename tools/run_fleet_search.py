"""Fleet planning: sequential beam searches, then leave-one-out improvement rounds and weakest-spacecraft elimination.
Usage: run_fleet_search.py outdir [--beam 300] [--wfuel 0.7] [--eta 0.6] [--max-sc 12] [--mmargin 100] [--dvmax 2.5]
       [--wrare 0] [--rarity-csv analysis/ballistic/out/junction_per_asteroid.csv] [--improve-rounds 2] [--eliminate]
       [--launch-max-days 400] [--launch-max-days-late 1000] [--nproc 6] [--resume]"""
import sys, json, time, pathlib, argparse, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ctoc14.kepler import Ephemeris
from ctoc14.search import beam_search, Params, save_tour, tour_json, UNREACHABLE
from ctoc14.constants import DAY, cost_sc

ap = argparse.ArgumentParser(); ap.add_argument('outdir'); ap.add_argument('--beam', type=int, default=300); ap.add_argument('--wfuel', type=float, default=0.7)
ap.add_argument('--eta', type=float, default=0.6); ap.add_argument('--max-sc', type=int, default=12); ap.add_argument('--min-gain', type=float, default=0.0)
ap.add_argument('--launch-max-days', type=float, default=400); ap.add_argument('--launch-max-days-late', type=float, default=1000); ap.add_argument('--launch-step', type=float, default=20)
ap.add_argument('--nproc', type=int, default=6); ap.add_argument('--resume', action='store_true'); ap.add_argument('--mmargin', type=float, default=100.0)
ap.add_argument('--dvmax', type=float, default=2.5); ap.add_argument('--wrare', type=float, default=0.0)
ap.add_argument('--rarity-csv', default='analysis/ballistic/out/junction_per_asteroid.csv'); ap.add_argument('--improve-rounds', type=int, default=2)
ap.add_argument('--eliminate', action='store_true'); ap.add_argument('--late-from', type=int, default=5)
ap.add_argument('--tail-from', type=int, default=99, help='from this spacecraft index on, use tail mode: launch anywhere, TOF up to --tail-tof-max')
ap.add_argument('--tail-tof-max', type=float, default=600); ap.add_argument('--tail-launch-step', type=float, default=30)
a = ap.parse_args()
outdir = pathlib.Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)
eph = Ephemeris()
rarity = None
if a.wrare > 0:
    import pandas as pd
    d = pd.read_csv(a.rarity_csv).sort_values('B'); f = d['frac_feas_plausible_tof180_a0.5'].to_numpy()
    rarity = np.nan_to_num(np.clip(1.0 - f / np.nanmax(f), 0, 1), nan=0.5)
    print(f'rarity weighting on: mean {rarity.mean():.2f}, hardest {np.argsort(-rarity)[:10] + 1}')
P = Params(beam=a.beam, eta=a.eta, w_fuel=a.wfuel, m_margin=a.mmargin, dv_max=a.dvmax, w_rare=a.wrare, rarity=rarity)

def plan(excluded, sc_index):
    if sc_index >= a.tail_from:
        PT = Params(beam=a.beam, eta=a.eta, w_fuel=a.wfuel, m_margin=a.mmargin, dv_max=a.dvmax, w_rare=a.wrare, rarity=rarity,
                    tofs=np.arange(15, a.tail_tof_max + 1e-9, 5) * DAY)
        grid = np.arange(0, 5400, a.tail_launch_step) * DAY
        best, _ = beam_search(eph, excluded=excluded, m0=2000.0, t_launch_grid=grid, P=PT, n_proc=a.nproc, verbose=False)
        return best
    lmax = a.launch_max_days if sc_index < a.late_from else a.launch_max_days_late
    best, _ = beam_search(eph, excluded=excluded, m0=2000.0, t_launch_grid=np.arange(0, lmax, a.launch_step) * DAY, P=P, n_proc=a.nproc, verbose=False)
    return best

def targets(t): return {l['ast'] for l in t['legs']}
def report(tours, tag):
    cov = set().union(*[targets(t) for t in tours]) if tours else set()
    J = sum(cost_sc(t['m0']) for t in tours) + (300 - len(cov))
    print(f'[{tag}] {len(tours)} spacecraft, flybys {[t["n_flybys"] for t in tours]}, covered {len(cov)}, planned J = {J:.2f}', flush=True)
    return cov, J

tours = []
if a.resume:
    for f in sorted(outdir.glob('tour_sc*.json'), key=lambda f: int(f.stem.split('sc')[1])):
        tours.append(json.load(open(f)))
    print(f'resumed {len(tours)} tours')
# ---- sequential pass
covered = set().union(*[targets(t) for t in tours]) if tours else set()
for sc in range(len(tours) + 1, a.max_sc + 1):
    tic = time.time(); best = plan(covered, sc)
    if best is None or best.n() - cost_sc(best.m0) < a.min_gain:
        print('stopping sequential pass'); break
    tours.append(tour_json(best)); covered |= targets(tours[-1])
    print(f'SC{sc}: {best.n()} flybys, fuel {best.fuel:.0f} kg, launch {best.t_launch/DAY:.0f} d, {time.time()-tic:.0f} s', flush=True)
cov, J = report(tours, 'sequential')
# ---- leave-one-out improvement rounds
for rnd in range(a.improve_rounds):
    improved = False
    order = sorted(range(len(tours)), key=lambda k: tours[k]['n_flybys'])
    for k in order:
        others = set().union(*[targets(t) for j, t in enumerate(tours) if j != k])
        tic = time.time(); best = plan(others, k + 1)
        if best is not None and (best.n() > tours[k]['n_flybys'] or (best.n() == tours[k]['n_flybys'] and best.fuel < tours[k]['fuel_est'])):
            print(f'  round {rnd+1}: SC{k+1} {tours[k]["n_flybys"]} -> {best.n()} ({time.time()-tic:.0f} s)', flush=True)
            tours[k] = tour_json(best); improved = True
    cov, J = report(tours, f'improve round {rnd+1}')
    if not improved: break
# ---- elimination of the weakest spacecraft
if a.eliminate:
    while len(tours) > 1:
        k = min(range(len(tours)), key=lambda i: tours[i]['n_flybys'])
        trial = [t for i, t in enumerate(tours) if i != k]
        # let the others absorb the freed targets (one leave-one-out round, weakest first)
        for j in sorted(range(len(trial)), key=lambda i: trial[i]['n_flybys']):
            others = set().union(*[targets(t) for i, t in enumerate(trial) if i != j])
            best = plan(others, j + 1)
            if best is not None and best.n() >= trial[j]['n_flybys']:
                trial[j] = tour_json(best)
        cov2, J2 = report(trial, f'elimination trial (removed SC{k+1} with {tours[k]["n_flybys"]})')
        if J2 < J - 1e-9:
            tours = trial; cov, J = cov2, J2; print('  -> accepted', flush=True)
        else:
            print('  -> rejected', flush=True); break
for f in outdir.glob('tour_sc*.json'): f.unlink()
for i, t in enumerate(tours, 1):
    json.dump(t, open(outdir / f'tour_sc{i}.json', 'w'), indent=1)
json.dump(sorted(cov), open(outdir / 'covered.json', 'w'))
json.dump(dict(J=J, n_sc=len(tours), covered=len(cov), flybys=[t['n_flybys'] for t in tours], args=vars(a)), open(outdir / 'plan_summary.json', 'w'), indent=1)
print(f'FLEET: {len(tours)} spacecraft, covered {len(cov)}, planned J = {J:.2f}; leftovers {sorted(set(range(1,301)) - cov - UNREACHABLE)}')
