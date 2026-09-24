"""Insert uncovered targets into existing planned tours (Lambert-junction model, see ctoc14/insertion.py).
Usage: run_insertion.py tourdir outdir (--targets 15,137,... | --targets-json list.json)
          [--tours tour_sc1,tour_sc2,...] [--full-tank-only] [--margin 60] [--eta 0.6] [--tstep 1] [--delta-max 10]
          [--delta-step 2.5] [--tof-min 15] [--tof-max-last 400] [--no-launch-gap] [--u-max 1.0] [--nproc 4]
          [--grow-m0] [--try-unreachable] [--tag name]
Writes every tour (modified or not) as outdir/<name>.json (a complete tour set for convert_fleet.py), plus
outdir/insertion_<tag>.json (history + per-target table) and prints the tables."""
import sys, json, time, pathlib, argparse, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ctoc14.kepler import Ephemeris
from ctoc14.insertion import InsParams, load_tours, check_reconstruction, insert_targets
from ctoc14.constants import DAY, M_DRY, M0_MAX, cost_sc

ap = argparse.ArgumentParser()
ap.add_argument('tourdir'); ap.add_argument('outdir')
ap.add_argument('--targets', default=None); ap.add_argument('--targets-json', default=None)
ap.add_argument('--tours', default=None, help='comma-separated tour names (file stems); default: all tour_sc*.json')
ap.add_argument('--full-tank-only', action='store_true', help='only tours with m0 = 2000 kg')
ap.add_argument('--margin', type=float, default=60.0); ap.add_argument('--eta', type=float, default=0.6)
ap.add_argument('--tstep', type=float, default=1.0); ap.add_argument('--delta-max', type=float, default=10.0)
ap.add_argument('--delta-step', type=float, default=2.5); ap.add_argument('--tof-min', type=float, default=15.0)
ap.add_argument('--tof-max-last', type=float, default=400.0); ap.add_argument('--no-launch-gap', action='store_true')
ap.add_argument('--u-max', type=float, default=1.0); ap.add_argument('--nproc', type=int, default=4)
ap.add_argument('--grow-m0', action='store_true'); ap.add_argument('--try-unreachable', action='store_true')
ap.add_argument('--tag', default='run')
ap.add_argument('--from-report', default=None, help='convert_fleet report JSON: use the CONVERTED flyby sequences/times and m0 '
                'of results[k] (tour_sc<k+1>) as the base tours instead of the planned legs')
ap.add_argument('--allow-infeasible-base', action='store_true', help='with --from-report: keep bases whose Lambert-chain model is infeasible')
a = ap.parse_args()

targets = []
if a.targets: targets += [int(x) for x in a.targets.split(',') if x.strip()]
if a.targets_json: targets += [int(x) for x in json.load(open(a.targets_json))]
targets = sorted(set(targets))
if not targets:
    sys.exit('no targets given')
td = pathlib.Path(a.tourdir); od = pathlib.Path(a.outdir); od.mkdir(parents=True, exist_ok=True)
if a.tours:
    files = [td / f'{nm}.json' for nm in a.tours.split(',')]
else:
    files = sorted([f for f in td.glob('tour_sc*.json') if f.stem.split('sc')[-1].isdigit()], key=lambda f: int(f.stem.split('sc')[-1]))
eph = Ephemeris()
tours = load_tours(eph, files, eta=a.eta)
if a.from_report:
    from ctoc14.insertion import tour_from_converted, n_infeasible
    rep = json.load(open(a.from_report))['results']
    conv_tours = {}
    for f in files:
        nm = f.stem; k = int(nm.split('sc')[-1]) - 1
        if k < len(rep) and rep[k].get('ok') and rep[k].get('legs'):
            conv_tours[nm] = tour_from_converted(eph, json.load(open(f)), rep[k], name=nm, eta=a.eta)
    print('converted bases (planner model of the converted flyby chain vs actual conversion):')
    for nm, T in conv_tours.items():
        print(f'  {nm:12s} m0={T.m0:7.1f} n={T.n:2d} (planned {tours[nm].n:2d}) fuel_model={T.fuel:7.1f} kg fuel_actual={T.extra["fuel_actual"]:7.1f} kg '
              f'infeasible_legs={n_infeasible(T)} min_slack={np.min(T.slack):+.3f} budget(margin)={T.budget(a.margin):+7.1f} kg')
    if not a.allow_infeasible_base:
        bad = [nm for nm, T in conv_tours.items() if n_infeasible(T) > 0]
        if bad:
            print(f'  WARNING: the Lambert-junction model of these converted chains is infeasible/overpriced -> excluded as bases: {bad} '
                  f'(use --allow-infeasible-base to force; extra-propellant numbers are then model artefacts)')
        conv_tours = {nm: T for nm, T in conv_tours.items() if nm not in bad}
    tours = conv_tours
if a.full_tank_only:
    tours = {nm: T for nm, T in tours.items() if T.m0 >= M0_MAX - 1e-6 or (a.from_report and json.load(open(td / f'{nm}.json'))['m0'] >= M0_MAX - 1e-6)}
P = InsParams(eta=a.eta, margin=a.margin, t_step=a.tstep * DAY, tof_min=a.tof_min * DAY, tof_max_last=a.tof_max_last * DAY,
              delta_max=a.delta_max * DAY, delta_step=a.delta_step * DAY, launch_gap=not a.no_launch_gap, u_max=a.u_max,
              try_unreachable=a.try_unreachable)
print(f'targets ({len(targets)}): {targets}')
print(f'params: margin={a.margin} kg eta={a.eta} tstep={a.tstep} d delta=+-{a.delta_max} d step {a.delta_step} d tof_min={a.tof_min} d '
      f'launch_gap={P.launch_gap} u_max={a.u_max} grow_m0={a.grow_m0}')
# --- reconstruction check
print('\nreconstruction check (dv vs stored dv_est):' if not a.from_report else '\nbase tours:')
for f in files:
    nm = f.stem
    if nm not in tours: continue
    if a.from_report:
        T = tours[nm]
        print(f'  {nm:12s} m0={T.m0:7.1f} n={T.n:2d} fuel_model={T.fuel:7.1f} kg budget(margin)={T.budget(a.margin):+7.1f} kg '
              f'already_has={sorted(set(targets) & set(T.asts.tolist()))}')
        continue
    c = check_reconstruction(tours[nm], json.load(open(f)))
    print(f'  {nm:12s} m0={tours[nm].m0:7.1f} n={c["n"]:2d} fuel={c["fuel_recomputed"]:7.1f} kg (est {c["fuel_est"]:7.1f}) max|ddv|={c["max_ddv"]:.1e} km/s '
          f'min_slack={c["min_slack"]:+.3f} feasible={c["feasible"]} budget(margin)={tours[nm].budget(a.margin):+7.1f} kg '
          f'already_has={sorted(set(targets) & set(tours[nm].asts.tolist()))}')
# --- greedy insertion
tic = time.time()
tours, history, table, remaining = insert_targets(eph, tours, targets, P, nproc=a.nproc, verbose=True, grow_m0=a.grow_m0)
print(f'\ninsertion finished in {time.time() - tic:.0f} s: {len(history)} inserted, {len(remaining)} remaining {remaining}')
print('\ninserted:')
print('  target tour          pos  after before   t_j[d]  delta[d]  extra[kg]   dvA   dvB   dvC   u   budget_after   m0    dJ')
for h in history:
    print(f'  {h["target"]:5d}  {h["tour"]:12s} {h["pos"]:3d}  {h["after_ast"]:5d} {h["before_ast"]:5d}  {h["t_j_d"]:8.1f}  {h["delta_d"]:+6.1f}  '
          f'{h["extra_kg"]:8.2f}  {h["dvA"]:5.2f} {h["dvB"]:5.2f} {h["dvC"] if np.isfinite(h["dvC"]) else 0:5.2f} {h["u_mod"]:5.2f}  '
          f'{h["budget_after"]:9.1f}  {h["m0"]:7.1f}  {h["dJ"]:.3f}')
print('\nremaining targets: cheapest leg-feasible candidate per tour ignoring the budget (extra kg needed vs budget):')
rem_rows = []
for j in remaining:
    row = []
    for nm, TM in tours.items():
        r = table.get((nm, j))
        if r is None:
            continue
        b = r['best_nobudget']
        if b is None:
            row.append((nm, None, r['reason'], r['budget']))
        else:
            row.append((nm, b['extra_kg'], f'after {b["after_ast"]} at {b["t_j_d"]:.0f} d d={b["delta_d"]:+.1f} dv {b["dvA"]:.2f}/{b["dvB"]:.2f}/{(b["dvC"] if np.isfinite(b["dvC"]) else 0):.2f} u={b["u_mod"]:.2f}', r['budget']))
    row.sort(key=lambda x: (x[1] is None, x[1] if x[1] is not None else 0))
    best = row[0] if row else None
    rem_rows.append(dict(target=j, options=[dict(tour=x[0], extra_kg=x[1], note=x[2], budget=x[3]) for x in row]))
    print(f'  {j:4d}: ' + '; '.join(f'{x[0]} {x[1]:.1f} kg (budget {x[3]:+.1f}) {x[2]}' if x[1] is not None else f'{x[0]} -- ({x[2]})' for x in row[:4]))
# --- write tours
for nm, TM in tours.items():
    json.dump(TM.to_json(), open(od / f'{nm}.json', 'w'), indent=1)
J = sum(cost_sc(T.m0) for T in tours.values())
summary = dict(targets=targets, inserted=[h['target'] for h in history], remaining=remaining, history=history, remaining_table=rem_rows,
               params=dict(margin=a.margin, eta=a.eta, tstep=a.tstep, delta_max=a.delta_max, delta_step=a.delta_step, tof_min=a.tof_min,
                           tof_max_last=a.tof_max_last, launch_gap=P.launch_gap, u_max=a.u_max, grow_m0=a.grow_m0),
               tours={nm: dict(m0=T.m0, n=T.n, fuel=T.fuel, budget=T.budget(a.margin), t_end_d=T.t_end / DAY, inserted=T.extra.get('inserted', []))
                      for nm, T in tours.items()}, sum_J_i=J)
json.dump(summary, open(od / f'insertion_{a.tag}.json', 'w'), indent=1, default=float)
print(f'\nwrote {len(tours)} tours to {od}; sum J_i of these tours = {J:.3f}; report {od / f"insertion_{a.tag}.json"}')
