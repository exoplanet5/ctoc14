"""Greedy set-cover campaign with a route pool: every round runs a single-craft beam search on the still-uncovered targets
(previous rounds' targets excluded), keeps the best route, dumps EVERY beam state of every depth into the pool, and finally
runs the set-cover MILP over the whole pool (tools/select_routes_milp.py) for the fleet.
Usage: run_pool_campaign.py outdir [--rounds 12] [--m0 1000] [--wfuel 1.5] [--beam 400] [--nproc 10] [--min-new 6] [--exclude ids.json]"""
import sys, json, time, pathlib, argparse, subprocess, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ctoc14.kepler import Ephemeris
from ctoc14.search import beam_search, Params, save_tour, tour_json, UNREACHABLE
from ctoc14.constants import DAY, AU, cost_sc
ap = argparse.ArgumentParser(); ap.add_argument('outdir'); ap.add_argument('--rounds', type=int, default=12); ap.add_argument('--m0', type=float, default=1000)
ap.add_argument('--wfuel', type=float, default=1.5); ap.add_argument('--beam', type=int, default=400); ap.add_argument('--nproc', type=int, default=10)
ap.add_argument('--dvmax', type=float, default=1.2); ap.add_argument('--lin-tofmax', type=float, default=600); ap.add_argument('--lin-drmax', type=float, default=0.15)
ap.add_argument('--vinf-cap', type=float, default=2.0); ap.add_argument('--mmargin', type=float, default=40); ap.add_argument('--launch-max', type=float, default=600)
ap.add_argument('--launch-step', type=float, default=20); ap.add_argument('--min-new', type=int, default=6); ap.add_argument('--exclude', default=None)
ap.add_argument('--tof-refine', action='store_true'); ap.add_argument('--wrare', type=float, default=0.0); ap.add_argument('--rare-ids', default=None, help='json list of target IDs given rarity 1 (others 0)'); ap.add_argument('--rare-weights', default=None, help='json {id: weight} (e.g. from price_targets.py)')
a = ap.parse_args(); out = pathlib.Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
eph = Ephemeris(); covered = set(json.load(open(a.exclude))) if a.exclude else set()
lin = np.arange(20, a.lin_tofmax + 1e-9, 10) * DAY
rarity = None
if a.wrare > 0 and a.rare_weights:
    rarity = np.zeros(300)
    for k, v in json.load(open(a.rare_weights)).items(): rarity[int(k) - 1] = float(v)
elif a.wrare > 0 and a.rare_ids:
    rarity = np.zeros(300); rarity[[i - 1 for i in json.load(open(a.rare_ids))]] = 1.0
elif a.wrare > 0:
    import pandas as pd
    d = pd.read_csv('analysis/population/population_stats.csv').sort_values('id')
    rarity = np.clip(3.0 / np.maximum(d['n_node_crossings_in_band'].to_numpy(float), 1.0), 0, 1) ** 2
pool = out / 'pool.jsonl'; log = open(out / 'campaign.log', 'a')
def say(s): print(s, flush=True); log.write(s + '\n'); log.flush()
for r in range(1, a.rounds + 1):
    P = Params(beam=a.beam, w_fuel=a.wfuel, m_margin=a.mmargin, dv_max=a.dvmax, lin_tofs=lin, lin_drmax=a.lin_drmax * AU, vinf_cap=a.vinf_cap,
               w_rare=a.wrare, rarity=rarity, tof_refine=a.tof_refine); P.collect = []
    tic = time.time()
    best, beam = beam_search(eph, excluded=covered, m0=a.m0, t_launch_grid=np.arange(0, a.launch_max, a.launch_step) * DAY, P=P, n_proc=a.nproc, verbose=False)
    if best is None: say(f'round {r}: nothing'); break
    seen = set(); n = 0
    with open(pool, 'a') as f:
        for s in P.collect + list(beam):
            key = (s.visited, round(s.t_launch))
            if key in seen: continue
            seen.add(key); f.write(json.dumps(tour_json(s)) + '\n'); n += 1
    new = {x[0] for x in best.seq} - covered
    save_tour(best, out / f'round{r}_best.json')
    say(f'round {r}: best {best.n()} flybys, fuel {best.fuel:.0f} kg, launch {best.t_launch/DAY:.0f} d, end {best.t/DAY/365.25:.1f} yr, new {len(new)}, pool +{n} ({time.time()-tic:.0f} s); covered so far {len(covered | new)}')
    covered |= new
    json.dump(sorted(covered), open(out / 'covered.json', 'w'))
    if len(new) < a.min_new: break
say(f'campaign done: covered {len(covered)} / 298; uncovered {sorted(set(range(1,301)) - covered - UNREACHABLE)}')
for n in (10, 11, 12):
    say(f'--- MILP max-craft {n}')
    res = subprocess.run([sys.executable, 'tools/select_routes_milp.py', str(out / f'milp_n{n}'), str(pool), '--fuel-margin', '0.3', '--reserve', '20', '--max-craft', str(n), '--time-limit', '600'], capture_output=True, text=True)
    say('\n'.join(l for l in res.stdout.splitlines() if l.startswith('SELECTED') or l.startswith('  SC') or 'routes' in l))
