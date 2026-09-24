"""Run the beam search for one spacecraft.
Usage: run_search.py out.json [--beam 300] [--eta 0.6] [--wfuel 1.0] [--exclude ids.json] [--m0 2000] [--launch-max-days 1096] [--launch-step 20]"""
import sys, json, time, pathlib, argparse, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ctoc14.kepler import Ephemeris
from ctoc14.search import beam_search, Params, save_tour
from ctoc14.constants import DAY
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--beam', type=int, default=300); ap.add_argument('--eta', type=float, default=0.6)
ap.add_argument('--wfuel', type=float, default=1.0); ap.add_argument('--exclude', default=None); ap.add_argument('--m0', type=float, default=2000.0)
ap.add_argument('--launch-max-days', type=float, default=1096); ap.add_argument('--launch-min-days', type=float, default=0); ap.add_argument('--launch-step', type=float, default=20); ap.add_argument('--nproc', type=int, default=8)
ap.add_argument('--tof-step', type=float, default=5); ap.add_argument('--tof-max', type=float, default=400); ap.add_argument('--tof-min', type=float, default=15)
ap.add_argument('--mmargin', type=float, default=100.0); ap.add_argument('--dvmax', type=float, default=2.5)
ap.add_argument('--lin-tofmax', type=float, default=0, help='enable the linearised leg model with TOF grid up to this many days (0 = off)')
ap.add_argument('--pool-out', default=None, help='write every beam state (all depths) as routes to this jsonl'); ap.add_argument('--wait', type=float, default=0, help='unused placeholder')
ap.add_argument('--tof-refine', action='store_true'); ap.add_argument('--vinf-cap', type=float, default=4.0, help='max |v_inf| of launch legs (km/s)'); ap.add_argument('--lin-step', type=float, default=10); ap.add_argument('--lin-ub', type=float, default=0.8); ap.add_argument('--lin-drmax', type=float, default=0.25, help='AU')
a = ap.parse_args()
from ctoc14.constants import AU
lin = np.arange(20, a.lin_tofmax + 1e-9, a.lin_step) * DAY if a.lin_tofmax > 0 else None
P = Params(beam=a.beam, eta=a.eta, w_fuel=a.wfuel, tofs=np.arange(a.tof_min, a.tof_max + 1e-9, a.tof_step) * DAY, m_margin=a.mmargin, dv_max=a.dvmax,
           lin_tofs=lin, lin_ub=a.lin_ub, lin_drmax=a.lin_drmax * AU, vinf_cap=a.vinf_cap, tof_refine=a.tof_refine)
excluded = set(json.load(open(a.exclude))) if a.exclude else set()
if a.pool_out: P.collect = []
eph = Ephemeris()
tic = time.time()
best, beam = beam_search(eph, excluded=excluded, m0=a.m0, t_launch_grid=np.arange(a.launch_min_days, a.launch_max_days, a.launch_step) * DAY, P=P, n_proc=a.nproc)
print(f'BEST: {best.n()} flybys, fuel {best.fuel:.0f} kg, launch day {best.t_launch/DAY:.0f}, end {best.t/DAY/365.25:.2f} yr, runtime {time.time()-tic:.0f} s')
print('sequence:', [x[0] for x in best.seq])
save_tour(best, a.out)
if a.pool_out:
    from ctoc14.search import tour_json
    seen = set(); n = 0
    with open(a.pool_out, 'w') as f:
        for s in P.collect + list(beam):
            key = (s.visited, round(s.t_launch)); 
            if key in seen: continue
            seen.add(key); f.write(json.dumps(tour_json(s)) + '\n'); n += 1
    print(f'pool: {n} routes -> {a.pool_out}')
