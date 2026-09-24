"""Tail planning: search cheap chains among a given target set with launch anywhere in the window.
Usage: run_tail_search.py outprefix targets.json --m0 700,850,1000 [--wfuel 0.3] [--dvmax 2.5] [--lin-tofmax 0] [--beam 200] [--nproc 8]"""
import sys, json, time, pathlib, argparse, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ctoc14.kepler import Ephemeris
from ctoc14.search import beam_search, Params, save_tour, tour_json
from ctoc14.constants import DAY, AU, cost_sc
ap = argparse.ArgumentParser(); ap.add_argument('prefix'); ap.add_argument('targets'); ap.add_argument('--m0', default='700,850,1000')
ap.add_argument('--wfuel', type=float, default=0.3); ap.add_argument('--dvmax', type=float, default=2.5); ap.add_argument('--eta', type=float, default=0.6)
ap.add_argument('--beam', type=int, default=200); ap.add_argument('--nproc', type=int, default=8); ap.add_argument('--mmargin', type=float, default=40)
ap.add_argument('--tofmax', type=float, default=600); ap.add_argument('--launch-max', type=float, default=5400); ap.add_argument('--launch-step', type=float, default=15)
ap.add_argument('--lin-tofmax', type=float, default=0); ap.add_argument('--tof-refine', action='store_true'); ap.add_argument('--vinf-cap', type=float, default=4.0, help='max |v_inf| of launch legs (km/s)'); ap.add_argument('--lin-step', type=float, default=10); ap.add_argument('--lin-drmax', type=float, default=0.3)
ap.add_argument('--pool-out', default=None, help='append every beam state of every run to this jsonl (route pool)')
a = ap.parse_args()
targets = set(json.load(open(a.targets))); m0s = [float(x) for x in a.m0.split(',')]
eph = Ephemeris(); excluded = set(range(1, 301)) - targets
lin = np.arange(20, a.lin_tofmax + 1e-9, a.lin_step) * DAY if a.lin_tofmax > 0 else None
for m0 in m0s:
    P = Params(beam=a.beam, eta=a.eta, w_fuel=a.wfuel, tofs=np.arange(10, a.tofmax + 1e-9, 5) * DAY, m_margin=a.mmargin, dv_max=a.dvmax,
               lin_tofs=lin, lin_drmax=a.lin_drmax * AU, vinf_cap=a.vinf_cap, tof_refine=a.tof_refine)
    if a.pool_out: P.collect = []
    tic = time.time()
    best, beam = beam_search(eph, excluded=excluded, m0=m0, t_launch_grid=np.arange(0, a.launch_max, a.launch_step) * DAY, P=P, n_proc=a.nproc, verbose=False)
    if best is None: print(f'm0={m0}: nothing'); continue
    print(f'm0={m0} (cost {cost_sc(m0):.3f}): {best.n()} targets {[x[0] for x in best.seq]} fuel {best.fuel:.0f} kg launch {best.t_launch/DAY:.0f} d end {best.t/DAY:.0f} d ({time.time()-tic:.0f}s)', flush=True)
    save_tour(best, f'{a.prefix}_m{int(m0)}.json')
    if a.pool_out:
        seen = set(); n = 0
        with open(a.pool_out, 'a') as f:
            for s in P.collect + list(beam):
                key = (s.visited, round(s.t_launch))
                if key in seen: continue
                seen.add(key); f.write(json.dumps(tour_json(s)) + '\n'); n += 1
        print(f'  pool: +{n} routes')
