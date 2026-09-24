"""Sequential fleet of PATIENT small-tank craft: high fuel weight, long TOF grid, all launched early.
Usage: run_swarm_fleet.py outdir --m0 900 --wfuel 3 --tofmax 600 [--lin-tofmax 900] [--max-sc 16] [--nproc 4] [--mmargin 20]"""
import sys, json, time, pathlib, argparse, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ctoc14.kepler import Ephemeris
from ctoc14.search import beam_search, Params, tour_json
from ctoc14.constants import DAY, cost_sc, AU
ap = argparse.ArgumentParser(); ap.add_argument('outdir'); ap.add_argument('--m0', type=float, default=900); ap.add_argument('--wfuel', type=float, default=3.0)
ap.add_argument('--tofmax', type=float, default=600); ap.add_argument('--tofstep', type=float, default=5); ap.add_argument('--lin-tofmax', type=float, default=0)
ap.add_argument('--max-sc', type=int, default=16); ap.add_argument('--nproc', type=int, default=4); ap.add_argument('--mmargin', type=float, default=20)
ap.add_argument('--beam', type=int, default=300); ap.add_argument('--dvmax', type=float, default=2.5); ap.add_argument('--launch-max', type=float, default=400)
ap.add_argument('--exclude', default=None)
a = ap.parse_args(); out = pathlib.Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
eph = Ephemeris()
lin = np.arange(20, a.lin_tofmax + 1e-9, 10) * DAY if a.lin_tofmax > 0 else None
P = Params(beam=a.beam, w_fuel=a.wfuel, m_margin=a.mmargin, dv_max=a.dvmax, tofs=np.arange(15, a.tofmax + 1e-9, a.tofstep) * DAY, lin_tofs=lin, lin_drmax=0.3 * AU)
grid = np.arange(0, a.launch_max, 20) * DAY
cov = set(json.load(open(a.exclude))) if a.exclude else set(); tours = []; cost = 0.0
for k in range(1, a.max_sc + 1):
    tic = time.time(); best, _ = beam_search(eph, excluded=cov, m0=a.m0, t_launch_grid=grid, P=P, n_proc=a.nproc, verbose=False)
    if best is None or best.n() == 0: break
    tj = tour_json(best); json.dump(tj, open(out / f'tour_sc{k}.json', 'w'), indent=1); tours.append(tj)
    new = {l['ast'] for l in tj['legs']} - cov; cov |= new; cost += cost_sc(a.m0)
    tof = [l['tof'] / DAY for l in tj['legs']]
    line = (f'SC{k:2d}: {best.n():2d} flybys ({len(new)} new), fuel {best.fuel:.0f}/{a.m0-600:.0f} kg, launch {best.t_launch/DAY:.0f} d, end {best.t/DAY/365.25:.1f} yr, '
            f'tof_med {np.median(tof):.0f} d, kg/flyby {best.fuel/max(best.n(),1):.1f} | fleet: covered {len(cov)}, cost {cost:.2f}, planned J {cost + 300 - len(cov):.2f} ({time.time()-tic:.0f}s)')
    print(line, flush=True); open(out / 'log.txt', 'a').write(line + '\n')
    if best.n() <= 2: break
