"""Rate vs pool size: beam search on RANDOM target subsets of size D (vs greedy leftovers of the same size)."""
import sys, json, time, pathlib, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ctoc14.kepler import Ephemeris
from ctoc14.search import beam_search, Params, tour_json
from ctoc14.constants import DAY
eph = Ephemeris(); rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
grid = np.arange(0, 400, 20) * DAY; P = Params(beam=300, w_fuel=0.7, m_margin=100.0, dv_max=2.5)
allids = [i for i in range(1, 301) if i not in (131, 144)]
for D in (250, 200, 150, 100, 60, 30):
    keep = set(rng.choice(allids, D, replace=False).tolist()); excl = set(range(1, 301)) - keep
    tic = time.time(); best, _ = beam_search(eph, excluded=excl, m0=2000.0, t_launch_grid=grid, P=P, n_proc=8, verbose=False)
    tj = tour_json(best); json.dump(tj, open(f'results/exp_j20/random_D{D}.json', 'w'))
    dv = [l['dv_est'] for l in tj['legs']]; tof = [l['tof'] / DAY for l in tj['legs']]
    print(f'random pool D={D:3d}: {best.n():2d} flybys, fuel {best.fuel:.0f} kg, end {best.t/DAY/365.25:.1f} yr, dv_med {np.median(dv):.2f}, tof_med {np.median(tof):.0f} d ({time.time()-tic:.0f}s)', flush=True)
