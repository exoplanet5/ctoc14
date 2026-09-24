"""Joint multi-craft (swarm) beam search. Usage: run_joint.py outdir --K 12 --m0 1000 --wfuel 3 --beam 200 --nproc 6 [--tofmax 400] [--launch-max 400] [--mmargin 20] [--exclude ids.json]"""
import sys, json, time, pathlib, argparse, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ctoc14.kepler import Ephemeris
from ctoc14.search import Params
from ctoc14.jointsearch import joint_beam_search, joint_tours
from ctoc14.constants import DAY, cost_sc
ap = argparse.ArgumentParser(); ap.add_argument('outdir'); ap.add_argument('--K', type=int, default=12); ap.add_argument('--m0', type=float, default=1000)
ap.add_argument('--m0list', default=None, help='comma-separated m0 per craft (overrides --m0/--K)')
ap.add_argument('--wfuel', type=float, default=3.0); ap.add_argument('--beam', type=int, default=200); ap.add_argument('--nproc', type=int, default=6)
ap.add_argument('--tofmax', type=float, default=400); ap.add_argument('--tofstep', type=float, default=5); ap.add_argument('--launch-max', type=float, default=400)
ap.add_argument('--launch-step', type=float, default=20); ap.add_argument('--mmargin', type=float, default=20); ap.add_argument('--dvmax', type=float, default=2.5)
ap.add_argument('--exclude', default=None); ap.add_argument('--max-children', type=int, default=60)
ap.add_argument('--wait', type=float, default=120, help='coast step (days) when no leg is feasible'); ap.add_argument('--parent-cap', type=int, default=0)
ap.add_argument('--wrare', type=float, default=0.0); ap.add_argument('--rarity-csv', default='analysis/ballistic/out/junction_per_asteroid.csv')
ap.add_argument('--lin-tofmax', type=float, default=0, help='linearised long-leg model TOF max (days), 0=off'); ap.add_argument('--tof-refine', action='store_true'); ap.add_argument('--vinf-cap', type=float, default=4.0, help='max |v_inf| of launch legs (km/s)'); ap.add_argument('--lin-step', type=float, default=10); ap.add_argument('--lin-drmax', type=float, default=0.3, help='AU')
ap.add_argument('--rarity-mode', default='junction', choices=['junction', 'crossings'], help='crossings: rarity = min(1, 4 / n_node_crossings_in_band) from population_stats.csv')
a = ap.parse_args(); out = pathlib.Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
m0s = [float(x) for x in a.m0list.split(',')] if a.m0list else [a.m0] * a.K
rarity = None
if a.wrare > 0 and a.rarity_mode == 'crossings':
    import pandas as pd
    d = pd.read_csv('analysis/population/population_stats.csv').sort_values('id')
    rarity = np.clip(3.0 / np.maximum(d["n_node_crossings_in_band"].to_numpy(float), 1.0), 0, 1) ** 2
elif a.wrare > 0:
    import pandas as pd
    d = pd.read_csv(a.rarity_csv).sort_values('B'); f = d['frac_feas_plausible_tof180_a0.5'].to_numpy()
    rarity = np.nan_to_num(np.clip(1.0 - f / np.nanmax(f), 0, 1), nan=0.5)
    print(f'rarity weighting on: mean {rarity.mean():.2f}, hardest {np.argsort(-rarity)[:12] + 1}')
from ctoc14.constants import AU
lin = np.arange(20, a.lin_tofmax + 1e-9, a.lin_step) * DAY if a.lin_tofmax > 0 else None
P = Params(beam=a.beam, w_fuel=a.wfuel, m_margin=a.mmargin, dv_max=a.dvmax, tofs=np.arange(15, a.tofmax + 1e-9, a.tofstep) * DAY, max_children=a.max_children,
           w_rare=a.wrare, rarity=rarity, lin_tofs=lin, lin_drmax=a.lin_drmax * AU, vinf_cap=a.vinf_cap, tof_refine=a.tof_refine)
P.wait = a.wait * DAY; P.parent_cap = a.parent_cap
excl = set(json.load(open(a.exclude))) if a.exclude else set()
eph = Ephemeris(); tic = time.time()
best = joint_beam_search(eph, m0s=m0s, launch_grid=np.arange(0, a.launch_max + 1e-9, a.launch_step) * DAY, P=P, excluded=excl, n_proc=a.nproc, verbose=True)
tours = joint_tours(best)
cov = set().union(*[{l['ast'] for l in t['legs']} for t in tours]) if tours else set()
J_fixed = sum(cost_sc(t['m0']) for t in tours) + 300 - len(cov) - len(excl)
# right-sized tanks: m0 = 600 + fuel_est + margin
J_rs = sum(cost_sc(min(2000.0, 600 + t['fuel_est'] + 30)) for t in tours) + 300 - len(cov) - len(excl)
for k, t in enumerate(tours, 1):
    json.dump(t, open(out / f'tour_sc{k}.json', 'w'), indent=1)
    dv = [l['dv_est'] for l in t['legs']]; tof = [l['tof'] / DAY for l in t['legs']]
    print(f'SC{k:2d}: {t["n_flybys"]:2d} flybys, fuel {t["fuel_est"]:.0f} kg, launch {t["t_launch"]/DAY:.0f} d, end {t["t_end"]/DAY/365.25:.1f} yr, dv_med {np.median(dv):.2f}, tof_med {np.median(tof):.0f} d, kg/flyby {t["fuel_est"]/t["n_flybys"]:.1f}')
print(f'FLEET: {len(tours)} craft, covered {len(cov)} (+{len(excl)} excluded), planned J (fixed tanks) = {J_fixed:.2f}, right-sized tanks = {J_rs:.2f}, runtime {time.time()-tic:.0f} s')
json.dump(dict(m0s=m0s, covered=sorted(cov), J_fixed=J_fixed, J_rs=J_rs), open(out / 'summary.json', 'w'))
