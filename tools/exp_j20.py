"""Feasibility experiments for J < 20 (docs/j20_research.md).
Phases: A tank-size sweep (single craft, all targets); B fuel-free ceiling (time-limited count); C beam width;
D eta; E sequential fleets with small tanks (lower bound on joint planning). Usage: exp_j20.py outdir [phases]"""
import sys, json, time, pathlib, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import ctoc14.search as S
from ctoc14.kepler import Ephemeris
from ctoc14.search import beam_search, Params, tour_json
from ctoc14.constants import DAY, cost_sc
out = pathlib.Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
phases = sys.argv[2] if len(sys.argv) > 2 else 'ABCDE'
eph = Ephemeris(); NP = 8
grid = np.arange(0, 1096, 20) * DAY
def P(beam=300, eta=0.6): return Params(beam=beam, eta=eta, w_fuel=0.7, m_margin=100.0, dv_max=2.5)
def run(tag, m0, excluded=(), beam=300, eta=0.6):
    tic = time.time(); best, _ = beam_search(eph, excluded=set(excluded), m0=m0, t_launch_grid=grid, P=P(beam, eta), n_proc=NP, verbose=False)
    tj = tour_json(best); json.dump(tj, open(out / f'{tag}.json', 'w'), indent=1)
    dv = [l['dv_est'] for l in tj['legs']]; tof = [l['tof'] / DAY for l in tj['legs']]
    line = (f'{tag}: m0={m0:.0f} n={best.n()} fuel={best.fuel:.0f} cost={cost_sc(m0):.3f} cost/flyby={cost_sc(m0)/max(best.n(),1):.4f} '
            f'dv_med={np.median(dv):.2f} tof_med={np.median(tof):.0f}d launch={best.t_launch/DAY:.0f}d end={best.t/DAY/365.25:.2f}yr ({time.time()-tic:.0f}s)')
    print(line, flush=True); open(out / 'log.txt', 'a').write(line + '\n'); return tj
if 'A' in phases:
    for m0 in (900, 1100, 1300, 1650, 2000): run(f'A_m{m0}', m0)
if 'B' in phases:
    VE0 = S.VE; S.VE = 1e12; run('B_fuelfree_m2000', 2000); S.VE = VE0
if 'C' in phases:
    run('C_beam1000', 2000, beam=1000)
if 'D' in phases:
    run('D_eta08', 2000, eta=0.8)
if 'E' in phases:
    for m0 in (1300, 1650):
        cov = set(); tours = []
        for k in range(1, 13):
            tj = run(f'E_m{m0}_sc{k}', m0, excluded=cov); tours.append(tj); cov |= {l['ast'] for l in tj['legs']}
            J = sum(cost_sc(t['m0']) for t in tours) + 300 - len(cov)
            line = f'   fleet m0={m0}: {k} craft covered {len(cov)} planned J={J:.2f}'; print(line, flush=True); open(out / 'log.txt', 'a').write(line + '\n')
            if tj['n_flybys'] <= 2: break
