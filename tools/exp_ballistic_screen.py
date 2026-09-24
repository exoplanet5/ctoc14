"""Monte Carlo ballistic screening (CTOC14FullSolver.md section 16): sample Earth departures (t_L, v_inf), propagate the
ballistic orbit through the mission and count asteroids that pass close to (a) the spacecraft itself (min distance over time)
and (b) the spacecraft's orbital PATH (phase-free distance, cheap to fix by along-track drift given lead time).
Usage: exp_ballistic_screen.py out.json [--n 400] [--seed 0] [--dt-days 2]"""
import sys, json, time, pathlib, argparse, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from scipy.spatial import cKDTree
from ctoc14.kepler import Ephemeris, propagate_twobody
from ctoc14.constants import DAY, T_MISSION, AU, VINF_MAX, MU
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--n', type=int, default=400); ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--dt-days', type=float, default=2.0); ap.add_argument('--launch-max-yr', type=float, default=3.0)
a = ap.parse_args()
eph = Ephemeris(); rng = np.random.default_rng(a.seed)
ts = np.arange(0, T_MISSION, a.dt_days * DAY)
tic = time.time()
A = np.stack([eph.all_ast_states(t)[0] for t in ts])            # (nT, 300, 3)
print(f'asteroid grid {A.shape} in {time.time()-tic:.1f} s', flush=True)
unreach = np.zeros(300, bool); unreach[[130, 143]] = True
rows = []
def sample_vinf(k):
    if k == 0: return np.zeros(3)
    u = rng.normal(size=3); u /= np.linalg.norm(u)
    return u * rng.uniform(0, VINF_MAX)
for k in range(a.n):
    tL = rng.uniform(0, a.launch_max_yr * 365.25) * DAY if k > 0 else 0.0
    rE, vE = eph.earth_state(tL); vinf = sample_vinf(k); v0 = vE + vinf
    sel = ts > tL; dts = ts[sel] - tL
    S, _ = propagate_twobody(rE, v0, dts)                          # (nT',3)
    D = np.linalg.norm(A[sel] - S[:, None, :], axis=-1)            # (nT',300)
    dmin = D.min(0); dmin[unreach] = np.inf
    # path (phase-free) distance: one orbital period of the ballistic orbit sampled at 1 d
    E = 0.5 * v0 @ v0 - MU / np.linalg.norm(rE); sma = -MU / (2 * E); per = 2 * np.pi * np.sqrt(sma ** 3 / MU)
    Pth, _ = propagate_twobody(rE, v0, np.arange(0, per, DAY))
    tree = cKDTree(Pth)
    dpath, _ = tree.query(A[sel].reshape(-1, 3), workers=1); dpath = dpath.reshape(D.shape)
    dpmin = dpath.min(0); dpmin[unreach] = np.inf
    rows.append(dict(tL_d=tL / DAY, vinf=list(map(float, vinf)), vinf_mag=float(np.linalg.norm(vinf)), sma_au=float(sma / AU),
                     n_d001=int((dmin < 0.01 * AU).sum()), n_d003=int((dmin < 0.03 * AU).sum()), n_d010=int((dmin < 0.1 * AU).sum()),
                     n_p0005=int((dpmin < 0.005 * AU).sum()), n_p001=int((dpmin < 0.01 * AU).sum()), n_p002=int((dpmin < 0.02 * AU).sum()),
                     n_p005=int((dpmin < 0.05 * AU).sum())))
    if k % 25 == 0:
        r = rows[-1]; print(f'{k:4d} tL={r["tL_d"]:5.0f} d |vinf|={r["vinf_mag"]:.2f} a={r["sma_au"]:.3f}: self<0.01/0.03/0.1 AU {r["n_d001"]:3d}/{r["n_d003"]:3d}/{r["n_d010"]:3d}; '
              f'path<0.005/0.01/0.02/0.05 AU {r["n_p0005"]:3d}/{r["n_p001"]:3d}/{r["n_p002"]:3d}/{r["n_p005"]:3d}  ({time.time()-tic:.0f} s)', flush=True)
json.dump(rows, open(a.out, 'w'), indent=1)
import statistics as st
for key in ['n_d001', 'n_d003', 'n_d010', 'n_p0005', 'n_p001', 'n_p002', 'n_p005']:
    v = [r[key] for r in rows]; print(f'{key}: mean {st.mean(v):6.1f} max {max(v):3d} (orbit tL={rows[int(np.argmax(v))]["tL_d"]:.0f} |vinf|={rows[int(np.argmax(v))]["vinf_mag"]:.2f})')
