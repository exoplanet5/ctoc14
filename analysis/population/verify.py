"""
Independent cross-checks of the population_analysis.py outputs.

1. Propagator vs the two sample lines of the problem PDF (sec. 6.2).
2. Propagator vs the lead's ctoc14/kepler.py Ephemeris (independent implementation).
3. Kepler-equation residual over the full mission grid.
4. Fine-grid (0.1 d) brute-force recount of perihelion passages, ecliptic crossings,
   r_min/r_max and Earth-distance minimum, compared with the analytic/refined values in population_stats.csv.
5. Nodal-crossing table: |z| at the analytic crossing epochs, r at the node.
6. MOID: independent brute-force grid (true-anomaly parametrisation, 2400 x 1200) + local refinement.
7. Class/count re-derivations from the CSVs (per-year in-band crossings etc.).
Run:  ~/.venvs/astro313/bin/python analysis/population/verify.py
"""
import os, sys, json, time
import numpy as np, pandas as pd
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE); sys.path.insert(0, ROOT)
from kepler import (load_mea, propagate, earth_state, mean_anomaly_at, solve_kepler, state_from_E, rotation_matrix,
                    AU, DAY, MJD_T0, MJD_END, T_MAX_S, MJD_EPH_AST, MJD_EPH_EARTH, EARTH_ELEMENTS, YEAR_D, mjd_to_date)
from ctoc14.kepler import Ephemeris

t_start = time.time()
ids, el = load_mea(os.path.join(ROOT, 'MEA.txt'))
N = len(ids)
stats = pd.read_csv(os.path.join(HERE, 'population_stats.csv'))
nodes = pd.read_csv(os.path.join(HERE, 'nodal_crossings.csv'))
ca = pd.read_csv(os.path.join(HERE, 'earth_close_approaches.csv'))
out = {}

# ---------- 1. PDF sample lines ----------
r1 = np.array([-1.9500328647e+07, 1.4581848347e+08, -7.6372048832e+03]); v1 = np.array([-2.9389477847e+01, -4.3166912970e+00, -3.9417153392e+00])
rE, vE = earth_state(MJD_T0)
t2 = 1.0008479152e+07
r2 = np.array([-1.1398611209e+08, -8.7213961221e+07, -1.6523817305e+07])
ra, va = propagate(el[174 - 1], MJD_EPH_AST, MJD_T0 + t2 / DAY)
out['pdf_check'] = dict(earth_pos_err_km=float(np.linalg.norm(rE - r1)), launch_vinf_kms=float(np.linalg.norm(v1 - vE)),
                        ast174_dist_to_sample_flyby_km=float(np.linalg.norm(ra - r2)))
print('1. PDF sample: Earth(t0) err = %.4f km, v_inf = %.6f km/s, ast174 dist to sample flyby pos = %.4f km'
      % (out['pdf_check']['earth_pos_err_km'], out['pdf_check']['launch_vinf_kms'], out['pdf_check']['ast174_dist_to_sample_flyby_km']))

# ---------- 2. vs lead's Ephemeris ----------
eph = Ephemeris(os.path.join(ROOT, 'MEA.txt'))
rng = np.random.default_rng(1)
ts = np.concatenate([[0.0, T_MAX_S], rng.uniform(0, T_MAX_S, 60)])
dp = 0.0; dv = 0.0; dpe = 0.0
for t in ts:
    ra_all, va_all = eph.all_ast_states(t)
    pa, va = propagate(el, MJD_EPH_AST, MJD_T0 + t / DAY)
    dp = max(dp, np.linalg.norm(pa - ra_all, axis=1).max()); dv = max(dv, np.linalg.norm(va - va_all, axis=1).max())
    re_, ve_ = eph.earth_state(t); pe, _ = earth_state(MJD_T0 + t / DAY)
    dpe = max(dpe, np.linalg.norm(pe - re_))
out['vs_lead_ephemeris'] = dict(n_times=len(ts), max_pos_diff_km=float(dp), max_vel_diff_kms=float(dv), max_earth_pos_diff_km=float(dpe))
print('2. vs ctoc14.kepler (62 epochs x 300): max |dr| = %.3e km, max |dv| = %.3e km/s, Earth max |dr| = %.3e km' % (dp, dv, dpe))

# ---------- 3. Kepler residual ----------
t_grid = np.arange(MJD_T0, MJD_END + 1e-9, 0.5)
M = mean_anomaly_at(el[:, None, :], MJD_EPH_AST, t_grid[None, :])
E = solve_kepler(M, el[:, None, 1])
res = np.abs(E - el[:, None, 1] * np.sin(E) - np.mod(M, 2 * np.pi)).max()
out['kepler_residual_max'] = float(res)
print('3. Kepler residual max over %d x %d = %.2e rad' % (N, len(t_grid), res))

# ---------- 4. fine-grid brute force ----------
DT = 0.1
tf = np.arange(MJD_T0, MJD_END + 1e-9, DT)
if tf[-1] < MJD_END: tf = np.append(tf, MJD_END)
pe_f, _ = earth_state(tf)
n_peri_g = np.zeros(N, int); n_cross_g = np.zeros(N, int); rmin_g = np.zeros(N); rmax_g = np.zeros(N)
dmin_g = np.zeros(N); tdmin_g = np.zeros(N); zmax_g = np.zeros(N)
CH = 30
for c0 in range(0, N, CH):
    sl = slice(c0, min(c0 + CH, N))
    pa, _ = propagate(el[sl, None, :], MJD_EPH_AST, tf[None, :])
    r = np.linalg.norm(pa, axis=-1) / AU; z = pa[:, :, 2] / AU
    d = np.linalg.norm(pa - pe_f[None], axis=-1) / AU
    loc = (r[:, 1:-1] < r[:, :-2]) & (r[:, 1:-1] <= r[:, 2:])
    n_peri_g[sl] = loc.sum(axis=1)
    n_cross_g[sl] = (np.sign(z[:, :-1]) * np.sign(z[:, 1:]) < 0).sum(axis=1)
    rmin_g[sl] = r.min(axis=1); rmax_g[sl] = r.max(axis=1); zmax_g[sl] = np.abs(z).max(axis=1)
    k = d.argmin(axis=1); dmin_g[sl] = d[np.arange(d.shape[0]), k]; tdmin_g[sl] = tf[k]
d_peri = n_peri_g - stats.n_perihelion_passages.values
d_cross = n_cross_g - stats.n_ecliptic_crossings.values
out['fine_grid'] = dict(
    dt_d=DT, n_samples=int(len(tf)),
    perihelion_count_mismatch_ids=[int(x) for x in ids[d_peri != 0]],
    ecliptic_crossing_count_mismatch_ids=[int(x) for x in ids[d_cross != 0]],
    rmin_grid_minus_analytic_au=dict(max=float((rmin_g - stats.r_min_window_au.values).max()), min=float((rmin_g - stats.r_min_window_au.values).min())),
    rmax_analytic_minus_grid_au=dict(max=float((stats.r_max_window_au.values - rmax_g).max()), min=float((stats.r_max_window_au.values - rmax_g).min())),
    dearth_grid_minus_refined_au=dict(max=float((dmin_g - stats.d_earth_min_window_au.values).max()), min=float((dmin_g - stats.d_earth_min_window_au.values).min())),
    dearth_time_diff_d_max=float(np.abs(tdmin_g - stats.t_d_earth_min_mjd.values).max()),
    zabs_max_diff_au_max=float(np.abs(zmax_g - stats.zabs_max_window_au.values).max()),
    total_perihelia_grid=int(n_peri_g.sum()), total_crossings_grid=int(n_cross_g.sum()))
print('4. fine grid (%.1f d, %d samples): perihelion-count mismatches: %s; crossing-count mismatches: %s'
      % (DT, len(tf), out['fine_grid']['perihelion_count_mismatch_ids'], out['fine_grid']['ecliptic_crossing_count_mismatch_ids']))
print('   r_min(grid)-r_min(analytic) in [%.2e, %.2e] AU (>=0 expected); r_max(analytic)-r_max(grid) in [%.2e, %.2e] AU'
      % (out['fine_grid']['rmin_grid_minus_analytic_au']['min'], out['fine_grid']['rmin_grid_minus_analytic_au']['max'],
         out['fine_grid']['rmax_analytic_minus_grid_au']['min'], out['fine_grid']['rmax_analytic_minus_grid_au']['max']))
print('   d_earth_min(grid)-d_earth_min(refined) in [%.2e, %.2e] AU; |t diff| max %.3f d; totals grid: %d perihelia, %d crossings'
      % (out['fine_grid']['dearth_grid_minus_refined_au']['min'], out['fine_grid']['dearth_grid_minus_refined_au']['max'],
         out['fine_grid']['dearth_time_diff_d_max'], n_peri_g.sum(), n_cross_g.sum()))
if len(out['fine_grid']['ecliptic_crossing_count_mismatch_ids']):
    for x in out['fine_grid']['ecliptic_crossing_count_mismatch_ids']:
        j = x - 1
        print('     id %d: grid %d vs analytic %d (asc %d, desc %d); i=%.2f' % (x, n_cross_g[j], stats.n_ecliptic_crossings[j], stats.n_asc_crossings[j], stats.n_desc_crossings[j], el[j, 2]))
if len(out['fine_grid']['perihelion_count_mismatch_ids']):
    for x in out['fine_grid']['perihelion_count_mismatch_ids']:
        j = x - 1
        print('     id %d: grid %d vs analytic %d; dates %s' % (x, n_peri_g[j], stats.n_perihelion_passages[j], stats.perihelion_dates_mjd[j]))

# ---------- 5. nodal table ----------
zk = np.empty(len(nodes)); rk = np.empty(len(nodes))
for k, row in enumerate(nodes.itertuples()):
    pa, _ = propagate(el[row.id - 1], MJD_EPH_AST, row.t_mjd)
    zk[k] = pa[2]; rk[k] = np.linalg.norm(pa) / AU
out['nodes'] = dict(n=int(len(nodes)), max_abs_z_km=float(np.abs(zk).max()), max_abs_r_err_au=float(np.abs(rk - nodes.r_au.values).max()),
                    in_window=bool(((nodes.t_mjd >= MJD_T0) & (nodes.t_mjd <= MJD_END)).all()),
                    asc_z_slope_positive=bool((nodes[nodes.node == 'asc'].vz_kms > 0).all()), desc_z_slope_negative=bool((nodes[nodes.node == 'desc'].vz_kms < 0).all()))
print('5. nodal table: %d rows, max |z| at crossing = %.3e km, max |r - r_node| = %.2e AU, all in window: %s, asc vz>0: %s, desc vz<0: %s'
      % (len(nodes), out['nodes']['max_abs_z_km'], out['nodes']['max_abs_r_err_au'], out['nodes']['in_window'], out['nodes']['asc_z_slope_positive'], out['nodes']['desc_z_slope_negative']))

# ---------- 6. MOID brute force (true anomaly parametrisation) ----------
def ellipse_points_nu(elm, nu):
    a, e, i, Om, w, _ = elm
    p = a * AU * (1 - e * e); r = p / (1 + e * np.cos(nu))
    xp = r * np.cos(nu); yp = r * np.sin(nu)
    R = rotation_matrix(i, Om, w)
    return R[:, 0][None, :] * xp[:, None] + R[:, 1][None, :] * yp[:, None]
nuA = np.linspace(0, 2 * np.pi, 2400, endpoint=False); nuE = np.linspace(0, 2 * np.pi, 1200, endpoint=False)
PE = ellipse_points_nu(EARTH_ELEMENTS, nuE)
moid_bf = np.empty(N)
for j in range(N):
    PA = ellipse_points_nu(el[j], nuA)
    D2 = ((PA[:, None, :] - PE[None, :, :]) ** 2).sum(-1)
    flat = np.argsort(D2, axis=None)[:6]
    best = np.sqrt(D2.min()) / AU
    fun = lambda v: np.linalg.norm(ellipse_points_nu(el[j], np.array([v[0]]))[0] - ellipse_points_nu(EARTH_ELEMENTS, np.array([v[1]]))[0]) / AU
    for f in flat:
        ia, ie = np.unravel_index(f, D2.shape)
        r_ = minimize(fun, x0=[nuA[ia], nuE[ie]], method='Nelder-Mead', options={'xatol': 1e-10, 'fatol': 1e-13, 'maxiter': 3000})
        best = min(best, r_.fun)
    moid_bf[j] = best
dm = moid_bf - stats.moid_au.values
out['moid'] = dict(max_abs_diff_au=float(np.abs(dm).max()), max_abs_diff_km=float(np.abs(dm).max() * AU),
                   ids_diff_gt_1em6_au=[int(x) for x in ids[np.abs(dm) > 1e-6]],
                   n_lt_0p05=int((moid_bf < 0.05).sum()), n_lt_0p02=int((moid_bf < 0.02).sum()), n_lt_0p01=int((moid_bf < 0.01).sum()),
                   ids_ge_0p05=[int(x) for x in ids[moid_bf >= 0.05]], min_moid_id=int(ids[moid_bf.argmin()]), min_moid_au=float(moid_bf.min()))
print('6. MOID brute force vs stored: max |diff| = %.2e AU (%.1f km); ids with |diff|>1e-6 AU: %s; <0.05: %d, <0.02: %d, <0.01: %d; >=0.05 AU: %s; min %.5f AU (id %d)'
      % (out['moid']['max_abs_diff_au'], out['moid']['max_abs_diff_km'], out['moid']['ids_diff_gt_1em6_au'], out['moid']['n_lt_0p05'],
         out['moid']['n_lt_0p02'], out['moid']['n_lt_0p01'], out['moid']['ids_ge_0p05'], out['moid']['min_moid_au'], out['moid']['min_moid_id']))

# ---------- 7. re-derivations from CSVs ----------
years = np.arange(2030, 2045); edges = np.append(MJD_T0 + (years - 2030) * YEAR_D, MJD_END + 1e-6)
nb = nodes[(nodes.r_au >= 0.8) & (nodes.r_au <= 1.6)]
per_year = np.histogram(nb.t_mjd, bins=edges)[0]
cls = stats.difficulty.values
out['rederived'] = dict(
    classes={c: int((cls == c).sum()) for c in ['easy', 'moderate', 'hard', 'very_hard', 'unreachable']},
    in_band_crossings=int(len(nb)), in_band_per_year=[int(x) for x in per_year], in_band_per_year_mean=float(per_year.mean()), in_band_per_year_std=float(per_year.std()),
    in_band_asteroids=int(nb.id.nunique()),
    in_band_crossings_by_class={c: int(nb.id.isin(ids[cls == c]).sum()) for c in ['easy', 'moderate', 'hard']},
    near_1au_0p9_1p1=int(((nodes.r_au >= 0.9) & (nodes.r_au <= 1.1)).sum()),
    earth_ca_lt_0p1=int((ca.d_au <= 0.1).sum()), earth_ca_lt_0p05=int((ca.d_au <= 0.05).sum()), earth_ca_lt_0p02=int((ca.d_au <= 0.02).sum()),
    atira_ids=[int(x) for x in ids[(el[:, 0] * (1 + el[:, 1])) < 0.983]],
    unreachable_ids=[int(x) for x in ids[cls == 'unreachable']],
    ids_rmin_window_gt_2p5=[int(x) for x in ids[stats.r_min_window_au.values > 2.5]],
    ids_retrograde=[int(x) for x in ids[el[:, 2] > 90]],
    all_q_le_1p3=bool((stats.q_au <= 1.3).all()), all_cross_annulus=bool(stats.crosses_annulus_0p7_1p5.all()),
    n_rad_s_col_ok=bool((stats.n_rad_s.values > 0).all()))
for k in (131, 144):
    j = k - 1
    tt = np.array([MJD_T0, MJD_T0 + 5 * YEAR_D, MJD_T0 + 10 * YEAR_D, MJD_END])
    pa, _ = propagate(el[j], MJD_EPH_AST, tt); pe, _ = earth_state(tt)
    out['rederived'][f'id{k}'] = dict(r_au_at=[float(x) for x in np.linalg.norm(pa, axis=1) / AU], dates=[mjd_to_date(t) for t in tt],
                                      d_earth_au_at=[float(x) for x in np.linalg.norm(pa - pe, axis=1) / AU],
                                      d_earth_min_finegrid_au=float(dmin_g[j]), r_min_finegrid_au=float(rmin_g[j]))
print('7. classes %s; in-band crossings %d (per-year mean %.1f, std %.1f) from %d asteroids; by class %s; near-1AU %d; Earth CA <0.1/<0.05/<0.02: %d/%d/%d'
      % (out['rederived']['classes'], len(nb), per_year.mean(), per_year.std(), nb.id.nunique(), out['rederived']['in_band_crossings_by_class'],
         out['rederived']['near_1au_0p9_1p1'], out['rederived']['earth_ca_lt_0p1'], out['rederived']['earth_ca_lt_0p05'], out['rederived']['earth_ca_lt_0p02']))
for k in (131, 144):
    print('   id %d: r(AU) at %s = %s; d_earth = %s; fine-grid min d_earth %.3f AU, min r %.3f AU'
          % (k, out['rederived'][f'id{k}']['dates'], np.round(out['rederived'][f'id{k}']['r_au_at'], 3), np.round(out['rederived'][f'id{k}']['d_earth_au_at'], 3),
             out['rederived'][f'id{k}']['d_earth_min_finegrid_au'], out['rederived'][f'id{k}']['r_min_finegrid_au']))
print('   atira ids %s, retrograde %s, r_min>2.5 %s, unreachable %s, all q<=1.3: %s, all cross annulus: %s, n_rad_s column ok: %s'
      % (out['rederived']['atira_ids'], out['rederived']['ids_retrograde'], out['rederived']['ids_rmin_window_gt_2p5'], out['rederived']['unreachable_ids'],
         out['rederived']['all_q_le_1p3'], out['rederived']['all_cross_annulus'], out['rederived']['n_rad_s_col_ok']))
out['elapsed_s'] = round(time.time() - t_start, 1)
with open(os.path.join(HERE, 'verify.json'), 'w') as f:
    json.dump(out, f, indent=1)
print('verify done in %.1f s' % out['elapsed_s'])
