"""
CTOC14 target-population analysis.

Computes per-asteroid orbital statistics, mission-window events (perihelion passages,
ecliptic crossings, min/max heliocentric distance, Earth close approaches, MOID),
accessibility/difficulty classification, and target opportunity density.

Outputs (all under analysis/population/):
    population_stats.csv        one row per asteroid
    nodal_crossings.csv         one row per ecliptic-plane crossing in the window
    earth_close_approaches.csv  one row per local minimum of Earth distance < 0.3 AU in the window
    access_windows.csv          one row per contiguous interval with 0.8<=r<=1.6 AU and |z|<=0.1 AU
    summary.json                aggregate counts used by docs/population_analysis.md
    *.png                       plots
Run:  ~/.venvs/astro313/bin/python analysis/population/population_analysis.py
"""
import os, sys, json, time
T_START = time.time()
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar, minimize

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
from kepler import (load_mea, propagate, earth_state, mean_motion, period_days, true_to_mean_anomaly,
                    state_from_E, solve_kepler, mjd_to_date, mjd_to_decimal_year,
                    AU, DAY, MJD_T0, MJD_END, T_MAX_D, MJD_EPH_AST, MJD_EPH_EARTH, EARTH_ELEMENTS, YEAR_D)

# ---------------- thresholds (documented in docs/population_analysis.md) ----------------
BAND_LO, BAND_HI = 0.8, 1.6        # heliocentric band a ~1 AU spacecraft can sweep cheaply
ANNULUS_LO, ANNULUS_HI = 0.7, 1.5  # annulus requested in the task
Z_ACCESS = 0.10                    # |z| <= 0.10 AU counts as "near the ecliptic" for access windows
R_UNREACH = 2.5                    # never closer than this to the Sun in window -> unreachable
I_EASY = 10.0                      # deg
I_MOD = 20.0                       # deg
Q_EASY = 1.3                       # AU
DT_GRID = 0.5                      # days, time grid step for window scans
CLOSE_APPROACH_MAX = 0.3           # AU, record Earth approaches below this

ids, el = load_mea(os.path.join(ROOT, 'MEA.txt'))
N = len(ids)
a, e, inc, Om, w, M0 = el.T
n_rad_s = mean_motion(a)
P_d = period_days(a)
q = a * (1 - e); Q = a * (1 + e)
p = a * (1 - e ** 2)
# ---------------- node geometry (constant for fixed Keplerian orbits) ----------------
nu_asc = np.radians(-w)              # true anomaly at ascending node (u = w + nu = 0)
nu_desc = np.radians(180.0 - w)      # descending node (u = 180 deg)
r_asc = p / (1 + e * np.cos(nu_asc))
r_desc = p / (1 + e * np.cos(nu_desc))
M_asc = true_to_mean_anomaly(nu_asc, e)
M_desc = true_to_mean_anomaly(nu_desc, e)
crosses_annulus = (q <= ANNULUS_HI) & (Q >= ANNULUS_LO)
crosses_band = (q <= BAND_HI) & (Q >= BAND_LO)

# ---------------- window events from mean anomaly ----------------
M_t0 = np.radians(M0) + n_rad_s * (MJD_T0 - MJD_EPH_AST) * DAY
M_tf = np.radians(M0) + n_rad_s * (MJD_END - MJD_EPH_AST) * DAY
twopi = 2 * np.pi

def event_times(M_event):
    """All MJDs in the window at which M(t) = M_event (mod 2pi), for each asteroid -> list of arrays."""
    out = []
    for j in range(N):
        k_lo = np.ceil((M_t0[j] - M_event[j]) / twopi)
        k_hi = np.floor((M_tf[j] - M_event[j]) / twopi)
        ks = np.arange(k_lo, k_hi + 1)
        t = MJD_EPH_AST + (M_event[j] + twopi * ks - np.radians(M0[j])) / n_rad_s[j] / DAY
        out.append(t)
    return out

peri_times = event_times(np.zeros(N))
apo_times = event_times(np.full(N, np.pi))
asc_times = event_times(M_asc)
desc_times = event_times(M_desc)
n_peri = np.array([len(t) for t in peri_times])
n_asc = np.array([len(t) for t in asc_times])
n_desc = np.array([len(t) for t in desc_times])

# min/max heliocentric distance in window: extremes only at perihelion/aphelion or endpoints
def r_at(j, t_mjd):
    pos, _ = propagate(el[j], MJD_EPH_AST, np.asarray(t_mjd))
    return np.linalg.norm(pos, axis=-1) / AU
r_min_w = np.empty(N); r_max_w = np.empty(N); t_rmin = np.empty(N); t_rmax = np.empty(N)
for j in range(N):
    cand_t = np.concatenate([[MJD_T0, MJD_END], peri_times[j], apo_times[j]])
    cand_r = r_at(j, cand_t)
    r_min_w[j] = cand_r.min(); t_rmin[j] = cand_t[np.argmin(cand_r)]
    r_max_w[j] = cand_r.max(); t_rmax[j] = cand_t[np.argmax(cand_r)]

# ---------------- time-grid scan over the window ----------------
t_grid = np.arange(MJD_T0, MJD_END + 1e-9, DT_GRID)
if t_grid[-1] < MJD_END: t_grid = np.append(t_grid, MJD_END)
T = len(t_grid)
pos_a, vel_a = propagate(el[:, None, :], MJD_EPH_AST, t_grid[None, :])    # (N,T,3)
pos_e, vel_e = earth_state(t_grid)                                           # (T,3)
r_grid = np.linalg.norm(pos_a, axis=-1) / AU
z_grid = pos_a[:, :, 2] / AU
d_earth_grid = np.linalg.norm(pos_a - pos_e[None, :, :], axis=-1) / AU
vrel_grid = np.linalg.norm(vel_a - vel_e[None, :, :], axis=-1)
in_band = (r_grid >= BAND_LO) & (r_grid <= BAND_HI)
near_ecl = np.abs(z_grid) <= Z_ACCESS
access = in_band & near_ecl
frac_in_band = in_band.mean(axis=1)
frac_access = access.mean(axis=1)
days_in_band = frac_in_band * T_MAX_D
days_access = frac_access * T_MAX_D
zmax_grid = np.abs(z_grid).max(axis=1)
r_at_t0 = r_grid[:, 0]; z_at_t0 = z_grid[:, 0]

def intervals(mask, t):
    """Contiguous True runs of mask along t -> list of (t_start, t_end)."""
    m = np.concatenate([[False], mask, [False]])
    d = np.diff(m.astype(int))
    starts = np.where(d == 1)[0]; ends = np.where(d == -1)[0] - 1
    return [(t[s], t[e]) for s, e in zip(starts, ends)]

access_rows = []
n_access_windows = np.zeros(N, int)
for j in range(N):
    for (ts, te) in intervals(access[j], t_grid):
        seg = slice(np.searchsorted(t_grid, ts), np.searchsorted(t_grid, te) + 1)
        access_rows.append(dict(id=ids[j], t_start_mjd=ts, t_end_mjd=te, date_start=mjd_to_date(ts), date_end=mjd_to_date(te),
                                duration_d=te - ts + DT_GRID, r_mean_au=r_grid[j, seg].mean(),
                                zabs_min_au=np.abs(z_grid[j, seg]).min(), d_earth_min_au=d_earth_grid[j, seg].min()))
        n_access_windows[j] += 1
access_df = pd.DataFrame(access_rows)

# ---------------- Earth close approaches (local minima of d_earth in the window) ----------------
def d_earth_fun(j):
    def f(t):
        pa, _ = propagate(el[j], MJD_EPH_AST, t); pe, _ = earth_state(t)
        return np.linalg.norm(pa - pe) / AU
    return f

ca_rows = []
d_earth_min = np.empty(N); t_d_earth_min = np.empty(N); vrel_at_min = np.empty(N)
n_ca_010 = np.zeros(N, int); n_ca_005 = np.zeros(N, int); n_ca_002 = np.zeros(N, int)
for j in range(N):
    d = d_earth_grid[j]
    f = d_earth_fun(j)
    # local minima (incl. endpoints)
    loc = np.where((d[1:-1] <= d[:-2]) & (d[1:-1] <= d[2:]))[0] + 1
    cands = list(loc)
    if d[0] < d[1]: cands.append(0)
    if d[-1] < d[-2]: cands.append(T - 1)
    best = (np.inf, None)
    for k in cands:
        if d[k] > CLOSE_APPROACH_MAX and k not in (0, T - 1) and d[k] > d.min() + 0.01:
            continue
        lo = max(t_grid[max(k - 1, 0)], MJD_T0); hi = min(t_grid[min(k + 1, T - 1)], MJD_END)
        if hi - lo < 1e-6:
            tm, dm = t_grid[k], d[k]
        else:
            res = minimize_scalar(f, bounds=(lo, hi), method='bounded', options={'xatol': 1e-6})
            tm, dm = res.x, res.fun
            if dm > d[k]: tm, dm = t_grid[k], d[k]
        pa, va = propagate(el[j], MJD_EPH_AST, tm); pe, ve = earth_state(tm)
        vr = np.linalg.norm(va - ve)
        if dm <= CLOSE_APPROACH_MAX:
            ca_rows.append(dict(id=ids[j], t_mjd=tm, date=mjd_to_date(tm), t_s=(tm - MJD_T0) * DAY, d_au=dm, d_km=dm * AU,
                                v_rel_kms=vr, r_helio_au=np.linalg.norm(pa) / AU, z_au=pa[2] / AU))
            if dm <= 0.10: n_ca_010[j] += 1
            if dm <= 0.05: n_ca_005[j] += 1
            if dm <= 0.02: n_ca_002[j] += 1
        if dm < best[0]: best = (dm, tm, vr)
    d_earth_min[j], t_d_earth_min[j], vrel_at_min[j] = best
ca_df = pd.DataFrame(ca_rows).sort_values('t_mjd').reset_index(drop=True)

# ---------------- MOID (geometric min distance between the two fixed ellipses) ----------------
def moid(j, ngrid=720):
    Ea = np.linspace(0, twopi, ngrid, endpoint=False)
    pa, _ = state_from_E(el[j], Ea)                       # (ngrid,3)
    pe, _ = state_from_E(EARTH_ELEMENTS, Ea)              # (ngrid,3)
    D = np.linalg.norm(pa[:, None, :] - pe[None, :, :], axis=-1)   # (ngrid,ngrid)
    flat = np.argsort(D, axis=None)[:40]
    ia, ie = np.unravel_index(flat, D.shape)
    # keep distinct starting cells
    starts = []
    for x, y in zip(ia, ie):
        if all(min(abs(x - sx), ngrid - abs(x - sx)) > 10 or min(abs(y - sy), ngrid - abs(y - sy)) > 10 for sx, sy in starts):
            starts.append((x, y))
        if len(starts) >= 4: break
    def fun(v):
        p1, _ = state_from_E(el[j], v[0]); p2, _ = state_from_E(EARTH_ELEMENTS, v[1])
        return np.linalg.norm(p1 - p2) / AU
    best = (D.min() / AU, Ea[ia[0]], Ea[ie[0]])
    for (x, y) in starts:
        res = minimize(fun, x0=[Ea[x], Ea[y]], method='Nelder-Mead', options={'xatol': 1e-10, 'fatol': 1e-12, 'maxiter': 2000})
        if res.fun < best[0]: best = (res.fun, res.x[0], res.x[1])
    dmin, Ea_b, Ee_b = best
    pa_b, _ = state_from_E(el[j], Ea_b)
    return dmin, np.linalg.norm(pa_b) / AU, np.degrees(np.arctan2(pa_b[1], pa_b[0])) % 360, pa_b[2] / AU

moid_au = np.empty(N); moid_r = np.empty(N); moid_lon = np.empty(N); moid_z = np.empty(N)
for j in range(N):
    moid_au[j], moid_r[j], moid_lon[j], moid_z[j] = moid(j)

# ---------------- nodal crossings table ----------------
node_rows = []
for j in range(N):
    for kind, tt, rr in (('asc', asc_times[j], r_asc[j]), ('desc', desc_times[j], r_desc[j])):
        lon = Om[j] if kind == 'asc' else (Om[j] + 180.0) % 360.0
        for t in tt:
            pa, va = propagate(el[j], MJD_EPH_AST, t); pe, ve = earth_state(t)
            node_rows.append(dict(id=ids[j], node=kind, t_mjd=t, date=mjd_to_date(t), t_s=(t - MJD_T0) * DAY,
                                  r_au=rr, lon_deg=lon, x_au=pa[0] / AU, y_au=pa[1] / AU, z_km=pa[2],
                                  d_earth_au=np.linalg.norm(pa - pe) / AU, v_rel_earth_kms=np.linalg.norm(va - ve),
                                  vz_kms=va[2], inc_deg=inc[j], in_band=(BAND_LO <= rr <= BAND_HI)))
node_df = pd.DataFrame(node_rows).sort_values('t_mjd').reset_index(drop=True)
n_node_band = node_df[node_df.in_band].groupby('id').size().reindex(ids, fill_value=0).values
n_node_wide = node_df[(node_df.r_au >= 0.6) & (node_df.r_au <= 2.0)].groupby('id').size().reindex(ids, fill_value=0).values
n_node_total = n_asc + n_desc

# ---------------- difficulty classification ----------------
retro = inc > 90.0
cls = np.empty(N, dtype=object)
for j in range(N):
    if retro[j] or r_min_w[j] > R_UNREACH:
        cls[j] = 'unreachable'
    elif inc[j] <= I_EASY and q[j] <= Q_EASY and days_access[j] > 0:
        cls[j] = 'easy'
    elif inc[j] <= I_MOD and days_access[j] > 0:
        cls[j] = 'moderate'
    elif inc[j] > I_MOD and (n_node_band[j] > 0 or days_access[j] > 0):
        cls[j] = 'hard'
    else:
        cls[j] = 'very_hard'
# secondary tag: crosses the band only through a node with r outside band, etc.
reason = []
for j in range(N):
    if cls[j] == 'unreachable':
        reason.append('retrograde' if retro[j] else f'r_min_window={r_min_w[j]:.2f}AU')
    elif cls[j] == 'very_hard':
        if days_in_band[j] == 0: reason.append(f'never in {BAND_LO}-{BAND_HI}AU band in window (r_min={r_min_w[j]:.2f})')
        elif days_access[j] == 0 and n_node_band[j] == 0: reason.append('in band but never near ecliptic; nodes at r=%.2f/%.2f AU' % (r_asc[j], r_desc[j]))
        else: reason.append('')
    else:
        reason.append('')

# ---------------- assemble per-asteroid table ----------------
def fmt_list(ts, rs=None):
    if rs is None:
        return ';'.join(f'{t:.2f}' for t in ts)
    return ';'.join(f'{t:.2f}:{r:.3f}' for t, r in zip(ts, rs))

stats = pd.DataFrame(dict(
    id=ids, a_au=a, e=e, i_deg=inc, Omega_deg=Om, omega_deg=w, M0_deg=M0,
    q_au=q, Q_au=Q, period_d=P_d, period_yr=P_d / YEAR_D, n_rad_s=[f"{x:.10e}" for x in n_rad_s], n_deg_per_day=np.degrees(n_rad_s) * DAY,
    r_asc_node_au=r_asc, r_desc_node_au=r_desc, lon_asc_node_deg=Om, lon_desc_node_deg=(Om + 180) % 360,
    crosses_annulus_0p7_1p5=crosses_annulus, crosses_band_0p8_1p6=crosses_band,
    retrograde=retro,
    r_at_t0_au=r_at_t0, z_at_t0_au=z_at_t0,
    r_min_window_au=r_min_w, t_r_min_mjd=t_rmin, date_r_min=[mjd_to_date(t) for t in t_rmin],
    r_max_window_au=r_max_w, t_r_max_mjd=t_rmax, date_r_max=[mjd_to_date(t) for t in t_rmax],
    zabs_max_window_au=zmax_grid,
    n_perihelion_passages=n_peri, perihelion_dates_mjd=[fmt_list(t) for t in peri_times],
    n_asc_crossings=n_asc, n_desc_crossings=n_desc, n_ecliptic_crossings=n_node_total,
    n_node_crossings_in_band=n_node_band, n_node_crossings_0p6_2p0=n_node_wide,
    asc_crossings_mjd_r=[fmt_list(asc_times[j], np.full(len(asc_times[j]), r_asc[j])) for j in range(N)],
    desc_crossings_mjd_r=[fmt_list(desc_times[j], np.full(len(desc_times[j]), r_desc[j])) for j in range(N)],
    moid_au=moid_au, moid_km=moid_au * AU, moid_r_helio_au=moid_r, moid_lon_deg=moid_lon, moid_z_au=moid_z,
    d_earth_min_window_au=d_earth_min, d_earth_min_window_km=d_earth_min * AU, t_d_earth_min_mjd=t_d_earth_min,
    date_d_earth_min=[mjd_to_date(t) for t in t_d_earth_min], v_rel_earth_at_min_kms=vrel_at_min,
    n_earth_approaches_lt_0p10au=n_ca_010, n_earth_approaches_lt_0p05au=n_ca_005, n_earth_approaches_lt_0p02au=n_ca_002,
    days_in_band=days_in_band, frac_in_band=frac_in_band,
    days_accessible=days_access, frac_accessible=frac_access, n_access_windows=n_access_windows,
    v_rel_earth_min_kms=vrel_grid.min(axis=1), v_rel_earth_max_kms=vrel_grid.max(axis=1),
    difficulty=cls, note=reason,
))
stats.to_csv(os.path.join(HERE, 'population_stats.csv'), index=False, float_format='%.12g')
node_df.to_csv(os.path.join(HERE, 'nodal_crossings.csv'), index=False, float_format='%.12g')
ca_df.to_csv(os.path.join(HERE, 'earth_close_approaches.csv'), index=False, float_format='%.12g')
access_df.to_csv(os.path.join(HERE, 'access_windows.csv'), index=False, float_format='%.12g')

# ---------------- aggregate summary ----------------
years = np.arange(2030, 2045)
yr_edges = MJD_T0 + (years - 2030) * YEAR_D
yr_edges = np.append(yr_edges, MJD_END + 1e-6)
nb = node_df[node_df.in_band]
per_year_band = np.histogram(nb.t_mjd, bins=yr_edges)[0]
per_year_all = np.histogram(node_df.t_mjd, bins=yr_edges)[0]
nb_hi = nb[nb.inc_deg > I_MOD]
per_year_band_hi = np.histogram(nb_hi.t_mjd, bins=yr_edges)[0]
per_year_access_starts = np.histogram(access_df.t_start_mjd, bins=yr_edges)[0] if len(access_df) else np.zeros(15, int)
ca01 = ca_df[ca_df.d_au <= 0.1]
per_year_ca01 = np.histogram(ca01.t_mjd, bins=yr_edges)[0]
per_year_ca005 = np.histogram(ca_df[ca_df.d_au <= 0.05].t_mjd, bins=yr_edges)[0]

def cnt(mask): return int(np.sum(mask))
summary = dict(
    N=N, window=dict(t0_mjd=MJD_T0, tend_mjd=MJD_END, days=T_MAX_D, grid_step_d=DT_GRID),
    thresholds=dict(BAND_LO=BAND_LO, BAND_HI=BAND_HI, ANNULUS_LO=ANNULUS_LO, ANNULUS_HI=ANNULUS_HI, Z_ACCESS=Z_ACCESS,
                    R_UNREACH=R_UNREACH, I_EASY=I_EASY, I_MOD=I_MOD, Q_EASY=Q_EASY, CLOSE_APPROACH_MAX=CLOSE_APPROACH_MAX),
    classes={c: cnt(cls == c) for c in ['easy', 'moderate', 'hard', 'very_hard', 'unreachable']},
    class_ids={c: [int(x) for x in ids[cls == c]] for c in ['very_hard', 'unreachable']},
    elements=dict(a=dict(min=float(a.min()), median=float(np.median(a)), max=float(a.max())),
                  e=dict(min=float(e.min()), median=float(np.median(e)), max=float(e.max())),
                  i=dict(min=float(inc.min()), median=float(np.median(inc)), max=float(inc.max())),
                  q=dict(min=float(q.min()), median=float(np.median(q)), max=float(q.max())),
                  Q=dict(min=float(Q.min()), median=float(np.median(Q)), max=float(Q.max())),
                  period_yr=dict(min=float((P_d / YEAR_D).min()), median=float(np.median(P_d / YEAR_D)), max=float((P_d / YEAR_D).max()))),
    counts=dict(
        q_lt_1=cnt(q < 1), q_lt_1p3=cnt(q < 1.3), Q_gt_1=cnt(Q > 1), a_lt_1=cnt(a < 1),
        atens=cnt((a < 1) & (Q > 0.983)), apollos=cnt((a >= 1) & (q < 1.017)), amors=cnt((a > 1) & (q > 1.017) & (q < 1.3)),
        atiras=cnt(Q < 0.983), other=cnt(q >= 1.3),
        i_le_5=cnt(inc <= 5), i_le_10=cnt(inc <= 10), i_le_15=cnt(inc <= 15), i_le_20=cnt(inc <= 20), i_le_30=cnt(inc <= 30),
        i_gt_30=cnt(inc > 30), i_gt_45=cnt(inc > 45), retrograde=cnt(retro), a_gt_3=cnt(a > 3), a_gt_5=cnt(a > 5),
        crosses_annulus_0p7_1p5=cnt(crosses_annulus), crosses_band_0p8_1p6=cnt(crosses_band),
        ever_in_band_window=cnt(days_in_band > 0), ever_accessible_window=cnt(days_access > 0),
        r_min_window_gt_2p5=cnt(r_min_w > 2.5), r_min_window_gt_1p6=cnt(r_min_w > 1.6),
        n_peri_0=cnt(n_peri == 0), n_peri_ge_5=cnt(n_peri >= 5), n_peri_ge_10=cnt(n_peri >= 10),
        moid_lt_0p05=cnt(moid_au < 0.05), moid_lt_0p02=cnt(moid_au < 0.02), moid_lt_0p01=cnt(moid_au < 0.01),
        d_earth_min_lt_0p1=cnt(d_earth_min < 0.1), d_earth_min_lt_0p05=cnt(d_earth_min < 0.05), d_earth_min_lt_0p02=cnt(d_earth_min < 0.02),
        d_earth_min_lt_0p01=cnt(d_earth_min < 0.01),
        node_band_ge_1=cnt(n_node_band >= 1), node_band_0=cnt(n_node_band == 0),
    ),
    totals=dict(
        perihelion_passages=int(n_peri.sum()), ecliptic_crossings=int(n_node_total.sum()),
        ecliptic_crossings_in_band=int(n_node_band.sum()), ecliptic_crossings_in_band_per_year=float(n_node_band.sum() / 15.0),
        ecliptic_crossings_in_band_highinc=int(len(nb_hi)),
        ecliptic_crossings_0p6_2p0=int(n_node_wide.sum()),
        access_windows=int(len(access_df)), access_window_median_days=float(access_df.duration_d.median()) if len(access_df) else 0.0,
        earth_approaches_lt_0p1=int(len(ca01)), earth_approaches_lt_0p05=int((ca_df.d_au <= 0.05).sum()),
        earth_approaches_lt_0p02=int((ca_df.d_au <= 0.02).sum()),
        asteroid_days_in_band=float(days_in_band.sum()), asteroid_days_accessible=float(days_access.sum()),
        mean_asteroids_in_band_at_any_time=float(in_band.sum(axis=0).mean()),
        mean_asteroids_accessible_at_any_time=float(access.sum(axis=0).mean()),
    ),
    per_year=dict(years=[int(y) for y in years], node_crossings_in_band=[int(x) for x in per_year_band],
                  node_crossings_in_band_inc_gt_20=[int(x) for x in per_year_band_hi],
                  node_crossings_all=[int(x) for x in per_year_all],
                  access_window_starts=[int(x) for x in per_year_access_starts],
                  earth_approaches_lt_0p1=[int(x) for x in per_year_ca01], earth_approaches_lt_0p05=[int(x) for x in per_year_ca005]),
    near_1au_nodes=dict(  # crossings with 0.9 <= r <= 1.1 AU
        total=int(((node_df.r_au >= 0.9) & (node_df.r_au <= 1.1)).sum()),
        per_year=float(((node_df.r_au >= 0.9) & (node_df.r_au <= 1.1)).sum() / 15.0),
        asteroids=int(node_df[(node_df.r_au >= 0.9) & (node_df.r_au <= 1.1)].id.nunique()),
        highinc_total=int(((node_df.r_au >= 0.9) & (node_df.r_au <= 1.1) & (node_df.inc_deg > I_MOD)).sum()),
        highinc_asteroids=int(node_df[(node_df.r_au >= 0.9) & (node_df.r_au <= 1.1) & (node_df.inc_deg > I_MOD)].id.nunique())),
    node_r_split=dict(r_lt_0p8=int((node_df.r_au < 0.8).sum()), in_band=int(node_df.in_band.sum()), r_gt_1p6=int((node_df.r_au > 1.6).sum())),
    in_band_node_vrel_earth_kms=dict(q10=float(nb.v_rel_earth_kms.quantile(.1)), median=float(nb.v_rel_earth_kms.median()), q90=float(nb.v_rel_earth_kms.quantile(.9))),
    in_band_node_abs_vz_kms=dict(q25=float(nb.vz_kms.abs().quantile(.25)), median=float(nb.vz_kms.abs().median()), q75=float(nb.vz_kms.abs().quantile(.75))),
    access_windows_by_class={c: dict(n=int((access_df.id.isin(ids[cls == c])).sum()),
                                     median_days=float(access_df[access_df.id.isin(ids[cls == c])].duration_d.median()) if (access_df.id.isin(ids[cls == c])).sum() else 0.0)
                             for c in ['easy', 'moderate', 'hard']},
    n_access_windows_quantiles={str(qq): float(np.quantile(n_access_windows, qq)) for qq in [0, .1, .25, .5, .75, .9, 1]},
    asteroids_with_ge3_windows=cnt(n_access_windows >= 3), asteroids_with_ge5_windows=cnt(n_access_windows >= 5),
    asteroids_with_earth_approach_lt_0p1=int(ca01.id.nunique()), asteroids_with_earth_approach_lt_0p05=int(ca_df[ca_df.d_au <= 0.05].id.nunique()),
    earth_approach_lt_0p1_vrel_kms=dict(min=float(ca01.v_rel_kms.min()), median=float(ca01.v_rel_kms.median()), max=float(ca01.v_rel_kms.max()),
                                        n_lt_10kms=int((ca01.v_rel_kms < 10).sum())),
    class_stats={c: dict(n=cnt(cls == c), i_median=float(np.median(inc[cls == c])) if cnt(cls == c) else None,
                         q_median=float(np.median(q[cls == c])) if cnt(cls == c) else None,
                         node_band_total=int(n_node_band[cls == c].sum()), days_access_median=float(np.median(days_access[cls == c])) if cnt(cls == c) else None)
                 for c in ['easy', 'moderate', 'hard', 'very_hard', 'unreachable']},
    special={int(k): dict(a=float(a[k - 1]), e=float(e[k - 1]), i=float(inc[k - 1]), q=float(q[k - 1]), Q=float(Q[k - 1]),
                          period_yr=float(P_d[k - 1] / YEAR_D), r_t0=float(r_at_t0[k - 1]), r_min_window=float(r_min_w[k - 1]),
                          date_r_min=mjd_to_date(t_rmin[k - 1]), r_max_window=float(r_max_w[k - 1]),
                          M_at_t0_deg=float(np.degrees(M_t0[k - 1]) % 360), M_at_tend_deg=float(np.degrees(M_tf[k - 1]) % 360),
                          last_peri_mjd=float(MJD_EPH_AST + (twopi * np.floor(M_t0[k - 1] / twopi) - np.radians(M0[k - 1])) / n_rad_s[k - 1] / DAY),
                          next_peri_mjd=float(MJD_EPH_AST + (twopi * np.ceil(M_t0[k - 1] / twopi) - np.radians(M0[k - 1])) / n_rad_s[k - 1] / DAY),
                          d_earth_min=float(d_earth_min[k - 1]), cls=str(cls[k - 1]))
             for k in [131, 144] + [int(x) for x in ids[(cls == 'unreachable') | (cls == 'very_hard')] if x not in (131, 144)]},
)
for k in summary['special']:
    summary['special'][k]['last_peri_date'] = mjd_to_date(summary['special'][k]['last_peri_mjd'])
    summary['special'][k]['next_peri_date'] = mjd_to_date(summary['special'][k]['next_peri_mjd'])
with open(os.path.join(HERE, 'summary.json'), 'w') as f:
    json.dump(summary, f, indent=1)
print(json.dumps(summary, indent=1)); print("elapsed_s", round(__import__("time").time() - T_START, 1))

# ---------------- plots ----------------
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
colors = {'easy': '#2a9d8f', 'moderate': '#457b9d', 'hard': '#e9a23b', 'very_hard': '#d1495b', 'unreachable': '#222222'}
order = ['easy', 'moderate', 'hard', 'very_hard', 'unreachable']

fig, ax = plt.subplots(1, 2, figsize=(13, 5.5))
for c in order:
    m = cls == c
    ax[0].scatter(a[m], e[m], s=18, c=colors[c], label=f'{c} ({m.sum()})', alpha=0.85, edgecolor='none')
    ax[1].scatter(a[m], inc[m], s=18, c=colors[c], label=f'{c} ({m.sum()})', alpha=0.85, edgecolor='none')
aa = np.linspace(0.3, 20, 500)
ax[0].plot(aa, np.clip(1 - 1.0 / aa, 0, 1), 'k--', lw=0.8, label='q = 1 AU')
ax[0].plot(aa, np.clip(1.0 / aa - 1, 0, 1), 'k:', lw=0.8, label='Q = 1 AU')
ax[0].plot(aa, np.clip(1 - 1.3 / aa, 0, 1), color='gray', ls='--', lw=0.8, label='q = 1.3 AU')
ax[0].set_xscale('log'); ax[0].set_xlabel('a [AU]'); ax[0].set_ylabel('e'); ax[0].set_ylim(0, 1); ax[0].set_title('a - e (300 targets)')
ax[0].legend(fontsize=8, loc='lower right'); ax[0].grid(alpha=0.3)
for j in [130, 143]: ax[0].annotate(str(ids[j]), (a[j], e[j]), fontsize=8, xytext=(4, 4), textcoords='offset points'); ax[1].annotate(str(ids[j]), (a[j], inc[j]), fontsize=8, xytext=(4, 4), textcoords='offset points')
ax[1].set_xscale('log'); ax[1].set_xlabel('a [AU]'); ax[1].set_ylabel('i [deg]'); ax[1].set_title('a - i'); ax[1].grid(alpha=0.3)
ax[1].axhline(I_EASY, color='gray', lw=0.7, ls='--'); ax[1].axhline(I_MOD, color='gray', lw=0.7, ls='--'); ax[1].axhline(90, color='k', lw=0.7, ls=':')
ax[1].legend(fontsize=8, loc='upper left')
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'scatter_a_e_a_i.png'), dpi=140); plt.close(fig)

fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
ax[0].hist(inc, bins=np.arange(0, 140, 5), color='#457b9d', edgecolor='white'); ax[0].set_xlabel('inclination [deg]'); ax[0].set_ylabel('count'); ax[0].set_title(f'i: median {np.median(inc):.1f} deg, {cnt(inc<=10)} <=10, {cnt(inc>30)} >30')
ax[1].hist(q, bins=np.arange(0, 1.4, 0.05), color='#2a9d8f', edgecolor='white'); ax[1].set_xlabel('perihelion q [AU]'); ax[1].set_title(f'q: median {np.median(q):.2f} AU, {cnt(q<1)} < 1 AU')
ax[2].hist(np.clip(Q, 0, 6), bins=np.arange(0.5, 6.1, 0.25), color='#e9a23b', edgecolor='white'); ax[2].set_xlabel('aphelion Q [AU] (clipped at 6)'); ax[2].set_title(f'Q: median {np.median(Q):.2f} AU')
for x in ax: x.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'hist_i_q_Q.png'), dpi=140); plt.close(fig)

# timeline of nodal crossings in band
fig, ax = plt.subplots(2, 1, figsize=(14, 8), gridspec_kw={'height_ratios': [3, 1.2]}, sharex=True)
nbp = nb.copy(); nbp['year'] = mjd_to_decimal_year(nbp.t_mjd.values)
sc = ax[0].scatter(nbp.year, nbp.r_au, c=nbp.inc_deg, cmap='viridis', s=14, vmin=0, vmax=60, alpha=0.85)
cb = fig.colorbar(sc, ax=ax[0], pad=0.01); cb.set_label('inclination [deg]')
ax[0].axhline(1.0, color='k', lw=0.6, ls='--'); ax[0].set_ylabel('heliocentric r at ecliptic crossing [AU]')
ax[0].set_title(f'Ecliptic-plane crossings with {BAND_LO} <= r <= {BAND_HI} AU during the mission window: {len(nb)} events ({len(nb)/15:.1f}/yr), {nb.id.nunique()} asteroids')
ax[0].grid(alpha=0.3)
ax[1].bar(years + 0.5, per_year_band, width=0.9, color='#457b9d', label='all i')
ax[1].bar(years + 0.5, per_year_band_hi, width=0.9, color='#e9a23b', label=f'i > {I_MOD:.0f} deg')
ax[1].set_ylabel('crossings / yr'); ax[1].set_xlabel('year'); ax[1].legend(); ax[1].grid(alpha=0.3, axis='y')
ax[1].set_xlim(2030, 2045)
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'timeline_nodal_crossings.png'), dpi=140); plt.close(fig)

# timeline of Earth close approaches + accessibility count
fig, ax = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
cay = mjd_to_decimal_year(ca_df.t_mjd.values)
ax[0].scatter(cay, ca_df.d_au, s=12, c='#d1495b', alpha=0.8)
for _, row in ca_df[ca_df.d_au < 0.02].iterrows():
    ax[0].annotate(str(row.id), (mjd_to_decimal_year(row.t_mjd), row.d_au), fontsize=6, xytext=(2, 2), textcoords='offset points')
ax[0].set_yscale('log'); ax[0].set_ylabel('closest approach to Earth [AU]'); ax[0].grid(alpha=0.3, which='both')
ax[0].set_title(f'Earth close approaches < {CLOSE_APPROACH_MAX} AU in window: {len(ca_df)} events, {len(ca01)} below 0.1 AU, {(ca_df.d_au<=0.05).sum()} below 0.05 AU')
ty = mjd_to_decimal_year(t_grid)
ax[1].plot(ty, in_band.sum(axis=0), color='#457b9d', lw=0.8, label=f'asteroids with {BAND_LO}<=r<={BAND_HI} AU')
ax[1].plot(ty, access.sum(axis=0), color='#2a9d8f', lw=0.8, label=f'... and |z| <= {Z_ACCESS} AU')
ax[1].set_ylabel('number of asteroids'); ax[1].set_xlabel('year'); ax[1].legend(); ax[1].grid(alpha=0.3); ax[1].set_xlim(2030, 2045)
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'timeline_earth_approaches_and_band_count.png'), dpi=140); plt.close(fig)

# polar map of nodal crossing positions in band (where in the ecliptic they happen)
fig = plt.figure(figsize=(7, 7)); axp = fig.add_subplot(111, projection='polar')
axp.scatter(np.radians(nb.lon_deg), nb.r_au, c=nb.inc_deg, cmap='viridis', s=12, vmin=0, vmax=60, alpha=0.8)
axp.plot(np.linspace(0, twopi, 200), np.ones(200), 'k--', lw=0.6)
axp.set_rmax(BAND_HI + 0.05); axp.set_title('Ecliptic longitude and r of in-band ecliptic crossings (window)')
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'polar_nodal_crossings.png'), dpi=140); plt.close(fig)
print('done')
