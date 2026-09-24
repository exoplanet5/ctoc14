"""
Task 1: Earth -> asteroid direct ballistic intercept scan.

Grid: launch date tL in [0, 15 yr) step DT_L days  x  time of flight in [30, 1500] d step DT_TOF.
For each asteroid and each cell: Lambert arc(s) from r_Earth(tL) to r_ast(tL+tof); launch
v_inf = |v1 - v_Earth(tL)|.  Arc types solved per cell:
    N=0 prograde, N=0 retrograde, and N = 1..NREV_MAX revolution arcs (both branches, prograde)
    for tof >= 250*N d.  The cell value is the minimum v_inf over all arc types.
Arrival must be inside the mission window.  A flyby only needs position matching, so the
arrival relative velocity is recorded but not constrained.

Outputs (analysis/ballistic/out/):
  earth_to_asteroid_summary.csv    one row per asteroid
  earth_to_asteroid_vinf_grid.npz  vinf[asteroid, tL, tof] (float16, km/s; min over arc types),
                                   vinf_n0[...] (single-rev only), nrev_best[...] (int8, -1 = none),
                                   vinf_min_tof[asteroid, tL] (float32)
  fig_min_vinf_hist.png, fig_launch_calendar.png
"""
import os, sys, time
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ephem import (MU_SUN, AU, DAY, T_MAX_DAYS, earth_state, asteroid_state, load_mea, elements_from_state)
from lambert import lambert, lambert_multirev

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
DT_L = 5.0          # launch-date step [d]
DT_TOF = 10.0       # TOF step [d]
TOF_MIN, TOF_MAX = 30.0, 1500.0
VINF_MAX = 4.0
NREV_MAX = 5
TOF_MIN_REV = {n: 250.0 * n for n in range(1, NREV_MAX + 1)}   # N-rev arcs only solved for tof >= this

TL = np.arange(0.0, T_MAX_DAYS, DT_L)
TOF = np.arange(TOF_MIN, TOF_MAX + 0.5 * DT_TOF, DT_TOF)

_IDS, _EL = load_mea()
_RE, _VE = earth_state(TL)                       # (nL,3)


def vinf_grids(el, TL_, TOF_, rE, vE):
    """Return dict of v_inf grids (nL, nT) per arc type, plus arrival-relative-speed and v1 for the best type."""
    t_arr = TL_[:, None] + TOF_[None, :]
    valid = t_arr <= T_MAX_DAYS
    rA, vA = asteroid_state(el, t_arr)                    # (nL, nT, 3)
    r1 = np.broadcast_to(rE[:, None, :], rA.shape)
    tof_s = np.broadcast_to(TOF_[None, :] * DAY, t_arr.shape)
    grids, v1s, vrels = {}, {}, {}
    for retro in (False, True):
        v1, v2 = lambert(r1, rA, tof_s, MU_SUN, retrograde=retro)
        vinf = np.linalg.norm(v1 - vE[:, None, :], axis=-1)
        key = 'N0retro' if retro else 'N0'
        grids[key] = np.where(valid & np.isfinite(vinf), vinf, np.inf)
        v1s[key] = v1
        vrels[key] = np.linalg.norm(v2 - vA, axis=-1)
    for n in range(1, NREV_MAX + 1):
        cols = TOF_ >= TOF_MIN_REV[n]
        if not cols.any():
            continue
        v1m, v2m = lambert_multirev(r1[:, cols], rA[:, cols], tof_s[:, cols], MU_SUN, N=n)   # (nL, nc, 2, 3)
        vinf_m = np.linalg.norm(v1m - vE[:, None, None, :], axis=-1)                         # (nL, nc, 2)
        vinf_m = np.where(np.isfinite(vinf_m), vinf_m, np.inf)
        ib = np.argmin(vinf_m, axis=-1)
        g = np.full(t_arr.shape, np.inf)
        g[:, cols] = np.where(valid[:, cols], np.take_along_axis(vinf_m, ib[..., None], -1)[..., 0], np.inf)
        v1full = np.full(rA.shape, np.nan)
        v1full[:, cols] = np.take_along_axis(v1m, ib[..., None, None], -2)[..., 0, :]
        vr = np.full(t_arr.shape, np.nan)
        vr[:, cols] = np.linalg.norm(np.take_along_axis(v2m, ib[..., None, None], -2)[..., 0, :] - vA[:, cols], axis=-1)
        grids[f'N{n}'] = g
        v1s[f'N{n}'] = v1full
        vrels[f'N{n}'] = vr
    return grids, v1s, vrels, valid


def scan_one(j):
    el = _EL[j]
    grids, v1s, vrels, valid = vinf_grids(el, TL, TOF, _RE, _VE)
    keys = list(grids.keys())
    stack = np.stack([grids[k] for k in keys])            # (nK, nL, nT)
    ibest = np.argmin(stack, axis=0)
    vinf = np.min(stack, axis=0)
    vinf_n0 = np.minimum(grids['N0'], grids['N0retro'])
    nrev_of_key = np.array([0 if k.startswith('N0') else int(k[1:]) for k in keys])
    nrev_best = np.where(np.isfinite(vinf), nrev_of_key[ibest], -1).astype(np.int8)
    vrel = np.choose(ibest, [vrels[k] for k in keys])
    k = np.unravel_index(np.argmin(vinf), vinf.shape)
    k0 = np.unravel_index(np.argmin(vinf_n0), vinf_n0.shape)
    feas = vinf <= VINF_MAX
    feas0 = vinf_n0 <= VINF_MAX
    tof_short = TOF <= 365.0
    res = dict(
        id=int(_IDS[j]), a=el[0], e=el[1], inc=el[2], q=el[0] * (1 - el[1]), Q=el[0] * (1 + el[1]),
        min_vinf=float(vinf[k]), best_nrev=int(nrev_best[k]),
        min_vinf_n0=float(vinf_n0[k0]), min_vinf_n0_pro=float(grids['N0'].min()), min_vinf_n0_retro=float(grids['N0retro'].min()),
        min_vinf_multirev=float(min(grids[k_].min() for k_ in keys if not k_.startswith('N0'))),
        best_tL_day=float(TL[k[0]]), best_tof_day=float(TOF[k[1]]), best_arr_day=float(TL[k[0]] + TOF[k[1]]),
        best_flyby_vrel=float(vrel[k]),
        best_n0_tL_day=float(TL[k0[0]]), best_n0_tof_day=float(TOF[k0[1]]),
        n_cells_le4=int(feas.sum()), n_cells_le3=int((vinf <= 3).sum()), n_cells_le2=int((vinf <= 2).sum()),
        n_cells_le1=int((vinf <= 1).sum()),
        n_cells_le4_n0=int(feas0.sum()), n_cells_le4_multirev=int((feas & (nrev_best > 0)).sum()),
        n_launch_dates_le4=int(feas.any(axis=1).sum()), n_launch_dates_le4_n0=int(feas0.any(axis=1).sum()),
        first_launch_le4=float(TL[feas.any(axis=1)].min()) if feas.any() else np.nan,
        last_launch_le4=float(TL[feas.any(axis=1)].max()) if feas.any() else np.nan,
        min_vinf_tof_le365=float(vinf[:, tof_short].min()),
        n_cells_le4_tof_le365=int(feas[:, tof_short].sum()),
        n_cells_retro_le4=int((grids['N0retro'] <= VINF_MAX).sum()),
        median_flyby_vrel_le4=float(np.median(vrel[feas])) if feas.any() else np.nan,
        n_valid_cells=int(valid.sum()),
    )
    v1b = v1s[keys[ibest[k]]][k]
    a_, e_, i_, q_, Q_ = elements_from_state(_RE[k[0]][None], v1b[None])
    res.update(best_arc_a=float(a_[0] / AU), best_arc_e=float(e_[0]), best_arc_inc=float(np.degrees(i_[0])))
    return j, res, vinf.astype(np.float16), vinf_n0.astype(np.float16), nrev_best, vinf.min(axis=1).astype(np.float32)


def main():
    import multiprocessing as mp
    import pandas as pd
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    nproc = int(os.environ.get('NPROC', 8))
    ctx = mp.get_context('fork')
    rows = [None] * len(_IDS)
    grid = np.zeros((len(_IDS), len(TL), len(TOF)), dtype=np.float16)
    grid0 = np.zeros_like(grid)
    nrev = np.zeros(grid.shape, dtype=np.int8)
    gmin = np.zeros((len(_IDS), len(TL)), dtype=np.float32)
    with ctx.Pool(nproc) as pool:
        for j, res, g, g0, nr, gm in pool.imap_unordered(scan_one, range(len(_IDS)), chunksize=2):
            rows[j] = res
            grid[j] = g
            grid0[j] = g0
            nrev[j] = nr
            gmin[j] = gm
            if j % 50 == 0:
                print(f"  asteroid {j} done, elapsed {time.time()-t0:.1f}s", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, 'earth_to_asteroid_summary.csv'), index=False, float_format='%.6g')
    np.savez_compressed(os.path.join(OUT, 'earth_to_asteroid_vinf_grid.npz'), vinf=grid, vinf_n0=grid0, nrev_best=nrev,
                        vinf_min_tof=gmin, tL_days=TL, tof_days=TOF, ids=_IDS)
    print(f"scan done in {time.time()-t0:.1f}s; grid {len(TL)} launch dates x {len(TOF)} TOFs = {len(TL)*len(TOF)} cells/asteroid; "
          f"arc types: N=0 pro/retro + N=1..{NREV_MAX} (both branches) for tof >= 250*N d")
    report(df, grid, nrev)
    plots(df, gmin)


def report(df, grid, nrev):
    print("\n=== Earth -> asteroid ballistic reachability (Lambert arcs, launch v_inf) ===")
    print(f"grid: launch step {DT_L} d ({len(TL)} dates), TOF {TOF_MIN}-{TOF_MAX} d step {DT_TOF} d ({len(TOF)} values)")
    print("  threshold | asteroids with min v_inf <= thr: all arcs | single-rev only")
    for thr in (0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10):
        print(f"  {thr:5.1f} km/s | {(df.min_vinf <= thr).sum():3d} | {(df.min_vinf_n0 <= thr).sum():3d}")
    print(f"  asteroids reachable (v_inf<=4) with TOF<=365 d: {(df.min_vinf_tof_le365 <= 4).sum()}")
    print(f"  asteroids whose best cell is a multi-rev arc: {(df.best_nrev > 0).sum()}; best_nrev histogram: "
          f"{dict(zip(*np.unique(df.best_nrev, return_counts=True)))}")
    print(f"  retrograde branch ever better than prograde at the optimum: {(df.min_vinf_n0_retro < df.min_vinf_n0_pro).sum()} asteroids; retro cells <=4 km/s total: {df.n_cells_retro_le4.sum()}")
    print("  min v_inf percentiles (5,25,50,75,95) all arcs:", np.percentile(df.min_vinf, [5, 25, 50, 75, 95]).round(2),
          " single-rev:", np.percentile(df.min_vinf_n0, [5, 25, 50, 75, 95]).round(2))
    r = df[df.min_vinf <= 4]
    tot_cells = int((grid <= 4).sum())
    tot_cells_multi = int(((grid <= 4) & (nrev > 0)).sum())
    print(f"  feasible cells (v_inf<=4) summed over all asteroids: {tot_cells} of {300*len(TL)*len(TOF)} "
          f"({tot_cells_multi} = {100*tot_cells_multi/max(tot_cells,1):.1f}% won by multi-rev arcs); by N: "
          f"{ {int(n): int(((grid <= 4) & (nrev == n)).sum()) for n in range(0, NREV_MAX+1)} }")
    print(f"  among reachable: feasible cells per asteroid median {r.n_cells_le4.median():.0f}, mean {r.n_cells_le4.mean():.0f}, "
          f"min {r.n_cells_le4.min()}, max {r.n_cells_le4.max()}; launch dates with a feasible TOF: median {r.n_launch_dates_le4.median():.0f} of {len(TL)} "
          f"(single-rev only: median cells {r.n_cells_le4_n0.median():.0f}, launch dates {r.n_launch_dates_le4_n0.median():.0f})")
    print(f"  flyby relative speed at best cell among reachable: median {r.best_flyby_vrel.median():.1f} km/s, "
          f"min {r.best_flyby_vrel.min():.1f}, max {r.best_flyby_vrel.max():.1f}")
    print(f"  best TOF among reachable: median {r.best_tof_day.median():.0f} d; 25/75%: {np.percentile(r.best_tof_day,[25,75])}")
    print("  unreachable asteroids (min v_inf > 4): ids", df[df.min_vinf > 4].id.tolist(),
          " with min v_inf", df[df.min_vinf > 4].min_vinf.round(2).tolist())
    print("  min v_inf vs orbit: corr with i:", np.corrcoef(df.min_vinf, df.inc)[0, 1].round(2),
          " with q:", np.corrcoef(df.min_vinf, df.q)[0, 1].round(2), " with e:", np.corrcoef(df.min_vinf, df.e)[0, 1].round(2))
    # per-year availability: how many asteroids have some feasible launch in each calendar year
    yrs = (TL / 365.25).astype(int)
    anyf = (grid <= 4).any(axis=2)                                    # (300, nL)
    per_year = [int(anyf[:, yrs == y].any(axis=1).sum()) for y in range(15)]
    print("  asteroids with >=1 feasible (v_inf<=4) launch date in year 2030..2044:", per_year)
    # how many asteroids could be *simultaneously* started from the same launch date (any TOF)
    n_same = anyf.sum(axis=0)
    print(f"  number of asteroids ballistically reachable from a given launch date: median {np.median(n_same):.0f}, min {n_same.min()}, max {n_same.max()}")


def plots(df, gmin):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    v = np.clip(df.min_vinf.values, 0, 20)
    ax[0].hist(v, bins=np.arange(0, 20.5, 0.5), color='#4C72B0', label='all arc types (N=0..5)')
    ax[0].hist(np.clip(df.min_vinf_n0.values, 0, 20), bins=np.arange(0, 20.5, 0.5), histtype='step', color='#DD8452', label='single-rev only')
    ax[0].axvline(4, color='k', ls='--', label='launch limit 4 km/s')
    ax[0].set_xlabel('min launch v_inf over grid [km/s] (clipped at 20)')
    ax[0].set_ylabel('asteroids')
    ax[0].legend(fontsize=8)
    ax[0].set_title(f'{(df.min_vinf<=4).sum()} / 300 reachable ballistically')
    vs = np.sort(df.min_vinf.values)
    ax[1].plot(vs, np.arange(1, 301), color='#4C72B0')
    ax[1].axvline(4, color='k', ls='--')
    ax[1].set_xlim(0, 15)
    ax[1].set_xlabel('min launch v_inf [km/s]')
    ax[1].set_ylabel('cumulative number of asteroids')
    ax[1].grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'fig_min_vinf_hist.png'), dpi=130)

    order = np.argsort(df.min_vinf.values)
    fig, ax = plt.subplots(figsize=(12, 7))
    img = np.where(gmin[order] <= 4.0, gmin[order], np.nan)
    im = ax.imshow(img, aspect='auto', origin='lower', extent=[0, TL[-1] / 365.25, 0, 300], cmap='viridis', vmin=0, vmax=4)
    ax.set_xlabel('launch date [years after 2030-01-01]')
    ax.set_ylabel('asteroid rank (sorted by min v_inf)')
    ax.set_title('Ballistic launch windows: min over TOF (30-1500 d) of launch v_inf (only cells <= 4 km/s shown)')
    plt.colorbar(im, ax=ax, label='v_inf [km/s]')
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'fig_launch_calendar.png'), dpi=130)
    plt.close('all')


if __name__ == '__main__':
    main()
