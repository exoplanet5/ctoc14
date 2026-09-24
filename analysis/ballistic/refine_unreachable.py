"""Finer Earth->asteroid grid for the asteroids that were not reachable with v_inf <= 4 km/s on the coarse
grid (launch step 2 d, TOF 20..5470 d step 4 d, arcs N=0 pro/retro and N=1..8 both branches), plus a local
1 d x 1 d refinement of the coarse optimum for all asteroids (all arc types)."""
import os, sys, time
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ephem import MU_SUN, AU, DAY, T_MAX_DAYS, earth_state, asteroid_state, load_mea
from lambert import lambert, lambert_multirev
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
ids, el = load_mea()
df = pd.read_csv(os.path.join(OUT, 'earth_to_asteroid_summary.csv'))
NREV_MAX_FINE = 8


def best_vinf(elj, TL, TOF, nrev_max):
    """min over the grid and arc types of launch v_inf; returns (vmin, tL, tof, nrev, vrel, n_cells_le4)"""
    rE, vE = earth_state(TL)
    t_arr = TL[:, None] + TOF[None, :]
    valid = t_arr <= T_MAX_DAYS
    rA, vA = asteroid_state(elj, t_arr)
    r1 = np.broadcast_to(rE[:, None, :], rA.shape)
    tofs = np.broadcast_to(TOF[None, :] * DAY, t_arr.shape)
    best = (np.inf, None, None, None, None)
    n4 = 0
    for retro in (False, True):
        v1, v2 = lambert(r1, rA, tofs, MU_SUN, retrograde=retro)
        vinf = np.linalg.norm(v1 - vE[:, None, :], axis=-1)
        vinf = np.where(valid & np.isfinite(vinf), vinf, np.inf)
        n4 += int((vinf <= 4).sum())
        k = np.unravel_index(np.argmin(vinf), vinf.shape)
        if vinf[k] < best[0]:
            best = (vinf[k], TL[k[0]], TOF[k[1]], 0, np.linalg.norm(v2[k] - vA[k]))
    for n in range(1, nrev_max + 1):
        cols = TOF >= 250.0 * n
        if not cols.any():
            continue
        v1m, v2m = lambert_multirev(r1[:, cols], rA[:, cols], tofs[:, cols], MU_SUN, N=n)
        vinf = np.linalg.norm(v1m - vE[:, None, None, :], axis=-1)
        vinf = np.where(valid[:, cols][..., None] & np.isfinite(vinf), vinf, np.inf)
        n4 += int((vinf <= 4).sum())
        k = np.unravel_index(np.argmin(vinf), vinf.shape)
        if vinf[k] < best[0]:
            best = (vinf[k], TL[k[0]], TOF[cols][k[1]], n, np.linalg.norm(v2m[k] - vA[:, cols][k[:2]]))
    return best + (n4,)


def refine_bad(aid):
    j = int(np.where(ids == aid)[0][0])
    t0 = time.time()
    TL = np.arange(0.0, T_MAX_DAYS, 2.0)
    best = (np.inf,)
    n4 = 0
    for tof_lo in range(20, 5480, 400):
        TOF = np.arange(tof_lo, min(tof_lo + 400, 5480), 4.0)
        if (TL[0] + TOF[0]) > T_MAX_DAYS:
            continue
        b = best_vinf(el[j], TL, TOF, NREV_MAX_FINE)
        n4 += b[-1]
        if b[0] < best[0]:
            best = b
    a, e, i = el[j][:3]
    t = np.arange(0, T_MAX_DAYS, 1.0)
    r, _ = asteroid_state(el[j], t)
    R = np.linalg.norm(r, axis=-1) / AU
    cross = np.where(np.diff(np.sign(R - 1.05)) != 0)[0]
    msg = (f"id {aid}: a={a:.3f} e={e:.3f} i={i:.2f} q={a*(1-e):.3f} Q={a*(1+e):.2f} P={a**1.5:.2f} yr | fine grid (2 d x 4 d, TOF 20-5476 d, N=0..{NREV_MAX_FINE}): "
           f"min v_inf = {best[0]:.3f} km/s at tL={best[1]:.0f} d, tof={best[2]:.0f} d, N={best[3]}, flyby vrel {best[4]:.1f} km/s; cells <= 4 km/s: {n4}\n"
           f"     |r| range in window {R.min():.3f}-{R.max():.2f} AU; crossings of 1.05 AU at days {t[cross].astype(int).tolist()}  [{time.time()-t0:.0f} s]")
    return msg


def refine_local(row):
    j = int(np.where(ids == row.id)[0][0])
    tl = np.arange(max(row.best_tL_day - 10, 0), min(row.best_tL_day + 10, T_MAX_DAYS) + 0.5, 1.0)
    tf = np.arange(max(row.best_tof_day - 20, 5), row.best_tof_day + 20 + 0.5, 1.0)
    b = best_vinf(el[j], tl, tf, 5)
    return b[0]


if __name__ == '__main__':
    import multiprocessing as mp
    bad = df[df.min_vinf > 4].id.values
    print("refining asteroids", bad, flush=True)
    ctx = mp.get_context('fork')
    with ctx.Pool(min(8, len(bad))) as pool:
        for msg in pool.imap_unordered(refine_bad, bad):
            print(msg, flush=True)
    print("\nlocal 1d x 1d refinement around the coarse optimum (all asteroids, all arc types):", flush=True)
    with ctx.Pool(8) as pool:
        ref = np.array(pool.map(refine_local, [row for _, row in df.iterrows()], chunksize=10))
    print(f"  coarse min v_inf vs refined: mean decrease {np.mean(df.min_vinf.values - ref):.3f} km/s, max decrease {np.max(df.min_vinf.values - ref):.3f} km/s; "
          f"reachable after refinement: {(ref <= 4).sum()}")
    df['min_vinf_refined'] = ref
    df.to_csv(os.path.join(OUT, 'earth_to_asteroid_summary.csv'), index=False, float_format='%.6g')
