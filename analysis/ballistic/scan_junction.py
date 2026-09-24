"""
Task 2: asteroid -> asteroid Lambert legs and "flyby-junction delta-v" statistics.

For every asteroid B and every grid epoch t_B (step DT_EP days), consider
  incoming arcs  A -> B : Lambert(r_A(t_B - tof_in), r_B(t_B), tof_in)  -> arc velocity at B, v_in
  outgoing arcs  B -> C : Lambert(r_B(t_B), r_C(t_B + tof_out), tof_out) -> arc velocity at B, v_out
for all A, C != B (all 300 catalogue asteroids) and tof_in, tof_out in TOF (60..600 d step 30 d plus 365 d).
Arc types: single-revolution prograde (all TOF) and, for TOF >= 250 d, the two 1-revolution
prograde solutions (a near-1-AU spacecraft doing 1.x revolutions).
A chain A -> B -> C with ballistic legs needs an impulsive velocity change at B
    dv_B = |v_out - v_in|.
This scan accumulates, per tof_out and per arc class:
  * the unconditional distribution of dv_B over (incoming arc, C) pairs (best arc type per C),
  * the "best-next-hop" distribution min_C dv_B for each incoming arc,
  * the number of C with dv_B below several thresholds (branching factor),
  * greedy "next flyby as soon as possible" statistics: for each incoming arc, the shortest tof_out
    for which min_C dv_B <= eta * a * tof_out, and the dv actually used on that leg
    (accumulated only for epochs t_B <= 15 yr - 600 d, so that the window end does not masquerade as 'no hop').
Arc classes (applied to BOTH the incoming and the outgoing arc):
  plausible : bound heliocentric arc with perihelion > 0.4 AU and aphelion < 3.0 AU
  earthlike : perihelion > 0.6 AU, aphelion < 1.8 AU, inclination < 10 deg
              (what an Earth launch with v_inf <= 4 km/s can produce: q >= 0.60, Q <= 1.80 AU, i <= 7.7 deg).
Incoming arcs are randomly subsampled (K_IN per (B, t_B)) to bound the cost; outgoing arcs are
exhaustive.  Results are accumulated as histograms (bins 0..DV_MAX km/s step DV_BIN, overflow bin).

Outputs (analysis/ballistic/out/):
  junction_hist.npz        histograms  H_all[class, tof_out, bin], H_min[class, tof_out, bin], branch counts,
                           tmin_hist[class, eta, accel, tof-index], dvused_hist[class, eta, accel, bin]
  junction_per_epoch.csv   per (B, t_B): number of valid/plausible/earthlike in/out arcs and the fraction of
                           sampled incoming arcs with a feasible next hop for several (tof_out, accel), eta = 0.8
  junction_per_asteroid.csv  the same averaged over epochs
  fig_junction_dv_cdf.png, fig_junction_feasible_vs_tof.png
"""
import os, sys, time
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ephem import MU_SUN, AU, DAY, T_MAX_DAYS, asteroid_state, load_mea
from lambert import lambert, lambert_multirev

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
DT_EP = 30.0
TOF = np.unique(np.concatenate([np.arange(60.0, 600.0 + 1, 30.0), [365.0]]))   # 20 values
TB = np.arange(0.0, T_MAX_DAYS + 1e-9, DT_EP)      # 183 epochs (0..5460)
K_IN = 200                                         # sampled incoming arcs per (B, t_B)
TOF_MIN_1REV = 250.0                               # 1-rev arcs solved for tof >= this
DV_BIN, DV_MAX = 0.1, 40.0
NB = int(DV_MAX / DV_BIN) + 1                      # last bin = overflow
CLASSES = ('plausible', 'earthlike')
ACCELS = (0.25, 0.5, 0.8)                          # mm/s^2
ETAS = (0.7, 0.8, 0.9)
ETA_REPORT = 0.8
TOF_REPORT = (60.0, 120.0, 180.0, 365.0)
THRESH = (0.5, 1.0, 2.0, 3.0, 5.0)                 # km/s for branching-factor counts

_IDS, _EL = load_mea()
_NAST = len(_IDS)


def arc_class(r, v):
    """Boolean masks (plausible, earthlike) for arcs given position r and velocity v at that point."""
    R = np.linalg.norm(r, axis=-1)
    V2 = np.sum(v * v, axis=-1)
    h = np.cross(r, v)
    hm = np.linalg.norm(h, axis=-1)
    with np.errstate(divide='ignore', invalid='ignore'):
        a = 1.0 / (2.0 / R - V2 / MU_SUN)
        p = hm ** 2 / MU_SUN
        e = np.sqrt(np.maximum(1.0 - p / a, 0.0))
        e = np.where(a > 0, e, 2.0)
        q = p / (1 + e)
        Q = np.where(e < 1, p / (1 - e), np.inf)
        cosi = h[..., 2] / hm
    ok = np.isfinite(v[..., 0])
    plaus = ok & (a > 0) & (q > 0.4 * AU) & (Q < 3.0 * AU)
    earth = plaus & (q > 0.6 * AU) & (Q < 1.8 * AU) & (cosi > np.cos(np.radians(10.0)))
    return plaus, earth


def arcs_at_B(rB, rX, tofs, incoming):
    """Velocity at B of all arc types between B (nE,3) and X (nE,nO,nT,3).  Returns (nE,nO,nT,3 types,3)."""
    nE, nO, nT = rX.shape[:3]
    out = np.full((nE, nO, nT, 3, 3), np.nan)
    rB4 = rB[:, None, None, :]
    if incoming:
        _, v = lambert(rX, rB4, tofs[None, None, :], MU_SUN)
    else:
        v, _ = lambert(rB4, rX, tofs[None, None, :], MU_SUN)
    out[..., 0, :] = v
    cols = TOF >= TOF_MIN_1REV
    if cols.any():
        if incoming:
            _, vm = lambert_multirev(rX[:, :, cols], rB4, tofs[cols][None, None, :], MU_SUN, N=1)
        else:
            vm, _ = lambert_multirev(rB4, rX[:, :, cols], tofs[cols][None, None, :], MU_SUN, N=1)
        out[:, :, cols, 1:, :] = vm
    return out


def junction_one(jB, seed=0):
    rng = np.random.default_rng(seed + jB)
    others = np.array([j for j in range(_NAST) if j != jB])
    elO = _EL[others]
    nE, nO, nT = len(TB), len(others), len(TOF)
    rB, vB = asteroid_state(_EL[jB], TB)                           # (nE,3)
    tofs = TOF * DAY
    # ---- incoming arcs
    tA = TB[:, None] - TOF[None, :]                                # (nE,nT)
    valid_in = tA >= 0.0
    rA, _ = asteroid_state(elO, np.maximum(tA, 0.0))               # (nO,nE,nT,3)
    rA = np.transpose(rA, (1, 0, 2, 3))                            # (nE,nO,nT,3)
    v_in = arcs_at_B(rB, rA, tofs, incoming=True)                  # (nE,nO,nT,3,3)
    pl_in, ea_in = arc_class(rB[:, None, None, None, :], v_in)
    pl_in &= valid_in[:, None, :, None]
    ea_in &= valid_in[:, None, :, None]
    del rA
    # ---- outgoing arcs
    tC = TB[:, None] + TOF[None, :]
    valid_out = tC <= T_MAX_DAYS
    rC, _ = asteroid_state(elO, np.minimum(tC, T_MAX_DAYS))
    rC = np.transpose(rC, (1, 0, 2, 3))
    v_out = arcs_at_B(rB, rC, tofs, incoming=False)
    pl_out, ea_out = arc_class(rB[:, None, None, None, :], v_out)
    pl_out &= valid_out[:, None, :, None]
    ea_out &= valid_out[:, None, :, None]
    del rC
    # ---- junction statistics
    nA, nEta = len(ACCELS), len(ETAS)
    H_all = np.zeros((2, nT, NB), dtype=np.int64)
    H_min = np.zeros((2, nT, NB), dtype=np.int64)
    n_in_sampled = np.zeros(2, dtype=np.int64)
    branch = np.zeros((2, nT, len(THRESH)), dtype=np.int64)
    tmin_hist = np.zeros((2, nEta, nA, nT + 1), dtype=np.int64)
    dvused_hist = np.zeros((2, nEta, nA, NB), dtype=np.int64)
    budget = np.array([[[eta * acc * 1e-6 * tf * DAY for tf in TOF] for acc in ACCELS] for eta in ETAS])   # (nEta,nA,nT) km/s
    ieta_rep = ETAS.index(ETA_REPORT)
    itof_rep = [int(np.argmin(np.abs(TOF - tf))) for tf in TOF_REPORT]
    itf180 = int(np.argmin(np.abs(TOF - 180.0)))
    rows = []
    for ie in range(nE):
        vin = v_in[ie].reshape(-1, 3)                              # (nO*nT*3, 3)
        vout = v_out[ie]                                           # (nO,nT,3,3)
        masks_in = (pl_in[ie].reshape(-1), ea_in[ie].reshape(-1))
        masks_out = (pl_out[ie], ea_out[ie])                       # (nO,nT,3)
        row = dict(B=int(_IDS[jB]), tB=float(TB[ie]),
                   n_in_valid=int(valid_in[ie].sum() * nO), n_in_plausible=int(masks_in[0].sum()), n_in_earthlike=int(masks_in[1].sum()),
                   n_out_valid=int(valid_out[ie].sum() * nO), n_out_plausible=int(masks_out[0].sum()), n_out_earthlike=int(masks_out[1].sum()),
                   n_in_plausible_1rev=int(pl_in[ie][..., 1:].sum()), n_out_plausible_1rev=int(pl_out[ie][..., 1:].sum()))
        for ic in range(2):
            cand = np.where(masks_in[ic])[0]
            if len(cand) == 0 or masks_out[ic].sum() == 0:
                for tf in TOF_REPORT:
                    for acc in ACCELS:
                        row[f'frac_feas_{CLASSES[ic]}_tof{int(tf)}_a{acc}'] = np.nan
                row[f'medmin_{CLASSES[ic]}_tof180'] = np.nan
                continue
            idx = cand if len(cand) <= K_IN else rng.choice(cand, K_IN, replace=False)
            n_in_sampled[ic] += len(idx)
            dv = np.linalg.norm(vout[None] - vin[idx][:, None, None, None, :], axis=-1)   # (K,nO,nT,3)
            dv = np.where(masks_out[ic][None], dv, np.inf)
            dvC = dv.min(axis=3)                                                          # best arc type per C  (K,nO,nT)
            fin = np.isfinite(dvC)
            b = np.minimum((dvC[fin] / DV_BIN).astype(np.int64), NB - 1)
            it = np.broadcast_to(np.arange(nT)[None, None, :], dvC.shape)[fin]
            H_all[ic] += np.bincount(it * NB + b, minlength=nT * NB).reshape(nT, NB)
            mn = dvC.min(axis=1)                                                          # (K,nT)
            finm = np.isfinite(mn)
            bm = np.minimum((mn[finm] / DV_BIN).astype(np.int64), NB - 1)
            itm = np.broadcast_to(np.arange(nT)[None, :], mn.shape)[finm]
            H_min[ic] += np.bincount(itm * NB + bm, minlength=nT * NB).reshape(nT, NB)
            for k, th in enumerate(THRESH):
                branch[ic, :, k] += (dvC <= th).sum(axis=(0, 1))
            fits = TB[ie] + TOF[-1] <= T_MAX_DAYS          # greedy-hop stats only where every TOF fits in the window
            for ieta in range(nEta if fits else 0):
                for ia in range(nA):
                    feas = mn <= budget[ieta, ia][None, :]                                # (K,nT)
                    anyf = feas.any(axis=1)
                    first = np.where(anyf, feas.argmax(axis=1), nT)
                    tmin_hist[ic, ieta, ia] += np.bincount(first, minlength=nT + 1)
                    dvu = mn[np.arange(len(first))[anyf], first[anyf]]
                    dvused_hist[ic, ieta, ia] += np.bincount(np.minimum((dvu / DV_BIN).astype(np.int64), NB - 1), minlength=NB)
            for tf, itf in zip(TOF_REPORT, itof_rep):
                for ia, acc in enumerate(ACCELS):
                    row[f'frac_feas_{CLASSES[ic]}_tof{int(tf)}_a{acc}'] = \
                        float((mn[:, itf] <= budget[ieta_rep, ia, itf]).mean()) if TB[ie] + tf <= T_MAX_DAYS else np.nan
            m180 = mn[:, itf180]
            row[f'medmin_{CLASSES[ic]}_tof180'] = float(np.median(m180[np.isfinite(m180)])) if np.isfinite(m180).any() else np.nan
        rows.append(row)
    return jB, H_all, H_min, branch, n_in_sampled, rows, tmin_hist, dvused_hist


def report_only():
    """Re-run the report/plots from the saved npz + per-epoch CSV (REPORT_ONLY=1)."""
    import pandas as pd
    d = np.load(os.path.join(OUT, 'junction_hist.npz'))
    df = pd.read_csv(os.path.join(OUT, 'junction_per_epoch.csv'))
    report(d['H_all'], d['H_min'], d['branch'], d['n_in_sampled'], df)
    report_greedy(d['tmin_hist'], d['dvused_hist'])
    plots(d['H_all'], d['H_min'])


def main():
    import multiprocessing as mp
    import pandas as pd
    if os.environ.get('REPORT_ONLY'):
        return report_only()
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    nproc = int(os.environ.get('NPROC', 8))
    nast = int(os.environ.get('NAST', _NAST))       # subset of junction asteroids B (A and C are always all 300)
    jB_list = [int(j) for j in np.round(np.linspace(0, _NAST - 1, nast))]   # evenly spaced over the catalogue
    H_all = np.zeros((2, len(TOF), NB), dtype=np.int64)
    H_min = np.zeros_like(H_all)
    branch = np.zeros((2, len(TOF), len(THRESH)), dtype=np.int64)
    n_in = np.zeros(2, dtype=np.int64)
    tmin_hist = np.zeros((2, len(ETAS), len(ACCELS), len(TOF) + 1), dtype=np.int64)
    dvused_hist = np.zeros((2, len(ETAS), len(ACCELS), NB), dtype=np.int64)
    rows = []
    ctx = mp.get_context('fork')
    with ctx.Pool(nproc) as pool:
        for k, (jB, ha, hm, br, ni, rw, th, dh) in enumerate(pool.imap_unordered(junction_one, jB_list, chunksize=1)):
            H_all += ha; H_min += hm; branch += br; n_in += ni; rows.extend(rw); tmin_hist += th; dvused_hist += dh
            if k % 25 == 0:
                print(f"  {k+1}/{nast} asteroids done, elapsed {time.time()-t0:.1f}s", flush=True)
    print(f"junction scan done in {time.time()-t0:.1f}s")
    np.savez_compressed(os.path.join(OUT, 'junction_hist.npz'), H_all=H_all, H_min=H_min, branch=branch, n_in_sampled=n_in,
                        tmin_hist=tmin_hist, dvused_hist=dvused_hist, accels=np.array(ACCELS), etas=np.array(ETAS),
                        tof_days=TOF, dv_bin=DV_BIN, dv_max=DV_MAX, classes=np.array(CLASSES), thresh=np.array(THRESH),
                        tB_days=TB, K_IN=K_IN, nast=nast, jB_ids=_IDS[jB_list])
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, 'junction_per_epoch.csv'), index=False, float_format='%.5g')
    per_ast = df.groupby('B').mean(numeric_only=True).drop(columns=['tB'])
    per_ast.to_csv(os.path.join(OUT, 'junction_per_asteroid.csv'), float_format='%.5g')
    report(H_all, H_min, branch, n_in, df)
    report_greedy(tmin_hist, dvused_hist)
    plots(H_all, H_min)


def _cdf_at(H, x):
    """fraction of samples with dv <= x from histogram H[..., NB]"""
    edges = np.arange(NB + 1) * DV_BIN
    n = H.sum(axis=-1)
    k = np.searchsorted(edges, x, side='right') - 1
    full = H[..., :k].sum(axis=-1)
    part = H[..., k] * ((x - edges[k]) / DV_BIN) if k < NB else 0
    with np.errstate(invalid='ignore', divide='ignore'):
        return (full + part) / n


def _pct(H, p):
    """percentiles p (scalar or vector) of the dv distribution in histogram H[..., NB]; result shape H.shape[:-1] + p.shape"""
    p = np.atleast_1d(np.asarray(p, float))
    c = np.cumsum(H, axis=-1) / np.maximum(H.sum(axis=-1, keepdims=True), 1)
    edges = (np.arange(NB) + 1) * DV_BIN
    out = np.array([edges[np.minimum(np.searchsorted(ci, p), NB - 1)] if ci[-1] > 0 else np.full(p.shape, np.nan)
                    for ci in c.reshape(-1, NB)])
    return out.reshape(H.shape[:-1] + p.shape)


def report(H_all, H_min, branch, n_in, df):
    print("\n=== Asteroid -> asteroid junction delta-v statistics ===")
    print(f"grid: epochs step {DT_EP} d ({len(TB)}), TOF {TOF[0]:.0f}-{TOF[-1]:.0f} d step 30 d + 365 d ({len(TOF)} values), {df.B.nunique()} junction asteroids B (evenly spaced ids), all {_NAST} asteroids as A and C; "
          f"arc types: 1 single-rev + 2 one-rev (tof >= {TOF_MIN_1REV:.0f} d); {K_IN} sampled incoming arcs per (B,t_B)")
    print(f"sampled incoming arcs: plausible {n_in[0]}, earthlike {n_in[1]}; "
          f"mean valid/plausible/earthlike incoming arcs per (B,t_B): {df.n_in_valid.mean():.0f}/{df.n_in_plausible.mean():.0f}/{df.n_in_earthlike.mean():.0f} "
          f"(of which 1-rev plausible {df.n_in_plausible_1rev.mean():.0f}); "
          f"outgoing: {df.n_out_valid.mean():.0f}/{df.n_out_plausible.mean():.0f}/{df.n_out_earthlike.mean():.0f} (1-rev plausible {df.n_out_plausible_1rev.mean():.0f})")
    for ic, cl in enumerate(CLASSES):
        print(f"\n--- class '{cl}' (both arcs) ---")
        print("  tof_out[d] | unconditional dv_B pct 10/50/90 [km/s] | best-next-hop min_C dv_B pct 10/50/90 | "
              f"P(min_C dv <= eta*a*tof), eta={ETA_REPORT}, a=0.25/0.5/0.8 mm/s2 | same for eta=0.7 (a=0.25/0.8) | eta=0.9 (a=0.25/0.8) | "
              "P(dv<=budget) unconditional a=0.25/0.8 | mean #C with dv<=0.5/1/2/3/5 km/s")
        for it, tf in enumerate(TOF):
            pa = _pct(H_all[ic, it], np.array([0.1, 0.5, 0.9]))
            pm = _pct(H_min[ic, it], np.array([0.1, 0.5, 0.9]))
            feas = [_cdf_at(H_min[ic, it], ETA_REPORT * a * 1e-6 * tf * DAY) for a in ACCELS]
            feas7 = [_cdf_at(H_min[ic, it], 0.7 * a * 1e-6 * tf * DAY) for a in (0.25, 0.8)]
            feas9 = [_cdf_at(H_min[ic, it], 0.9 * a * 1e-6 * tf * DAY) for a in (0.25, 0.8)]
            feas_u = [_cdf_at(H_all[ic, it], ETA_REPORT * a * 1e-6 * tf * DAY) for a in (0.25, 0.8)]
            nin = H_min[ic, it].sum()
            bf = branch[ic, it] / max(nin, 1)
            print(f"  {tf:5.0f} | {pa[0]:5.2f} {pa[1]:5.2f} {pa[2]:5.2f} | {pm[0]:5.2f} {pm[1]:5.2f} {pm[2]:5.2f} | "
                  f"{feas[0]*100:5.1f}% {feas[1]*100:5.1f}% {feas[2]*100:5.1f}% | {feas7[0]*100:5.1f}% {feas7[1]*100:5.1f}% | {feas9[0]*100:5.1f}% {feas9[1]*100:5.1f}% | "
                  f"{feas_u[0]*100:6.3f}% {feas_u[1]*100:6.3f}% | "
                  f"{bf[0]:5.2f} {bf[1]:5.2f} {bf[2]:5.2f} {bf[3]:5.2f} {bf[4]:5.2f}   (n_in={nin})")
    print("\nper-asteroid: fraction of sampled incoming arcs with a feasible next hop (plausible, tof 180 d, a=0.5, eta=0.8):")
    col = 'frac_feas_plausible_tof180_a0.5'
    pa = df.groupby('B')[col].mean()
    print("  percentiles 5/25/50/75/95:", np.nanpercentile(pa, [5, 25, 50, 75, 95]).round(3))
    print("  asteroids with mean frac < 0.05:", pa[pa < 0.05].index.tolist())
    print("  asteroids with mean frac < 0.20:", pa[pa < 0.20].index.tolist())
    col = 'frac_feas_plausible_tof365_a0.5'
    pa = df.groupby('B')[col].mean()
    print("  same for tof 365 d: percentiles 5/25/50/75/95:", np.nanpercentile(pa, [5, 25, 50, 75, 95]).round(3),
          "; asteroids < 0.2:", pa[pa < 0.2].index.tolist())
    # inclination dependence of the per-asteroid feasibility
    inc = dict(zip(_IDS, _EL[:, 2]))
    pa180 = df.groupby('B')['frac_feas_plausible_tof180_a0.5'].mean()
    incs = np.array([inc[b] for b in pa180.index])
    for lo, hi in ((0, 10), (10, 20), (20, 30), (30, 180)):
        m = (incs >= lo) & (incs < hi)
        print(f"  inclination {lo:3d}-{hi:3d} deg ({m.sum():3d} asteroids): mean feasible fraction (tof180, a=0.5) {np.nanmean(pa180.values[m]):.3f}")


def report_greedy(tmin_hist, dvused_hist):
    print("\n=== Greedy 'next flyby as soon as possible' statistics (shortest grid TOF with a feasible next hop) ===")
    print("  criterion: min_C |v_out - v_in| <= eta * a * TOF_out")
    print("  class | eta | a[mm/s2] | P(no feasible hop within 600 d) | shortest feasible TOF: mean [d], pct 25/50/75 | implied flybys/yr = 365/mean | "
          "junction dv on that leg: mean, pct 25/50/75 [km/s] | dv per year of flying [km/s/yr]")
    edges_dv = (np.arange(NB) + 0.5) * DV_BIN
    for ic, cl in enumerate(CLASSES):
        for ieta, eta in enumerate(ETAS):
            for ia, acc in enumerate(ACCELS):
                h = tmin_hist[ic, ieta, ia]
                n = h.sum()
                pnone = h[-1] / max(n, 1)
                hf = h[:-1]
                if hf.sum() == 0:
                    continue
                mean_t = np.sum(hf * TOF) / hf.sum()
                c = np.cumsum(hf) / hf.sum()
                pct_t = [TOF[np.searchsorted(c, p)] for p in (0.25, 0.5, 0.75)]
                hd = dvused_hist[ic, ieta, ia]
                mean_dv = np.sum(hd * edges_dv) / max(hd.sum(), 1)
                cd = np.cumsum(hd) / max(hd.sum(), 1)
                pct_dv = [edges_dv[np.searchsorted(cd, p)] for p in (0.25, 0.5, 0.75)]
                print(f"  {cl:9s} | {eta:.1f} | {acc:4.2f} | {pnone*100:5.1f}% | {mean_t:6.1f}  {pct_t[0]:4.0f}/{pct_t[1]:4.0f}/{pct_t[2]:4.0f} | {365/mean_t:5.2f} | "
                      f"{mean_dv:5.2f}  {pct_dv[0]:4.2f}/{pct_dv[1]:4.2f}/{pct_dv[2]:4.2f} | {mean_dv*365/mean_t:5.2f}   (n={n})")


def plots(H_all, H_min):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    edges = np.arange(NB) * DV_BIN
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for ic, cl in enumerate(CLASSES):
        ax = axes[ic]
        for it, tf in enumerate(TOF):
            if tf not in (60, 120, 180, 365, 600):
                continue
            c = np.cumsum(H_min[ic, it]) / max(H_min[ic, it].sum(), 1)
            ax.plot(edges + DV_BIN, c, label=f'best next hop, tof_out={tf:.0f} d')
            c2 = np.cumsum(H_all[ic, it]) / max(H_all[ic, it].sum(), 1)
            ax.plot(edges + DV_BIN, c2, ls=':', color=ax.lines[-1].get_color())
        for a in ACCELS:
            ax.axvline(ETA_REPORT * a * 1e-6 * 180 * DAY, color='grey', ls='--', lw=0.8)
        ax.set_xlim(0, 15); ax.set_ylim(0, 1); ax.grid(alpha=.3)
        ax.set_xlabel('junction delta-v at B [km/s]'); ax.set_ylabel('CDF')
        ax.set_title(f"class '{cl}' (solid: min over C; dotted: all C)\nvertical: eta*a*180d for a=0.25/0.5/0.8 mm/s2")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'fig_junction_dv_cdf.png'), dpi=130)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
    for ic, cl in enumerate(CLASSES):
        ax = axes[ic]
        for a in ACCELS:
            for eta, ls in zip(ETAS, (':', '-', '--')):
                f = [_cdf_at(H_min[ic, it], eta * a * 1e-6 * tf * DAY) for it, tf in enumerate(TOF)]
                ax.plot(TOF, f, ls=ls, marker='o' if eta == ETA_REPORT else None, ms=3,
                        label=f'a={a} mm/s2, eta={eta}' if eta == ETA_REPORT else None,
                        color={0.25: '#4C72B0', 0.5: '#55A868', 0.8: '#C44E52'}[a])
        ax.set_xlabel('next-leg TOF [d]'); ax.set_ylim(0, 1); ax.grid(alpha=.3)
        ax.set_title(f"class '{cl}': P(best next hop feasible)\n(solid eta=0.8, dotted 0.7, dashed 0.9)")
        ax.legend(fontsize=8)
    axes[0].set_ylabel('fraction of incoming arcs with a feasible next hop')
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'fig_junction_feasible_vs_tof.png'), dpi=130)
    plt.close('all')


if __name__ == '__main__':
    main()
