"""
Zero-fuel DOUBLE flybys (Earth -> A -> B on one Keplerian arc, launch v_inf <= 4 km/s).

A ballistic spacecraft has 4 free parameters (launch date + v_inf vector); a flyby costs 2 net constraints
(3 position components minus the free flyby time), so a ballistic spacecraft can generically hit 2 asteroids
at isolated parameter points.  This script finds those points:

  1. coarse: for every asteroid A and every feasible cell (tL, tof) of the Earth scan (v_inf <= 4 km/s, arc type
     = the best one recorded in earth_to_asteroid_vinf_grid.npz), the arc is coasted from the flyby of A to the
     end of the window (step DT_COAST d) and the minimum distance to every other asteroid B is recorded.
     (A -> B pairs with miss < D_CAND are candidates; they are clustered per (A, B, arc type, ~month of t_B).)
  2. refine: for the best MAX_PER_PAIR clusters of each ordered (A, B) pair solve r_sc(t_B; tL, t_A) - r_B(t_B) = 0
     for (tL, t_A, t_B) with a vectorised damped Gauss-Newton (central-difference Jacobian; the spacecraft arc is the
     Lambert arc Earth(tL) -> A(t_A) with the same N/branch as the coarse cell).  refine_one() is the scalar
     scipy.optimize.least_squares version of the same solve, kept for spot checks.
     A pair is accepted when the residual miss at B is <= 1000 km (the contest flyby criterion), the miss at A is
     exact by construction, tL >= 0, t_B <= 15 yr and the launch v_inf <= 4 km/s.
  3. greedy maximum matching of the accepted pairs -> how many asteroids could be covered by cost-1 spacecraft
     flying 2 asteroids each.

Outputs (analysis/ballistic/out/): ballistic_pairs.csv (accepted pairs), ballistic_pair_candidates.csv (clusters),
ballistic_pairs_summary.txt.
"""
import os, sys, time
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ephem import MU_SUN, AU, DAY, T_MAX_DAYS, earth_state, asteroid_state, load_mea, propagate_kepler
from lambert import lambert, lambert_multirev

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
DT_COAST = 2.0            # d, coast sampling after the flyby of A
D_CAND = 0.06 * AU        # km, candidate threshold on the coarse grid
MAX_CELLS = 500           # feasible cells per asteroid used (random subsample above this)
VINF_MAX = 4.0
MISS_OK = 1000.0          # km, contest flyby criterion
MAX_PER_PAIR = 1          # clusters refined per ordered (A,B) pair

_IDS, _EL = load_mea()
with np.load(os.path.join(OUT, 'earth_to_asteroid_vinf_grid.npz')) as _G:   # eager loads: an open NpzFile is not fork-safe
    TL, TOF = _G['tL_days'], _G['tof_days']
    _VINF = _G['vinf']                                              # (300, nL, nT) float16
    _NREV = _G['nrev_best']                                         # (300, nL, nT) int8
_TC = np.arange(0.0, T_MAX_DAYS + 1e-9, DT_COAST)                 # coast time grid
_RAST, _ = asteroid_state(_EL, _TC)                                 # (300, nTc, 3)  ~ 20 MB


def arc_v1(tL, tA, N, branch, elA):
    """Earth(tL) -> A(tA) Lambert arc velocity at Earth (km/s), Earth state; arc type (N, branch)."""
    rE, vE = earth_state(tL)
    rA, _ = asteroid_state(elA, tA)
    tof = (tA - tL) * DAY
    if N == 0:
        v1, _ = lambert(rE, rA, tof, MU_SUN)
    else:
        v1m, _ = lambert_multirev(rE, rA, tof, MU_SUN, N=N)
        v1 = v1m[branch]
    return rE, vE, v1


def coarse_one(jA, seed=0):
    rng = np.random.default_rng(seed + jA)
    vinf = _VINF[jA].astype(np.float32)
    nrev = _NREV[jA]
    cells = np.argwhere(vinf <= VINF_MAX)
    if len(cells) == 0:
        return jA, []
    if len(cells) > MAX_CELLS:
        cells = cells[rng.choice(len(cells), MAX_CELLS, replace=False)]
    out = []
    for iL, iT in cells:
        tL, tof = TL[iL], TOF[iT]
        tA = tL + tof
        N = int(nrev[iL, iT])
        rE, vE = earth_state(tL)
        rA, _ = asteroid_state(_EL[jA], tA)
        if N == 0:
            v1, _ = lambert(rE, rA, tof * DAY, MU_SUN)
            branch = 0
            vi = np.linalg.norm(v1 - vE)
        else:
            v1m, _ = lambert_multirev(rE, rA, tof * DAY, MU_SUN, N=N)
            vv = np.linalg.norm(v1m - vE, axis=-1)
            branch = int(np.nanargmin(vv))
            v1 = v1m[branch]
            vi = vv[branch]
        if not np.isfinite(vi) or vi > VINF_MAX + 0.05:
            continue
        # coast from the flyby of A to the end of the window
        k0 = int(np.searchsorted(_TC, tA))
        tc = _TC[k0:]
        if len(tc) < 2:
            continue
        rsc, _ = propagate_kepler(np.repeat(rE[None], len(tc), 0), np.repeat(v1[None], len(tc), 0), (tc - tL) * DAY)
        d = np.linalg.norm(_RAST[:, k0:, :] - rsc[None], axis=-1)          # (300, ntc)
        d[jA] = np.inf
        # exclude the first ~10 d after the flyby of A (same place) -- B must be a different asteroid anyway
        kmin = np.argmin(d, axis=1)
        dmin = d[np.arange(300), kmin]
        for jB in np.where(dmin < D_CAND)[0]:
            out.append((jA, int(jB), float(tL), float(tA), N, branch, float(vi), float(tc[kmin[jB]]), float(dmin[jB] / AU)))
    return jA, out


def refine_one(c):
    """c = (jA, jB, tL, tA, N, branch, vinf, tB, miss_au) -> refined dict"""
    from scipy.optimize import least_squares
    jA, jB, tL0, tA0, N, branch, vi0, tB0, miss0 = c
    elA, elB = _EL[jA], _EL[jB]

    def resid(p):
        tL, tA, tB = p
        if tA - tL < 2.0 or tB - tA < 1.0:
            return np.full(3, 1e6)
        rE, vE, v1 = arc_v1(tL, tA, N, branch, elA)
        if not np.isfinite(v1[0]):
            return np.full(3, 1e6)
        rsc, _ = propagate_kepler(rE[None], v1[None], np.array([(tB - tL) * DAY]))
        rB, _ = asteroid_state(elB, tB)
        return (rsc[0] - rB) / 1e3                                          # units of 1000 km

    try:
        sol = least_squares(resid, [tL0, tA0, tB0], bounds=([0.0, 0.0, 0.0], [T_MAX_DAYS] * 3),
                            x_scale=[5.0, 5.0, 5.0], xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=400)
        tL, tA, tB = sol.x
        miss = np.linalg.norm(sol.fun) * 1e3
        rE, vE, v1 = arc_v1(tL, tA, N, branch, elA)
        vi = float(np.linalg.norm(v1 - vE))
        rsc, vsc = propagate_kepler(rE[None], v1[None], np.array([(tB - tL) * DAY]))
        _, vB = asteroid_state(elB, tB)
        _, vA = asteroid_state(elA, tA)
        rscA, vscA = propagate_kepler(rE[None], v1[None], np.array([(tA - tL) * DAY]))
        ok = bool(miss <= MISS_OK and vi <= VINF_MAX and tL >= 0 and tB <= T_MAX_DAYS and tA > tL and tB > tA)
    except Exception as ex:                                                  # pragma: no cover
        return dict(A=int(_IDS[jA]), B=int(_IDS[jB]), ok=False, err=str(ex))
    return dict(A=int(_IDS[jA]), B=int(_IDS[jB]), ok=ok, tL=float(tL), tA=float(tA), tB=float(tB), N=N, branch=branch,
                vinf=vi, miss_km=float(miss), vrelA=float(np.linalg.norm(vscA[0] - vA)), vrelB=float(np.linalg.norm(vsc[0] - vB)),
                tL0=tL0, tA0=tA0, tB0=tB0, miss0_au=miss0, nfev=int(sol.nfev))


def _ast_states(jidx, t):
    """Elementwise asteroid positions: asteroid index jidx[k] at time t[k] (days). Returns (K,3)."""
    r = np.empty((len(t), 3))
    for j in np.unique(jidx):
        m = jidx == j
        r[m] = asteroid_state(_EL[j], t[m])[0]
    return r


def resid_vec(p, jA, jB, N, branch):
    """Miss vector at B (km) for candidates p = (tL, tA, tB) rows; also returns v_inf and the arc velocity."""
    K = len(p)
    tL, tA, tB = p[:, 0], p[:, 1], p[:, 2]
    rE, vE = earth_state(tL)
    rA = _ast_states(jA, tA)
    tof = (tA - tL) * DAY
    v1 = np.full((K, 3), np.nan)
    for n in np.unique(N):
        for b in (0, 1) if n > 0 else (0,):
            m = (N == n) & (branch == b) if n > 0 else (N == n)
            if not m.any():
                continue
            if n == 0:
                v1[m] = lambert(rE[m], rA[m], tof[m], MU_SUN)[0]
            else:
                v1[m] = lambert_multirev(rE[m], rA[m], tof[m], MU_SUN, N=int(n))[0][:, b]
    rsc, vsc = propagate_kepler(rE, v1, (tB - tL) * DAY)
    rB = _ast_states(jB, tB)
    F = rsc - rB
    bad = (tA - tL < 2.0) | (tB - tA < 1.0) | ~np.isfinite(v1[:, 0]) | ~np.isfinite(F[:, 0])
    F[bad] = 1e9
    return F, np.linalg.norm(v1 - vE, axis=1), vsc


def refine_vec(cl, itmax=40, h=0.02, lam0=1e-3):
    """Vectorised Levenberg-Marquardt on (tL, tA, tB) for all candidate rows of DataFrame cl.
    Residual F = r_sc(tB) - r_B(tB) [km]; dF/dtB = (v_sc - v_B) analytically, dF/dtL and dF/dtA by central differences.
    Per-row damping lambda: step accepted (lambda/3) if |F| decreases, else rejected (lambda*10)."""
    jA = cl.jA.values.astype(int); jB = cl.jB.values.astype(int)
    N = cl.N.values.astype(int); br = cl.branch.values.astype(int)
    p = cl[['tL', 'tA', 'tB']].values.astype(float).copy()
    K = len(p)
    F, _, vsc = resid_vec(p, jA, jB, N, br)
    _, vB = _ast_states_v(jB, p[:, 2])
    miss = np.linalg.norm(F, axis=1)
    lam = np.full(K, lam0)
    active = miss > 1.0
    for it in range(itmax):
        idx = np.where(active)[0]
        if len(idx) == 0:
            break
        pa = p[idx]
        J = np.empty((len(idx), 3, 3))
        for k in range(2):
            dp = np.zeros(3); dp[k] = h
            Fp, _, _ = resid_vec(pa + dp, jA[idx], jB[idx], N[idx], br[idx])
            Fm, _, _ = resid_vec(pa - dp, jA[idx], jB[idx], N[idx], br[idx])
            J[:, :, k] = (Fp - Fm) / (2 * h)
        J[:, :, 2] = (vsc[idx] - vB[idx]) * DAY
        Fa = F[idx]
        JtJ = np.einsum('kij,kil->kjl', J, J)
        JtF = np.einsum('kij,ki->kj', J, Fa)
        diag = np.einsum('kjj->kj', JtJ)
        A = JtJ + lam[idx][:, None, None] * np.eye(3)[None] * np.maximum(diag, 1e-12)[:, :, None]
        good = np.isfinite(A).all(axis=(1, 2)) & np.isfinite(JtF).all(axis=1)
        step = np.zeros((len(idx), 3))
        if good.any():
            try:
                step[good] = -np.linalg.solve(A[good], JtF[good][..., None])[..., 0]
            except np.linalg.LinAlgError:
                for i in np.where(good)[0]:
                    try:
                        step[i] = -np.linalg.solve(A[i], JtF[i])
                    except np.linalg.LinAlgError:
                        pass
        step = np.clip(step, -200.0, 200.0)
        pn = pa + step
        Fn, _, vscn = resid_vec(pn, jA[idx], jB[idx], N[idx], br[idx])
        missn = np.linalg.norm(Fn, axis=1)
        acc = missn < miss[idx]
        ia = idx[acc]
        p[ia] = pn[acc]; F[ia] = Fn[acc]; miss[ia] = missn[acc]; vsc[ia] = vscn[acc]
        _, vB[ia] = _ast_states_v(jB[ia], p[ia, 2])
        lam[ia] = np.maximum(lam[ia] / 3.0, 1e-9)
        ir = idx[~acc]
        lam[ir] = lam[ir] * 10.0
        active[idx] = (miss[idx] > 1.0) & (lam[idx] < 1e6)
    F, vinf, vsc = resid_vec(p, jA, jB, N, br)
    miss = np.linalg.norm(F, axis=1)
    rE, vE = earth_state(p[:, 0])
    _, vA = _ast_states_v(jA, p[:, 1])
    _, vB = _ast_states_v(jB, p[:, 2])
    rA = _ast_states(jA, p[:, 1])
    v1 = np.full((K, 3), np.nan)
    tof = (p[:, 1] - p[:, 0]) * DAY
    for n in np.unique(N):
        for b in (0, 1) if n > 0 else (0,):
            m = (N == n) & (br == b) if n > 0 else (N == n)
            if not m.any():
                continue
            v1[m] = lambert(rE[m], rA[m], tof[m], MU_SUN)[0] if n == 0 else lambert_multirev(rE[m], rA[m], tof[m], MU_SUN, N=int(n))[0][:, b]
    _, vscA = propagate_kepler(rE, v1, tof)
    out = cl[['jA', 'jB', 'N', 'branch', 'tL', 'tA', 'tB', 'miss_au', 'vinf']].copy()
    out = out.rename(columns=dict(tL='tL0', tA='tA0', tB='tB0', miss_au='miss0_au', vinf='vinf0'))
    out['A'] = _IDS[jA]; out['B'] = _IDS[jB]
    out['tL'] = p[:, 0]; out['tA'] = p[:, 1]; out['tB'] = p[:, 2]
    out['miss_km'] = miss; out['vinf'] = vinf
    out['vrelA'] = np.linalg.norm(vscA - vA, axis=1); out['vrelB'] = np.linalg.norm(vsc - vB, axis=1)
    out['ok'] = (miss <= MISS_OK) & (vinf <= VINF_MAX) & (out.tL >= 0) & (out.tB <= T_MAX_DAYS) & (out.tA > out.tL) & (out.tB > out.tA)
    return out


def _ast_states_v(jidx, t):
    r = np.empty((len(t), 3)); v = np.empty((len(t), 3))
    for j in np.unique(jidx):
        m = jidx == j
        r[m], v[m] = asteroid_state(_EL[j], t[m])
    return r, v


def greedy_matching(pairs):
    """pairs: list of (A,B) ids -> greedy maximum matching (lowest-degree first)."""
    from collections import defaultdict
    adj = defaultdict(set)
    for a, b in pairs:
        adj[a].add(b); adj[b].add(a)
    matched = {}
    for a in sorted(adj, key=lambda k: len(adj[k])):
        if a in matched:
            continue
        cands = [b for b in adj[a] if b not in matched]
        if cands:
            b = min(cands, key=lambda k: len(adj[k]))
            matched[a] = b; matched[b] = a
    return matched


def main():
    import multiprocessing as mp
    import pandas as pd
    t0 = time.time()
    nproc = int(os.environ.get('NPROC', 8))
    ctx = mp.get_context('fork')
    cands = []
    n_cells = 0
    with ctx.Pool(nproc) as pool:
        for k, (jA, out) in enumerate(pool.imap_unordered(coarse_one, range(300), chunksize=1)):
            cands.extend(out)
            if k % 50 == 0:
                print(f"  coarse {k+1}/300, {len(cands)} candidate (cell, B) hits, elapsed {time.time()-t0:.0f}s", flush=True)
    print(f"coarse scan done in {time.time()-t0:.0f}s: {len(cands)} (A-cell, B) hits with miss < {D_CAND/AU:.2f} AU")
    if not cands:
        return
    cd = pd.DataFrame(cands, columns=['jA', 'jB', 'tL', 'tA', 'N', 'branch', 'vinf', 'tB', 'miss_au'])
    cd['A'] = _IDS[cd.jA.values]; cd['B'] = _IDS[cd.jB.values]
    cd['tB_month'] = (cd.tB / 30.0).round().astype(int)
    cl = cd.sort_values('miss_au').groupby(['jA', 'jB', 'N', 'tB_month'], as_index=False).first()
    cl.to_csv(os.path.join(OUT, 'ballistic_pair_candidates.csv'), index=False, float_format='%.6g')
    print(f"  {len(cl)} candidate clusters (A, B, N, month of t_B); miss < 0.01 AU: {(cl.miss_au < 0.01).sum()}, < 0.02: {(cl.miss_au < 0.02).sum()}, "
          f"< 0.04: {(cl.miss_au < 0.04).sum()}; distinct A: {cl.jA.nunique()}, distinct B: {cl.jB.nunique()}", flush=True)
    t1 = time.time()
    # refine: best MAX_PER_PAIR clusters per ordered (A, B) pair, vectorised Gauss-Newton in chunks (parallel over chunks)
    cl = cl.sort_values('miss_au')
    cl = cl.groupby(['jA', 'jB'], as_index=False).head(MAX_PER_PAIR).reset_index(drop=True)
    chunks = [cl.iloc[i:i + 2000] for i in range(0, len(cl), 2000)]
    with ctx.Pool(nproc) as pool:
        parts = pool.map(refine_vec, chunks, chunksize=1)
    rf = pd.concat(parts, ignore_index=True)
    rf.to_csv(os.path.join(OUT, 'ballistic_pairs_refined_all.csv'), index=False, float_format='%.12g')
    print(f"refinement of {len(rf)} clusters (<= {MAX_PER_PAIR} per ordered pair) done in {time.time()-t1:.0f}s; converged to a flyby pair "
          f"(miss <= {MISS_OK:.0f} km, v_inf <= 4): {rf.ok.sum()}; miss <= 1000 km but v_inf > 4: {((rf.miss_km <= MISS_OK) & (rf.vinf > VINF_MAX)).sum()}; "
          f"not converged: {(rf.miss_km > MISS_OK).sum()}")
    good = rf[rf.ok].copy()
    good = good.sort_values('vinf')
    good.to_csv(os.path.join(OUT, 'ballistic_pairs.csv'), index=False, float_format='%.12g')
    lines = []
    lines.append(f"coarse grid: launch cells from earth_to_asteroid_vinf_grid.npz (step {TL[1]-TL[0]:.0f} d x {TOF[1]-TOF[0]:.0f} d), max {MAX_CELLS} cells per A, coast step {DT_COAST:.0f} d, candidate miss < {D_CAND/AU:.2f} AU")
    lines.append(f"candidate clusters: {len(cl)}; refined OK pairs: {len(good)} (ordered A->B; {len(set(map(tuple, np.sort(good[['A','B']].values, axis=1))))} unordered)")
    if len(good):
        pairs = list(zip(good.A, good.B))
        m = greedy_matching(pairs)
        cov = len(m)
        deg = pd.Series([a for a, b in pairs] + [b for a, b in pairs]).value_counts()
        lines.append(f"asteroids appearing in >= 1 feasible pair: {deg.size}; greedy matching covers {cov} asteroids with {cov//2} cost-1 spacecraft (0.5 per asteroid)")
        lines.append(f"v_inf of accepted pairs: min {good.vinf.min():.2f}, median {good.vinf.median():.2f} km/s; arc type N=0: {(good.N==0).sum()}, N>=1: {(good.N>=1).sum()}")
        lines.append(f"pairs with v_inf <= 3.5 km/s (margin for a real solution): {(good.vinf <= 3.5).sum()}; <= 3.0: {(good.vinf <= 3.0).sum()}")
        lines.append(f"time between the two flybys: median {np.median(good.tB - good.tA):.0f} d, min {np.min(good.tB - good.tA):.0f}, max {np.max(good.tB - good.tA):.0f}")
        lines.append(f"asteroids with the most pair partners: {deg.head(10).to_dict()}")
        # how many asteroids can be covered by 2-flyby ballistic spacecraft that also have no other option, etc. is left to the planner
    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(OUT, 'ballistic_pairs_summary.txt'), 'w') as f:
        f.write(txt + "\n")
    print(f"total {time.time()-t0:.0f}s")


if __name__ == '__main__':
    main()
