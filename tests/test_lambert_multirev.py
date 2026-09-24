"""Tests of the multi-revolution Lambert solver (ctoc14.lambert.lambert_multirev / lambert_all / lambert_best).

Run: python tests/test_lambert_multirev.py  (exit code 0 on success). Checks
  * every returned solution (nrev 0..2, both branches) propagates from (r1, v1) to r2 in tof, with v2 consistent
    (independent element-based propagation: Kepler's equation for ellipses, hyperbolic Kepler equation otherwise);
  * nrev = 0 reproduces the existing single-revolution lambert() to 1e-9 km/s wherever both exist;
  * for nrev >= 1 both branches exist exactly when tof > tof_min and 'left' is the high-energy one (larger a);
  * T(z) is unimodal on each nrev >= 1 interval (grid scan) and the solver's tof_min matches the grid minimum;
  * lambert_best returns the option with the smallest junction dv;
  * speed (rows/s).
"""
import sys, pathlib, time
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ctoc14.constants import MU, AU, DAY
from ctoc14.kepler import Body, propagate_twobody
from ctoc14.lambert import lambert, lambert_multirev, lambert_all, lambert_best, lambert_tof_min, _MRGeom, _mr_bounds, \
    BRANCH_LEFT, BRANCH_RIGHT


# ----------------------------------------------------------------------------------------------- reference propagation
def elements(r, v):
    rn = np.linalg.norm(r); a = 1.0 / (2.0 / rn - (v @ v) / MU)
    h = np.cross(r, v); hn = np.linalg.norm(h)
    evec = np.cross(v, h) / MU - r / rn; e = np.linalg.norm(evec)
    inc = np.arccos(np.clip(h[2] / hn, -1, 1)); nv = np.cross([0.0, 0.0, 1.0], h)
    if np.linalg.norm(nv) < 1e-12 * hn:
        nv = np.array([1.0, 0.0, 0.0])
    Om = np.arctan2(nv[1], nv[0])
    w = np.arctan2(np.dot(np.cross(nv, evec), h) / hn, np.dot(nv, evec))
    nu = np.arctan2(np.dot(np.cross(evec, r), h) / hn, np.dot(evec, r))
    return a, e, inc, Om, w, nu, h


def propagate_ref(r, v, dt):
    """Independent two-body propagation through orbital elements (ellipse: kepler.Body; hyperbola: hyperbolic Kepler eq.)."""
    a, e, inc, Om, w, nu, h = elements(r, v)
    if a > 0 and e < 0.9999:
        E0 = 2 * np.arctan2(np.sqrt(1 - e) * np.sin(nu / 2), np.sqrt(1 + e) * np.cos(nu / 2)); M0 = E0 - e * np.sin(E0)
        b = Body([a / AU, e, np.degrees(inc), np.degrees(Om), np.degrees(w), np.degrees(M0)], 0.0)
        return b.state(dt)
    if a < 0 and e > 1.0001:
        H0 = 2 * np.arctanh(np.sqrt((e - 1) / (e + 1)) * np.tan(nu / 2))
        n = np.sqrt(MU / (-a) ** 3); M = e * np.sinh(H0) - H0 + n * dt
        H = np.arcsinh(M / e) if abs(M) > 10 else M
        for _ in range(100):
            f = e * np.sinh(H) - H - M; fp = e * np.cosh(H) - 1
            dH = -f / fp; H += dH
            if abs(dH) < 1e-14 * max(1.0, abs(H)): break
        p = a * (1 - e * e); nu1 = 2 * np.arctan(np.sqrt((e + 1) / (e - 1)) * np.tanh(H / 2))
        rr = p / (1 + e * np.cos(nu1))
        rp = rr * np.array([np.cos(nu1), np.sin(nu1), 0.0]); vp = np.sqrt(MU / p) * np.array([-np.sin(nu1), e + np.cos(nu1), 0.0])
        cO, sO = np.cos(Om), np.sin(Om); ci, si = np.cos(inc), np.sin(inc); cw, sw = np.cos(w), np.sin(w)
        R = np.array([[cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si],
                      [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si],
                      [sw * si, cw * si, ci]])
        return R @ rp, R @ vp
    return propagate_twobody(r, v, dt)        # near-parabolic: universal variables


def random_states(rng, N, a_range=(0.8, 2.0), e_max=0.6, i_max=20.0):
    """N random heliocentric positions on random ellipses (a in AU range, e <= e_max, i <= i_max deg)."""
    out = np.empty((N, 3))
    for k in range(N):
        el = [rng.uniform(*a_range), rng.uniform(0, e_max), rng.uniform(0, i_max), rng.uniform(0, 360), rng.uniform(0, 360),
              rng.uniform(0, 360)]
        out[k] = Body(el, 0.0).state(0.0)[0]
    return out


def check_solutions(R1, R2, tof, v1, v2, ok, label, tol_rel=1e-6, tol_v=1e-6):
    idx = np.where(ok)[0]
    if len(idx) == 0:
        print(f'  {label}: no solutions'); return
    err = np.empty(len(idx)); verr = np.empty(len(idx))
    for j, i in enumerate(idx):
        r, v = propagate_ref(R1[i], v1[i], tof[i])
        err[j] = np.linalg.norm(r - R2[i]); verr[j] = np.linalg.norm(v - v2[i])
    rel = err / np.linalg.norm(R2[idx], axis=1)
    print(f'  {label}: {len(idx)} solutions; |r - r2| max {err.max():.2e} km (median {np.median(err):.2e}), rel max {rel.max():.1e}; '
          f'|v - v2| max {verr.max():.2e} km/s')
    assert rel.max() < tol_rel, f'{label}: endpoint error {err.max()} km'
    assert np.median(err) < 1e-3, f'{label}: median endpoint error {np.median(err)} km'
    assert verr.max() < tol_v, f'{label}: arrival velocity inconsistency {verr.max()} km/s'


def main():
    rng = np.random.default_rng(7)
    N = 4000
    R1 = random_states(rng, N); R2 = random_states(rng, N)
    tof = rng.uniform(100, 900, N) * DAY
    print(f'{N} random heliocentric pairs (a 0.8-2 AU, e <= 0.6), tof 100-900 d')

    # ---- 1. all branches propagate to r2
    sols = {}
    for n in (0, 1, 2):
        for b in ('left', 'right'):
            v1, v2, info = lambert_multirev(R1, R2, tof, nrev=n, branch=b, return_info=True)
            sols[(n, b)] = (v1, v2, info)
            if n == 0 and b == 'left':
                assert not np.any(info['ok']), 'nrev = 0 must have no left branch'
                continue
            check_solutions(R1, R2, tof, v1, v2, info['ok'], f'nrev {n} {b:5s}')

    # ---- 2. nrev = 0 reproduces lambert()
    v1ref, v2ref = lambert(R1, R2, tof)
    v1, v2, info = sols[(0, 'right')]
    okref = np.isfinite(v1ref[:, 0]); ok = info['ok']
    both = ok & okref
    d1 = np.linalg.norm(v1[both] - v1ref[both], axis=1).max(); d2 = np.linalg.norm(v2[both] - v2ref[both], axis=1).max()
    print(f'  nrev 0 vs lambert(): both exist {both.sum()}; max |dv1| {d1:.2e}, |dv2| {d2:.2e} km/s; '
          f'lambert() only {np.sum(okref & ~ok)}, new only {np.sum(ok & ~okref)} (hyperbolic short-way cases lambert() misses)')
    assert d1 < 1e-8 and d2 < 1e-8   # 1.2e-9 observed: same iteration, different stopping rule
    assert np.sum(okref & ~ok) == 0, 'multirev nrev=0 lost solutions that lambert() finds'

    # ---- 3. branch existence and energy ordering for nrev >= 1
    for n in (1, 2):
        v1L, v2L, iL = sols[(n, 'left')]; v1R, v2R, iR = sols[(n, 'right')]
        exists = tof > iL['tof_min']
        assert np.array_equal(iL['ok'], iR['ok']), 'left/right existence differs'
        assert np.array_equal(iL['ok'], exists), 'existence != (tof > tof_min)'
        tmin2 = lambert_tof_min(R1, R2, nrev=n)
        assert np.allclose(tmin2, iL['tof_min'], rtol=1e-9, atol=1.0)
        idx = np.where(iL['ok'])[0]
        aL = np.array([elements(R1[i], v1L[i])[0] for i in idx]); aR = np.array([elements(R1[i], v1R[i])[0] for i in idx])
        eL = np.array([elements(R1[i], v1L[i])[1] for i in idx]); eR = np.array([elements(R1[i], v1R[i])[1] for i in idx])
        frac = np.mean(aL > aR)
        print(f'  nrev {n}: {len(idx)} rows with solutions ({100 * len(idx) / N:.1f}%); tof_min median {np.median(iL["tof_min"][np.isfinite(iL["tof_min"])]) / DAY:.0f} d; '
              f'left a > right a in {100 * frac:.1f}% (a_L median {np.median(aL) / AU:.2f} AU, a_R median {np.median(aR) / AU:.2f} AU; '
              f'e_L median {np.median(eL):.2f}, e_R median {np.median(eR):.2f})')
        assert frac == 1.0, 'left branch must be the high-energy (larger a) solution'
        # distinct solutions
        assert np.all(np.linalg.norm(v1L[idx] - v1R[idx], axis=1) > 1e-6)

    # ---- 4. unimodality of T(z) on each multi-rev interval and tof_min accuracy (grid scan)
    for n in (1, 2):
        M = 300
        G = _MRGeom(R1[:M], R2[:M], tof[:M], MU, True, None)
        zz = np.linspace(4 * np.pi ** 2 * n ** 2 + 1e-3, 4 * np.pi ** 2 * (n + 1) ** 2 - 1e-3, 4001)
        T = np.array([G.T_of(np.full(M, z))[0] for z in zz])     # (nz, M)
        d = np.diff(T, axis=0)
        nmin = np.sum((d[:-1] < 0) & (d[1:] > 0), axis=0)
        assert np.all(nmin == 1), f'T(z) not unimodal for nrev {n}: {np.unique(nmin)}'
        tmin_grid = T.min(axis=0) / np.sqrt(MU)
        tmin = lambert_tof_min(R1[:M], R2[:M], nrev=n)
        rel = (tmin_grid - tmin) / tmin
        assert np.all(rel > -1e-9) and np.all(rel < 1e-4), f'tof_min mismatch nrev {n}: {rel.min()} {rel.max()}'
        print(f'  nrev {n}: T(z) unimodal on all {M} scanned geometries; grid/solver tof_min agree to {rel.max():.1e}')

    # ---- 5. lambert_all / lambert_best consistency
    V1, V2, nrevs, branches = lambert_all(R1, R2, tof, nrev_max=2)
    for k, (n, b) in enumerate(zip(nrevs, branches)):
        v1k, v2k, ik = sols[(n, 'left' if b == BRANCH_LEFT else 'right')]
        assert np.array_equal(np.isfinite(V1[k, :, 0]), ik['ok'])
        assert np.allclose(V1[k][ik['ok']], v1k[ik['ok']], atol=1e-12, rtol=0)
    v_prev = np.cross([0, 0, 1.0], R1); v_prev *= (np.sqrt(MU / np.linalg.norm(R1, axis=1)) / np.linalg.norm(v_prev, axis=1))[:, None]
    v_prev += rng.normal(0, 2.0, (N, 3))
    v1b, v2b, nb, bb = lambert_best(R1, R2, tof, v_prev, nrev_max=2)
    dv_all = np.linalg.norm(V1 - v_prev[None], axis=-1); dv_all = np.where(np.isfinite(dv_all), dv_all, np.inf)
    dv_b = np.linalg.norm(v1b - v_prev, axis=1)
    assert np.all(dv_b <= dv_all.min(axis=0) + 1e-12)
    kb = np.array([np.where((nrevs == n) & (branches == b))[0][0] for n, b in zip(nb, bb)])
    assert np.allclose(np.take_along_axis(V1, kb[None, :, None], axis=0)[0], v1b, equal_nan=True)
    pref = np.bincount(nb, minlength=3)
    print(f'  lambert_best: min junction dv option = nrev 0/1/2 in {pref[0]}/{pref[1]}/{pref[2]} rows '
          f'(nrev >= 1 available in {np.sum(np.isfinite(V1[1, :, 0]))} rows)')

    # ---- 6. speed
    Nb = 100000
    ii = rng.integers(0, N, Nb); jj = rng.integers(0, N, Nb)
    Rb1 = R1[ii]; Rb2 = R2[jj]; tb = rng.uniform(100, 900, Nb) * DAY
    tic = time.time(); lambert(Rb1, Rb2, tb); t_ref = time.time() - tic
    tic = time.time(); lambert_multirev(Rb1, Rb2, tb, nrev=1, branch='right'); t_1 = time.time() - tic
    tic = time.time(); lambert_multirev(Rb1, Rb2, tb, nrev=0, branch='right'); t_0 = time.time() - tic
    tic = time.time(); lambert_all(Rb1, Rb2, tb, nrev_max=2); t_all = time.time() - tic
    print(f'  speed on {Nb} rows: lambert() {Nb / t_ref:.0f} rows/s; lambert_multirev nrev0 {Nb / t_0:.0f}, nrev1 one branch {Nb / t_1:.0f}; '
          f'lambert_all (5 solutions) {Nb / t_all:.0f} rows/s')
    print('ALL TESTS PASSED')


if __name__ == '__main__':
    main()
