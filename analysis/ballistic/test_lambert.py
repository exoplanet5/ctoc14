"""Validation of ephem.py and lambert.py.

Checks: (1) Earth/asteroid ephemeris sanity and the PDF sample lines; (2) Lambert (N=0) consistency:
orbit elements from (r1,v1) and (r2,v2) agree, and the Keplerian time between the two points equals tof;
(3) for elliptic arcs propagate (r1,v1) with the universal-variable propagator and hit r2;
(4) multi-revolution solutions (N=1..3, both branches): Kepler time between the points equals
tof - N*period, propagation hits r2, and every solution that exists (T >= T_min) converges.
"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ephem import (MU_SUN, AU, DAY, earth_state, asteroid_state, load_mea, propagate_universal, propagate_kepler,
                   elements_from_state, T_MAX_DAYS)
from lambert import lambert, lambert_multirev

rng = np.random.default_rng(1)


def kepler_tof(r1, v1, r2, v2, mu=MU_SUN):
    """Time from r1 to r2 along the conic defined by (r1,v1), less than one revolution, along motion."""
    a, e, inc, q, Q = elements_from_state(r1, v1)
    h = np.cross(r1, v1)
    hn = np.linalg.norm(h, axis=1)
    R1 = np.linalg.norm(r1, axis=1); R2 = np.linalg.norm(r2, axis=1)
    p = hn ** 2 / mu

    def nu(r, v, R):
        c = np.clip((p / R - 1) / e, -1, 1)
        n_ = np.arccos(c)
        rv = np.sum(r * v, axis=1)
        return np.where(rv < 0, 2 * np.pi - n_, n_)
    nu1 = nu(r1, v1, R1); nu2 = nu(r2, v2, R2)
    ell = e < 1
    out = np.full(len(e), np.nan)
    ee = e[ell]; aa = a[ell]
    E1 = 2 * np.arctan2(np.sqrt(1 - ee) * np.sin(nu1[ell] / 2), np.sqrt(1 + ee) * np.cos(nu1[ell] / 2))
    E2 = 2 * np.arctan2(np.sqrt(1 - ee) * np.sin(nu2[ell] / 2), np.sqrt(1 + ee) * np.cos(nu2[ell] / 2))
    M1 = E1 - ee * np.sin(E1); M2 = E2 - ee * np.sin(E2)
    dM = np.mod(M2 - M1, 2 * np.pi)
    out[ell] = dM / np.sqrt(mu / aa ** 3)
    hy = ~ell
    eh = e[hy]; ah = -a[hy]
    H1 = 2 * np.arctanh(np.clip(np.sqrt((eh - 1) / (eh + 1)) * np.tan(nu1[hy] / 2), -1 + 1e-15, 1 - 1e-15))
    H2 = 2 * np.arctanh(np.clip(np.sqrt((eh - 1) / (eh + 1)) * np.tan(nu2[hy] / 2), -1 + 1e-15, 1 - 1e-15))
    N1 = eh * np.sinh(H1) - H1; N2 = eh * np.sinh(H2) - H2
    out[hy] = (N2 - N1) / np.sqrt(mu / ah ** 3)
    return out


def period(r, v, mu=MU_SUN):
    a = elements_from_state(r, v)[0]
    return 2 * np.pi * np.sqrt(np.abs(a) ** 3 / mu)


# --- Earth ephemeris sanity + PDF sample lines (problem statement sec. 6.2)
t = np.arange(0, 5479, 1.0)
rE, vE = earth_state(t)
R = np.linalg.norm(rE, axis=-1) / AU
V = np.linalg.norm(vE, axis=-1)
print(f"Earth |r| AU: min {R.min():.5f} max {R.max():.5f}; |v| km/s: min {V.min():.4f} max {V.max():.4f}")
n = np.sqrt(MU_SUN / (1.0009175020 * AU) ** 3)
print(f"Earth period {2*np.pi/n/DAY:.4f} d")
r_s1 = np.array([-1.9500328647e+07, 1.4581848347e+08, -7.6372048832e+03]); v_s1 = np.array([-2.9389477847e+01, -4.3166912970e+00, -3.9417153392e+00])
r0, v0 = earth_state(0.0)
print(f"PDF sample line 1: Earth r(t0) error {np.linalg.norm(r0 - r_s1)*1e3:.2f} m; sample launch v_inf {np.linalg.norm(v_s1 - v0):.6f} km/s")
t2 = 1.0008479152e+07
r_s2 = np.array([-1.1398611209e+08, -8.7213961221e+07, -1.6523817305e+07])
ids, el = load_mea()
ra2, _ = asteroid_state(el[173], t2 / DAY)
print(f"PDF sample line 2: asteroid 174 at t={t2/DAY:.2f} d is {np.linalg.norm(ra2 - r_s2)*1e3:.1f} m from the sample flyby position")
r1, v1 = propagate_universal(r0[None], v0[None], np.array([1000 * DAY]))
r1a, v1a = earth_state(1000.0)
print(f"propagate_universal vs ephemeris @1000d: dr {np.linalg.norm(r1[0]-r1a):.3e} km, dv {np.linalg.norm(v1[0]-v1a):.3e} km/s")
ra, va = asteroid_state(el, np.array([0.0, 2000.0]))
a_, e_, i_, q_, Q_ = elements_from_state(ra[:, 0], va[:, 0])
print(f"asteroid element round-trip: max |a-a0|/a0 {np.max(np.abs(a_/AU-el[:,0])/el[:,0]):.2e}, max |e-e0| {np.max(np.abs(e_-el[:,1])):.2e}, max |i-i0| deg {np.max(np.abs(np.degrees(i_)-el[:,2])):.2e}")

# --- Lambert N=0 validation on random Earth->asteroid geometries
N = 300000
tL = rng.uniform(0, T_MAX_DAYS, N)
tof = np.exp(rng.uniform(np.log(5), np.log(3000), N))
k = rng.integers(0, 300, N)
r1, v1e = earth_state(tL)
r2 = np.empty((N, 3))
for j in range(300):
    m = k == j
    if m.any():
        r2[m] = asteroid_state(el[j], tL[m] + tof[m])[0]
allok = True
for retro in (False, True):
    t0_ = time.time()
    v1, v2, ex = lambert(r1, r2, tof * DAY, MU_SUN, retrograde=retro, return_extra=True)
    dt_ = time.time() - t0_
    ok = np.isfinite(v1[:, 0])
    a1, e1, i1, q1, Q1 = elements_from_state(r1[ok], v1[ok])
    a2, e2, i2, q2, Q2 = elements_from_state(r2[ok], v2[ok])
    h1 = np.cross(r1[ok], v1[ok]); h2 = np.cross(r2[ok], v2[ok])
    dh = np.linalg.norm(h1 - h2, axis=1) / np.linalg.norm(h1, axis=1)
    de = np.abs(e1 - e2)
    tk = kepler_tof(r1[ok], v1[ok], r2[ok], v2[ok])
    dtof = np.abs(tk - tof[ok] * DAY) / (tof[ok] * DAY)
    prog = h1[:, 2] > 0
    print(f"N=0 retrograde={retro}: {dt_:.2f}s for {N} solves ({dt_/N*1e6:.2f} us/solve); converged {ok.mean()*100:.4f}% ({(~ok).sum()} fail); "
          f"iters mean {ex['iters'].mean():.2f} max {ex['iters'].max()}")
    print(f"   h-vector consistency rel err: max {dh.max():.2e}; |e1-e2| max {de.max():.2e}; "
          f"Kepler TOF rel err: median {np.median(dtof):.1e}, 99.99% {np.percentile(dtof,99.99):.1e}, max {np.nanmax(dtof):.1e}; "
          f"h_z>0 fraction {prog.mean()*100:.2f}% (expect {'100' if not retro else '0'}%)")
    sel = np.abs(e1 - 1) > 1e-6               # element propagator is exact for all conics except the parabola
    rp, vp = propagate_kepler(r1[ok][sel], v1[ok][sel], tof[ok][sel] * DAY)
    err_r = np.linalg.norm(rp - r2[ok][sel], axis=1)
    err_v = np.linalg.norm(vp - v2[ok][sel], axis=1)
    print(f"   propagation check (element propagator, {sel.sum()} arcs): pos err median {np.median(err_r):.1e} km, 99.9% {np.percentile(err_r,99.9):.1e} km, max {np.nanmax(err_r):.1e} km; vel err max {np.nanmax(err_v):.1e} km/s")
    print(f"   arcs: hyperbolic fraction {(e1>=1).mean()*100:.2f}%, |x-1|<0.01 fraction {(np.abs(ex['x'][ok]-1)<0.01).mean()*100:.3f}%; iters>=10: {(ex['iters']>=10).sum()}")
    if (~ok).any():
        print("   non-converged sample (T, lam):", ex['T'][~ok][:5], ex['lam'][~ok][:5])
        allok = False

# --- Multi-revolution validation
print("\n--- multi-revolution (prograde) ---")
Nm = 200000
tLm = rng.uniform(0, T_MAX_DAYS, Nm)
tofm = rng.uniform(200, 3000, Nm)
km = rng.integers(0, 300, Nm)
r1m, _ = earth_state(tLm)
r2m = np.empty((Nm, 3))
for j in range(300):
    m = km == j
    if m.any():
        r2m[m] = asteroid_state(el[j], tLm[m] + tofm[m])[0]
for Nrev in (1, 2, 3):
    t0_ = time.time()
    v1, v2, ex = lambert_multirev(r1m, r2m, tofm * DAY, MU_SUN, N=Nrev, return_extra=True)
    dt_ = time.time() - t0_
    exists = ex['T'] >= ex['TM']
    # T_min check by brute force on a subset: T(x) on a fine x grid must be >= TM
    sub = rng.choice(Nm, 2000, replace=False)
    xg = np.linspace(-0.999, 0.999, 4001)
    from lambert import _x2tof
    Tg = _x2tof(xg[None, :], ex['lam'][sub][:, None], Nrev)
    tm_err = (Tg.min(axis=1) - ex['TM'][sub]) / ex['TM'][sub]
    print(f"N={Nrev}: {dt_:.2f}s for {Nm} cells ({dt_/Nm*1e6:.2f} us/cell, both branches); solutions exist in {exists.mean()*100:.1f}% of cells; "
          f"T_min check: min(T_grid)-T_M rel = [{tm_err.min():.1e}, {tm_err.max():.1e}] (must be >= ~0)")
    for ib, name in enumerate(('left', 'right')):
        ok = np.isfinite(v1[:, ib, 0])
        conv_of_exist = ok[exists].mean()
        v1b, v2b = v1[ok, ib], v2[ok, ib]
        a1, e1, i1, q1, Q1 = elements_from_state(r1m[ok], v1b)
        h1 = np.cross(r1m[ok], v1b); h2 = np.cross(r2m[ok], v2b)
        dh = np.linalg.norm(h1 - h2, axis=1) / np.linalg.norm(h1, axis=1)
        P = period(r1m[ok], v1b)
        tk = kepler_tof(r1m[ok], v1b, r2m[ok], v2b)
        nrev_meas = (tofm[ok] * DAY - tk) / P
        rp, vp = propagate_kepler(r1m[ok], v1b, tofm[ok] * DAY)
        err_r = np.linalg.norm(rp - r2m[ok], axis=1)
        side_ok = (ex['x'][ok, ib] <= ex['xM'][ok]) if ib == 0 else (ex['x'][ok, ib] >= ex['xM'][ok])
        print(f"   {name:5s}: converged {conv_of_exist*100:.4f}% of existing ({(exists & ~ok).sum()} fail), iters mean {ex['iters'][exists, ib].mean():.1f} max {ex['iters'][:, ib].max()}; "
              f"all elliptic: {(e1 < 1).all()}; h consistency max {dh.max():.1e}; (tof - t_kepler)/P = {Nrev} to within max |dev| {np.abs(nrev_meas - Nrev).max():.1e}; "
              f"propagation pos err median {np.median(err_r):.1e} km, 99.9% {np.percentile(err_r, 99.9):.1e} km, max {err_r.max():.1e} km; "
              f"x on correct side of x_M: {side_ok.mean()*100:.3f}%; a median {np.median(a1)/AU:.3f} AU, h_z>0: {(h1[:,2]>0).mean()*100:.1f}%")
        if (exists & ~ok).any():
            allok = False
            bad = np.where(exists & ~ok)[0][:5]
            print("      failures (T, TM, lam, x):", ex['T'][bad], ex['TM'][bad], ex['lam'][bad], ex['x'][bad, ib])

# --- stress: near-parabolic and short-TOF cases (N=0)
r1s = np.array([[AU, 0, 0]] * 5)
r2s = np.array([[0, AU, 0], [-AU, 0.1 * AU, 0], [0.5 * AU, 0.5 * AU, 0.1 * AU], [2 * AU, 0.3 * AU, 0.0], [AU, 0.01 * AU, 0]])
for tofd in [5, 30, 60, 100, 200, 400, 800, 3000]:
    v1, v2, ex = lambert(r1s, r2s, np.full(5, tofd * DAY), MU_SUN, return_extra=True)
    rp, vp = propagate_kepler(r1s, v1, np.full(5, tofd * DAY))
    err = np.linalg.norm(rp - r2s, axis=1)
    print(f"tof {tofd:5d} d: x = {np.array2string(ex['x'], precision=4)}  pos err km = {np.array2string(err, precision=2)}  iters {ex['iters']}")
# multi-rev stress: 1 AU circular-ish geometry, TOF from 300 to 1500 d, N=1..3
for tofd in [300, 400, 600, 800, 1200, 1500]:
    line = f"multirev tof {tofd:5d} d: "
    for Nrev in (1, 2, 3):
        v1, v2, ex = lambert_multirev(r1s, r2s, np.full(5, tofd * DAY), MU_SUN, N=Nrev, return_extra=True)
        rp, vp = propagate_kepler(np.repeat(r1s, 2, axis=0), v1.reshape(-1, 3), np.full(10, tofd * DAY))
        err = np.linalg.norm(rp - np.repeat(r2s, 2, axis=0), axis=1).reshape(5, 2)
        nsol = np.isfinite(v1[:, :, 0]).sum()
        line += f"N={Nrev}: {nsol}/10 sols, max pos err {np.nanmax(err) if nsol else 0:.1e} km | "
    print(line)
print("ALL CONVERGED" if allok else "SOME FAILURES")
