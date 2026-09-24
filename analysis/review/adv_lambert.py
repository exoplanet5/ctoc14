"""Adversarial probe of ctoc14.lambert.lambert (independent of the earlier lambert_stress.py).

Reference: an independent Lambert solver based on a robust bisection of the same universal-variable time equation over a
wide z bracket (z in (-1e4, 4 pi^2)), plus verification of every returned v1 by an accurate two-body propagation
(analytic f/g through ctoc14.propagate_twobody cross-checked with solve_ivp when the two disagree).
Probes: transfer angle near 0 / 180 deg, long way, tof 10 d and 1500 d, retrograde flag, hyperbolic transfers, and a
random sweep over asteroid-to-asteroid legs to count (a) NaN where a solution exists, (b) finite solutions that are WRONG.
"""
import sys, pathlib, numpy as np, warnings
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from scipy.integrate import solve_ivp
from ctoc14.lambert import lambert
from ctoc14.kepler import Ephemeris, propagate_twobody
from ctoc14.constants import MU, AU, DAY

np.set_printoptions(precision=6, linewidth=160)
FAILS = []
def check(name, cond, detail=''):
    print(('ok   ' if cond else 'FAIL ') + name + (': ' + detail if detail else ''))
    if not cond: FAILS.append(name)

def stumpff(z):
    if z > 1e-6: s = np.sqrt(z); return (1 - np.cos(s)) / z, (s - np.sin(s)) / s ** 3
    if z < -1e-6: s = np.sqrt(-z); return (np.cosh(s) - 1) / (-z), (np.sinh(s) - s) / s ** 3
    return 0.5 - z / 24 + z * z / 720, 1 / 6 - z / 120 + z * z / 5040

def lambert_ref(r1, r2, tof, long_way=False, z_min=-1e4):
    """Independent single-revolution Lambert: bisection on z of Curtis' F(z) with a wide bracket. Returns (v1, v2, z) or None."""
    r1n = np.linalg.norm(r1); r2n = np.linalg.norm(r2)
    cosd = np.clip(r1 @ r2 / (r1n * r2n), -1, 1); dnu = np.arccos(cosd)
    if long_way: dnu = 2 * np.pi - dnu
    if abs(np.sin(dnu)) < 1e-14: return None
    A = np.sin(dnu) * np.sqrt(r1n * r2n / (1 - np.cos(dnu)))
    def yF(z):
        C, S = stumpff(z); y = r1n + r2n + A * (z * S - 1) / np.sqrt(C)
        if y <= 0: return y, None
        return y, (y / C) ** 1.5 * S + A * np.sqrt(y) - np.sqrt(MU) * tof
    lo, hi = z_min, 4 * np.pi ** 2 - 1e-9
    # raise lo until y > 0 (y increases with z when A > 0)
    y, F = yF(lo)
    while F is None:
        lo = 0.5 * (lo + hi) if lo < -1e-3 else lo + 0.5 * (hi - lo)
        y, F = yF(lo)
        if hi - lo < 1e-12: return None
    if F > 0: return None            # tof shorter than the minimum reachable with this bracket
    for _ in range(300):
        mid = 0.5 * (lo + hi); y, F = yF(mid)
        if F is None or F < 0: lo = mid
        else: hi = mid
        if hi - lo < 1e-13 * max(1, abs(mid)): break
    z = 0.5 * (lo + hi); y, F = yF(z)
    if F is None: return None
    f = 1 - y / r1n; g = A * np.sqrt(y / MU); gd = 1 - y / r2n
    return (r2 - f * r1) / g, (gd * r2 - r1) / g, z

def endpoint_err(r1, v1, tof, r2):
    rp, _ = propagate_twobody(r1, v1, tof)
    return np.linalg.norm(rp - r2)

def ivp_err(r1, v1, tof, r2):
    f = lambda t, y: np.concatenate([y[3:], -MU * y[:3] / np.linalg.norm(y[:3]) ** 3])
    s = solve_ivp(f, (0, tof), np.concatenate([r1, v1]), method='DOP853', rtol=1e-12, atol=1e-6)
    return np.linalg.norm(s.y[:3, -1] - r2)

eph = Ephemeris()
# ------------------------------------------------------------------ 1. transfer-angle extremes on a 1 AU -> 1.3 AU geometry
print('=== 1. transfer angle extremes (r1 = 1 AU on x, r2 = 1.3 AU rotated by dnu about z, tof 200 d)')
r1 = np.array([AU, 0.0, 0.0]); tof = 200 * DAY
for dnu_deg in (1e-6, 1e-4, 1e-2, 0.1, 1.0, 10.0, 90.0, 170.0, 179.0, 179.9, 179.99, 179.999, 180.0, 180.001, 180.1, 181.0, 190.0, 270.0, 350.0, 359.9, 359.999):
    d = np.radians(dnu_deg); r2 = 1.3 * AU * np.array([np.cos(d), np.sin(d), 0.0])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore'); v1, v2 = lambert(r1, r2, tof)
    fin = np.all(np.isfinite(v1))
    err = endpoint_err(r1, v1, tof, r2) if fin else np.nan
    ref = lambert_ref(r1, r2, tof, long_way=(dnu_deg > 180))
    ref_err = endpoint_err(r1, ref[0], tof, r2) if ref else np.nan
    hz = np.cross(r1, v1)[2] if fin else np.nan
    tag = 'ok' if (fin and err < 1e-3) or (not fin and ref is None) else 'FAIL'
    if tag == 'FAIL': FAILS.append(f'dnu={dnu_deg}')
    print(f'  {tag:4} dnu={dnu_deg:9.4f} deg: finite={fin!s:5} |v1|={np.linalg.norm(v1) if fin else np.nan:8.3f} km/s hz={hz:+.2e} endpoint err={err:.3e} km | ref: {"z=%.4f err=%.2e" % (ref[2], ref_err) if ref else "none"}')

# ------------------------------------------------------------------ 2. tof extremes (10 d, 1500 d) on the same geometry at 60 deg
print('=== 2. tof extremes, dnu = 60 deg, 1 AU -> 1.3 AU')
d = np.radians(60.0); r2 = 1.3 * AU * np.array([np.cos(d), np.sin(d), 0.0])
for tofd in (1.0, 5.0, 10.0, 30.0, 60.0, 80.0, 100.0, 150.0, 300.0, 700.0, 1000.0, 1500.0, 3000.0, 5478.0):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore'); v1, v2 = lambert(r1, r2, tofd * DAY)
    fin = np.all(np.isfinite(v1)); err = endpoint_err(r1, v1, tofd * DAY, r2) if fin else np.nan
    ref = lambert_ref(r1, r2, tofd * DAY)
    ref_err = endpoint_err(r1, ref[0], tofd * DAY, r2) if ref else np.nan
    a = None
    if fin:
        a = 1 / (2 / AU - v1 @ v1 / MU)
    tag = 'ok' if (fin and err < 1e-3) or (not fin and ref is None) else 'FAIL'
    if tag == 'FAIL': FAILS.append(f'tof={tofd}')
    print(f'  {tag:4} tof={tofd:7.1f} d: finite={fin!s:5} |v1|={np.linalg.norm(v1) if fin else np.nan:8.3f} a={(a/AU if a else np.nan):+9.3f} AU err={err:.3e} km | ref: {"z=%.5f |v1|=%.3f err=%.2e" % (ref[2], np.linalg.norm(ref[0]), ref_err) if ref else "none"}')

# ------------------------------------------------------------------ 3. retrograde / long_way flags
print('=== 3. prograde / long_way flag semantics (dnu = 60 deg short way, z-component of r1 x r2 > 0)')
tof = 200 * DAY
for pro, lw in ((True, None), (False, None), (True, True), (True, False), (False, True), (False, False)):
    v1, v2 = lambert(r1, r2, tof, prograde=pro, long_way=lw)
    hz = np.cross(r1, v1)[2]; err = endpoint_err(r1, v1, tof, r2)
    print(f'  prograde={pro!s:5} long_way={lw!s:5}: hz={hz:+.3e} (retrograde={hz < 0}) err={err:.2e} km')
v1p, _ = lambert(r1, r2, tof, prograde=True); v1r, _ = lambert(r1, r2, tof, prograde=False)
check('prograde=False gives negative h_z and a valid transfer', np.cross(r1, v1r)[2] < 0 and endpoint_err(r1, v1r, tof, r2) < 1e-3, '')
check('prograde=True gives positive h_z', np.cross(r1, v1p)[2] > 0, '')
# inclined geometry: r2 below the ecliptic, prograde flag decided by the z-component only
r2i = 1.3 * AU * np.array([np.cos(d) * np.cos(0.6), np.sin(d) * np.cos(0.6), -np.sin(0.6)])
v1i, _ = lambert(r1, r2i, tof, prograde=True)
check('inclined geometry prograde solution valid', endpoint_err(r1, v1i, tof, r2i) < 1e-3 and np.cross(r1, v1i)[2] > 0, f'err={endpoint_err(r1, v1i, tof, r2i):.2e}')

# ------------------------------------------------------------------ 4. random asteroid-to-asteroid sweep: wrong-but-finite and spurious-NaN counts
print('=== 4. random asteroid->asteroid legs')
rng = np.random.default_rng(7)
for tofd in (10.0, 60.0, 200.0, 500.0, 1000.0, 1500.0):
    N = 3000
    t0 = rng.uniform(0, 4e8, N); k1 = rng.integers(0, 300, N); k2 = rng.integers(0, 300, N)
    R1 = np.array([eph.ast_state(a, t)[0] for a, t in zip(k1, t0)]); R2 = np.array([eph.ast_state(b, t + tofd * DAY)[0] for b, t in zip(k2, t0)])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore'); V1, V2 = lambert(R1, R2, np.full(N, tofd * DAY))
    fin = np.isfinite(V1[:, 0])
    errs = np.full(N, np.nan)
    for i in np.where(fin)[0]:
        errs[i] = endpoint_err(R1[i], V1[i], tofd * DAY, R2[i])
    wrong = np.where(fin & ~(errs < 1.0))[0]
    # among wrong ones, is it the Lambert v1 or the propagate check that is wrong? re-check the worst 5 with solve_ivp
    wrong_confirmed = []
    for i in wrong[:5]:
        e_ivp = ivp_err(R1[i], V1[i], tofd * DAY, R2[i])
        q = None
        h = np.cross(R1[i], V1[i]); en = 0.5 * V1[i] @ V1[i] - MU / np.linalg.norm(R1[i]); a = -MU / (2 * en)
        ecc = np.sqrt(max(0.0, 1 - (h @ h) / (MU * a))) if a > 0 else np.nan; q = a * (1 - ecc) if a > 0 else np.nan
        wrong_confirmed.append((i, errs[i], e_ivp, a / AU, ecc, q / AU if q == q else np.nan))
    # spurious NaN: reference finds a solution (any z) -> classify by z sign
    nan_idx = np.where(~fin)[0][:200]
    ref_found = [];
    for i in nan_idx:
        lw = np.cross(R1[i], R2[i])[2] < 0
        ref = lambert_ref(R1[i], R2[i], tofd * DAY, long_way=lw)
        if ref is not None and endpoint_err(R1[i], ref[0], tofd * DAY, R2[i]) < 1.0:
            ref_found.append((i, ref[2], np.linalg.norm(ref[0])))
    zs = np.array([z for _, z, _ in ref_found]) if ref_found else np.array([])
    print(f'  tof={tofd:6.0f} d: finite {fin.sum()}/{N}; finite-but-wrong (endpoint err >= 1 km): {len(wrong)}; '
          f'NaN with a reference solution: {len(ref_found)}/{len(nan_idx)} checked (z range {zs.min() if len(zs) else np.nan:.2f}..{zs.max() if len(zs) else np.nan:.2f}; {int((zs > 0).sum())} elliptic)')
    for (i, e1, e2, a, ecc, q) in wrong_confirmed:
        print(f'      wrong #{i}: err(f/g)={e1:.2e} km err(solve_ivp)={e2:.2e} km |v1|={np.linalg.norm(V1[i]):.2f} km/s  orbit a={a:+.3f} AU e={ecc:.6f} q={q:.2e} AU')
    if len(wrong): FAILS.append(f'finite-but-wrong tof={tofd}')

# ------------------------------------------------------------------ 5. hyperbolic case explicit
print('=== 5. hyperbolic transfer (1 AU -> 2 AU, 90 deg, tof 60 d): F(0) > 0 so the root is at z < 0')
r2h = 2 * AU * np.array([0.0, 1.0, 0.0])
with warnings.catch_warnings():
    warnings.simplefilter('ignore'); v1h, _ = lambert(r1, r2h, 60 * DAY)
ref = lambert_ref(r1, r2h, 60 * DAY)
check('hyperbolic single-rev transfer found', np.all(np.isfinite(v1h)), f'lambert v1={v1h}; reference z={ref[2] if ref else None:.4f} |v1|={np.linalg.norm(ref[0]) if ref else np.nan:.3f} km/s err={endpoint_err(r1, ref[0], 60*DAY, r2h) if ref else np.nan:.2e} km')

# ------------------------------------------------------------------ 6. accepted-solution tolerance: how large an endpoint error can pass the `ok` test?
print('=== 6. ok-test slack: |F| < 1e-6 sqrt(mu) tof  <=>  time-of-flight error up to 1e-6 tof')
print(f'  at tof = 900 d that is {1e-6 * 900 * DAY:.1f} s of tof error, i.e. up to ~{1e-6 * 900 * DAY * 30:.0f} km endpoint error if the iteration stops at maxit without converging')
print('\nFAILURES:', FAILS if FAILS else 'none')
