"""Adversarial tests of ctoc14.lambert: transfer angles near 0 / 180 deg, long way, tof extremes, retrograde flag,
hyperbolic transfers. Every returned solution is verified by propagating (r1, v1) for tof and comparing with r2.
An independent bisection Lambert (same universal-variable formulation but a full bracket search) is used to decide whether
NaN results are genuine 'no single-rev solution' cases or solver failures."""
import sys, pathlib, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from ctoc14.lambert import lambert, _stumpff
from ctoc14.kepler import Ephemeris, propagate_twobody
from ctoc14.constants import MU, AU, DAY

fails = []
def check(name, cond, detail=''):
    (print(f'ok   {name}: {detail}') if cond else (fails.append(name), print(f'FAIL {name}: {detail}')))

def ref_lambert(r1, r2, tof, long_way):
    """Independent scalar single-rev Lambert: bisection on z over the full admissible range (y>0), Curtis formulation."""
    r1 = np.asarray(r1, float); r2 = np.asarray(r2, float)
    r1n = np.linalg.norm(r1); r2n = np.linalg.norm(r2)
    cosd = np.clip(r1 @ r2 / (r1n * r2n), -1, 1); dnu = np.arccos(cosd)
    if long_way: dnu = 2 * np.pi - dnu
    A = np.sin(dnu) * np.sqrt(r1n * r2n / (1 - cosd))
    if abs(A) < 1e-12: return None
    def y(z):
        C, S = _stumpff(np.array([z])); return float(r1n + r2n + A * (z * S[0] - 1) / np.sqrt(C[0]))
    def F(z):
        C, S = _stumpff(np.array([z])); yy = y(z)
        if yy <= 0: return np.nan
        return float((yy / C[0]) ** 1.5 * S[0] + A * np.sqrt(yy) - np.sqrt(MU) * tof)
    # find z_lo with y>0 and F<0 by scanning downward
    zhi = 4 * np.pi ** 2 - 1e-9
    zlo = zhi
    step = 1.0
    while True:
        zlo -= step; step *= 1.5
        f = F(zlo)
        if np.isnan(f):          # y<=0 region reached without a sign change
            return None
        if f < 0: break
        if zlo < -1e4: return None
    for _ in range(300):
        zm = 0.5 * (zlo + zhi); fm = F(zm)
        if np.isnan(fm) or fm > 0: zhi = zm
        else: zlo = zm
    z = 0.5 * (zlo + zhi); C, S = _stumpff(np.array([z])); yy = y(z)
    f = 1 - yy / r1n; g = A * np.sqrt(yy / MU); gd = 1 - yy / r2n
    return (r2 - f * r1) / g, (gd * r2 - r1) / g, z

def endpoint_err(r1, v1, tof, r2):
    rr, _ = propagate_twobody(r1, v1, tof)
    return np.linalg.norm(rr - r2)

eph = Ephemeris()
rE, vE = eph.earth_state(0.0)

# ---------------------------------------------------------------- 1. PDF example (scalar inputs)
t2 = 1.0008479152e+07; ra, _ = eph.ast_state(173, t2)
v1, v2 = lambert(rE, ra, t2)
check('PDF example scalar', np.linalg.norm(v1 - np.array([-29.389477847, -4.3166912970, -3.9417153392])) < 1e-6, f'v1={v1}')

# ---------------------------------------------------------------- 2. transfer angles near 0 and 180 deg, coplanar circular 1 AU -> 1.5 AU
def state_on_circle(R, ang, inc=0.0):
    return R * np.array([np.cos(ang), np.sin(ang) * np.cos(inc), np.sin(ang) * np.sin(inc)])
print('\n--- transfer angle sweep (r1 = 1 AU, r2 = 1.5 AU, tof = 200 d) ---')
for deg in [1e-6, 1e-4, 1e-2, 0.1, 1, 10, 90, 170, 179, 179.9, 179.99, 179.999, 179.9999, 180.0, 180.0001, 180.01, 181, 190, 270, 350, 359.9, 359.999]:
    ang = np.radians(deg)
    r1 = state_on_circle(AU, 0.0); r2 = state_on_circle(1.5 * AU, ang)
    tof = 200 * DAY
    v1, v2 = lambert(r1, r2, tof)
    ok = np.all(np.isfinite(v1))
    err = endpoint_err(r1, v1, tof, r2) if ok else np.nan
    v1b, v2b = lambert(r1, r2, tof, long_way=True)
    errb = endpoint_err(r1, v1b, tof, r2) if np.all(np.isfinite(v1b)) else np.nan
    hz = np.cross(r1, v1)[2] if ok else np.nan
    ref = ref_lambert(r1, r2, tof, long_way=(deg > 180))
    referr = endpoint_err(r1, ref[0], tof, r2) if ref is not None else np.nan
    print(f'dnu={deg:10.5f}: auto-prograde err={err:10.3e} km (hz sign {np.sign(hz):+.0f})   long_way=True err={errb:10.3e} km   ref err={referr:10.3e}')
    if 0.01 <= deg <= 179.9 or 180.1 <= deg <= 359.99:
        check(f'lambert dnu={deg}', ok and err < 1.0, f'err={err:.3e} km')

# ---------------------------------------------------------------- 3. tof extremes: 10 d and 1500 d, random asteroid pairs
print('\n--- tof extremes ---')
rng = np.random.default_rng(3)
for tof_d in [1.0, 10.0, 30.0, 100.0, 365.0, 800.0, 1500.0, 3000.0]:
    tof = tof_d * DAY
    N = 3000
    t0 = rng.uniform(0, 4e8, N); k1 = rng.integers(0, 300, N); k2 = rng.integers(0, 300, N)
    R1 = np.array([eph.ast_state(k, t)[0] for k, t in zip(k1, t0)]); R2 = np.array([eph.ast_state(k, t + tof)[0] for k, t in zip(k2, t0)])
    V1, V2 = lambert(R1, R2, np.full(N, tof))
    ok = np.isfinite(V1[:, 0])
    errs = np.array([endpoint_err(R1[i], V1[i], tof, R2[i]) for i in np.where(ok)[0]])
    # check NaN cases against the reference
    nan_idx = np.where(~ok)[0]
    ref_solvable = 0; ref_err = []; ref_z = []
    for i in nan_idx[:300]:
        cr = np.cross(R1[i], R2[i])
        ref = ref_lambert(R1[i], R2[i], tof, long_way=(cr[2] < 0))
        if ref is not None:
            e = endpoint_err(R1[i], ref[0], tof, R2[i])
            if e < 1.0:
                ref_solvable += 1; ref_err.append(e); ref_z.append(ref[2])
    print(f'tof={tof_d:7.1f} d: solved {ok.sum()}/{N}; endpoint err max={errs.max() if len(errs) else np.nan:.3e} km, '
          f'99%={np.percentile(errs, 99) if len(errs) else np.nan:.3e}; of first {min(300, len(nan_idx))} NaN cases the reference solves {ref_solvable} '
          f'(z range {min(ref_z) if ref_z else np.nan:.2f}..{max(ref_z) if ref_z else np.nan:.2f})')
    check(f'tof={tof_d} d endpoint accuracy', len(errs) and errs.max() < 1.0, f'max err {errs.max():.3e} km')
    check(f'tof={tof_d} d no spurious NaN', ref_solvable == 0, f'{ref_solvable} NaN results have a valid single-rev solution (hyperbolic z<0 or otherwise)')

# ---------------------------------------------------------------- 4. explicit hyperbolic transfer (z < 0): 1 AU -> 2 AU, 90 deg, 40 d
print('\n--- hyperbolic transfer ---')
r1 = state_on_circle(AU, 0.0); r2 = state_on_circle(2 * AU, np.pi / 2); tof = 40 * DAY
v1, v2 = lambert(r1, r2, tof)
ref = ref_lambert(r1, r2, tof, long_way=False)
print('lambert v1 =', v1, ' ref v1 =', ref[0] if ref else None, ' ref z =', ref[2] if ref else None, ' ref err =', endpoint_err(r1, ref[0], tof, r2) if ref else None)
check('hyperbolic transfer solved', np.all(np.isfinite(v1)) and endpoint_err(r1, v1, tof, r2) < 1.0, f'v1={v1}')
# mildly hyperbolic (z slightly negative): scan tof downward from the parabolic tof
for tof_d in [120, 100, 90, 85, 80, 75, 70, 60, 50]:
    tof = tof_d * DAY
    v1, _ = lambert(r1, r2, tof); ref = ref_lambert(r1, r2, tof, False)
    print(f'  tof={tof_d} d: lambert finite={np.all(np.isfinite(v1))}, ref z={ref[2] if ref else None}')

# ---------------------------------------------------------------- 5. retrograde flag and long_way flag semantics
print('\n--- retrograde / long-way flags ---')
r1 = state_on_circle(AU, 0.0); r2 = state_on_circle(1.2 * AU, np.radians(60), inc=np.radians(10)); tof = 150 * DAY
for pro, lw in [(True, None), (False, None), (True, True), (True, False), (False, True), (False, False)]:
    v1, v2 = lambert(r1, r2, tof, prograde=pro, long_way=lw)
    ok = np.all(np.isfinite(v1)); err = endpoint_err(r1, v1, tof, r2) if ok else np.nan
    hz = np.cross(r1, v1)[2] if ok else np.nan
    print(f'prograde={pro} long_way={lw}: finite={ok} err={err:.3e} km hz={hz:+.3e}')
    if ok:
        if lw is None:
            check(f'prograde={pro} hz sign', (hz >= 0) == pro, f'hz={hz:+.3e}')
        check(f'prograde={pro} long_way={lw} endpoint', err < 1.0, f'{err:.3e}')

# ---------------------------------------------------------------- 6. broadcasting shapes
V1, V2 = lambert(r1[None, None, :], np.stack([r2, r2 * 1.1])[None, :, :], np.array([[150 * DAY, 200 * DAY]]))
check('broadcast shape', V1.shape == (1, 2, 3), f'{V1.shape}')
V1, V2 = lambert(r1, r2, np.array([150 * DAY, 200 * DAY]))
check('vector tof shape', V1.shape == (2, 3), f'{V1.shape}')

# ---------------------------------------------------------------- 7. tof <= 0 and degenerate geometry
v1, _ = lambert(r1, r2, 0.0); check('tof=0 -> NaN', np.all(np.isnan(v1)), f'{v1}')
v1, _ = lambert(r1, r1 * 1.5, 100 * DAY); check('collinear same direction -> NaN', np.all(np.isnan(v1)), f'{v1}')
v1, _ = lambert(r1, -r1 * 1.5, 100 * DAY); print('exactly 180 deg:', v1)

print('\nFAILURES:', fails if fails else 'none')
