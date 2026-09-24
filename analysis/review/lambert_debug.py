"""Debug finite-but-wrong Lambert solutions and confirm that hyperbolic (z<0) transfers are never found."""
import sys, pathlib, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from ctoc14.lambert import lambert, _stumpff
from ctoc14.kepler import Ephemeris, propagate_twobody
from ctoc14.constants import MU, AU, DAY

def y_F(r1, r2, tof, z, lw):
    r1n = np.linalg.norm(r1); r2n = np.linalg.norm(r2)
    cosd = np.clip(r1 @ r2 / (r1n * r2n), -1, 1); dnu = np.arccos(cosd)
    if lw: dnu = 2 * np.pi - dnu
    A = np.sin(dnu) * np.sqrt(r1n * r2n / (1 - cosd))
    C, S = _stumpff(np.array([z])); C = C[0]; S = S[0]
    y = r1n + r2n + A * (z * S - 1) / np.sqrt(C)
    F = (y / C) ** 1.5 * S + A * np.sqrt(y) - np.sqrt(MU) * tof if y > 0 else np.nan
    return A, y, F, np.degrees(dnu)

def ref_lambert(r1, r2, tof, lw):
    """Robust scalar reference: locate the admissible z range (y>0) by bisection, then bisection on F."""
    r1n = np.linalg.norm(r1); r2n = np.linalg.norm(r2)
    A = y_F(r1, r2, tof, 0.0, lw)[0]
    def y(z): return y_F(r1, r2, tof, z, lw)[1]
    def F(z): return y_F(r1, r2, tof, z, lw)[2]
    zmax = 4 * np.pi ** 2 * (1 - 1e-12)
    if A > 0:
        # y increasing in z; find z0 with y(z0)=0 below zmax
        lo = -1.0
        while y(lo) > 0: lo *= 2
        if lo < -1e6: return None
        hi = zmax
        for _ in range(200):
            m = 0.5 * (lo + hi)
            if y(m) > 0: hi = m
            else: lo = m
        zlo = hi + 1e-12; zhi = zmax
    else:
        zhi = zmax
        while np.isnan(F(zhi)) or y(zhi) <= 0:
            zhi = 4 * np.pi ** 2 - 2 * (4 * np.pi ** 2 - zhi) if zhi < 4 * np.pi ** 2 - 1e-9 else zhi - 1.0
            if zhi < -1e6: return None
        zlo = zhi - 1.0
        while F(zlo) > 0:
            zlo -= 2 * (zhi - zlo)
            if zlo < -1e6: return None
    flo = F(zlo); fhi = F(zhi)
    if not (flo < 0 < fhi): return None
    for _ in range(300):
        m = 0.5 * (zlo + zhi); fm = F(m)
        if fm > 0: zhi = m
        else: zlo = m
    z = 0.5 * (zlo + zhi); yy = y(z)
    f = 1 - yy / r1n; g = A * np.sqrt(yy / MU); gd = 1 - yy / r2n
    return (r2 - f * r1) / g, (gd * r2 - r1) / g, z

def endpoint_err(r1, v1, tof, r2):
    rr, _ = propagate_twobody(r1, v1, tof); return np.linalg.norm(rr - r2)

eph = Ephemeris()
rng = np.random.default_rng(3)
# reproduce the same draws as lambert_stress.py (tof list order) and collect bad solutions
bad_cases = []
for tof_d in [1.0, 10.0, 30.0, 100.0, 365.0, 800.0, 1500.0, 3000.0]:
    tof = tof_d * DAY; N = 3000
    t0 = rng.uniform(0, 4e8, N); k1 = rng.integers(0, 300, N); k2 = rng.integers(0, 300, N)
    R1 = np.array([eph.ast_state(k, t)[0] for k, t in zip(k1, t0)]); R2 = np.array([eph.ast_state(k, t + tof)[0] for k, t in zip(k2, t0)])
    V1, V2 = lambert(R1, R2, np.full(N, tof))
    for i in range(N):
        if np.all(np.isfinite(V1[i])):
            e = endpoint_err(R1[i], V1[i], tof, R2[i])
            if not (e < 1.0):
                bad_cases.append((tof, R1[i], R2[i], V1[i], e))
print(f'{len(bad_cases)} finite solutions with endpoint error >= 1 km (or NaN on propagation)')
for (tof, r1, r2, v1, e) in bad_cases[:12]:
    cr = np.cross(r1, r2); lw = cr[2] < 0
    A, y0, F0, dnu = y_F(r1, r2, tof, 0.0, lw)
    ref = ref_lambert(r1, r2, tof, lw)
    # recover z from the returned solution: y = r1n(1 - f) where f from v1... instead scan F for sign structure
    zz = np.linspace(-40, 4 * np.pi ** 2 - 1e-6, 9)
    Fs = [y_F(r1, r2, tof, z, lw)[2] for z in zz]
    print(f'\n tof={tof/DAY:.0f} d dnu={dnu:.4f} deg A={A:.4e} |r1|={np.linalg.norm(r1)/AU:.3f} AU |r2|={np.linalg.norm(r2)/AU:.3f} AU')
    print(f'   lambert v1={v1} |v1|={np.linalg.norm(v1):.3e} km/s endpoint err={e:.3e} km')
    print('   F(z) at z=', np.round(zz, 2), ':', np.array([f'{f:.2e}' for f in Fs]))
    if ref is not None:
        print(f'   reference: z={ref[2]:.6f} v1={ref[0]} endpoint err={endpoint_err(r1, ref[0], tof, r2):.3e} km')
    else:
        print('   reference: none')

# ---------------------------------------------------------------- hyperbolic explicit: 1 AU -> 2 AU, 90 deg, parabolic tof = 110 d
r1 = np.array([AU, 0, 0]); r2 = np.array([0, 2 * AU, 0])
print('\n1 AU -> 2 AU, 90 deg (parabolic tof ~ 110 d):')
for tof_d in [200, 150, 120, 111, 110, 109, 100, 60, 30, 10]:
    tof = tof_d * DAY
    v1, _ = lambert(r1, r2, tof); ref = ref_lambert(r1, r2, tof, False)
    e1 = endpoint_err(r1, v1, tof, r2) if np.all(np.isfinite(v1)) else np.nan
    print(f'  tof={tof_d:4d} d: lambert finite={np.all(np.isfinite(v1))!s:5} err={e1:9.2e} km | ref z={ref[2] if ref else np.nan:9.4f} err={endpoint_err(r1, ref[0], tof, r2) if ref else np.nan:9.2e} km')

# why: y(z_lo) <= 0 at z_lo = -4pi^2 makes the pre-loop raise z_lo to 0, discarding every z<0 root
z_lo = -4 * np.pi ** 2
A, y, F, d = y_F(r1, r2, 60 * DAY, z_lo, False); print(f'\ny(z=-4pi^2) = {y:.3e} (<=0 -> z_lo jumps to 0);  y(0) = {y_F(r1, r2, 60*DAY, 0.0, False)[1]:.3e}; F(0) = {y_F(r1, r2, 60*DAY, 0.0, False)[2]:.3e} (>0: root is at z<0)')
