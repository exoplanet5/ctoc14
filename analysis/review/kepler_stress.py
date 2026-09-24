"""Adversarial tests of ctoc14.kepler: solve_kepler at high e / M near 0, 2pi; Body.state consistency; propagate_twobody
for very long dt, negative dt, near-parabolic and hyperbolic orbits (compared against closed-form Kepler propagation)."""
import sys, pathlib, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from ctoc14.kepler import solve_kepler, Ephemeris, Body, propagate_twobody
from ctoc14.constants import MU, AU, T_MISSION, DAY, T_OFF_AST

np.set_printoptions(precision=6, linewidth=150)
fails = []
def check(name, cond, detail=''):
    if not cond:
        fails.append(name)
        print(f'FAIL {name}: {detail}')
    else:
        print(f'ok   {name}: {detail}')

# ---------------------------------------------------------------- 1. solve_kepler residuals at e = 0.947, 0.99, 0.999
for e in [0.0, 0.3, 0.7999, 0.8, 0.947, 0.9469937, 0.99, 0.999, 0.9999]:
    Ms = np.concatenate([np.array([0.0, 1e-300, 1e-16, 1e-12, 1e-8, 1e-4, np.pi, 2 * np.pi - 1e-4, 2 * np.pi - 1e-8, 2 * np.pi - 1e-12,
                                   2 * np.pi - 1e-15, 2 * np.pi, -1e-12, -1e-8, 2 * np.pi + 1e-12, 1000.0, -1000.0, 1e6, -1e6]),
                        np.linspace(0, 2 * np.pi, 20001), np.random.default_rng(0).uniform(-500, 500, 20000)])
    E = solve_kepler(Ms, e)
    res = np.abs(E - e * np.sin(E) - np.mod(Ms, 2 * np.pi))
    # scale-aware residual: the true error in E is res/(1 - e cos E)
    dE = res / (1 - e * np.cos(E))
    check(f'kepler e={e} residual', np.max(res) < 1e-11 and np.max(dE) < 1e-9,
          f'max |E - e sinE - M| = {res.max():.2e}, max |dE| = {dE.max():.2e} at M={Ms[np.argmax(dE)]:.6e}')
    # iteration count / convergence: rerun with maxit=1000 and compare
    E2 = solve_kepler(Ms, e, maxit=1000)
    check(f'kepler e={e} converged within 60 iters', np.max(np.abs(E2 - E)) < 1e-12, f'max |E60 - E1000| = {np.max(np.abs(E2 - E)):.2e}')

# scalar input and 0-d array input
E0 = solve_kepler(0.5, 0.947); check('kepler scalar input', np.ndim(E0) == 0 and abs(E0 - 0.947 * np.sin(E0) - 0.5) < 1e-12, f'E={E0}')

# ---------------------------------------------------------------- 2. ephemeris ids / consistency of Body.state vs all_ast_states
eph = Ephemeris()
check('MEA ids are 1..300 in order', np.array_equal(eph.ids, np.arange(1, 301)), f'ids[:5]={eph.ids[:5]}, n={len(eph.ids)}')
for t in [0.0, 1234567.0, T_MISSION, -T_OFF_AST]:
    ra, va = eph.all_ast_states(t)
    rb = np.array([eph.ast_state(k, t)[0] for k in range(300)]); vb = np.array([eph.ast_state(k, t)[1] for k in range(300)])
    check(f'all_ast_states == ast_state t={t:.3e}', np.max(np.abs(ra - rb)) < 1e-6 and np.max(np.abs(va - vb)) < 1e-12,
          f'dr={np.max(np.abs(ra - rb)):.2e} km dv={np.max(np.abs(va - vb)):.2e}')
# energy / angular momentum of Body.state over time
for k in [143, 22, 295, 130]:
    b = eph.ast[k]
    tt = np.linspace(0, T_MISSION, 5001)
    r, v = b.state(tt)
    rn = np.linalg.norm(r, axis=-1)
    eng = 0.5 * np.sum(v * v, -1) - MU / rn
    h = np.cross(r, v)
    check(f'ast {k+1} (e={b.e:.4f}, a={b.a/AU:.3f} AU) energy const', np.ptp(eng) / abs(eng.mean()) < 1e-11, f'rel ptp energy {np.ptp(eng)/abs(eng.mean()):.2e}, a_from_energy/a = {-MU/(2*eng.mean())/b.a - 1:.2e}')
    check(f'ast {k+1} h const', np.ptp(np.linalg.norm(h, -1)) / np.linalg.norm(h, -1).mean() < 1e-12, f'rel ptp h {np.ptp(np.linalg.norm(h, -1)) / np.linalg.norm(h, -1).mean():.2e}')
    # r == a(1 - e cos E) and periapsis distance
    check(f'ast {k+1} rmin ~ q', abs(rn.min() - b.a * (1 - b.e)) / b.a < 5e-3, f'rmin/a={rn.min()/b.a:.5f} q/a={1-b.e:.5f}')
    # velocity is d r / dt (finite difference check, 5-point)
    hfd = 10.0
    tt2 = np.array([5e7 + m * hfd for m in (-2, -1, 1, 2)])
    r2, _ = b.state(tt2)
    vfd = (-r2[3] + 8 * r2[2] - 8 * r2[1] + r2[0]) / (12 * hfd)
    _, v5 = b.state(5e7)
    check(f'ast {k+1} v == dr/dt', np.linalg.norm(vfd - v5) < 1e-6, f'|v_fd - v| = {np.linalg.norm(vfd - v5):.2e} km/s')

# Earth: position at t=0 vs PDF, velocity direction check
rE, vE = eph.earth_state(0.0)
check('Earth r(0) vs PDF', np.linalg.norm(rE - np.array([-1.9500328647e+07, 1.4581848347e+08, -7.6372048832e+03])) < 0.01, f'{np.linalg.norm(rE - np.array([-1.9500328647e+07, 1.4581848347e+08, -7.6372048832e+03])):.2e} km')

# ---------------------------------------------------------------- 3. propagate_twobody vs analytic Kepler for asteroid orbits
def prop_check(name, k, t0, dt, tol_r=1e-4, tol_v=1e-9):
    b = eph.ast[k]
    r0, v0 = b.state(t0)
    r1, v1 = propagate_twobody(r0, v0, dt)
    ra, va = b.state(t0 + dt)
    er = np.linalg.norm(r1 - ra); ev = np.linalg.norm(v1 - va)
    check(name, np.isfinite(er) and er < tol_r and ev < tol_v, f'ast {k+1} e={b.e:.4f} a={b.a/AU:.2f} AU dt={dt/DAY:.1f} d: dr={er:.3e} km dv={ev:.3e} km/s')
    return er, ev

worst = 0
for k in [143, 22, 277, 295, 130, 0, 2, 3]:
    for dt in [T_MISSION, -T_MISSION, 2 * T_MISSION, 10.0, -10.0, 0.0, 1e-3]:
        er, ev = prop_check(f'prop ast{k+1} dt={dt:.3e}', k, 1.0e7, dt)
        worst = max(worst, er)
# all asteroids, 15 yr forward and backward, from several start epochs
for t0 in [0.0, 1.2e8, T_MISSION]:
    ers = []
    for k in range(300):
        b = eph.ast[k]; r0, v0 = b.state(t0)
        for dt in [T_MISSION - t0, -t0 - 1.0]:
            r1, v1 = propagate_twobody(r0, v0, dt); ra, va = b.state(t0 + dt)
            ers.append((np.linalg.norm(r1 - ra), np.linalg.norm(v1 - va), k, dt))
    ers = np.array(ers)
    kk = int(np.argmax(ers[:, 0]))
    check(f'prop all asteroids from t0={t0:.2e}', np.all(np.isfinite(ers[:, :2])) and ers[:, 0].max() < 1e-3 and ers[:, 1].max() < 1e-8,
          f'max dr={ers[:,0].max():.3e} km (ast {int(ers[kk,2])+1}, dt={ers[kk,3]:.3e}), max dv={ers[:,1].max():.3e} km/s')

# ---------------------------------------------------------------- 4. near-parabolic and hyperbolic (compare with analytic via Body for ellipse; via
# energy/angular momentum + reversibility for hyperbola)
def make_state(a_km, e, nu):
    p = a_km * (1 - e * e)
    r = p / (1 + e * np.cos(nu))
    rp = np.array([r * np.cos(nu), r * np.sin(nu), 0.0])
    vp = np.sqrt(MU / p) * np.array([-np.sin(nu), e + np.cos(nu), 0.0])
    return rp, vp

for e in [0.999, 0.999999, 1 - 1e-9, 1 - 1e-12, 1.0, 1 + 1e-12, 1 + 1e-9, 1.000001, 1.001, 1.5, 5.0]:
    q = 0.3 * AU
    if e == 1.0:
        a_km = np.inf
        p = 2 * q
        for nu in [-2.0, -0.5, 0.0, 0.5, 2.0]:
            r = p / (1 + np.cos(nu)); rp = np.array([r * np.cos(nu), r * np.sin(nu), 0]); vp = np.sqrt(MU / p) * np.array([-np.sin(nu), 1 + np.cos(nu), 0])
            for dt in [30 * DAY, -30 * DAY, 1000 * DAY, -1000 * DAY]:
                r1, v1 = propagate_twobody(rp, vp, dt)
                r2, v2 = propagate_twobody(r1, v1, -dt)
                eng0 = 0.5 * vp @ vp - MU / np.linalg.norm(rp); eng1 = 0.5 * v1 @ v1 - MU / np.linalg.norm(r1)
                check(f'parabola nu={nu} dt={dt/DAY:+.0f}d', np.all(np.isfinite(r1)) and np.linalg.norm(r2 - rp) < 1e-3 and abs(eng1 - eng0) < 1e-9,
                      f'round-trip dr={np.linalg.norm(r2 - rp):.2e} km, energy {eng0:.2e}->{eng1:.2e}')
        continue
    a_km = q / (1 - e)
    for nu in [-2.0, -0.5, 0.0, 0.5, 2.0]:
        rp, vp = make_state(a_km, e, nu)
        for dt in [30 * DAY, -30 * DAY, 1000 * DAY, -1000 * DAY]:
            r1, v1 = propagate_twobody(rp, vp, dt)
            r2, v2 = propagate_twobody(r1, v1, -dt)
            eng0 = 0.5 * vp @ vp - MU / np.linalg.norm(rp); eng1 = 0.5 * v1 @ v1 - MU / np.linalg.norm(r1)
            h0 = np.linalg.norm(np.cross(rp, vp)); h1 = np.linalg.norm(np.cross(r1, v1))
            ok = np.all(np.isfinite(r1)) and np.linalg.norm(r2 - rp) < 1e-2 and abs(eng1 - eng0) < 1e-8 * max(1, abs(eng0)) and abs(h1 - h0) < 1e-6 * h0
            if e < 1:
                # compare against exact Kepler propagation through Body
                nn = np.sqrt(MU / a_km ** 3)
                E0 = 2 * np.arctan2(np.sqrt(1 - e) * np.sin(nu / 2), np.sqrt(1 + e) * np.cos(nu / 2))
                M0 = E0 - e * np.sin(E0)
                b = Body([a_km / AU, e, 0.0, 0.0, 0.0, np.degrees(M0)], 0.0)
                ra, va = b.state(dt)
                ok = ok and np.linalg.norm(r1 - ra) < 1e-2 and np.linalg.norm(v1 - va) < 1e-9
                extra = f' vs Kepler dr={np.linalg.norm(r1 - ra):.2e} km dv={np.linalg.norm(v1 - va):.2e}'
            else:
                extra = ''
            check(f'e={e} nu={nu} dt={dt/DAY:+.0f}d', ok, f'round-trip dr={np.linalg.norm(r2 - rp):.2e} km, dE/E={abs(eng1 - eng0)/max(1e-30,abs(eng0)):.1e}, dh/h={abs(h1-h0)/h0:.1e}' + extra)

# ---------------------------------------------------------------- 5. array dt and scalar dt shapes; dt=0
r0, v0 = eph.earth_state(0.0)
rr, vv = propagate_twobody(r0, v0, np.array([0.0, 1.0, DAY, -DAY]))
check('array dt shape', rr.shape == (4, 3) and vv.shape == (4, 3) and np.allclose(rr[0], r0), f'{rr.shape}')
r1, v1 = propagate_twobody(r0, v0, 0.0)
check('dt=0 identity', np.array_equal(r1, r0) and np.array_equal(v1, v0), '')

print('\nFAILURES:', fails if fails else 'none')
