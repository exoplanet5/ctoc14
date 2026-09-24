"""Adversarial probe of ctoc14.kepler (independent of the earlier kepler_stress.py).

1. solve_kepler at e = 0.947 with M near 0 and 2pi (residual, iteration count, monotone convergence), e = 0, e = 0.99.
2. Body.state / all_ast_states self-consistency and against an independent element->state routine.
3. propagate_twobody: 15-year dt (all 300 asteroids + Earth vs the ephemeris itself), negative dt, dt = 0, mixed-sign arrays,
   near-parabolic (|alpha| around the 1e-12 branch threshold) both signs of dt, hyperbolic both signs, near-rectilinear
   ellipse (the orbit type that appears in long-tof / small-dnu Lambert solutions).
Reference for propagation: high-accuracy solve_ivp (DOP853, rtol 1e-13) of the two-body EOM.
"""
import sys, pathlib, numpy as np, warnings
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from scipy.integrate import solve_ivp
from ctoc14.kepler import solve_kepler, Ephemeris, Body, propagate_twobody, DEG
from ctoc14.constants import MU, AU, DAY, T_MISSION, T_OFF_AST, T_OFF_EARTH

np.set_printoptions(precision=5, linewidth=160)
FAILS = []
def check(name, cond, detail=''):
    print(('ok   ' if cond else 'FAIL ') + name + (': ' + detail if detail else ''))
    if not cond: FAILS.append(name)

def ivp(r0, v0, dt):
    f = lambda t, y: np.concatenate([y[3:], -MU * y[:3] / np.linalg.norm(y[:3]) ** 3])
    s = solve_ivp(f, (0, dt), np.concatenate([r0, v0]), method='DOP853', rtol=1e-13, atol=1e-9)
    return s.y[:3, -1], s.y[3:, -1]

# ------------------------------------------------------------------ 1. Kepler equation
print('=== 1. solve_kepler')
def kepler_iters(M, e, tol=1e-13, maxit=60):
    """Re-implementation of the solver loop that also returns the iteration count."""
    M = np.mod(np.asarray(M, float), 2 * np.pi); E = np.where(e < 0.8, M + e * np.sin(M), np.pi * np.ones_like(M))
    for it in range(1, maxit + 1):
        f = E - e * np.sin(E) - M; dE = np.clip(-f / (1 - e * np.cos(E)), -1, 1); E = E + dE
        if np.all(np.abs(dE) < tol): return E, it
    return E, maxit + 1
for e in (0.0, 0.5, 0.8, 0.947, 0.99):
    Ms = np.array([0.0, 1e-15, 1e-12, 1e-9, 1e-6, 1e-3, 0.05, np.pi - 1e-9, np.pi, np.pi + 1e-9, 2 * np.pi - 0.05, 2 * np.pi - 1e-3,
                   2 * np.pi - 1e-6, 2 * np.pi - 1e-9, 2 * np.pi - 1e-12, 2 * np.pi - 1e-15, 2 * np.pi, 2 * np.pi + 1e-12, -1e-12, -1e-9, 4 * np.pi + 1e-9, 1e4, -1e4])
    E = solve_kepler(Ms, e)
    res = np.abs(E - e * np.sin(E) - np.mod(Ms, 2 * np.pi))
    # residual in mean anomaly -> position error ~ a * res / (1-e); at a = 1 AU that is 1.5e8 * res / (1-e) km
    E2, its = kepler_iters(Ms, e)
    worst = np.argmax(res)
    check(f'solve_kepler e={e} residual', res.max() < 1e-12, f'max |E - e sinE - M| = {res.max():.2e} at M={Ms[worst]:.3e}; iterations {its}; E in [{E.min():.3e},{E.max():.6f}]')
    check(f'solve_kepler e={e} E range', np.all(E >= -1e-12) and np.all(E <= 2 * np.pi + 1e-12), '')
# iteration count per M (scalar) at e = 0.947 to see whether any single value stalls near maxit
its_list = []
for M in np.linspace(0, 2 * np.pi, 20001):
    _, it = kepler_iters(M, 0.947); its_list.append(it)
its_list = np.array(its_list)
check('e=0.947 iteration count over 20001 M values', its_list.max() < 60, f'max iterations {its_list.max()} (mean {its_list.mean():.1f}); count >= 20: {(its_list >= 20).sum()}')
# residual sweep over all 300 asteroids at 2000 random mission times
eph = Ephemeris()
rng = np.random.default_rng(0)
tt = rng.uniform(0, T_MISSION, 2000)
worst = 0.0
for t in tt[:300]:
    M = eph._M0 + eph._n * (t + T_OFF_AST); E = solve_kepler(M, eph._e)
    worst = max(worst, np.abs(E - eph._e * np.sin(E) - np.mod(M, 2 * np.pi)).max())
check('all-asteroid Kepler residual over 300 random times', worst < 1e-12, f'max residual {worst:.2e}')

# ------------------------------------------------------------------ 2. element -> state, independent implementation
print('=== 2. Body.state vs independent implementation')
def state_ref(elem, t_off, t):
    a, e, i, Om, w, M0 = elem; a *= AU; i *= DEG; Om *= DEG; w *= DEG; M0 *= DEG
    n = np.sqrt(MU / a ** 3); M = np.mod(M0 + n * (t + t_off), 2 * np.pi)
    E = M
    for _ in range(200):
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
    p = a * (1 - e * e); r = p / (1 + e * np.cos(nu))
    rp = r * np.array([np.cos(nu), np.sin(nu), 0.0]); vp = np.sqrt(MU / p) * np.array([-np.sin(nu), e + np.cos(nu), 0.0])
    def R3(th): c, s = np.cos(th), np.sin(th); return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    def R1(th): c, s = np.cos(th), np.sin(th); return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    Q = R3(Om) @ R1(i) @ R3(w)
    return Q @ rp, Q @ vp
worst_r = worst_v = 0.0
for k in range(300):
    for t in (0.0, 1.234e8, T_MISSION):
        r1, v1 = eph.ast_state(k, t); r2, v2 = state_ref(eph.elem[k], T_OFF_AST, t)
        worst_r = max(worst_r, np.linalg.norm(r1 - r2)); worst_v = max(worst_v, np.linalg.norm(v1 - v2))
rE1, vE1 = eph.earth_state(T_MISSION); rE2, vE2 = state_ref(eph.earth.a / AU * np.ones(1)[0] if False else np.array([1.0009175020, 0.017566762041, 0.002976847126, 189.953211282428, 273.196254000254, 357.4135031077]), T_OFF_EARTH, T_MISSION)
check('asteroid states vs independent routine (300 x 3 epochs)', worst_r < 1e-6 and worst_v < 1e-12, f'max dr={worst_r:.2e} km, dv={worst_v:.2e} km/s')
check('Earth state at T_MISSION vs independent routine', np.linalg.norm(rE1 - rE2) < 1e-6, f'dr={np.linalg.norm(rE1 - rE2):.2e} km')
# vectorised all_ast_states vs per-body
ra, va = eph.all_ast_states(5.5e7); rb = np.array([eph.ast_state(k, 5.5e7)[0] for k in range(300)]); vb = np.array([eph.ast_state(k, 5.5e7)[1] for k in range(300)])
check('all_ast_states == ast_state', np.abs(ra - rb).max() < 1e-7 and np.abs(va - vb).max() < 1e-13, f'{np.abs(ra - rb).max():.2e} km')
# array-time evaluation shape
r_arr, v_arr = eph.ast_state(0, np.array([0.0, 1e6, 2e6]))
check('Body.state array time shape', r_arr.shape == (3, 3) and np.allclose(r_arr[1], eph.ast_state(0, 1e6)[0]), f'{r_arr.shape}')
# energy / angular momentum conservation along the ephemeris (catches rotation-matrix or velocity-scaling errors)
worst = 0.0
for k in (0, 50, 100, 143, 200, 299):
    b = eph.ast[k]; ts = np.linspace(0, T_MISSION, 50); r, v = b.state(ts)
    En = 0.5 * np.sum(v * v, 1) - MU / np.linalg.norm(r, axis=1); worst = max(worst, np.ptp(En) / abs(En.mean()))
check('ephemeris energy constant along orbit', worst < 1e-12, f'max relative spread {worst:.2e}')

# ------------------------------------------------------------------ 3. propagate_twobody
print('=== 3. propagate_twobody')
# 3a. 15-yr dt for every asteroid and Earth: ephemeris state at t0 propagated by T_MISSION must equal the ephemeris at T_MISSION
worst_r = worst_v = 0.0; worst_k = -1
for k in range(300):
    r0, v0 = eph.ast_state(k, 0.0); r1, v1 = eph.ast_state(k, T_MISSION)
    rp, vp = propagate_twobody(r0, v0, T_MISSION)
    er = np.linalg.norm(rp - r1)
    if er > worst_r: worst_r, worst_k = er, k
    worst_v = max(worst_v, np.linalg.norm(vp - v1))
check('15-yr propagation of all asteroids vs ephemeris', worst_r < 1e-3, f'max dr={worst_r:.2e} km (asteroid index {worst_k}, e={eph.elem[worst_k, 1]:.3f}, a={eph.elem[worst_k, 0]:.3f} AU), dv={worst_v:.2e} km/s')
r0, v0 = eph.earth_state(0.0); rp, vp = propagate_twobody(r0, v0, T_MISSION); r1, v1 = eph.earth_state(T_MISSION)
check('15-yr propagation of Earth', np.linalg.norm(rp - r1) < 1e-3, f'dr={np.linalg.norm(rp - r1):.2e} km')
# 3b. negative dt (elliptic): propagate to T_MISSION and back
worst = 0.0
for k in range(0, 300, 7):
    r0, v0 = eph.ast_state(k, T_MISSION); r1, v1 = eph.ast_state(k, 0.0)
    rp, vp = propagate_twobody(r0, v0, -T_MISSION); worst = max(worst, np.linalg.norm(rp - r1))
check('negative 15-yr dt (elliptic)', worst < 1e-3, f'max dr={worst:.2e} km')
# 3c. dt = 0 and mixed-sign arrays
r0, v0 = eph.ast_state(3, 0.0)
rp, vp = propagate_twobody(r0, v0, 0.0)
check('dt = 0 returns the input', np.allclose(rp, r0) and np.allclose(vp, v0), '')
dts = np.array([-3e8, -1e5, 0.0, 1e5, 3e8]); rps, vps = propagate_twobody(r0, v0, dts)
ref = [eph.ast_state(3, d)[0] for d in dts]
check('mixed-sign dt array', max(np.linalg.norm(rps[i] - ref[i]) for i in range(5)) < 1e-4, f'max dr={max(np.linalg.norm(rps[i] - ref[i]) for i in range(5)):.2e} km')
# 3d. near-parabolic: build states with alpha = 1/a spanning the 1e-12 1/km branch threshold, q = 0.7 AU, at true anomaly nu
def state_from_qe(q, e, nu):
    p = q * (1 + e); r = p / (1 + e * np.cos(nu))
    return r * np.array([np.cos(nu), np.sin(nu), 0.0]), np.sqrt(MU / p) * np.array([-np.sin(nu), e + np.cos(nu), 0.0])
print('  near-parabolic / hyperbolic: e, alpha[1/km], nu, dt[d], |dr| km vs solve_ivp')
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    for e in (0.99, 0.9999, 0.999999, 1.0 - 1e-9, 1.0, 1.0 + 1e-9, 1.000001, 1.0001, 1.01, 1.5, 3.0):
        for nu in (-1.0, 0.0, 1.0):
            for dtd in (-200.0, -30.0, 30.0, 200.0):
                r0, v0 = state_from_qe(0.7 * AU, e, nu); alpha = 2 / np.linalg.norm(r0) - v0 @ v0 / MU
                rp, vp = propagate_twobody(r0, v0, dtd * DAY); rr, vr = ivp(r0, v0, dtd * DAY)
                err = np.linalg.norm(rp - rr)
                ok = np.isfinite(err) and err < 1.0
                if not ok:
                    print(f'  FAIL e={e:<14} alpha={alpha:+.2e} nu={nu:+.1f} dt={dtd:+6.0f} d: |dr|={err:.3e} km (ref |r|={np.linalg.norm(rr)/AU:.3f} AU)')
                    FAILS.append(f'near-parabolic e={e} nu={nu} dt={dtd}')
# 3e. near-rectilinear ellipse (q tiny, e -> 1, a = 2 AU): the kind of orbit small-dnu Lambert solutions produce
print('  near-rectilinear ellipses (a = 2 AU):')
for q_au in (1e-2, 1e-3, 1e-4, 1e-5):
    e = 1 - q_au / 2.0
    for nu in (0.5, 2.0, 3.0):
        for dtd in (30.0, 800.0):
            r0, v0 = state_from_qe(q_au * AU, e, nu)
            rp, vp = propagate_twobody(r0, v0, dtd * DAY); rr, vr = ivp(r0, v0, dtd * DAY)
            err = np.linalg.norm(rp - rr); rel = err / np.linalg.norm(rr)
            if not (np.isfinite(err) and rel < 1e-6):
                print(f'  FAIL q={q_au:.0e} AU e={e:.6f} nu={nu} dt={dtd:.0f} d: |dr|={err:.3e} km rel={rel:.2e}')
                FAILS.append(f'rectilinear q={q_au} nu={nu} dt={dtd}')
# 3f. Newton iteration count / residual instrumentation for the 15-yr case of the fastest asteroid (most revolutions)
kfast = int(np.argmin(eph._a)); r0, v0 = eph.ast_state(kfast, 0.0)
rp, vp = propagate_twobody(r0, v0, T_MISSION, maxit=100); r1, v1 = eph.ast_state(kfast, T_MISSION)
n_rev = T_MISSION * eph._n[kfast] / (2 * np.pi)
check(f'fastest asteroid (a={eph._a[kfast]/AU:.3f} AU, {n_rev:.1f} revs in 15 yr) default maxit', np.linalg.norm(rp - r1) < 1e-3, f'dr={np.linalg.norm(rp - r1):.2e} km')
rp, vp = propagate_twobody(r0, v0, T_MISSION, maxit=5)
print(f'  same with maxit=5: dr={np.linalg.norm(rp - r1):.2e} km (no convergence warning is ever emitted by propagate_twobody)')

print('\nFAILURES:', FAILS if FAILS else 'none')
