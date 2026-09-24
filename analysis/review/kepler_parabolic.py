"""propagate_twobody near the parabolic threshold (|alpha| ~ 1e-12 1/km): forward/backward accuracy against a tight
numerical integration, and the negative-dt NaN in the parabolic branch."""
import sys, pathlib, numpy as np, warnings
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from scipy.integrate import solve_ivp
from ctoc14.kepler import propagate_twobody
from ctoc14.constants import MU, AU, DAY

def rhs(t, y):
    r = y[:3]; rn = np.linalg.norm(r)
    return np.concatenate([y[3:], -MU * r / rn ** 3])

def ref_prop(r0, v0, dt):
    sol = solve_ivp(rhs, (0, dt), np.concatenate([r0, v0]), method='DOP853', rtol=1e-13, atol=1e-6)
    return sol.y[:3, -1], sol.y[3:, -1]

def state(q, e, nu):
    p = q * (1 + e)
    r = p / (1 + e * np.cos(nu))
    return np.array([r * np.cos(nu), r * np.sin(nu), 0.0]), np.sqrt(MU / p) * np.array([-np.sin(nu), e + np.cos(nu), 0.0])

print(f'{"e":>16s} {"alpha[1/km]":>12s} {"nu":>5s} {"dt[d]":>7s} {"dr[km]":>10s} {"dv[km/s]":>10s} branch')
q = 0.5 * AU
for e in [0.99, 0.9999, 0.999999, 1 - 3e-7, 1 - 1e-7, 1 - 3e-8, 1 - 1e-9, 1.0, 1 + 1e-9, 1 + 3e-8, 1 + 1e-7, 1 + 1e-6, 1.001, 1.5]:
    p = q * (1 + e); a_inv = (1 - e) / q      # alpha = 1/a = (1-e)/q
    for nu in [-1.5, 0.0, 1.5]:
        r0, v0 = state(q, e, nu)
        for dt in [30 * DAY, -30 * DAY, 400 * DAY, -400 * DAY]:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                r1, v1 = propagate_twobody(r0, v0, dt)
            rr, vr = ref_prop(r0, v0, dt)
            branch = 'par' if abs(a_inv) <= 1e-12 else ('hyp' if a_inv < 0 else 'ell')
            print(f'{e:16.12f} {a_inv:12.2e} {nu:5.1f} {dt/DAY:7.0f} {np.linalg.norm(r1 - rr):10.2e} {np.linalg.norm(v1 - vr):10.2e} {branch}')

# root cause of the NaN: negative base to fractional power
s = np.array([-0.3]); print('\nnp.tan(s)**(1/3) for s<0 ->', np.tan(s) ** (1 / 3), '   np.cbrt ->', np.cbrt(np.tan(s)))
