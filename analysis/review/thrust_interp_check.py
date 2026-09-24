"""Compare ctoc14.thrust.thrust_interp with a direct scalar implementation of the PDF rule (Sec 6.1 rule 3, 0-based centred
window) for n = 1..7 samples, t exactly at nodes, arc ends, just outside, non-uniform spacing; and check arc_max_thrust."""
import sys, pathlib, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from ctoc14.thrust import thrust_interp, arc_max_thrust, window_start, lagrange_basis

fails = []
def check(name, cond, detail=''):
    (print(f'ok   {name}: {detail}') if cond else (fails.append(name), print(f'FAIL {name}: {detail}')))

def pdf_thrust(ts, Ts, t):
    """Direct transcription of the PDF rule for scalar t."""
    n = len(ts)
    if t < ts[0] or t > ts[-1]:
        return np.zeros(3)
    if n < 4:
        idx = list(range(n))
    else:
        j = max(k for k in range(0, n - 1) if ts[k] <= t)   # largest k in [0, n-2] with t_k <= t
        l = min(max(j - 1, 0), n - 4)
        idx = [l, l + 1, l + 2, l + 3]
    T = np.zeros(3)
    for i in idx:
        L = 1.0
        for k in idx:
            if k != i:
                L *= (t - ts[k]) / (ts[i] - ts[k])
        T += L * Ts[i]
    return T

rng = np.random.default_rng(7)
for n in range(1, 8):
    for spacing in ['uniform', 'nonuniform']:
        if spacing == 'uniform':
            ts = 1e7 + 86400.0 * np.arange(n)
        else:
            ts = 1e7 + np.concatenate([[0.0], np.cumsum(rng.uniform(8640.0, 5 * 86400.0, n - 1))])
        Ts = rng.uniform(-0.5, 0.5, (n, 3))
        # evaluation times: nodes, arc ends, midpoints, random, just outside, node +- ulp
        tt = list(ts) + [ts[0] - 1e-6, ts[-1] + 1e-6, ts[0] - 1.0, ts[-1] + 1.0]
        if n > 1:
            tt += list(0.5 * (ts[:-1] + ts[1:])) + list(rng.uniform(ts[0], ts[-1], 200))
            tt += list(np.nextafter(ts, np.inf)) + list(np.nextafter(ts, -np.inf))
        tt = np.array(tt)
        A = thrust_interp(ts, Ts, tt)
        B = np.array([pdf_thrust(ts, Ts, t) for t in tt])
        d = np.abs(A - B).max()
        check(f'n={n} {spacing} vectorised == PDF rule', d < 1e-12, f'max diff {d:.2e} N')
        # scalar call
        ds = max(np.abs(thrust_interp(ts, Ts, t) - pdf_thrust(ts, Ts, t)).max() for t in tt)
        check(f'n={n} {spacing} scalar == PDF rule', ds < 1e-12, f'max diff {ds:.2e} N')
        # nodes reproduced exactly
        dn = np.abs(thrust_interp(ts, Ts, ts) - Ts).max()
        check(f'n={n} {spacing} nodes exact', dn < 1e-12, f'{dn:.2e}')
        # outside zero
        z = np.abs(thrust_interp(ts, Ts, np.array([ts[0] - 1e-9, ts[-1] + 1e-9, -1e9, 1e12]))).max()
        check(f'n={n} {spacing} zero outside', z == 0.0, f'{z}')

# window_start table for n=5,6
print('window_start n=5:', [window_start(j, 5) for j in range(4)], ' n=6:', [window_start(j, 6) for j in range(5)], ' n=4:', [window_start(j, 4) for j in range(3)], ' n=3:', [window_start(j, 3) for j in range(2)])
check('window_start n=5', [window_start(j, 5) for j in range(4)] == [0, 0, 1, 1], '')
check('window_start n=6', [window_start(j, 6) for j in range(5)] == [0, 0, 1, 2, 2], '')

# continuity at interior nodes and discontinuity of derivative (informational)
ts = 1e7 + 86400.0 * np.arange(8); Ts = rng.uniform(-0.5, 0.5, (8, 3))
for k in range(1, 7):
    lft = thrust_interp(ts, Ts, np.nextafter(ts[k], -np.inf)); rgt = thrust_interp(ts, Ts, np.nextafter(ts[k], np.inf))
    check(f'continuity at node {k}', np.abs(lft - rgt).max() < 1e-9, f'{np.abs(lft - rgt).max():.2e}')

# arc_max_thrust accuracy vs dense grid / exact polynomial maximisation
def exact_max(ts, Ts):
    """Exact max |T| over the arc: per interval, roots of d/dt |p(t)|^2 for the cubic window polynomial."""
    n = len(ts); best = 0.0; tb = ts[0]
    for j in range(n - 1):
        if n < 4: idx = list(range(n))
        else:
            l = min(max(j - 1, 0), n - 4); idx = [l, l + 1, l + 2, l + 3]
        tw = ts[idx]; tau = lambda t: (t - ts[j]) / (ts[j + 1] - ts[j])
        # fit polynomial coefficients per component in tau via Lagrange on the window nodes
        xs = np.array([tau(t) for t in tw])
        polys = [np.poly1d(np.polyfit(xs, Ts[idx, c], len(idx) - 1)) for c in range(3)]
        S = sum(p * p for p in polys)
        cands = [0.0, 1.0] + [float(r.real) for r in S.deriv().roots if abs(r.imag) < 1e-12 and 0 < r.real < 1]
        for c in cands:
            val = np.sqrt(S(c))
            if val > best: best = val; tb = ts[j] + c * (ts[j + 1] - ts[j])
    return best, tb

for name, Ts in [('step', np.array([[0, 0, 0], [0, 0, 0], [.5, 0, 0], [.5, 0, 0], [.5, 0, 0], [0, 0, 0]], float)),
                 ('plateau2', np.array([[0, 0, 0], [.5, 0, 0], [.5, 0, 0], [0, 0, 0]], float)),
                 ('turn90', np.array([[.5, 0, 0], [.5, 0, 0], [0, .5, 0], [0, .5, 0]], float)),
                 ('random', rng.uniform(-0.5, 0.5, (9, 3)))]:
    ts = 1e7 + 86400.0 * np.arange(len(Ts))
    g, tg = arc_max_thrust(ts, Ts); ex, tx = exact_max(ts, Ts)
    print(f'{name:8s}: grid max {g:.9f} at {tg-1e7:.1f}, exact max {ex:.9f} at {tx-1e7:.1f}, grid underestimate {ex - g:.2e} N')
    check(f'arc_max_thrust {name} within 1e-6 of exact', abs(ex - g) < 1e-6, f'{ex - g:.2e}')

# single-sample arc
check('n=1 arc_max_thrust', arc_max_thrust(np.array([5.0]), np.array([[0.3, 0.4, 0.0]])) == (0.5, 5.0), str(arc_max_thrust(np.array([5.0]), np.array([[0.3, 0.4, 0.0]]))))

print('\nFAILURES:', fails if fails else 'none')
