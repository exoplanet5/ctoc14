"""Checks of ctoc14/colgen.py: cost-model fixed point, master LP duals / reduced costs (complementary slackness),
whole-store pricing makes the LP exact, integer cover on a toy instance."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.optimize import linprog
from ctoc14.colgen import CostModel, ColumnStore, RMP, TARGETS
from ctoc14.constants import cost_sc

def main():
    cm = CostModel(s=1.15, reserve=20.0)
    # (a) tank fixed point: m0 = 600 + R + s * fuel * m0 / m0_search
    for F, ms in ((360.0, 1000.0), (450.0, 1100.0), (100.0, 700.0)):
        m0 = float(cm.tank(F, ms)); assert abs(m0 - (600 + 20 + 1.15 * F * m0 / ms)) < 1e-9, (F, ms, m0)
    assert np.isinf(cm.tank(900.0, 1000.0)) and float(cm.cost(900.0, 1000.0)) == 10.0
    print(f'  tank(360 kg @ 1000) = {float(cm.tank(360, 1000)):.1f} kg, J_i {float(cm.cost(360, 1000)):.4f}, w_fuel {cm.w_fuel():.3f}')

    # (b) random store: 3000 routes of 8-30 targets over the 298 coverable targets
    rng = np.random.default_rng(7); store = ColumnStore(cm)
    for k in range(3000):
        n = int(rng.integers(8, 31)); ids = rng.choice(TARGETS, n, replace=False)
        store.add(ids, fuel_est=float(rng.uniform(150, 420)), m0_search=1000.0, t_launch=0.0, src=('inline', {}), tag='rand')
    # duplicate target set with a higher cost must not replace the cheaper one
    j0 = 0; ids0 = store.targets(j0); c0 = store.cost[j0]
    j, new = store.add(ids0, fuel_est=419.0, m0_search=1000.0, t_launch=0.0, src=('inline', {}), tag='dup')
    assert j == j0 and not new and store.cost[j0] == c0
    rows = TARGETS - 1
    for N in (12, 10):
        rmp = RMP(store); res = rmp.solve(rows, N)
        C = store.csr(); cost = store.costs()
        rc = cost - C @ res.pi + res.mu
        assert rc.min() > -1e-6, f'LP not exact over the store: min rc {rc.min()}'
        frac = (res.y > 1e-6) & (res.y < 1 - 1e-6)
        assert np.all(np.abs(rc[res.cols[frac]]) < 1e-6), 'basic fractional columns must have rc = 0'
        # reference LP over the whole store
        n = len(store); nR = len(rows); A = np.zeros((nR + 1, n + nR)); A[:nR, :n] = -C[:, rows].T.toarray(); A[:nR, n:] = -np.eye(nR); A[nR, :n] = 1
        ref = linprog(np.concatenate([cost, np.ones(nR)]), A_ub=A, b_ub=np.concatenate([-np.ones(nR), [N]]), bounds=(0, 1), method='highs')
        assert abs(ref.fun + 2 - res.value) < 1e-6, (ref.fun + 2, res.value)
        assert np.all(res.pi >= -1e-9) and np.all(res.pi <= 1 + 1e-9) and res.mu >= -1e-9
        J, sel, missed = rmp.mip(rows, N, time_limit=30)
        assert len(sel) <= N and J >= res.value - 1e-6
        print(f'  N={N}: LP {res.value:.4f} (= reference), working set {len(rmp.active)} of {len(store)}, mu {res.mu:.3f}, MIP {J:.3f} with {len(sel)} routes, {len(missed)} missed')
    print('ALL TESTS PASSED')

if __name__ == '__main__':
    main()
