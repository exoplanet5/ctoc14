"""Checks for ctoc14/insertion.py: (1) the Lambert reconstruction of every planned tour reproduces the stored dv_est,
tof and fuel_est; (2) the vectorised insertion evaluation agrees with the rebuilt tour; (3) the greedy driver inserts
#137 into tour_sc4 with margin 0 (the known result) and the written JSON round-trips."""
import sys, json, pathlib, glob, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ctoc14.kepler import Ephemeris
from ctoc14.insertion import (tour_from_json, check_reconstruction, evaluate_insertions, best_insertion, apply_insertion,
                              insert_targets, InsParams, mass_chain)
from ctoc14.constants import DAY, VE, TMAX
root = pathlib.Path(__file__).resolve().parents[1]
eph = Ephemeris()
files = sorted(glob.glob(str(root / 'results/fleet_b300/tour_sc*.json')))
assert files, 'no tours found'
tours = {}
worst = 0.0
for f in files:
    tour = json.load(open(f)); nm = pathlib.Path(f).stem
    TM = tour_from_json(eph, tour, name=nm); tours[nm] = TM
    c = check_reconstruction(TM, tour)
    worst = max(worst, c['max_ddv'])
    assert c['max_ddv'] < 1e-6, (nm, c)
    assert abs(c['dfuel']) < 1e-6, (nm, c)
    assert c['max_dtof'] < 1e-6, (nm, c)
    assert c['dvinf'] < 1e-9, (nm, c)
    assert TM.feasible and c['n_nan'] == 0, (nm, c)
print(f'reconstruction: {len(files)} tours, max |dv - dv_est| = {worst:.2e} km/s (required < 1e-3)')
# planner mass model spot check against search.expand's formula
m = 1500.0; dv = 0.8; tof = 60 * DAY
a = TMAX / m * 1e-3; kappa = 1 + dv / (2 * a * tof)
dm, ms, slack, m_end = mass_chain(m, [dv], [tof], 0.6)
assert abs(dm[0] - m * (1 - np.exp(-kappa * dv / VE))) < 1e-12 and abs(slack[0] - (0.6 * a * tof - dv)) < 1e-12
# vectorised evaluation vs rebuilt tour (several candidates, including shifted ones)
P = InsParams(margin=0.0)
TM = tours['tour_sc4']
res = evaluate_insertions(eph, TM, 137, P)
ok = np.where(res['ok_nobudget'])[0]
assert len(ok) > 0
rng = np.random.default_rng(0)
for c in rng.choice(ok, size=min(12, len(ok)), replace=False):
    from ctoc14.insertion import _summ
    cand = _summ(TM, 137, res, int(c))
    new = apply_insertion(eph, TM, cand, P)
    assert abs(new.fuel - cand['fuel_new']) < 1e-6, (cand['fuel_new'], new.fuel)
    assert new.n == TM.n + 1 and 137 in new.asts.tolist()
    k = cand['k']
    assert np.all(new.slack[k + 1:k + 5] >= -1e-7)
print(f'evaluation vs rebuild: {min(12, len(ok))} candidates agree to < 1e-6 kg')
b = best_insertion(eph, TM, 137, P)
assert b['ok'] and b['best']['after_ast'] == 300 and abs(b['best']['extra_kg'] - 30.26) < 0.05, b['best']
# greedy driver on the known case
sub = {nm: tours[nm] for nm in ('tour_sc4', 'tour_sc9')}
sub, hist, table, remaining = insert_targets(eph, sub, [137, 138, 131], P, nproc=1, verbose=False)
assert [h['target'] for h in hist] == [137] and hist[0]['tour'] == 'tour_sc4' and remaining == [138], (hist, remaining)
js = sub['tour_sc4'].to_json()
back = tour_from_json(eph, js, name='rt')
assert abs(back.fuel - sub['tour_sc4'].fuel) < 1e-9 and js['inserted'] == [137] and js['n_flybys'] == 40
print('greedy driver + JSON round trip OK')
print('ALL INSERTION TESTS PASSED')
