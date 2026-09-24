"""Insertion-cost surrogate: what does it cost a FINISHED route to also fly one more target?

WHAT WAS MEASURED
-----------------
Hosts: real SETTLED impulsive routes taken from several independent fleets
(results/s7/deep1/fleet, results/s7/F1/fleet, results/s8/J9base, results/s7/base8,
results/n8s6/H2, results/s7/farmO), spanning 10-48 flybys and ~610-1150 kg of tank.

For each host, tools/insert_sample.py samples the host's Kepler-arc trajectory, finds the local
minima of the distance to every asteroid NOT already on the route (run_ialns.w_cands, dmax 0.35 AU),
keeps the closest approach per target, stratifies over 8 distance buckets, and then calls
run_ialns.w_insert on each candidate: aim-point homotopy (globalopt.insert_homotopy, 10 stages)
followed by a full restore / optimise / restore settle of the WHOLE route. The measured quantity is
therefore the REAL post-settle tank increase, not an estimate:

    kg_cost = tank_after - tank_before

Failures are recorded too. Every observed failure was a HOMOTOPY divergence (final miss > 1e4 km);
no insertion that reached the settle stage then failed to settle. So `predict_ok` answers
"can this route reach that target at all", and it falls off with host depth as well as with distance.

Features (all available BEFORE any settle, so a fleet designer can price a move for free):
    d_au     closest approach of the target to the host trajectory [AU]  (run_ialns.w_cands 'dist')
    depth    flybys already on the host route
    tank     host tank [kg] before the insertion
    margin_d days from the chosen insertion epoch to the host's nearest existing flyby

MEASURED (1225 trials, 572 settled, 55 distinct hosts; results/s10/insert/samples.jsonl)
----------------------------------------------------------------------------------------
    kg   = (1202.67 * d**1.128) * (tank/900)^-0.124 * (depth/30)^0.954 * (max(margin_d,3)/60)^-0.744
    p_ok = sigma(-2.433 - 14.288 d - 0.0541 depth - 1.792 tank/1000 + 4.298 log10 max(margin_d,1))

  R2 in log kg 0.651, median |error| 18.4 kg (49 % relative), 18.7 kg leave-one-HOST-out, against
  27.8 kg for distance alone and 35.2 kg for a constant. P(ok) is well calibrated (Brier 0.121,
  quintiles 0.03/0.17/0.44/0.74/0.96 predicted vs 0.03/0.16/0.42/0.78/0.96 observed).

THE HEADLINE -- DEPTH
---------------------
depth and tank are 0.91-correlated in real routes (a route gains ~14.8 kg of tank per flyby), so the
operationally meaningful depth effect is the one measured ALONG that ridge, not the partial exponent
at fixed tank. The fit is therefore done in the orthogonal pair (size exponent G, off-ridge tank
exponent P) and G is the number to quote:

    G = 0.91   ->  kg cost ~ depth^0.91,  i.e. DOUBLING host depth almost DOUBLES the price (x1.87),
                   +10 flybys multiplies it by ~1.37

and feasibility falls at the same time: at 0.05 AU and a 60 d margin, P(ok) runs
0.94 (depth 12) -> 0.89 (20) -> 0.78 (30) -> 0.61 (40) -> 0.49 (46).

Combining the two, in contest units (breakeven is 0.0369 J per flyby):

    host depth   median dJ of a forced insertion   P(settle)   fraction of ATTEMPTS under breakeven
      10-19            0.0160                        0.74               53 %
      20-29            0.0387                        0.60               29 %
      30-37            0.0971                        0.31                6 %
      38-48            0.0974                        0.21                1 %

That is the mechanism behind the measured 3x "freedom premium": a 35-flyby route is not merely
expensive to extend, it is ~9x less likely than a 15-flyby route to absorb an extra target at a
price the fleet can afford. (The last column is over the stratified sample, which is uniform in
distance bucket, not over a natural candidate distribution -- read it as a relative comparison
across depths.)

CAVEATS
-------
- No constant floor: the fitted A pinned at 0, so predict_kg -> 0 as d -> 0. Real costs at
  d < 0.01 AU are 2-6 kg on shallow hosts; the model is clipped at 0.05 kg, not at that floor.
- depth 38-48 is thin (23 settled samples out of 107 trials) because deep hosts mostly fail.
- margin_d helps only above ~20 d: measured medians at d < 0.06 AU are 45 kg (margin 0-20 d),
  60 kg (20-50), 20 kg (50-100), 7 kg (>100). Insert into a long empty gap, not into a crowded leg.
- 2 % of settled insertions came out FREE or negative (-9 to 0 kg), all on hosts of depth 18-28:
  the re-settle occasionally finds a better whole-route solution than the one it started from.

USE
---
    import sys; sys.path.insert(0, 'tools')
    from insert_model import predict_kg, predict_ok, predict_dJ, predict_expected_dJ
    kg = predict_kg(0.05, depth=34, tank=980, margin_d=60)

Model JSON: results/s10/insert_model.json (coefficients, bucketed tables, fit quality, n_samples).
Override with the INSERT_MODEL environment variable or model(path=...).
"""
import os, json, math, pathlib

_HERE = pathlib.Path(__file__).resolve().parent
_JSON = _HERE.parents[0] / 'results' / 's10' / 'insert_model.json'
_M = None

M_DRY = 600.0


def model(path=None):
    """The fitted model dict (lazily loaded from results/s10/insert_model.json)."""
    global _M
    if _M is None or path is not None:
        _M = json.load(open(path or os.environ.get('INSERT_MODEL', _JSON)))
    return _M


def predict_kg(d_au, depth, tank=900.0, margin_d=60.0):
    """Expected extra propellant [kg] for adding a target whose closest approach to the settled host
    trajectory is `d_au`, to a host carrying `depth` flybys and `tank` kg, inserted `margin_d` days
    from the host's nearest existing flyby.

    kg = (A + B*d**N) * (tank/900)^P * (depth/30)^Q * (max(margin_d,3)/60)^-S

    Fit quality is in model()['cost'] (r2_log, median_abs_error_kg). Valid over the measured support
    d 0.004-0.35 AU, depth 10-48, tank 610-1150 kg; extrapolate with care."""
    c = model()['cost']
    d = max(0.0, float(d_au))
    base = max(0.05, c['A'] + c['B'] * d ** c['N'])
    return (base * (max(float(tank), 300.0) / 900.0) ** c['P_tank']
            * (max(float(depth), 1.0) / 30.0) ** c['Q_depth']
            * (max(float(margin_d), 3.0) / 60.0) ** (-c['S_margin']))


def predict_kg_local(d_au, local_n, tank=900.0, margin_d=60.0):
    """Variant that prices LOCAL congestion (host flybys within +-200 d of the insertion epoch)
    instead of the route's global depth. Fits marginally better than predict_kg; use it when the
    epoch is already chosen, and predict_kg when it is not."""
    c = model()['cost_local_variant']
    d = max(0.0, float(d_au))
    base = max(0.05, c['A'] + c['B'] * d ** c['N'])
    return (base * (max(float(tank), 300.0) / 900.0) ** c['P_tank']
            * (max(float(local_n), 1.0) / c['local_ref']) ** c['Q_local']
            * (max(float(margin_d), 3.0) / 60.0) ** (-c['S_margin']))


def predict_ok(d_au, depth, tank=900.0, margin_d=60.0):
    """Probability that the insertion is feasible at all (aim-point homotopy converges, route
    re-settles). Below ~0.35 an insertion is a long shot; below ~0.15 treat the target as
    unreachable from that route. Accuracy / Brier in model()['success']."""
    w = model()['success']['w']
    z = (w[0] + w[1] * float(d_au) + w[2] * float(depth) / 10.0 + w[3] * float(tank) / 1000.0
         + w[4] * math.log10(max(float(margin_d), 1.0)))
    return 1.0 / (1.0 + math.exp(-max(-40.0, min(40.0, z))))


def cost_J(tank):
    """The contest per-craft objective J_i = 1 + x + x^2, x = (tank - 600)/1400."""
    x = (tank - M_DRY) / 1400.0
    return 1.0 + x + x * x


def predict_dJ(d_au, depth, tank=900.0, margin_d=60.0):
    """Predicted insertion cost in the contest objective: J_i(tank + kg) - J_i(tank)."""
    return cost_J(tank + predict_kg(d_au, depth, tank, margin_d)) - cost_J(tank)


def predict_expected_dJ(d_au, depth, tank=900.0, margin_d=60.0, miss_penalty=1.0):
    """Risk-weighted price of ATTEMPTING the insertion: p*dJ + (1-p)*miss_penalty.
    With miss_penalty = 1.0 (an uncovered target costs +1 in J) this is the number to compare
    against handing the target to a different craft."""
    p = predict_ok(d_au, depth, tank, margin_d)
    return p * predict_dJ(d_au, depth, tank, margin_d) + (1.0 - p) * miss_penalty


def table_kg(d_au, depth):
    """Empirical median kg from the measured (distance x depth) table, or None outside its support.
    A sanity check on predict_kg; it carries no tank or margin correction."""
    for k, v in model()['table_cost']['by_distance_depth'].items():
        db, qb = k.split('|')
        lo, hi = (float(x) for x in db.split('-')); qlo, qhi = (int(x) for x in qb.split('-'))
        if lo <= d_au < hi and qlo <= depth < qhi:
            return v['median_kg']
    return None


def table_ok(d_au, depth):
    """Empirical P(settle) from the measured (distance x depth) table, or None outside its support."""
    for k, v in model()['table_success']['by_distance_depth'].items():
        db, qb = k.split('|')
        lo, hi = (float(x) for x in db.split('-')); qlo, qhi = (int(x) for x in qb.split('-'))
        if lo <= d_au < hi and qlo <= depth < qhi:
            return v['p_ok']
    return None


def ridge_tank(depth):
    """Tank a real route of this depth typically carries (the depth<->tank ridge the fit rides)."""
    de = model()['depth_effect']
    return de['ridge_tank_per_flyby'] * float(depth) + de['ridge_tank_intercept']


def formula():
    m = model()
    return m['cost']['formula'], m['success']['formula']


if __name__ == '__main__':
    m = model()
    print(f"insertion surrogate: {m['n_samples']} trials, {m['n_settled']} settled "
          f"({m['settle_rate']:.0%} feasible); median |err| {m['cost']['median_abs_error_kg']} kg, "
          f"R2(log) {m['cost']['r2_log']} (distance alone {m['cost']['distance_only_r2_log']})")
    print('kg:   ' + m['cost']['formula'])
    print('p_ok: ' + m['success']['formula'])
    de = m['depth_effect']
    print(f"\ndepth: partial exponent Q = {de['exponent_Q_at_fixed_tank']:.3f} at fixed tank, "
          f"tank exponent P = {de['exponent_P_tank']:.3f}; real routes ride "
          f"tank ~ {de['ridge_tank_per_flyby']:.1f}*depth + {de['ridge_tank_intercept']:.0f} kg")
    print(f"  along that ridge: " + ', '.join(f'{k} x{v}' for k, v in de['ridge_cost_multiplier'].items()))
    print(f"  => +10 flybys multiplies the kg price by {de['per_10_flybys']:.2f}")
    print('\n   d_au  depth  tank  margin ->     kg   p_ok      dJ   E[dJ]')
    for d in (0.01, 0.03, 0.06, 0.10, 0.20, 0.30):
        for dep in (15, 30, 42):
            tk = ridge_tank(dep)
            print(f'  {d:5.2f} {dep:6d} {tk:5.0f} {60:7d} -> {predict_kg(d, dep, tk, 60):6.1f} '
                  f'{predict_ok(d, dep, tk, 60):6.2f} {predict_dJ(d, dep, tk, 60):7.4f} '
                  f'{predict_expected_dJ(d, dep, tk, 60):7.4f}')
