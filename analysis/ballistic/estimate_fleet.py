"""
Task 3: first-order fleet estimate from the ballistic scan statistics.

Inputs (analysis/ballistic/out/):
  earth_to_asteroid_summary.csv   (scan_earth_to_ast.py + refine_unreachable.py)
  junction_hist.npz               (scan_junction.py)

Leg model.  After a flyby the spacecraft is on some heliocentric arc; scan_junction.py measured, for each
sampled incoming arc, the shortest grid TOF (60..600 d) for which SOME catalogue asteroid C can be reached with a
junction delta-v  min_C |v_out - v_in| <= eta * a * TOF  (a = 0.25 / 0.5 / 0.8 mm/s^2), and the delta-v actually
used on that leg.  From the histograms tmin_hist / dvused_hist we take, per (class, eta, a):
    E[TOF_leg], E[dv_leg], P(no hop within 600 d)
and interpolate them linearly in log(a) to the spacecraft's current acceleration a = Tmax / m.
Arcs with no feasible hop within 600 d are charged TOF_NOHOP days and the mean dv (they are rare for a >= 0.5).

Spacecraft model.  m0 = 600 + m_fuel kg, T = 0.5 N, Isp = 4000 s.  First leg: ballistic Earth -> asteroid
(v_inf <= 4 km/s), duration = median over launch dates of the shortest TOF that reaches SOME asteroid with
v_inf <= 4 km/s (from the Earth-scan v_inf grid).  Each subsequent leg costs
    dv_lt = LT_PENALTY * dv_impulsive         (low-thrust arc vs impulse at the junction)
    dm    = m (1 - exp(-dv_lt / v_e)),  v_e = Isp g0 = 39.227 km/s
and takes E[TOF_leg](a).  The chain stops when the window (15 yr from launch at t = 0) or the fuel runs out.

Outputs: fleet_estimate.csv (one row per fuel load / class / eta / penalty), fig_fleet_estimate.png, and a text report.
"""
import os, sys
import numpy as np
import pandas as pd

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
DAY = 86400.0
YEAR = 365.25
T_WINDOW_D = 5478.75
ISP, G0 = 4000.0, 9.80665e-3
VE = ISP * G0                     # km/s
TMAX = 0.5                        # N
M_DRY = 600.0
TOF_NOHOP = 730.0                 # d charged to incoming arcs with no feasible hop within 600 d
FUEL_LOADS = [0, 100, 200, 350, 500, 700, 1000, 1400]
LT_PENALTIES = [1.0, 1.25, 1.5]
N_TARGETS_REACHABLE = None        # filled from the Earth scan


def cost(mfuel):
    x = mfuel / 1400.0
    return 1 + x + x * x


def leg_stats(d):
    """Return arrays E_tof[class, eta, acc], E_dv[class, eta, acc], P_nohop[class, eta, acc] from the histograms."""
    tmin, dvu = d['tmin_hist'], d['dvused_hist']            # (2, nEta, nA, nT+1), (2, nEta, nA, NB)
    tof = d['tof_days']
    dvb = float(d['dv_bin'])
    NB = dvu.shape[-1]
    edges_dv = (np.arange(NB) + 0.5) * dvb
    n = tmin.sum(axis=-1)
    p_none = tmin[..., -1] / np.maximum(n, 1)
    hf = tmin[..., :-1]
    e_tof_feas = (hf * tof).sum(-1) / np.maximum(hf.sum(-1), 1)
    e_dv = (dvu * edges_dv).sum(-1) / np.maximum(dvu.sum(-1), 1)
    e_tof = (1 - p_none) * e_tof_feas + p_none * TOF_NOHOP
    return e_tof, e_tof_feas, e_dv, p_none


def interp_loga(vals, accels, a):
    """vals[..., nA] interpolated (linear in log a) to acceleration a, clamped to the tabulated range."""
    la = np.log(np.clip(a, accels[0], accels[-1]))
    return np.array([np.interp(la, np.log(accels), v) for v in vals.reshape(-1, len(accels))]).reshape(vals.shape[:-1])


def simulate(mfuel, e_tof, e_dv, p_none, accels, first_leg_d, lt_penalty, t_launch_d=0.0):
    """Expected-value chain: returns (n_flybys, t_end_d, m_end, dv_total, legs list)."""
    m = M_DRY + mfuel
    t = t_launch_d + first_leg_d
    n = 1 if t <= T_WINDOW_D else 0
    dv_tot = 0.0
    legs = []
    while True:
        a = TMAX / m * 1000.0                          # N/kg = m/s^2 -> mm/s^2 (0.25 at 2000 kg, 0.83 at 600 kg)
        tof = float(interp_loga(e_tof, accels, a))
        dv = float(interp_loga(e_dv, accels, a)) * lt_penalty
        if mfuel <= 0:
            break                                      # a ballistic spacecraft cannot make a junction (0 dv budget)
        dm = m * (1 - np.exp(-dv / VE))
        if m - dm < M_DRY or t + tof > T_WINDOW_D:
            break
        m -= dm
        t += tof
        dv_tot += dv
        n += 1
        legs.append((t, m, a, tof, dv))
    return n, t, m, dv_tot, legs


def main():
    d = np.load(os.path.join(OUT, 'junction_hist.npz'))
    df = pd.read_csv(os.path.join(OUT, 'earth_to_asteroid_summary.csv'))
    vcol = 'min_vinf_refined' if 'min_vinf_refined' in df else 'min_vinf'
    reach = df[vcol] <= 4.0
    n_reach = int(reach.sum())
    first_leg_cheapest = float(df.loc[reach, 'best_tof_day'].median())
    # first leg of a fuelled spacecraft: the shortest TOF (over all asteroids) with v_inf <= 4 km/s from a given launch date
    with np.load(os.path.join(OUT, 'earth_to_asteroid_vinf_grid.npz')) as g:
        feas = g['vinf'].astype(np.float32) <= 4.0                     # (300, nL, nT)
        tofg = g['tof_days']
    anyA = feas.any(axis=0)                                             # (nL, nT): some asteroid reachable in this cell
    tmin_L = np.where(anyA.any(axis=1), tofg[np.argmax(anyA, axis=1)], np.nan)
    anyL = feas.any(axis=1)                                             # (300, nT)
    tmin_A = np.where(anyL.any(axis=1), tofg[np.argmax(anyL, axis=1)], np.nan)
    first_leg = float(np.nanmedian(tmin_L))
    print(f"first Earth->asteroid leg: shortest feasible TOF from a launch date (any asteroid): median {first_leg:.0f} d "
          f"(10/90 %: {np.nanpercentile(tmin_L, 10):.0f}/{np.nanpercentile(tmin_L, 90):.0f} d); shortest feasible TOF per asteroid: median "
          f"{np.nanmedian(tmin_A):.0f} d (10/90 %: {np.nanpercentile(tmin_A, 10):.0f}/{np.nanpercentile(tmin_A, 90):.0f} d); "
          f"cheapest-v_inf cell TOF: median {first_leg_cheapest:.0f} d (used only for the ballistic, m_fuel = 0 spacecraft)")
    classes = [str(c) for c in d['classes']]
    etas = d['etas']; accels = d['accels']
    e_tof, e_tof_feas, e_dv, p_none = leg_stats(d)
    print("=== Task 3: first-order fleet estimate ===")
    print(f"reachable asteroids (min v_inf <= 4 km/s, column {vcol}): {n_reach}/300")
    print("leg statistics from scan_junction.py (shortest feasible hop, TOF grid 60..600 d):")
    print("  class | eta | a[mm/s2] | E[TOF_leg] incl. no-hop [d] | E[TOF | hop] [d] | E[dv_junction] [km/s] | P(no hop <=600 d) | flybys/yr | dv/yr [km/s]")
    for ic, cl in enumerate(classes):
        for ie, eta in enumerate(etas):
            for ia, a in enumerate(accels):
                print(f"  {cl:9s} | {eta:.1f} | {a:4.2f} | {e_tof[ic,ie,ia]:6.1f} | {e_tof_feas[ic,ie,ia]:6.1f} | {e_dv[ic,ie,ia]:5.2f} | {p_none[ic,ie,ia]*100:5.1f}% | "
                      f"{YEAR/e_tof[ic,ie,ia]:5.2f} | {e_dv[ic,ie,ia]*YEAR/e_tof[ic,ie,ia]:5.2f}")
    rows = []
    for ic, cl in enumerate(classes):
        for ie, eta in enumerate(etas):
            for pen in LT_PENALTIES:
                for mf in FUEL_LOADS:
                    n, t, m, dvt, legs = simulate(mf, e_tof[ic, ie], e_dv[ic, ie], p_none[ic, ie], accels,
                                                  first_leg if mf > 0 else first_leg_cheapest, pen)
                    limit = 'ballistic'
                    if mf > 0:                         # which constraint stopped the chain?
                        a_end = TMAX / m * 1000
                        tof_next = float(interp_loga(e_tof[ic, ie], accels, a_end))
                        limit = 'time' if t + tof_next > T_WINDOW_D else 'fuel'
                    J = cost(mf)
                    rows.append(dict(cls=cl, eta=float(eta), lt_penalty=pen, m_fuel=mf, J=J, n_flybys=n, cost_per_flyby=J / max(n, 1),
                                     t_end_yr=t / YEAR, m_end=m, dv_used=dvt, limit=limit,
                                     n_sc_for_catalog=np.ceil(n_reach / max(n, 1)), J_total=np.ceil(n_reach / max(n, 1)) * J + (300 - n_reach)))
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUT, 'fleet_estimate.csv'), index=False, float_format='%.4g')
    print("\nexpected-value chain per spacecraft (launch at t=0, first leg ballistic), class 'plausible':")
    print("  eta | LT pen. | m_fuel [kg] |  J_i  | flybys | limited by | t_end [yr] | dv used [km/s] | J_i/flyby | N_sc for catalog | J_total (N_sc*J_i + misses)")
    for _, r in res[res.cls == 'plausible'].iterrows():
        print(f"  {r.eta:.1f} | {r.lt_penalty:4.2f} | {r.m_fuel:6.0f} | {r.J:5.2f} | {r.n_flybys:4d} | {r.limit:9s} | {r.t_end_yr:5.1f} | {r.dv_used:5.1f} | "
              f"{r.cost_per_flyby:6.3f} | {r.n_sc_for_catalog:4.0f} | {r.J_total:6.1f}")
    print("\nsame for class 'earthlike' (Earth-launch-like arcs only), eta=0.8, LT penalty 1.25:")
    for _, r in res[(res.cls == 'earthlike') & (np.isclose(res.eta, 0.8)) & (np.isclose(res.lt_penalty, 1.25))].iterrows():
        print(f"  m_fuel {r.m_fuel:6.0f} kg: J_i {r.J:5.2f}, flybys {r.n_flybys:3d} ({r.limit}), t_end {r.t_end_yr:4.1f} yr, dv {r.dv_used:5.1f} km/s, "
              f"J_i/flyby {r.cost_per_flyby:.3f}, N_sc {r.n_sc_for_catalog:.0f}, J_total {r.J_total:.0f}")
    # detail of the full-tank chain (plausible, eta 0.8, penalty 1.25)
    ic = classes.index('plausible'); ie = int(np.argmin(np.abs(etas - 0.8)))
    n, t, m, dvt, legs = simulate(1400, e_tof[ic, ie], e_dv[ic, ie], p_none[ic, ie], accels, first_leg, 1.25)
    print(f"\nfull-tank chain detail (plausible, eta 0.8, penalty 1.25): {n} flybys, ends t={t/YEAR:.1f} yr with m={m:.0f} kg, dv={dvt:.1f} km/s")
    yrs = np.array([l[0] for l in legs]) / YEAR
    for y in range(1, 16):
        k = int((yrs <= y).sum())
        if k == 0 or y > 15:
            continue
        if y in (1, 2, 3, 5, 8, 10, 15) or k == n - 1:
            l = legs[k - 1]
            print(f"   by year {y:2d}: {k+1:3d} flybys, m={l[1]:.0f} kg, a={l[2]:.2f} mm/s2, current leg TOF {l[3]:.0f} d, dv/leg {l[4]:.2f} km/s")
    budget_table(d, classes, first_leg)
    plots(res, classes, etas)


def budget_table(d, classes, first_leg):
    """Complementary estimate: a planner that fixes the leg TOF T and pays the typical best-next-hop junction dv(T)
    (median / 25th percentile of min_C dv over the sampled incoming arcs).  Flybys per full-tank spacecraft
    n = min( dv_budget / (pen * dv(T)),  (15 yr - first leg) / T ),  dv_budget = 47.2 km/s; the thrust
    capability constraint dv(T) <= eta * a * T is checked with a = 0.25 mm/s^2 (worst case, full tank)."""
    H = d['H_min']; tof = d['tof_days']; dvb = float(d['dv_bin']); NB = H.shape[-1]
    edges = (np.arange(NB) + 1) * dvb
    dv_budget = VE * np.log(2000.0 / 600.0)
    print(f"\n=== Complementary estimate: fixed leg TOF T, junction dv = typical best-next-hop dv(T) (dv budget {dv_budget:.1f} km/s, 15 yr window, first leg {first_leg:.0f} d) ===")
    print("  class | T [d] | dv(T) median / 25% [km/s] | eta*a*T at a=0.25 [km/s] | flybys per full tank: dv-limited (pen 1.0 / 1.25, median dv) | time-limited | n = min (pen 1.25, median) | n (pen 1.25, 25% dv) | J_i/flyby | N_sc for 298")
    rows = []
    for ic, cl in enumerate(classes):
        for it, T in enumerate(tof):
            h = H[ic, it]
            if h.sum() == 0:
                continue
            c = np.cumsum(h) / h.sum()
            med = edges[min(np.searchsorted(c, 0.5), NB - 1)]
            q25 = edges[min(np.searchsorted(c, 0.25), NB - 1)]
            cap = 0.8 * 0.25e-6 * T * DAY
            n_dv10 = dv_budget / med
            n_dv125 = dv_budget / (1.25 * med)
            n_time = (T_WINDOW_D - first_leg) / T
            n = int(min(n_dv125, n_time)) + 1            # +1: the first (ballistic) flyby
            n25 = int(min(dv_budget / (1.25 * q25), n_time)) + 1
            rows.append(dict(cls=cl, T=T, dv_med=med, dv_q25=q25, cap=cap, n=n, n25=n25))
            if T in (120, 180, 240, 300, 365, 450, 600):
                print(f"  {cl:9s} | {T:4.0f} | {med:5.2f} / {q25:5.2f} | {cap:5.2f}{' (!)' if med > cap else '    '} | {n_dv10:5.1f} / {n_dv125:5.1f} | {n_time:5.1f} | {n:3d} | {n25:3d} | {3.0/n:.3f} | {np.ceil(298/n):3.0f}")
    b = pd.DataFrame(rows)
    b.to_csv(os.path.join(OUT, 'fleet_budget_table.csv'), index=False, float_format='%.4g')
    for cl in classes:
        bb = b[b.cls == cl]
        k = bb.n.idxmax()
        print(f"  best fixed-TOF chain, class {cl}: T = {bb['T'][k]:.0f} d, {bb.n[k]} flybys per full-tank spacecraft (median dv), "
              f"{bb.n25[bb.n25.idxmax()]} with 25%-ile dv at T = {bb['T'][bb.n25.idxmax()]:.0f} d")


def plots(res, classes, etas):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for cl, ls in zip(classes, ('-', '--')):
        for pen, col in zip(LT_PENALTIES, ('#4C72B0', '#55A868', '#C44E52')):
            r = res[(res.cls == cl) & np.isclose(res.eta, 0.8) & np.isclose(res.lt_penalty, pen)].sort_values('m_fuel')
            ax[0].plot(r.m_fuel, r.n_flybys, ls=ls, marker='o', ms=3, color=col, label=f'{cl}, LT penalty {pen}')
            ax[1].plot(r.m_fuel, r.cost_per_flyby, ls=ls, marker='o', ms=3, color=col, label=f'{cl}, LT penalty {pen}')
    ax[0].set_xlabel('fuel load [kg]'); ax[0].set_ylabel('expected flybys per spacecraft (15 yr)'); ax[0].grid(alpha=.3); ax[0].legend(fontsize=7)
    ax[0].set_title('eta = 0.8; leg model from junction scan')
    ax[1].set_xlabel('fuel load [kg]'); ax[1].set_ylabel('cost J_i per flyby'); ax[1].grid(alpha=.3); ax[1].set_ylim(0, 1.0)
    ax[1].set_title('J_i = 1 + x + x^2, x = m_fuel/1400')
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'fig_fleet_estimate.png'), dpi=130)
    plt.close('all')


if __name__ == '__main__':
    main()
