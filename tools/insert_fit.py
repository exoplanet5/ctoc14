"""Fit the insertion-cost surrogate from tools/insert_sample.py rows.

Cost model (multiplicative, fitted in log space; the distance base is a free power law):

    kg = (A + B*d**N) * (tank/900)^P * (depth/30)^Q * (max(margin,3)/60)^-S

A second variant replaces `depth` by LOCAL congestion (flybys of the host within +-200 d of the
insertion epoch), which measures better.

Success model (logistic):

    logit P(ok) = w0 + w1*d + w2*depth/10 + w3*tank/1000 + w4*log10(max(margin,1))

Writes results/s10/insert_model.json.
"""
import os, sys, json, argparse, pathlib
import numpy as np
from scipy.optimize import least_squares, minimize

DAYs = 86400.0
DBUCK = [(0.0, 0.03), (0.03, 0.06), (0.06, 0.10), (0.10, 0.15), (0.15, 0.22), (0.22, 0.35)]
QBUCK = [(0, 20), (20, 30), (30, 38), (38, 60)]
MIN_KG = 1.0
_CACHE = {}


def load(paths):
    rows = []
    for p in paths:
        for line in open(p):
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    seen = {}
    for r in rows:
        seen[(r['dir'], r['host'], r['ast'], round(r['t'], 1))] = r
    return list(seen.values())


def host_tf(r):
    k = (r['dir'], r['host'])
    if k not in _CACHE:
        z = np.load(os.path.join(r['dir'], f"route_{r['host']}.npz"))
        _CACHE[k] = np.sort(np.asarray(z['tf'], float))
    return _CACHE[k]


def feats(R):
    f = {}
    f['d'] = np.array([r['dist'] for r in R])
    f['dep'] = np.array([r['depth'] for r in R], float)
    f['tk'] = np.array([r['tank0'] for r in R])
    f['mg'] = np.clip(np.array([r['margin_d'] for r in R]), 3.0, None)
    f['loc'] = np.array([float(np.sum(np.abs(host_tf(r) - r['t']) < 200 * DAYs)) for r in R])
    gap = []
    for r in R:
        tf = host_tf(r); nxt = tf[tf > r['t']]; prv = tf[tf < r['t']]
        gap.append(((nxt[0] if len(nxt) else r['t'] + 400 * DAYs) -
                    (prv[-1] if len(prv) else r['t'] - 400 * DAYs)) / DAYs)
    f['gap'] = np.array(gap)
    return f


def fit_cost(f, ylog, cong='dep', cong0=30.0, lam_P=0.7):
    """kg = (A + B*d**N) * (tk/900)^P * (cong/cong0)^Q * (mg/60)^-S, least squares in log kg.

    depth and tank are ~0.93 correlated in real routes, so (P, Q) individually are near-unidentified.
    The fit is therefore done in the orthogonal pair (G, P):
        G = Q + P*beta   the SIZE exponent along the real depth<->tank ridge (well determined)
        P                the residual tank effect off the ridge (weakly determined, mildly shrunk)
    where beta = d log(tank) / d log(depth) measured on the sample itself. Q = G - P*beta is recovered
    for the published formula, and G is the number to quote for "how fast does cost grow with depth".
    Returns (params [A,B,N,P,Q,S], model_fn, G, beta)."""
    d, tk, mg, c = f['d'], f['tk'], f['mg'], np.clip(f[cong], 1.0, None)
    u = np.log(c / cong0); v = np.log(tk / 900.0)
    beta = float(np.polyfit(u, v, 1)[0])
    w = v - beta * u

    def fmodel(q):
        A, B, N, P, G, S = q
        base = np.clip(A + B * d ** N, 0.05, None)
        return np.log(base) + P * w + G * u - S * np.log(mg / 60.0)

    best = None
    for n0 in (0.4, 0.7, 1.0, 1.5):
        q0 = np.array([2.0, 200.0, n0, 0.5, 0.8, 0.2])
        lo = np.array([0.0, 1.0, 0.15, -4.0, -3.0, -2.0])
        hi = np.array([60.0, 5000.0, 3.0, 6.0, 6.0, 2.0])
        s = least_squares(lambda q: np.concatenate([fmodel(q) - ylog, [np.sqrt(lam_P) * q[3]]]),
                          q0, bounds=(lo, hi), max_nfev=40000)
        if best is None or s.cost < best.cost:
            best = s
    A, B, N, P, G, S = best.x
    Q = G - P * beta
    par = np.array([A, B, N, P, Q, S])

    def model(p):
        A2, B2, N2, P2, Q2, S2 = p
        base = np.clip(A2 + B2 * d ** N2, 0.05, None)
        return np.log(base) + P2 * np.log(tk / 900.0) + Q2 * np.log(c / cong0) - S2 * np.log(mg / 60.0)

    return par, model, float(G), beta


def fit_ok(f, ok):
    X = np.column_stack([np.ones(len(ok)), f['d'], f['dep'] / 10.0, f['tk'] / 1000.0,
                         np.log10(np.clip(f['mg'], 1.0, None))])

    def nll(w):
        z = X @ w
        return float(np.sum(np.logaddexp(0, z) - ok * z) + 1e-3 * np.sum(w[1:] ** 2))

    w = minimize(nll, np.zeros(X.shape[1]), method='L-BFGS-B').x
    return w, 1.0 / (1.0 + np.exp(-(X @ w)))


def r2(pred, y):
    return float(1 - np.var(pred - y) / np.var(y))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('samples', nargs='+')
    ap.add_argument('--out', default='results/s10/insert_model.json')
    a = ap.parse_args()

    rows = load(a.samples)
    OKR = [r for r in rows if r['ok'] and 'dkg' in r]
    fa = feats(rows); fo = feats(OKR)
    y = np.array([r['dkg'] for r in OKR])
    ylog = np.log(np.clip(y, MIN_KG, None))
    okv = np.array([1.0 if r['ok'] else 0.0 for r in rows])
    print(f'{len(rows)} trials, {len(OKR)} settled ({len(OKR)/len(rows):.1%}); '
          f'dist {fa["d"].min():.3f}-{fa["d"].max():.3f}, depth {int(fa["dep"].min())}-{int(fa["dep"].max())}, '
          f'tank {fa["tk"].min():.0f}-{fa["tk"].max():.0f}, {len({(r["dir"],r["host"]) for r in rows})} distinct hosts')

    # ---------------- cost fit
    p, model, G, beta = fit_cost(fo, ylog)
    A, B, N, P, Q, S = p
    yhat = np.exp(model(p)); ae = np.abs(yhat - y)
    rel = ae / np.clip(y, MIN_KG, None)
    R2 = r2(model(p), ylog)
    print(f'\nCOST FIT   kg = ({A:.3f} + {B:.2f}*d^{N:.3f}) * (tank/900)^{P:.3f} * (depth/30)^{Q:.3f} * (margin/60)^{-S:.3f}')
    print(f'  R2(log) {R2:.3f}   median|err| {np.median(ae):.1f} kg   mean|err| {ae.mean():.1f} kg   '
          f'median rel {np.median(rel):.0%}   p90|err| {np.percentile(ae,90):.0f} kg')
    m0 = float(np.median(y))
    print(f'  baseline constant {m0:.0f} kg      : median|err| {np.median(np.abs(y-m0)):.1f} kg, R2(log) 0.000')
    sd = least_squares(lambda q: np.log(np.clip(q[0] + q[1] * fo['d'] ** q[2], 0.05, None)) - ylog,
                       [2.0, 200.0, 1.0], bounds=([0, 1, 0.15], [60, 5000, 3]), max_nfev=20000).x
    yd = np.clip(sd[0] + sd[1] * fo['d'] ** sd[2], 0.05, None)
    print(f'  distance only                : median|err| {np.median(np.abs(yd-y)):.1f} kg, R2(log) {r2(np.log(yd), ylog):.3f}')

    # variant on local congestion
    p2, model2, G2, beta2 = fit_cost(fo, ylog, cong='loc', cong0=6.0)
    y2 = np.exp(model2(p2)); ae2 = np.abs(y2 - y)
    print(f'  LOCAL-congestion variant     : median|err| {np.median(ae2):.1f} kg, R2(log) {r2(model2(p2), ylog):.3f}, '
          f'Q_local {p2[4]:.3f}  (n within +-200 d, ref 6)')

    # leave-one-HOST-out cross validation (honest error on a route the fit never saw)
    hk = np.array([f'{r["dir"]}:{r["host"]}' for r in OKR])
    cvp = np.zeros(len(OKR))
    for h in np.unique(hk):
        tr = hk != h; te = ~tr
        if tr.sum() < 20:
            cvp[te] = yhat[te]; continue
        ftr = {k: v[tr] for k, v in fo.items()}
        pt, _, _, _ = fit_cost(ftr, ylog[tr])
        At, Bt, Nt, Pt, Qt, St = pt
        base = np.clip(At + Bt * fo['d'][te] ** Nt, 0.05, None)
        cvp[te] = base * (fo['tk'][te] / 900.0) ** Pt * (fo['dep'][te] / 30.0) ** Qt * (fo['mg'][te] / 60.0) ** (-St)
    cvae = np.abs(cvp - y)
    print(f'  leave-one-host-out CV        : median|err| {np.median(cvae):.1f} kg, R2(log) {r2(np.log(np.clip(cvp,1e-3,None)), ylog):.3f}')

    # ---------------- depth effect
    # (a) partial exponent; (b) along the REAL depth<->tank ridge, which is what a designer sees
    rg = np.polyfit(fo['dep'], fo['tk'], 1)
    ridge = lambda dep: float(np.polyval(rg, dep))
    print(f'\nDEPTH EFFECT')
    print(f'  SIZE exponent along the real ridge  G = {G:.3f}  -> cost x (dep2/dep1)^{G:.2f}   (beta = dlog tank/dlog depth = {beta:.3f})')
    print(f'  partial exponent Q = {Q:.3f} at fixed tank (weakly identified: corr(depth,tank) is high)')
    print(f'  tank exponent    P = {P:.3f}; real fleets ride tank ~ {rg[0]:.1f}*depth + {rg[1]:.0f} kg, so depth and tank move together')
    ridge_mult = {}
    for a_, b_ in ((10, 20), (15, 30), (20, 40), (30, 45), (10, 45)):
        m = (b_ / a_) ** G
        ridge_mult[f'{a_}->{b_}'] = round(float(m), 3)
        print(f'  along the ridge depth {a_} (tank {ridge(a_):.0f}) -> {b_} (tank {ridge(b_):.0f}): cost x{m:.2f}')
    per10 = float(np.exp(np.log(ridge_mult['20->40']) / 2.0))
    print(f'  => +10 flybys of host depth multiplies the kg price by about {per10:.2f}')

    # ---------------- the freedom premium, in the contest objective
    BREAK = 0.0369          # J per flyby a fleet must average to reach J < 13 with full coverage
    Jc = lambda t: 1 + (t - 600) / 1400.0 + ((t - 600) / 1400.0) ** 2
    dJ = np.array([Jc(r['tank1']) - Jc(r['tank0']) for r in OKR])
    depa = fa['dep']
    prem = {}
    print(f'\nFREEDOM PREMIUM  (dJ of a FORCED insertion vs the {BREAK} J/flyby a chosen target costs)')
    print(f'  all settled insertions: median dJ {np.median(dJ):.4f} (p25 {np.percentile(dJ,25):.4f}, p75 {np.percentile(dJ,75):.4f})')
    for qlo, qhi in QBUCK:
        m = (fo['dep'] >= qlo) & (fo['dep'] < qhi); ma = (depa >= qlo) & (depa < qhi)
        if m.sum() < 3:
            continue
        under = float(np.mean(dJ[m] <= BREAK))
        # per ATTEMPT, folding in infeasibility: fraction of all trials that both settle and beat breakeven
        good = sum(1 for r in rows if qlo <= r['depth'] < qhi and r['ok'] and Jc(r['tank1']) - Jc(r['tank0']) <= BREAK)
        prem[f'{qlo}-{qhi}'] = dict(n=int(m.sum()), median_dJ=round(float(np.median(dJ[m])), 4),
                                    frac_settled_under_breakeven=round(under, 3),
                                    frac_attempts_under_breakeven=round(good / max(int(ma.sum()), 1), 3),
                                    p_settle=round(float(okv[ma].mean()), 3))
        print(f'  depth {qlo:2d}-{qhi:2d}: n {m.sum():3d}  median dJ {np.median(dJ[m]):.4f}  '
              f'{under:5.0%} of settled beat breakeven; only {good/max(int(ma.sum()),1):5.1%} of ATTEMPTS do '
              f'(P(settle) {okv[ma].mean():.2f})')

    # empirical check inside one distance band
    print('  empirical check (same distance band, shallow vs deep hosts):')
    for lo, hi in DBUCK:
        m = (fo['d'] >= lo) & (fo['d'] < hi)
        sh = m & (fo['dep'] < 25); de = m & (fo['dep'] >= 32)
        if sh.sum() >= 3 and de.sum() >= 3:
            print(f'    d {lo:.2f}-{hi:.2f}: depth<25 median {np.median(y[sh]):5.1f} kg (n{sh.sum():2d})  '
                  f'depth>=32 median {np.median(y[de]):6.1f} kg (n{de.sum():2d})  x{np.median(y[de])/max(np.median(y[sh]),1e-9):.2f}')

    # ---------------- bucketed tables
    tab, stab = {}, {}
    for lo, hi in DBUCK:
        for qlo, qhi in QBUCK:
            m = (fo['d'] >= lo) & (fo['d'] < hi) & (fo['dep'] >= qlo) & (fo['dep'] < qhi)
            if m.sum() >= 2:
                tab[f'{lo:.2f}-{hi:.2f}|{qlo}-{qhi}'] = dict(n=int(m.sum()), median_kg=round(float(np.median(y[m])), 1),
                                                             p25=round(float(np.percentile(y[m], 25)), 1),
                                                             p75=round(float(np.percentile(y[m], 75)), 1))
            ma = (fa['d'] >= lo) & (fa['d'] < hi) & (fa['dep'] >= qlo) & (fa['dep'] < qhi)
            if ma.sum() >= 2:
                stab[f'{lo:.2f}-{hi:.2f}|{qlo}-{qhi}'] = dict(n=int(ma.sum()), p_ok=round(float(okv[ma].mean()), 3))
    dtab, sdtab = {}, {}
    for lo, hi in DBUCK:
        m = (fo['d'] >= lo) & (fo['d'] < hi); ma = (fa['d'] >= lo) & (fa['d'] < hi)
        if m.sum():
            dtab[f'{lo:.2f}-{hi:.2f}'] = dict(n=int(m.sum()), median_kg=round(float(np.median(y[m])), 1),
                                              p25=round(float(np.percentile(y[m], 25)), 1),
                                              p75=round(float(np.percentile(y[m], 75)), 1),
                                              p90=round(float(np.percentile(y[m], 90)), 1))
        if ma.sum():
            sdtab[f'{lo:.2f}-{hi:.2f}'] = dict(n=int(ma.sum()), p_ok=round(float(okv[ma].mean()), 3))
    qtab, sqtab = {}, {}
    for qlo, qhi in QBUCK:
        m = (fo['dep'] >= qlo) & (fo['dep'] < qhi); ma = (fa['dep'] >= qlo) & (fa['dep'] < qhi)
        if m.sum():
            qtab[f'{qlo}-{qhi}'] = dict(n=int(m.sum()), median_kg=round(float(np.median(y[m])), 1))
        if ma.sum():
            sqtab[f'{qlo}-{qhi}'] = dict(n=int(ma.sum()), p_ok=round(float(okv[ma].mean()), 3))

    hdr = '  dist       ' + ''.join(f'{q[0]}-{q[1]:<7}' for q in QBUCK)
    print('\nMEDIAN kg  (distance x host depth)'); print(hdr)
    for lo, hi in DBUCK:
        line = f'  {lo:.2f}-{hi:.2f} '
        for qlo, qhi in QBUCK:
            k = f'{lo:.2f}-{hi:.2f}|{qlo}-{qhi}'
            line += (f'{tab[k]["median_kg"]:>5.0f}({tab[k]["n"]:<3d})' if k in tab else '    -     ')
        print(line)

    # ---------------- success model
    w, pok = fit_ok(fa, okv)
    acc = float(((pok > 0.5) == (okv > 0.5)).mean()); brier = float(np.mean((pok - okv) ** 2))
    print(f'\nSUCCESS FIT  w {np.round(w,3).tolist()}  accuracy {acc:.1%}  Brier {brier:.3f} '
          f'(base rate {okv.mean():.1%})')
    print('  calibration (predicted vs observed, quintiles):')
    o = np.argsort(pok)
    for k in range(5):
        idx = o[k * len(o) // 5:(k + 1) * len(o) // 5]
        print(f'    pred {pok[idx].mean():.2f} -> obs {okv[idx].mean():.2f}  (n {len(idx)})')
    print('\nP(settle)  (distance x host depth)'); print(hdr)
    for lo, hi in DBUCK:
        line = f'  {lo:.2f}-{hi:.2f} '
        for qlo, qhi in QBUCK:
            k = f'{lo:.2f}-{hi:.2f}|{qlo}-{qhi}'
            line += (f'{stab[k]["p_ok"]:>5.2f}({stab[k]["n"]:<3d})' if k in stab else '    -     ')
        print(line)
    print('  P(ok) along the real ridge at d = 0.05 AU, margin 60 d:')
    ridge_ok = {}
    for dep in (12, 20, 30, 40, 46):
        z = w[0] + w[1] * 0.05 + w[2] * dep / 10 + w[3] * ridge(dep) / 1000 + w[4] * np.log10(60.0)
        ridge_ok[dep] = round(float(1 / (1 + np.exp(-z))), 3)
        print(f'    depth {dep:2d} (tank {ridge(dep):4.0f}): {ridge_ok[dep]:.2f}')

    # ---------------- diagnostics
    diag = {}
    for nm in ('d', 'dep', 'tk', 'mg', 'loc', 'gap'):
        v = fo[nm]
        diag[nm] = dict(pearson_log_kg=round(float(np.corrcoef(np.log(np.clip(v, 1e-3, None)), ylog)[0, 1]), 3),
                        spearman_log_kg=round(float(np.corrcoef(np.argsort(np.argsort(v)), np.argsort(np.argsort(ylog)))[0, 1]), 3))
    diag['corr_depth_tank'] = round(float(np.corrcoef(fo['dep'], fo['tk'])[0, 1]), 3)
    diag['corr_depth_local'] = round(float(np.corrcoef(fo['dep'], fo['loc'])[0, 1]), 3)
    print('\nDIAGNOSTICS (correlation with log kg)')
    for k, v in diag.items():
        print(f'  {k}: {v}')

    formula = (f'kg = ({A:.4f} + {B:.4f}*d**{N:.4f}) * (tank/900.0)**{P:.4f} * (depth/30.0)**{Q:.4f} '
               f'* (max(margin_d,3.0)/60.0)**{-S:.4f}')
    formula2 = (f'kg = ({p2[0]:.4f} + {p2[1]:.4f}*d**{p2[2]:.4f}) * (tank/900.0)**{p2[3]:.4f} '
                f'* (local_n/6.0)**{p2[4]:.4f} * (max(margin_d,3.0)/60.0)**{-p2[5]:.4f}')
    formula_ok = ('p_ok = 1/(1+exp(-(' + ' + '.join(
        f'{c:.4f}*{n}' for c, n in zip(w, ['1', 'd', 'depth/10', 'tank/1000', 'log10(max(margin_d,1))'])) + ')))')

    out = dict(
        n_samples=len(rows), n_settled=len(OKR), settle_rate=round(len(OKR) / len(rows), 4),
        what='real post-settle tank increase from inserting ONE target into a finished impulsive route '
             '(run_ialns.w_insert = globalopt.insert_homotopy + full restore/optimise/restore)',
        features=dict(d_au='closest approach of the target to the host trajectory [AU] (run_ialns.w_cands)',
                      depth='flybys already on the host route',
                      tank='host tank [kg] before the insertion',
                      margin_d="days from the insertion epoch to the host's nearest existing flyby",
                      local_n='(variant) host flybys within +-200 d of the insertion epoch'),
        cost=dict(A=float(A), B=float(B), N=float(N), P_tank=float(P), Q_depth=float(Q), S_margin=float(S),
                  formula=formula, r2_log=round(R2, 4),
                  median_abs_error_kg=round(float(np.median(ae)), 2),
                  mean_abs_error_kg=round(float(ae.mean()), 2),
                  p90_abs_error_kg=round(float(np.percentile(ae, 90)), 2),
                  median_rel_error=round(float(np.median(rel)), 4), min_kg_floor=MIN_KG,
                  baseline_constant_median_abs_error_kg=round(float(np.median(np.abs(y - m0))), 2),
                  distance_only_r2_log=round(r2(np.log(yd), ylog), 4),
                  distance_only_median_abs_error_kg=round(float(np.median(np.abs(yd - y))), 2),
                  cv_median_abs_error_kg=round(float(np.median(cvae)), 2),
                  cv_r2_log=round(r2(np.log(np.clip(cvp, 1e-3, None)), ylog), 4)),
        cost_local_variant=dict(A=float(p2[0]), B=float(p2[1]), N=float(p2[2]), P_tank=float(p2[3]),
                                Q_local=float(p2[4]), S_margin=float(p2[5]), formula=formula2,
                                r2_log=round(r2(model2(p2), ylog), 4),
                                median_abs_error_kg=round(float(np.median(ae2)), 2),
                                local_ref=6.0),
        success=dict(w=[float(x) for x in w], formula=formula_ok, accuracy=round(acc, 4), brier=round(brier, 4),
                     base_rate=round(float(okv.mean()), 4),
                     terms=['1', 'd', 'depth/10', 'tank/1000', 'log10(max(margin_d,1))'],
                     note='every observed failure was an aim-point HOMOTOPY divergence (miss > 1e4 km); '
                          'no insertion that reached the settle stage failed to settle',
                     p_ok_along_ridge_d005={str(k): v for k, v in ridge_ok.items()}),
        depth_effect=dict(size_exponent_G=float(G), beta_dlogtank_dlogdepth=float(beta),
                          exponent_Q_at_fixed_tank=float(Q), exponent_P_tank=float(P),
                          ridge_tank_per_flyby=float(rg[0]), ridge_tank_intercept=float(rg[1]),
                          ridge_cost_multiplier=ridge_mult, per_10_flybys=round(per10, 3),
                          median_kg_by_depth=qtab, p_ok_by_depth=sqtab),
        freedom_premium=dict(breakeven_dJ_per_flyby=BREAK, by_depth=prem,
                             median_dJ_all=round(float(np.median(dJ)), 4),
                             note='a target a route CHOOSES costs ~0.033-0.037 J; this table is what a '
                                  'target FORCED onto a finished route costs, by host depth'),
        table_cost=dict(by_distance=dtab, by_distance_depth=tab),
        table_success=dict(by_distance=sdtab, by_distance_depth=stab),
        diagnostics=diag,
        data=dict(dist_range=[round(float(fa['d'].min()), 4), round(float(fa['d'].max()), 4)],
                  depth_range=[int(fa['dep'].min()), int(fa['dep'].max())],
                  tank_range=[round(float(fa['tk'].min()), 1), round(float(fa['tk'].max()), 1)],
                  margin_range=[round(float(fa['mg'].min()), 1), round(float(fa['mg'].max()), 1)],
                  n_hosts=len({(r['dir'], r['host']) for r in rows}),
                  host_dirs=sorted({r['dir'] for r in rows}), sources=list(a.samples)),
    )
    pathlib.Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(a.out, 'w'), indent=1)
    print(f'\nwrote {a.out}\n  {formula}\n  {formula_ok}')


if __name__ == '__main__':
    main()
