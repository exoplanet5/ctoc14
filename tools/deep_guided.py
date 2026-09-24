"""SURROGATE-GUIDED deepening: spend the insertion budget on the candidates the fitted model says will pay.

Why this exists.  `tools/fleet_deepen.py` ranks the untaken targets by CLOSEST APPROACH alone and tries the ten
nearest.  The insertion surrogate built this run (results/s10/insert_model.json, 1225 real settles) says that is the
wrong ranking on two counts:

  * only 47 % of stratified candidates converge at all, and P(converge) collapses with HOST DEPTH and with
    closeness of the epoch to an existing flyby -- exactly the regime a deepener works in.  Every diverged
    homotopy costs the same wall time as a successful one, so a deepener that ignores p_ok throws away half its
    budget.  `predict_ok` is calibrated (Brier 0.121, quintiles 0.03/0.17/0.44/0.74/0.96 vs 0.03/0.16/0.42/0.78/0.96).
  * the price is NOT a function of distance alone (R2 0.32) but of distance x depth x margin (R2 0.65); the
    measured medians at d < 0.06 AU are 45 kg into a crowded leg (margin < 20 d) against 7 kg into an empty one
    (margin > 100 d).  Distance-first ranking therefore picks the crowded, expensive insertion first.

So this deepener ranks by `predict_expected_dJ` = p*dJ + (1-p)*miss_penalty, drops anything below --pmin, and accepts
a settle only if the REAL tank step is under the contest breakeven (an added flyby is worth 0.0369 J fleet-wide, which
at tank T is 0.0369*1400/(1+2x) kg, x = (T-600)/1400 -- 34 kg at 950 kg, 30 kg at 1100 kg).  The candidate list is
re-derived from the settled trajectory every --recand acceptances, not every acceptance, because `w_cands` (the
trajectory sample against all 298 catalogue orbits) is a third of the wall time.

The verdict is still the twin: every accepted step is a real `run_ialns.w_insert` (aim-point homotopy + whole-route
re-settle) and the reported tank is the settled tank.

Usage: deep_guided.py fleet_dir out_dir [--tmax 1500] [--extra 20] [--pmin 0.3] [--nproc 4] [--breakeven 0.0369]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from insert_model import predict_kg, predict_ok, predict_dJ, predict_expected_dJ
from ctoc14.search import UNREACHABLE
from ctoc14.constants import DAY

ALL = [t for t in range(1, 301) if t not in UNREACHABLE]
_A = None


def breakeven_kg(tank, dJ_target):
    """kg step whose dJ equals dJ_target (J_i = 1 + x + x^2, x = (T-600)/1400)."""
    x = (tank - 600.0) / 1400.0
    return dJ_target * 1400.0 / (1.0 + 2.0 * x)


def rank(cands, st, tank, pmin, dJ_target, cap, avoid=()):
    """Best (ast, t) first by risk-weighted price; one entry per target, the cheapest epoch for it."""
    tf = np.asarray(st['tf'], float)
    have = set(int(x) for x in st['asts'])          # the candidate list may be stale (--recand > 1)
    depth = len(have)
    best = {}
    for c in cands:
        if c['ast'] in avoid or c['ast'] in have: continue
        m = float(np.abs(tf - c['t']).min() / DAY)
        if m < _A.gap: continue
        p = predict_ok(c['dist'], depth, tank, m)
        if p < pmin: continue
        kg = predict_kg(c['dist'], depth, tank, m)
        if kg > cap: continue
        e = p * (predict_dJ(c['dist'], depth, tank, m)) + (1.0 - p) * _A.miss_penalty
        cur = best.get(c['ast'])
        if cur is None or e < cur[0]: best[c['ast']] = (e, c, p, kg, m)
    return sorted(best.values(), key=lambda x: x[0])


def w_deep(job):
    n, st, tank, avoid = job; a = _A
    t0 = time.time(); added = []; tried = set(); nfail = 0; cands = None; since = 99
    trace = []
    while len(added) < a.extra and time.time() - t0 < a.tmax:
        if since >= a.recand:
            have = set(int(x) for x in st['asts'])
            left = [t for t in ALL if t not in have]
            if not left: break
            cands = RI.w_cands((n, st, left, a.dmax)); since = 0
        cap = min(a.max_step, breakeven_kg(tank, a.breakeven))
        R = [x for x in rank(cands, st, tank, a.pmin, a.breakeven, cap, avoid) if (x[1]['ast'], round(x[1]['t'] / DAY)) not in tried]
        if not R:
            if since == 0: break
            since = 99; continue
        got = False
        for e, c, p, kg, m in R[:a.ntrial]:
            tried.add((c['ast'], round(c['t'] / DAY)))
            r = RI.w_insert((n, st, c['ast'], c['t']))
            ok = bool(r.get('ok')) and r['tank'] <= a.max_tank
            real = (r['tank'] - tank) if ok else None
            trace.append(dict(ast=int(c['ast']), d=float(c['dist']), depth=len(st['asts']), tank=float(tank),
                              margin=float(m), p=float(p), pred_kg=float(kg), ok=bool(ok),
                              real_kg=(float(real) if ok else None)))
            if ok and real <= cap:
                st, tank = r['st'], r['tank']; added.append(int(c['ast'])); got = True; since += 1; nfail = 0
                break
            if ok:                       # settled but too dear: undo (keep the cheaper route)
                continue
        if not got:
            nfail += 1
            if nfail >= a.nfail: break
            since = 99
    return n, st, float(tank), added, time.time() - t0, trace


def main():
    global _A
    ap = argparse.ArgumentParser(); ap.add_argument('src'); ap.add_argument('out')
    ap.add_argument('--max-step', type=float, default=45.0); ap.add_argument('--max-tank', type=float, default=1250.0)
    ap.add_argument('--breakeven', type=float, default=0.0369, help='J a flyby is worth fleet-wide')
    ap.add_argument('--extra', type=int, default=24); ap.add_argument('--dmax', type=float, default=0.25)
    ap.add_argument('--ntrial', type=int, default=6); ap.add_argument('--gap', type=float, default=10.0)
    ap.add_argument('--pmin', type=float, default=0.30); ap.add_argument('--miss-penalty', type=float, default=1.0)
    ap.add_argument('--nfail', type=int, default=5); ap.add_argument('--recand', type=int, default=3)
    ap.add_argument('--tmax', type=float, default=1500.0); ap.add_argument('--nproc', type=int, default=4)
    ap.add_argument('--only', default='', help='comma separated route names to deepen')
    ap.add_argument('--avoid', default='', help='JSON file: {route_name: [targets to skip]} or [targets to skip]')
    a = ap.parse_args(); _A = a
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt')
    fl = RI.IFleet(a.src)
    names = [x for x in sorted(fl.routes) if not a.only or x in a.only.split(',')]
    say(f'guided deepening of {len(names)} routes of {a.src} (covered {len(fl.coverage())}), pmin {a.pmin}, '
        f'breakeven {a.breakeven} J = {breakeven_kg(950, a.breakeven):.0f} kg at 950 kg')
    AV = json.load(open(a.avoid)) if a.avoid else {}
    av = (lambda x: set(AV)) if isinstance(AV, list) else (lambda x: set(AV.get(x, [])))
    jobs = [(x, fl.routes[x]['st'], fl.routes[x]['tank'], av(x)) for x in names]
    ref = {j[0]: (len(j[1]['asts']), j[2]) for j in jobs}
    if AV: say(f'avoid list active ({len(AV)} entries)')
    new = RI.IFleet(); rows = []; TR = []
    with mp.get_context('fork').Pool(a.nproc) as pl:
        for n, st, tank, added, dt, trace in pl.imap_unordered(w_deep, jobs):
            n0, t0 = ref[n]; nf = len(set(int(x) for x in st['asts']))
            new.routes[n] = dict(st=st, tank=tank)
            np.savez(out / f'route_{n}.npz', **st)
            say(f'  {n}: {n0} -> {nf} fb, {t0:.0f} -> {tank:.0f} kg, {RI.cost(t0)/n0:.4f} -> {RI.cost(tank)/nf:.4f} J/fb'
                f'  (+{len(added)}: {added})  ({dt:.0f} s, {len(trace)} trials, {sum(1 for x in trace if x["ok"])} settled)')
            rows.append(dict(name=n, n0=n0, n=nf, tank0=t0, tank=tank, Jfb=RI.cost(tank) / nf, added=added))
            TR += [dict(host=n, **x) for x in trace]
    new.save(out / 'fleet', note='surrogate-guided deepening')
    tot = sum(len(set(int(x) for x in r['st']['asts'])) for r in new.routes.values())
    say(f'deepened fleet: {len(new.routes)} routes, covered {len(new.coverage())}, sum J_i '
        f'{sum(RI.cost(r["tank"]) for r in new.routes.values()):.3f}, flybys {tot} (surplus {tot - len(new.coverage())})')
    if rows:
        J = np.array([r['Jfb'] for r in rows])
        say(f'J/fb best {J.min():.4f} median {np.median(J):.4f}; bar deep-fixed 0.0339, fleet target 0.0369')
    json.dump(rows, open(out / 'deep.json', 'w'), indent=1)
    json.dump(TR, open(out / 'trials.jsonl', 'w'))
    ok = [x for x in TR if x['ok']]
    if ok:
        err = np.array([abs(x['real_kg'] - x['pred_kg']) for x in ok])
        say(f'surrogate check: {len(ok)}/{len(TR)} settled (model said {np.mean([x["p"] for x in TR]):.2f}), '
            f'median |pred - real| {np.median(err):.1f} kg')


if __name__ == '__main__':
    main()
