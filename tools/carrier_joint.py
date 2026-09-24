"""Joint design of k routes by Lagrangian pricing (stage 8).

Stage 7's invariant: a target FORCED onto a route costs 0.115-0.148 J, one the route CHOOSES costs 0.032-0.042.
Every fleet construction so far either lets routes choose (deep, overlapping) or forces an assignment (disjoint,
shallow).  Prices are the way out: relax "each target is flown at most once" with multipliers mu_t >= 0, and the
k-route problem separates into k INDEPENDENT single-route DPs whose node gain becomes (1 - mu_t).  No route is ever
forced -- each still takes its own optimum -- but a target two routes want gets dearer until one of them walks away.

    L(mu) = sum_r max_route [ sum_{t in r} (1 - mu_t) - lam dv_r ] + sum_t mu_t        (>= the true optimum)
    subgradient: mu_t <- max(0, mu_t + step (n_t - 1)),  n_t = how many routes took t

This is NOT the set-cover column generation that failed in stage 7 section 10.2: there the duals came from a fixed
pool of columns and were degenerate; here the subproblem is solved EXACTLY for any mu (0.03 s per route), which is the
condition under which subgradient ascent on a Lagrangian is well behaved.  Carriers are re-optimised every
--reopt iterations against the current prices, so the ORBITS also co-adapt (4-parameter pattern search, carrier_opt).

Usage: carrier_joint.py out_dir [--k 9] [--iters 400] [--lib results/s7/libopt.pkl]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pickle, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import carrier_dp as CD
from ctoc14.constants import DAY, VINF_MAX
from ctoc14.kepler import load_mea, Ephemeris
from ctoc14.search import UNREACHABLE
from exp_carrier import rv2frame
_S = {}


def car_of(x):
    tL = float(np.clip(x[0], 0.0, 730.0) * DAY); vi = np.array(x[1:4], float)
    nv = np.linalg.norm(vi)
    if nv > VINF_MAX: vi = vi * (VINF_MAX / nv)
    rE, vE = _S['eph'].earth_state(tL)
    return tL, vi, np.array(rE, float), np.array(vE, float) + vi


def route_of(x, w):
    """best route for carrier x under the price-adjusted weights w -> (value, targets, epochs, lags, tL, vinf)."""
    a = _S['a']
    tL, vi, r, v = car_of(x)
    try:
        t, lag, gap, ast, Pc = CD.events_for(dict(tL=tL, r=r, v=v), _S['ids'], _S['elem'], a.eps)
        if len(t) < 3: return None
        ac_ = float(rv2frame(r[None], v[None])[3][0])
        val, seq = CD.dp_route(t, lag, gap, ast, tL, a.lam, a.kgap, ac=ac_, w=w)
    except Exception:
        return None
    seen = set(); seq = [i for i in seq if not (ast[i] in seen or seen.add(ast[i]))]
    if not seq: return None
    return dict(val=float(val), tL=tL, vinf=vi.tolist(),
                base=[(int(ast[i]), float(t[i]), float(lag[i])) for i in seq])


def w_route(job):
    return route_of(job[0], job[1])


def w_reopt(job):
    """pattern search on (tL, vinf) maximising the PRICED DP value."""
    x0, w, seed = job; a = _S['a']
    rng = np.random.default_rng(seed)
    f0 = route_of(x0, w); x = np.array(x0, float)
    f = f0['val'] if f0 else -1e9
    step = np.array([25.0, 0.5, 0.5, 0.5])
    for it in range(a.reopt_iters):
        best = None
        for y in x + step * rng.normal(size=(6, 4)):
            r = route_of(y, w)
            if r and (best is None or r['val'] > best[0]): best = (r['val'], y)
        if best and best[0] > f: f, x = best; step *= 1.15
        else: step *= 0.75
        if step[0] < 0.5: break
    return x.tolist(), float(f)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--k', type=int, default=9)
    ap.add_argument('--iters', type=int, default=400); ap.add_argument('--lib', default='results/s7/libopt.pkl')
    ap.add_argument('--eps', type=float, default=0.012); ap.add_argument('--lam', type=float, default=1.0)
    ap.add_argument('--kgap', type=float, default=10.0); ap.add_argument('--kphase', type=float, default=1.5)
    ap.add_argument('--step', type=float, default=0.06); ap.add_argument('--step-min', type=float, default=0.004)
    ap.add_argument('--reopt', type=int, default=40); ap.add_argument('--reopt-iters', type=int, default=50)
    ap.add_argument('--nproc', type=int, default=8); ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--mu-max', type=float, default=0.95)
    a = ap.parse_args()
    CD.K_SLOPE = CD.K_PHYS * a.kphase
    eph = Ephemeris(); ids, elem = load_mea()
    keep = np.array([int(x) not in UNREACHABLE for x in ids]); ids = ids[keep]; elem = elem[keep]
    _S.update(a=a, eph=eph, ids=ids, elem=elem)
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    import run_ialns as RI; say = RI.logger(out / 'log.txt')
    lib = pickle.load(open(a.lib, 'rb')); rng = np.random.default_rng(a.seed)

    # start from k carriers that already prefer different targets (greedy max-coverage over the DP bases)
    left = set(int(t) for t in ids); X = []
    for _ in range(a.k):
        j = max(range(len(lib)), key=lambda j: len(left & set(t for t, _, _ in lib[j]['base'])))
        X.append([lib[j]['tL'] / DAY] + list(lib[j]['vinf'])); left -= set(t for t, _, _ in lib[j]['base'])
        lib = lib[:j] + lib[j + 1:]
    mu = np.zeros(301); best = None; step = a.step
    say(f'k={a.k} carriers, {a.iters} subgradient iterations, step {step}')
    with mp.get_context('fork').Pool(a.nproc) as pl:
        for it in range(1, a.iters + 1):
            w = np.clip(1.0 - mu, -0.2, 1.0)
            R = pl.map(w_route, [(x, w) for x in X])
            R = [r for r in R if r]
            cnt = np.zeros(301, int)
            for r in R:
                for t, _, _ in r['base']: cnt[t] += 1
            distinct = int((cnt > 0).sum()); dup = int(np.maximum(cnt - 1, 0).sum())
            tot = sum(len(r['base']) for r in R)
            if best is None or distinct > best[0]:
                best = (distinct, [list(x) for x in X], [dict(r) for r in R])
                json.dump(dict(distinct=distinct, X=best[1],
                               routes=[dict(tL=r['tL'], vinf=r['vinf'], base=r['base']) for r in best[2]]),
                          open(out / 'best.json', 'w'), indent=1)
                say(f'[{it:4d}] distinct {distinct:3d} (total {tot}, dup {dup}), per route '
                    f'{sorted((len(r["base"]) for r in R), reverse=True)}  <-- BEST')
            elif it % 25 == 0:
                say(f'[{it:4d}] distinct {distinct:3d} (total {tot}, dup {dup}), step {step:.4f}, '
                    f'mu>0 {int((mu > 1e-9).sum())}, max mu {mu.max():.2f}')
            # correct subgradient of "each target at most once": a target nobody took must get CHEAPER again
            # (masking this with cnt!=0 lets a priced-out target stay dead forever -- it cost 131 -> 83 coverage)
            mu = np.clip(mu + step * (cnt - 1.0), 0.0, a.mu_max)
            step = max(a.step_min, step * 0.995)
            if a.reopt and it % a.reopt == 0:
                res = pl.map(w_reopt, [(X[c], w, a.seed * 71 + it * 13 + c) for c in range(len(X))])
                X = [r[0] for r in res]
                say(f'       carriers re-optimised against prices (values {[round(r[1],1) for r in res]})')
    say(f'done: best distinct {best[0]} of {len(ids)}')


if __name__ == '__main__':
    main()
