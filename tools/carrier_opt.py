"""Optimise the CARRIER ORBIT itself (stage 7 step 1, the one lever never tested).

Every carrier so far was a random sample: launch epoch + v_inf drawn uniformly, then ranked.  The DP route value is a
0.03 s function of 4 numbers (t_launch, v_inf vector), so it can simply be MAXIMISED.  Deeper/cheaper BASES matter far
more than they look: a base of 21 flybys grew to 40 @ 992 kg (0.0339 J/fb, the pool's best route), while bases of 10-13
grow to 28-32 at 0.044-0.050 -- growth inherits the base's efficiency.

Local search: Gaussian pattern search on (tL, vinf) with a shrinking step, restarted from the best random samples.
Objective = the DP value (flybys - lam x dv), which is exactly what a route is worth before insertion growth.

Usage: carrier_opt.py out.pkl [--starts 64] [--iters 120] [--n 60000]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, time, pickle, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import carrier_dp as CD
from ctoc14.constants import DAY, VINF_MAX
from ctoc14.kepler import load_mea, Ephemeris
from ctoc14.search import UNREACHABLE
from exp_carrier import rv2frame
_S = {}


def evaluate(x):
    """x = (tL [d], vinf [km/s] x3) -> (value, nflybys, route)"""
    g = _S; a = g['a']
    tL = float(np.clip(x[0], 0.0, 730.0) * DAY); vi = np.array(x[1:4], float)
    nv = np.linalg.norm(vi)
    if nv > VINF_MAX: vi = vi * (VINF_MAX / nv)
    rE, vE = g['eph'].earth_state(tL); r = np.array(rE, float); v = np.array(vE, float) + vi
    car = dict(tL=tL, r=r, v=v)
    try:
        t, lag, gap, ast, Pc = CD.events_for(car, g['ids'], g['elem'], a.eps)
        if len(t) < 3: return -1e9, 0, None
        ac_ = float(rv2frame(r[None], v[None])[3][0])
        val, seq = CD.dp_route(t, lag, gap, ast, tL, a.lam, a.kgap, ac=ac_)
    except Exception:
        return -1e9, 0, None
    seen = set(); seq = [i for i in seq if not (ast[i] in seen or seen.add(ast[i]))]
    rec = dict(tL=tL, vinf=vi.tolist(), val=float(val),
               base=[(int(ast[i]), float(t[i]), float(lag[i])) for i in seq], neigh={})
    return float(val), len(seq), rec


def w_opt(job):
    x0, seed = job; a = _S['a']
    rng = np.random.default_rng(seed)
    x = np.array(x0, float); f, n, rec = evaluate(x)
    step = np.array([30.0, 0.6, 0.6, 0.6])
    for it in range(a.iters):
        cand = x + step * rng.normal(size=(a.batch, 4))
        best = None
        for y in cand:
            fy, ny, ry = evaluate(y)
            if best is None or fy > best[0]: best = (fy, ny, ry, y)
        if best[0] > f: f, n, rec, x = best[0], best[1], best[2], best[3]; step *= 1.15
        else: step *= 0.75
        if step[0] < 0.4: break
    return f, n, rec


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out')
    ap.add_argument('--n', type=int, default=60000); ap.add_argument('--starts', type=int, default=64)
    ap.add_argument('--iters', type=int, default=120); ap.add_argument('--batch', type=int, default=8)
    ap.add_argument('--eps', type=float, default=0.012); ap.add_argument('--lam', type=float, default=1.0)
    ap.add_argument('--kgap', type=float, default=10.0); ap.add_argument('--kphase', type=float, default=1.5)
    ap.add_argument('--nproc', type=int, default=8); ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--dmax', type=float, default=0.15)
    a = ap.parse_args()
    CD.K_SLOPE = CD.K_PHYS * a.kphase
    eph = Ephemeris(); ids, elem = load_mea()
    keep = np.array([int(x) not in UNREACHABLE for x in ids]); ids = ids[keep]; elem = elem[keep]
    _S.update(a=a, eph=eph, ids=ids, elem=elem)
    rng = np.random.default_rng(a.seed)
    tLs = rng.uniform(0, 730.0, a.n)
    d = rng.normal(size=(a.n, 3)); d /= np.linalg.norm(d, axis=1)[:, None]
    vinf = d * (VINF_MAX * rng.uniform(0, 1, a.n) ** (1 / 3))[:, None]
    X0 = np.column_stack([tLs, vinf])
    tic = time.time()
    with mp.get_context('fork').Pool(a.nproc) as pl:
        vals = np.array(pl.map(_score, [X0[i::a.nproc] for i in range(a.nproc)]), dtype=object)
        sc = np.concatenate([np.asarray(v) for v in vals]); idx = np.concatenate([np.arange(a.n)[i::a.nproc] for i in range(a.nproc)])
        order = idx[np.argsort(-sc)][:a.starts]
        print(f'{a.n} samples scored in {time.time()-tic:.0f} s: best {sc.max():.2f}, median {np.median(sc):.2f}', flush=True)
        res = pl.map(w_opt, [(X0[j], a.seed * 977 + i) for i, j in enumerate(order)])
    res = [r for r in res if r[2]]
    res.sort(key=lambda r: -r[0])
    print(f'optimised {len(res)} carriers in {time.time()-tic:.0f} s: value {res[0][0]:.2f} (was {sc.max():.2f}), '
          f'base flybys best {max(r[1] for r in res)}, median {np.median([r[1] for r in res]):.0f}')
    lib = [r[2] for r in res]
    with mp.get_context('fork').Pool(a.nproc) as pl:
        lib = [c for c in pl.map(_neigh, lib) if c]
    pickle.dump(lib, open(a.out, 'wb'))
    print(f'-> {a.out} ({len(lib)} carriers, base median {np.median([len(c["base"]) for c in lib]):.0f})')


def _neigh(c):
    """closest approach of every target to the carrier's seed trajectory (the farm ranks carriers with this)."""
    import run_ialns as RI
    from ctoc14.globalopt import extend_arc
    from ctoc14.constants import AU, T_MISSION
    g = _S; a = g['a']; eph = g['eph']; ids = g['ids']
    base = c['base']
    try:
        ip = CD.seed_ip(eph, c['tL'], np.array(c['vinf']), np.array([b[1] for b in base]), np.array([b[2] for b in base]))
        extend_arc(ip, T_MISSION - 2 * DAY)
        tt, r = RI.trajectory(ip, step=2 * DAY)
    except Exception:
        return None
    tb = np.array([b[1] for b in base]); free = np.abs(tt[:, None] - tb[None, :]).min(1) > 12 * DAY
    out = {}
    for X in ids:
        ra, _ = eph.ast_states_at(np.full(len(tt), int(X) - 1), tt)
        dd = np.linalg.norm(r - ra, axis=1) / AU; dd[~free] = 9.0
        k = int(dd.argmin())
        if dd[k] < a.dmax: out[int(X)] = (float(dd[k]), float(tt[k]))
    c['neigh'] = out
    return c


def _score(X):
    return [evaluate(x)[0] for x in X]


if __name__ == '__main__':
    main()
