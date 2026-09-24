"""Co-design ALL carriers against a PARTITION of the catalogue (stage 7, section 11).

Everything measured so far says the same thing three times over (stage 6 forcing cost, the 09-16 redistribution
campaign, tools/fleet_dissolve.py): a target FORCED onto a finished route costs 0.115-0.148 J, while a target a route
is DESIGNED AROUND costs 0.032-0.042.  Freedom is worth ~3x, so the fleet has to be designed as one object.

Lloyd's algorithm makes that tractable here, because a carrier is only 4 numbers and the DP scores a route in 0.03 s:

  assign : every target goes to the carrier whose orbit it is cheapest for (node gap + phase lag feasibility)
  update : each carrier is RE-OPTIMISED (carrier_opt pattern search) on its own targets ONLY -- the DP value counts
           only assigned targets, so the orbit migrates towards the cluster it serves
  repeat : until the assignment stops changing

The output is k carriers whose preferred target sets are disjoint BY CONSTRUCTION, which is exactly what every greedy
fleet failed to produce.  Growth (carrier_farm/fleet_deepen) then runs on each carrier's own partition.

Usage: carrier_lloyd.py out_dir [--k 8] [--iters 12] [--restarts 4]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pickle, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import carrier_dp as CD, carrier_opt as CO
from ctoc14.constants import DAY, VINF_MAX, AU
from ctoc14.kepler import load_mea, Ephemeris
from ctoc14.search import UNREACHABLE
from exp_carrier import frame, rv2frame, gaps
_S = {}


def car_of(x):
    tL = float(np.clip(x[0], 0.0, 730.0) * DAY); vi = np.array(x[1:4], float)
    nv = np.linalg.norm(vi)
    if nv > VINF_MAX: vi = vi * (VINF_MAX / nv)
    rE, vE = _S['eph'].earth_state(tL)
    return tL, vi, np.array(rE, float), np.array(vE, float) + vi


def gap_row(x):
    """|radial gap| at the mutual node for every target [AU] -- the carrier's geometric affinity."""
    tL, vi, r, v = car_of(x)
    pc, nc, Ec, ac = rv2frame(r[None], v[None])
    g, _ = gaps(pc, nc, Ec, _S['pa'], _S['na'], _S['Ea'])
    return np.abs(g[0]).min(1)


def w_fit(job):
    """re-optimise one carrier on its own targets: DP value counted over `mine` only."""
    x0, mine, seed = job
    a = _S['a']; ids = _S['ids']; elem = _S['elem']
    keep = np.array([int(t) in mine for t in ids])
    if keep.sum() < 5: return x0, -1e9, 0
    sub_ids = ids[keep]; sub_el = elem[keep]
    def val(x):
        tL, vi, r, v = car_of(x)
        car = dict(tL=tL, r=r, v=v)
        try:
            t, lag, gap, ast, Pc = CD.events_for(car, sub_ids, sub_el, a.eps)
            if len(t) < 3: return -1e9, 0
            ac_ = float(rv2frame(r[None], v[None])[3][0])
            f, seq = CD.dp_route(t, lag, gap, ast, tL, a.lam, a.kgap, ac=ac_)
        except Exception:
            return -1e9, 0
        return float(f), len(set(int(ast[i]) for i in seq))
    rng = np.random.default_rng(seed)
    x = np.array(x0, float); f, n = val(x); step = np.array([25.0, 0.5, 0.5, 0.5])
    for it in range(a.fit_iters):
        best = None
        for y in x + step * rng.normal(size=(6, 4)):
            fy, ny = val(y)
            if best is None or fy > best[0]: best = (fy, ny, y)
        if best[0] > f: f, n, x = best; step *= 1.15
        else: step *= 0.75
        if step[0] < 0.5: break
    return x.tolist(), float(f), int(n)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--k', type=int, default=8)
    ap.add_argument('--iters', type=int, default=12); ap.add_argument('--restarts', type=int, default=4)
    ap.add_argument('--fit-iters', type=int, default=60); ap.add_argument('--eps', type=float, default=0.012)
    ap.add_argument('--lam', type=float, default=1.0); ap.add_argument('--kgap', type=float, default=10.0)
    ap.add_argument('--kphase', type=float, default=1.5); ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--seed', type=int, default=0); ap.add_argument('--cap', type=int, default=45)
    ap.add_argument('--lib', default='results/s7/libopt.pkl')
    a = ap.parse_args()
    CD.K_SLOPE = CD.K_PHYS * a.kphase
    eph = Ephemeris(); ids, elem = load_mea()
    keep = np.array([int(x) not in UNREACHABLE for x in ids]); ids = ids[keep]; elem = elem[keep]
    pa, na, Ea = frame(elem[:, :5])
    _S.update(a=a, eph=eph, ids=ids, elem=elem, pa=pa, na=na, Ea=Ea)
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    import run_ialns as RI; say = RI.logger(out / 'log.txt')
    lib = pickle.load(open(a.lib, 'rb'))
    rng = np.random.default_rng(a.seed)
    best_all = None
    with mp.get_context('fork').Pool(a.nproc) as pl:
        for rs in range(a.restarts):
            pick = rng.choice(len(lib), size=a.k, replace=False)
            X = [[lib[int(j)]['tL'] / DAY] + list(lib[int(j)]['vinf']) for j in pick]
            prev = None
            for it in range(a.iters):
                G = np.array(pl.map(gap_row, X))                       # (k, ntarget)
                who = G.argmin(0); dmin = G.min(0)
                part = {c: set(int(ids[j]) for j in np.nonzero(who == c)[0]
                               if dmin[j] < 0.05) for c in range(a.k)}
                for c in range(a.k):                                   # cap: keep the closest --cap targets
                    if len(part[c]) > a.cap:
                        mine = sorted(part[c], key=lambda t: G[c, list(ids).index(t)])[:a.cap]
                        part[c] = set(mine)
                sig = tuple(sorted(tuple(sorted(part[c])) for c in range(a.k)))
                res = pl.map(w_fit, [(X[c], part[c], a.seed * 31 + 7 * it + c) for c in range(a.k)])
                X = [r[0] for r in res]; V = [r[1] for r in res]; N = [r[2] for r in res]
                asg = sum(len(part[c]) for c in range(a.k)); uni = len(set().union(*part.values()))
                say(f'restart {rs} iter {it}: assigned {asg} (distinct {uni}), DP flybys {sorted(N, reverse=True)} '
                    f'= {sum(N)}, value {sum(v for v in V if v > -1e8):.1f}')
                if sig == prev: break
                prev = sig
            tot = sum(N)
            if best_all is None or tot > best_all[0]:
                best_all = (tot, [list(x) for x in X], {c: sorted(part[c]) for c in part}, N)
                json.dump(dict(total=tot, X=best_all[1], part={str(k): v for k, v in best_all[2].items()}, N=N),
                          open(out / 'best.json', 'w'), indent=1)
                say(f'  NEW BEST: {tot} DP flybys over {a.k} co-designed carriers')
    say(f'done: best {best_all[0]} DP flybys, per carrier {sorted(best_all[3], reverse=True)}')


if __name__ == '__main__':
    main()
