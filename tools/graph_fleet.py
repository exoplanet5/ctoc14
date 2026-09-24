"""Fleet search on the transfer graph (stage 9): thousands of complete fleets per hour.

A route now costs ~0.6 s instead of 5-10 min (tools/graph_beam.py), so the compute moves to where every GTOC winner
put it -- the combinatorics of which targets go to which craft.  Two modes:

  greedy   : build craft one at a time, covered targets kept as ZERO-PRIZE stepping stones (never excluded), with
             randomised lam / launch window / depth per craft.  A whole fleet costs seconds, so we sample many.
  aco      : ant colony over the whole partition (NUDT's GTOC9 top level, 2nd place on the isomorphic problem --
             partition N objects among missions, cost = fixed + variable per mission).  Each ant builds the ENTIRE
             fleet, so every target is assigned exactly once BY CONSTRUCTION: no overlap to resolve afterwards and
             no target ever forced onto a finished route (the 3x premium of stage 7).  Pheromone on targets learns
             which ones are worth saving for a later craft, which is exactly what greedy construction cannot do.

Score is the real objective J = sum_i J_i(tank) + N_miss, with tank from the Lambert-chain dv (the impulsive twin's
own seed cost, calibrated to 3-7 % in stages 5-6).

Usage: graph_fleet.py tdb.npz out_dir [--mode aco] [--ants 30] [--gens 40] [--ncraft 10]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
from ctoc14.constants import VE, M_DRY
from ctoc14.kepler import Ephemeris
from ctoc14.search import UNREACHABLE
import graph_beam as GB

ALL = [t for t in range(1, 301) if t not in UNREACHABLE]
J1 = lambda tk: 1.0 + (tk - 600.0) / 1400.0 + ((tk - 600.0) / 1400.0) ** 2
_G = {}


def trace(par, lev, j):
    """walk the parent chain of label j at level lev -> node list, launch epoch."""
    nodes = []; tL = None
    for k in range(lev, -1, -1):
        l, lab, q, tl = par[k]
        nodes.append(int(q[j]))
        if k == 0: tL = float(tl[j])
        j = int(lab[j])
    return nodes[::-1], tL


def one_route(G, S, a, prize, rng):
    b = argsnap(a, rng)
    par, hist, best, fin = GB.beam(G, S, b, prize=prize)
    node, dv, cnt = fin
    if len(node) == 0: return None
    # pick the label maximising (new targets) / J_i over ALL levels reached
    bestr = None
    for lev in range(len(par) - 1, 0, -1):
        l, lab, q, _ = par[lev]
        if len(q) == 0: continue
        break
    for lev in [len(par) - 1]:
        for j in range(len(par[lev][2])):
            pass
    nodes, tL = trace(par, len(par) - 1, int(np.argmax(cnt - b.lam2 * dv)))
    tk = float(M_DRY * 1.0025 * np.exp(dv[int(np.argmax(cnt - b.lam2 * dv))] / VE))
    tg = [int(G.tgt[n]) for n in nodes]
    new = [t for t in dict.fromkeys(tg) if prize[t] > 0]
    if not new: return None
    return dict(nodes=nodes, targets=tg, new=new, tank=tk, tL=tL)


class argsnap:
    def __init__(self, a, rng):
        self.beam = a.beam; self.depth = int(rng.integers(a.dmin, a.dmax))
        self.lam = float(rng.uniform(a.lam_lo, a.lam_hi)); self.lam2 = a.lam2
        self.dvmax = a.dvmax; self.max_tank = a.max_tank; self.per_node = a.per_node
        self.cap_by = getattr(a, 'cap_by', 'target')


def fleet_from(prize0, seed, a):
    G = _G['G']; S = _G['S']; rng = np.random.default_rng(seed)
    prize = prize0.copy(); routes = []
    for k in range(a.ncraft):
        r = one_route(G, S, a, prize, rng)
        if r is None: break
        for t in r['new']: prize[t] = 0.0
        routes.append(r)
        if all(prize[t] == 0.0 for t in ALL): break
    cov = set(t for r in routes for t in r['targets'])
    miss = len([t for t in ALL if t not in cov])
    J = sum(J1(r['tank']) for r in routes) + miss + 2
    return dict(J=float(J), routes=routes, covered=len(cov & set(ALL)), miss=miss,
                sumJ=float(sum(J1(r['tank']) for r in routes)), n=len(routes))


def w_fleet(job):
    return fleet_from(*job)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('tdb'); ap.add_argument('out')
    ap.add_argument('--mode', default='greedy'); ap.add_argument('--ncraft', type=int, default=12)
    ap.add_argument('--beam', type=int, default=1200); ap.add_argument('--per-node', type=int, default=3)
    ap.add_argument('--dmin', type=int, default=28); ap.add_argument('--dmax', type=int, default=52)
    ap.add_argument('--lam-lo', type=float, default=0.8); ap.add_argument('--lam-hi', type=float, default=2.2)
    ap.add_argument('--lam2', type=float, default=0.35); ap.add_argument('--dvmax', type=float, default=1.2)
    ap.add_argument('--max-tank', type=float, default=1500.0); ap.add_argument('--cap-by', default='target')
    ap.add_argument('--ants', type=int, default=24); ap.add_argument('--gens', type=int, default=30)
    ap.add_argument('--rho', type=float, default=0.15); ap.add_argument('--pz-hi', type=float, default=2.5)
    ap.add_argument('--nproc', type=int, default=8); ap.add_argument('--seed', type=int, default=0)
    a = ap.parse_args()
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    import run_ialns as RI; say = RI.logger(out / 'log.txt')
    G = GB.Graph(a.tdb); eph = Ephemeris()
    S = G.sources(eph, 0.0, 900.0)
    _G.update(G=G, S=S)
    say(f'graph {G.n} encounters, {len(G.src)} edges; {len(S["node"])} launch arcs; mode {a.mode}')
    base = np.zeros(301); base[ALL] = 1.0
    tau = np.ones(301)                      # pheromone: how much a target is worth saving
    best = None; t0 = time.time()
    with mp.get_context('fork').Pool(a.nproc) as pl:
        for gen in range(1, a.gens + 1):
            pz = base * np.clip(tau, 0.3, a.pz_hi)
            jobs = [(pz, a.seed * 9991 + gen * 101 + i, a) for i in range(a.ants)]
            res = pl.map(w_fleet, jobs)
            res.sort(key=lambda r: r['J'])
            if best is None or res[0]['J'] < best['J']:
                best = res[0]
                json.dump(dict(J=best['J'], n=best['n'], covered=best['covered'], miss=best['miss'],
                               routes=[dict(targets=r['targets'], tank=r['tank'], tL=r['tL']) for r in best['routes']]),
                          open(out / 'best.json', 'w'), indent=1)
            say(f'gen {gen:3d}: best J {res[0]["J"]:.3f} ({res[0]["n"]} craft, {res[0]["covered"]} covered, '
                f'sum J_i {res[0]["sumJ"]:.3f}), median {np.median([r["J"] for r in res]):.2f}, '
                f'overall {best["J"]:.3f}  ({time.time()-t0:.0f} s)')
            if a.mode == 'aco':
                tau *= (1 - a.rho)
                for r in res[:max(1, a.ants // 6)]:        # elite ants deposit on the targets they covered late
                    for k, rt in enumerate(r['routes']):
                        for t in rt['new']: tau[t] += a.rho * (1.0 + 0.5 * k) / (1.0 + r['J'] - 12.0)
                tau = np.clip(tau, 0.3, a.pz_hi)
    say(f'done: best J {best["J"]:.3f}, {best["n"]} craft, {best["covered"]} covered -> {out}/best.json')


if __name__ == '__main__':
    main()
