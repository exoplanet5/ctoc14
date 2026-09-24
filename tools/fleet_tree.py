"""Branch A of docs/stage4_n8_search_tree.md: a beam search whose STATES are partial fleets.

A node at depth d = d disjoint, settled routes + the set R of still-uncovered targets.  Expanding a node runs one
prize beam over R (hardness prizes H, route length capped at D_d so easy connectors are saved for later routes),
settles the diverse top candidates (greedy_cover.plan_route), and turns up to `--branch` genuinely different ones
into children.  We keep the `--width` best partial fleets per depth by an optimistic node value and prune children
that can no longer reach the coverage floor or bring the hard leftovers to <=4, or whose new route is too heavy.

Prizes make early routes shallow-but-hard-dense and later routes deep-and-easy (the hard targets are absorbed first,
then only easy ones remain), so there is NO per-route depth floor -- a 26-flyby hard route is fine if the fleet can
still reach `--cover-goal`.  Node value (J units, lower=better): sum J_i(done) + n_left*[1 + c(nbar)] + forced_misses
+ 0.35*sum_{t in R}(H_t - 1), with nbar = min(|R|, 44*n_left)/n_left and forced_misses = max(0, |R| - 44*n_left);
c(n) = x + x^2 at tank 601.5*exp(0.45 n / VE).  A TERMINAL node (d routes = N, or R empty) is scored on the truth,
sum J_i(done) + |R| (each leftover is a 1 J miss).

Usage: fleet_tree.py out_dir [--n 8] [--width 3] [--branch 3] [--prizes f.json] [--m0 1600] [--vinf 4.0] [--nproc 8]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, math, pathlib, argparse
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from greedy_cover import plan_route, BP, ALL, CM
from ctoc14.constants import DAY, AU

VE = 39.2266


def est_Ji(n):
    """1 + c(n): the estimated J_i of a route of n flybys at 0.45 km/s each."""
    tank = 601.5 * math.exp(0.45 * n / VE); x = (tank - 600.0) / 1400.0
    return 1.0 + x + x * x


class Node:
    __slots__ = ('routes', 'R', 'term', 'value')

    def __init__(self, routes, R, H, N, term=None):
        self.routes = routes                 # list of (st, tank, targets)
        self.R = frozenset(R)
        self.term = (len(routes) >= N or not self.R) if term is None else term
        d = len(routes); rem = len(self.R)
        doneJ = sum(RI.cost(t) for _, t, _ in routes)
        if self.term:
            self.value = doneJ + rem         # each leftover is a full miss
        else:
            nleft = N - d
            coverable = min(rem, 44 * nleft)                       # remaining routes visit <=44 targets each
            nbar = coverable / nleft
            forced_miss = rem - coverable                          # targets no remaining route can reach -> misses
            self.value = (doneJ + nleft * est_Ji(nbar) + forced_miss
                          + 0.35 * sum(H[t - 1] - 1.0 for t in self.R))


def diverse(got, thr, b):
    picks = []; sets = []
    for ip, tank, tg, value in got:
        ts = set(tg)
        if any(len(ts & o) > thr * max(len(ts), len(o)) for o in sets): continue
        sets.append(ts); picks.append((ip, tank, tg, value))
        if len(picks) >= b: break
    return picks


def to_fleet(node):
    fl = RI.IFleet()
    for i, (st, tank, tg) in enumerate(node.routes, 1):
        fl.routes[f'{i:02d}'] = dict(st=st, tank=tank)
    return fl


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out')
    ap.add_argument('--n', type=int, default=8); ap.add_argument('--width', type=int, default=3)
    ap.add_argument('--branch', type=int, default=3); ap.add_argument('--prizes', default='')
    ap.add_argument('--m0', type=float, default=1600.0); ap.add_argument('--vinf', type=float, default=4.0)
    ap.add_argument('--beam', type=int, default=300); ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--dvmax', type=float, default=1.2); ap.add_argument('--tofmax', type=float, default=400.0)
    ap.add_argument('--lintofmax', type=float, default=600.0); ap.add_argument('--grid', default='0,800,20')
    ap.add_argument('--min-keep', type=int, default=20); ap.add_argument('--scan', type=int, default=30)
    ap.add_argument('--try-best', type=int, default=3); ap.add_argument('--iters', type=int, default=100)
    ap.add_argument('--settle-timeout', type=float, default=150.0); ap.add_argument('--settle-budget', type=float, default=700.0)
    ap.add_argument('--max-tank', type=float, default=1150.0); ap.add_argument('--overlap', type=float, default=0.80)
    ap.add_argument('--dvpf-max', type=float, default=0.65, help='reject a route above this km/s per flyby (prune iv); '
                    'hard-dense routes measure 0.59-0.60, so 0.60 is too tight -- missing a hard target costs a full 1 J')
    ap.add_argument('--cover-goal', type=int, default=282, help='coverage floor (gate G1); prune states that cannot reach it')
    ap.add_argument('--prize-gamma', type=float, default=1.0,
                    help='soften the prizes to H_soft = 1 + gamma*(H-1); <1 keeps routes deeper (prizes cost depth)')
    ap.add_argument('--columns', default='')
    a = ap.parse_args()
    N = a.n
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt')
    RI.eph()
    H = np.ones(300)
    if a.prizes:
        H = np.asarray(json.load(open(a.prizes)), float)
    Hfull = H.copy()                                  # the un-softened prizes define the "hard" set (prune/gate)
    if a.prize_gamma != 1.0:                           # soften the beam's prize gradient without changing what is "hard"
        H = np.where(Hfull > 0, 1.0 + a.prize_gamma * (Hfull - 1.0), 0.0)
    hard = set(t for t in ALL if Hfull[t - 1] >= 1.5)
    g = [float(x) for x in a.grid.split(',')]; grid = np.arange(g[0], g[1] + 1e-9, g[2]) * DAY
    wf = CM.w_fuel(0.5 * (a.m0 - 640.0), a.m0)
    P = BP(beam=a.beam, w_fuel=wf, m_margin=40.0, dv_max=a.dvmax,
           lin_tofs=np.arange(20, a.lintofmax + 1e-9, 10) * DAY, tofs=np.arange(15, a.tofmax + 1e-9, 5) * DAY,
           lin_drmax=0.15 * AU, vinf_cap=a.vinf, tof_refine=True)
    cf = open(a.columns, 'a') if a.columns else None
    say(f'fleet tree: N {N}, width {a.width}, branch {a.branch}, m0 {a.m0:.0f}, vinf {a.vinf}, '
        f'{len(hard)} hard targets (prize>=1.5), w_fuel {wf:.3f}')

    def prune(routes, R):
        # Prizes make early routes shallow-but-hard-dense and late routes deep-easy (the hard targets get absorbed
        # first), so a fixed per-route depth floor is wrong: what matters is that the fleet can still reach the
        # coverage floor and bring the hard leftovers to <=4.  Both are OPTIMISTIC feasibility checks (44/route max).
        d = len(routes); rem = len(R); nleft = N - d
        if nleft <= 0:
            return None
        covered = len(ALL) - rem
        if covered + 44 * nleft < a.cover_goal:                       # cannot reach the coverage floor
            return f'cov {covered}+44*{nleft}<{a.cover_goal}'
        if sum(1 for t in R if t in hard) > 18 * nleft + 4:           # cannot bring hard-left to <=4
            return 'hard'
        return None

    root = Node([], set(ALL), Hfull, N)
    beam = [root]; cache = {}; best_term = None; tic = time.time()
    for it in range(N):
        nxt = []; expanded = 0
        for node in beam:
            if node.term:
                nxt.append(node); continue
            expanded += 1
            R = node.R; d = len(node.routes); nleft = N - d
            Dd = min(44, math.ceil(len(R) / nleft) + 4)
            P.max_depth = Dd
            key = (R, Dd)
            got = cache.get(key)
            if got is None:
                got = plan_route(sorted(R), H, a, P, grid, a.nproc, say, tag=f'd{d} |R|{len(R)} D{Dd}', cf=cf)
                cache[key] = got
            made = 0
            for ip, tank, tg, value in diverse(got, a.overlap, a.branch):
                dvpf = ip.dv() / max(len(tg), 1)
                if tank > a.max_tank or dvpf > a.dvpf_max:
                    continue                                             # prune (iv)
                newroutes = node.routes + [(RI.ist(ip), tank, sorted(tg))]
                newR = R - set(tg)
                if prune(newroutes, newR):
                    continue
                nxt.append(Node(newroutes, newR, Hfull, N)); made += 1
            if made == 0:                                                # cannot extend: bank it as a terminal
                nxt.append(Node(node.routes, node.R, Hfull, N, term=True))
        if not expanded:
            break
        for nd in nxt:                                                   # track the best banked / complete fleet
            if nd.term and (best_term is None or nd.value < best_term.value):
                best_term = nd
        beam = sorted(nxt, key=lambda n: n.value)[:a.width]
        for k, nd in enumerate(beam, 1):
            to_fleet(nd).save(out / f'd{it+1}_{k}', note=f'depth {len(nd.routes)} value {nd.value:.3f}')
        top = beam[0]
        cov = len(ALL) - len(top.R)
        say(f'[depth {it+1}] beam {len(beam)}: best value {top.value:.3f} '
            f'({len(top.routes)} routes, covered {cov}, left {len(top.R)}, hard-left {sum(1 for t in top.R if t in hard)})'
            f'{" TERMINAL" if top.term else ""}  ({time.time()-tic:.0f} s)')
        if not top.term and len(top.routes) >= 5:
            nl = N - len(top.routes); cv = len(ALL) - len(top.R)
            if cv + 44 * nl < a.cover_goal + 6:
                say(f'DEAD-END: best node covers {cv}, {nl} routes left, cannot reach {a.cover_goal} with margin; seed branch B'); break
        if all(n.term for n in beam):
            say('all beam nodes terminal'); break

    best = best_term or beam[0]
    fl = to_fleet(best); fl.save(out)
    cov = fl.coverage(); left = sorted(set(ALL) - set(cov)); hleft = [t for t in left if t in hard]
    sumJi = sum(RI.cost(r['tank']) for r in fl.routes.values())
    json.dump(dict(n=len(fl.routes), covered=len(cov), covered_list=sorted(cov), leftovers=left,
                   hard_leftovers=hleft, sum_Ji=sumJi, J_impulsive=sumJi + 2 + len(left),
                   tanks={n: r['tank'] for n, r in fl.routes.items()},
                   routes={n: sorted(int(x) for x in r['st']['asts']) for n, r in fl.routes.items()}),
              open(out / 'result.json', 'w'))
    say(f'BEST: {fl.summary()}')
    say(f'covered {len(cov)}, leftovers {len(left)} ({len(hleft)} hard): {left}')
    say(f'sum J_i {sumJi:.4f}, impulsive J ~ {sumJi + 2 + len(left):.4f}  (gate G1: >=282 covered, <=16 left, <=4 hard)')
    if cf: cf.close()


if __name__ == '__main__':
    main()
