"""Partition local search: given a fleet of N impulsive routes, repeatedly RE-PLAN one route against the targets
nobody covers.  For route i: neighbourhood = its own targets + the K uncovered targets nearest to its trajectory;
run the beam inside that neighbourhood with a prize of 1 J on every target that no OTHER route covers (a target
lost by route i and picked up by nobody is a miss, worth exactly 1 J), so the beam optimises the true objective
dJ = cost(new) - cost(old) - (targets gained - targets lost).  Accept the best candidate, rebuild its impulsive
twin (settle_tour) for the real cost and the next distance ranking, and move on.

Usage: partition_ls.py src_dir out_dir [--rounds 6] [--knn 30] [--m0 1200] [--nproc 8]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from ctoc14.search import Params, beam_search, tour_json, UNREACHABLE
from ctoc14.colgen import CostModel
from ctoc14.impulsive import settle_tour
from ctoc14.constants import DAY, AU

ALL = [t for t in range(1, 301) if t not in UNREACHABLE]
CM = CostModel(s=0.70, reserve=2.0)


class DeepCollect(list):
    def __init__(self, minn): super().__init__(); self.minn = minn
    def extend(self, states):
        for s in states:
            if len(s.seq) >= self.minn: list.append(self, s)


class BP(Params):
    def __getstate__(self):
        d = dict(self.__dict__); d.pop('collect', None); return d


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('src'); ap.add_argument('dst')
    ap.add_argument('--rounds', type=int, default=6); ap.add_argument('--knn', type=int, default=30)
    ap.add_argument('--m0', type=float, default=1600.0); ap.add_argument('--beam', type=int, default=400)
    ap.add_argument('--prizes', default='', help='JSON file with 300 per-target prizes; biases the beam toward hard uncovered targets (the accept test stays at 1 J/target)')
    ap.add_argument('--nproc', type=int, default=8); ap.add_argument('--dmax', type=float, default=1.5)
    ap.add_argument('--min-keep', type=int, default=18); ap.add_argument('--grid', default='0,800,20')
    ap.add_argument('--max-tank', type=float, default=1250.0); ap.add_argument('--iters', type=int, default=100)
    ap.add_argument('--vinf', type=float, default=4.0); ap.add_argument('--dvmax', type=float, default=1.2)
    ap.add_argument('--tofmax', type=float, default=400.0); ap.add_argument('--lintofmax', type=float, default=600.0)
    a = ap.parse_args()
    out = pathlib.Path(a.dst); out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt')
    fleet = RI.IFleet(a.src); RI.eph()
    PZ = np.ones(300)
    if a.prizes:
        PZ = np.asarray(json.load(open(a.prizes)), float)
        say(f'prizes from {a.prizes}: median {np.median(PZ):.2f} max {PZ.max():.2f}')
    say(f'start: {fleet.summary()}')
    g = [float(x) for x in a.grid.split(',')]; grid = np.arange(g[0], g[1] + 1e-9, g[2]) * DAY
    wf = CM.w_fuel(0.5 * (a.m0 - 640.0), a.m0)
    P = BP(beam=a.beam, w_fuel=wf, m_margin=40.0, dv_max=a.dvmax, lin_tofs=np.arange(20, a.lintofmax + 1e-9, 10) * DAY,
           tofs=np.arange(15, a.tofmax + 1e-9, 5) * DAY,
           lin_drmax=0.15 * AU, vinf_cap=a.vinf, tof_refine=True)
    say(f'beam w_fuel {wf:.3f} (J per full-tank fraction at m0 {a.m0:.0f})')
    with mp.get_context('fork').Pool(a.nproc) as pool:
        for rnd in range(a.rounds):
            changed = 0
            order = sorted(fleet.routes, key=lambda n: len(fleet.routes[n]['st']['asts']))
            for name in order:
                t0 = time.time()
                cov = fleet.coverage()
                own = set(int(x) for x in fleet.routes[name]['st']['asts'])
                other = set(t for t, v in cov.items() if v - {name})
                U = [t for t in ALL if t not in cov]
                if not U:
                    continue
                cand = pool.map(RI.w_cands, [(name, fleet.routes[name]['st'], U, a.dmax)])[0]
                d = {}
                for c in cand: d[c['ast']] = min(d.get(c['ast'], 9e9), c['dist'])
                near = [t for _, t in sorted((v, k) for k, v in d.items())][:a.knn]
                sub = sorted(own | set(near))
                prize = np.zeros(300)
                for t in sub:
                    if t not in other: prize[t - 1] = PZ[t - 1]
                P.prize = prize; P.collect = DeepCollect(a.min_keep)
                best, beam = beam_search(RI.eph(), excluded=sorted(set(ALL) - set(sub)), m0=a.m0,
                                         t_launch_grid=grid, P=P, n_proc=a.nproc, verbose=False)
                states = list(P.collect) + list(beam)
                c_old = RI.cost(fleet.routes[name]['tank']); own_val = len(own - other)
                scored = []
                for s in states:
                    tg = set(x[0] for x in s.seq)
                    c_new = float(CM.cost(s.fuel, s.m0))
                    dJ = c_new - c_old - (len(tg - other) - own_val)
                    scored.append((dJ, c_new, s, tg))
                scored.sort(key=lambda x: x[0])
                deep = max(scored, key=lambda x: len(x[3])) if scored else None
                say(f'   {name}: subset {len(sub)} ({len(own)} own + {len(near)} uncovered, '
                    f'nearest {min(d.values()) if d else -1:.2f}-{(sorted(d.values())[len(near)-1] if near else -1):.2f} AU), '
                    f'{len(states)} states, deepest {len(deep[3]) if deep else 0} targets '
                    f'(cost {deep[1] if deep else 0:.3f}, dJ {deep[0] if deep else 0:+.3f}), best dJ {scored[0][0]:+.3f}')
                took = None
                for dJ, c_new, s, tg in scored[:6]:
                    if dJ > -0.02: break
                    tour = tour_json(s)
                    try:
                        ip, miss, lag = settle_tour(RI.eph(), tour, lambda ip: RI.settle(ip, a.iters))
                    except Exception as e:
                        say(f'   {name}: candidate {len(tg)} targets raised {type(e).__name__}'); continue
                    if ip is None or miss > 150.0: 
                        say(f'   {name}: candidate {len(tg)} targets failed to settle (miss {miss:.0f})'); continue
                    tank = float(ip.tank())
                    if tank > a.max_tank:
                        say(f'   {name}: candidate {len(tg)} targets tank {tank:.0f} > cap'); continue
                    dJ_true = RI.cost(tank) - c_old - (len(tg - other) - own_val)
                    if dJ_true < -0.02:
                        took = (ip, tank, tg, dJ_true, dJ, c_new); break
                    say(f'   {name}: candidate {len(tg)} targets planner dJ {dJ:+.3f} -> twin {dJ_true:+.3f}, rejected')
                if took is None:
                    say(f'[r{rnd}] {name}: no improvement ({len(own)} targets, {len(U)} uncovered, '
                        f'best planner dJ {scored[0][0] if scored else 0:+.3f})  ({time.time()-t0:.0f} s)')
                    continue
                ip, tank, tg, dJ_true, dJ, c_new = took
                fleet.routes[name] = dict(st=RI.ist(ip), tank=tank)
                fleet.save(out, note=f'round {rnd} {name}')
                changed += 1
                say(f'[r{rnd}] {name}: {len(own)} -> {len(tg)} targets ({len(tg-other)-own_val:+d} net), tank '
                    f'{fleet.routes[name]["tank"]:.0f}, dJ {dJ_true:+.3f}; fleet J {fleet.J():.4f}, '
                    f'covered {len(fleet.coverage())}  ({time.time()-t0:.0f} s)')
            say(f'--- round {rnd}: {changed} routes replaced; {fleet.summary()}')
            if not changed: break
    fleet.save(out, note='final'); say(f'final: {fleet.summary()}')


if __name__ == '__main__':
    main()
