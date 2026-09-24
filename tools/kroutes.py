"""k-routes: Lloyd's algorithm where the centroids are beam-planned trajectories.

J_i is convex in the route length, so the best fleet of N craft is N BALANCED routes of ~298/N targets each --
a greedy cover (46, 43, 38, 33, ... then 67 unchainable orphans) is exactly the wrong shape, and merging the routes
of an existing fleet cannot balance them either. So iterate:
  assign  : every target goes to the route whose trajectory it passes closest to, under a capacity cap (~298/N),
            so the groups stay balanced; unreachable ones go to the emptiest group;
  re-plan : beam inside each group (plus a small pad of the nearest outsiders), prize 1 J per target, and rebuild
            the impulsive twin of the best candidate;
  repeat  : the new routes are the new centroids.
Every round is scored with the true objective (sum cost + misses) and the best fleet is saved.
Usage: kroutes.py seed_fleet out_dir [--n 7] [--rounds 4] [--m0 1600] [--vinf 4.0] [--nproc 8]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, signal, shutil, pathlib, argparse, multiprocessing as mp
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
FAR = 3.0


class DeepCollect(list):
    def __init__(self, minn): super().__init__(); self.minn = minn
    def extend(self, states):
        for s in states:
            if len(s.seq) >= self.minn: list.append(self, s)


class BP(Params):
    def __getstate__(self):
        d = dict(self.__dict__); d.pop('collect', None); return d


def _timeout(sig, frm):
    raise TimeoutError('settle timeout')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('src'); ap.add_argument('out')
    ap.add_argument('--n', type=int, default=7); ap.add_argument('--rounds', type=int, default=4)
    ap.add_argument('--m0', type=float, default=1600.0); ap.add_argument('--vinf', type=float, default=4.0)
    ap.add_argument('--beam', type=int, default=300); ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--cap', type=float, default=1.12, help='capacity cap as a multiple of 298/n')
    ap.add_argument('--pad', type=int, default=14, help='nearest outsiders added to each group before the beam')
    ap.add_argument('--min-keep', type=int, default=20); ap.add_argument('--dmax', type=float, default=2.0)
    ap.add_argument('--try-best', type=int, default=3); ap.add_argument('--scan', type=int, default=20)
    ap.add_argument('--settle-timeout', type=float, default=150.0); ap.add_argument('--iters', type=int, default=100)
    ap.add_argument('--columns', default=''); ap.add_argument('--max-tank', type=float, default=1300.0)
    ap.add_argument('--prizes', default='', help='JSON file with 300 per-target prizes; own-group targets use these')
    ap.add_argument('--pad-prize', type=float, default=0.3, help='prize on padded outsiders (another route owns them)')
    ap.add_argument('--dual-up', type=float, default=1.6, help='multiply an uncovered target prize by this each round')
    ap.add_argument('--dual-cap', type=float, default=4.0, help='cap on the dual-scaled prize')
    a = ap.parse_args()
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt')
    seed = RI.IFleet(a.src); RI.eph()
    names = sorted(seed.routes, key=lambda n: -len(seed.routes[n]['st']['asts']))[:a.n]
    fleet = RI.IFleet(); fleet.routes = {f'{i:02d}': seed.routes[n] for i, n in enumerate(names, 1)}
    say(f'seed ({a.src}): {fleet.summary()}')
    wf = CM.w_fuel(0.5 * (a.m0 - 640.0), a.m0)
    P = BP(beam=a.beam, w_fuel=wf, m_margin=40.0, dv_max=1.2, lin_tofs=np.arange(20, 600 + 1e-9, 10) * DAY,
           lin_drmax=0.15 * AU, vinf_cap=a.vinf, tof_refine=True)
    cap = int(np.ceil(a.cap * len(ALL) / a.n))
    say(f'n {a.n}, m0 {a.m0:.0f}, vinf {a.vinf}, capacity {cap} targets per group, pad {a.pad}')
    cf = open(a.columns, 'a') if a.columns else None
    PZ = np.ones(300)
    if a.prizes:
        PZ = np.asarray(json.load(open(a.prizes)), float)
        say(f'prizes from {a.prizes}: median {np.median(PZ):.2f} max {PZ.max():.2f}; '
            f'pad-prize {a.pad_prize}, dual x{a.dual_up} cap {a.dual_cap}')
    pz = PZ.copy()                                   # dual-scaled working prizes
    unc = set(ALL) - set(fleet.coverage())           # targets the seed fleet misses (assigned first each round)
    best = (1e9, None)
    with mp.get_context('fork').Pool(a.nproc) as pool:
        for rnd in range(a.rounds):
            tic = time.time()
            keys = sorted(fleet.routes)
            res = pool.map(RI.w_cands, [(n, fleet.routes[n]['st'], ALL, a.dmax) for n in keys])
            D = {n: {} for n in keys}
            for n, lst in zip(keys, res):
                for c in lst: D[n][c['ast']] = min(D[n].get(c['ast'], 9e9), c['dist'])
            for n in keys:                       # a route's own targets are at distance 0
                for t in fleet.routes[n]['st']['asts']: D[n][int(t)] = 0.0
            # uncovered targets first: each to its nearest route with room (they are the ones the fleet keeps missing)
            grp = {n: set() for n in keys}; done = set()
            for t in sorted(unc, key=lambda t: min(D[n].get(t, FAR) for n in keys)):
                n = min(keys, key=lambda n: (D[n].get(t, FAR), len(grp[n])))
                if len(grp[n]) < cap: grp[n].add(t); done.add(t)
            # then balanced fill: cheapest (target, group) first, respecting the capacity
            opts = sorted((D[n].get(t, FAR), t, n) for t in ALL for n in keys if t not in done)
            for d, t, n in opts:
                if t in done or len(grp[n]) >= cap: continue
                grp[n].add(t); done.add(t)
            for t in ALL:                        # leftovers to the emptiest group
                if t not in done:
                    n = min(keys, key=lambda n: len(grp[n])); grp[n].add(t); done.add(t)
            say(f'[r{rnd}] groups: ' + ' '.join(f'{n}:{len(grp[n])}' for n in keys))
            newf = RI.IFleet()
            for n in keys:
                sub = set(grp[n])
                for t in sorted(D[n], key=lambda t: D[n][t]):
                    if len(sub) >= len(grp[n]) + a.pad: break
                    sub.add(t)
                prize = np.zeros(300)
                for t in grp[n]: prize[t - 1] = pz[t - 1]           # own group at the (dual-scaled) prize
                for t in sub - grp[n]: prize[t - 1] = a.pad_prize    # padded outsiders are cheap: another route owns them
                P.prize = prize; P.collect = DeepCollect(a.min_keep); t0 = time.time()
                bs, beam = beam_search(RI.eph(), excluded=sorted(set(ALL) - sub), m0=a.m0,
                                       t_launch_grid=np.arange(0, 800 + 1e-9, 20) * DAY, P=P,
                                       n_proc=a.nproc, verbose=False)
                uniq = {}
                for s in list(P.collect) + list(beam):
                    k = frozenset(x[0] for x in s.seq)
                    if k not in uniq or s.fuel < uniq[k].fuel: uniq[k] = s
                cand = sorted(uniq.values(), key=lambda s: float(CM.cost(s.fuel, s.m0)) - len(s.seq))
                if cf:
                    for s in cand[:250]:
                        r = tour_json(s); r['tag'] = f'kr{a.n}r{rnd}g{n}'; r['targets'] = sorted(x[0] for x in s.seq)
                        r['cost_planner'] = float(CM.cost(s.fuel, s.m0)); cf.write(json.dumps(r) + '\n')
                    cf.flush()
                trial = []; seen = []
                for s in cand:
                    if len(trial) >= a.scan: break
                    ts = set(x[0] for x in s.seq)
                    if any(len(ts & o) > 0.92 * max(len(ts), len(o)) for o in seen): continue
                    seen.append(ts); trial.append(s)
                took = None; bv = 0.0; nok = 0
                for s in trial:
                    if nok >= a.try_best: break
                    tg = sorted(x[0] for x in s.seq)
                    try:
                        signal.signal(signal.SIGALRM, _timeout); signal.setitimer(signal.ITIMER_REAL, a.settle_timeout)
                        ip, miss, lag = settle_tour(RI.eph(), tour_json(s), lambda ip: RI.settle(ip, a.iters))
                    except Exception as e:
                        say(f'    {n}: {len(tg)} flybys settle raised {type(e).__name__}'); continue
                    finally:
                        signal.setitimer(signal.ITIMER_REAL, 0)
                    if ip is None or miss > 150.0:
                        say(f'    {n}: {len(tg)} flybys settle failed (miss {miss:.0f} km)'); continue
                    nok += 1; tank = float(ip.tank()); v = RI.cost(tank) - len(tg)
                    if tank <= a.max_tank and (took is None or v < bv): took = (ip, tank, tg); bv = v
                if took is None:
                    say(f'    {n}: no settle, keeping the old route'); newf.routes[n] = fleet.routes[n]; continue
                ip, tank, tg = took
                newf.routes[n] = dict(st=RI.ist(ip), tank=tank)
                say(f'    {n}: group {len(grp[n])} (+{len(sub)-len(grp[n])} pad) -> {len(tg)} flybys, tank {tank:.0f}, '
                    f'J_i {RI.cost(tank):.3f} ({ip.dv()/len(tg):.3f} km/s per flyby)  ({time.time()-t0:.0f} s)')
            fleet = newf
            cov = fleet.coverage(); J = fleet.J()
            unc = set(ALL) - set(cov)                # next round assigns these first, and raises their prize
            for t in unc: pz[t - 1] = min(a.dual_cap, pz[t - 1] * a.dual_up)
            fleet.save(out / f'r{rnd}', note=f'round {rnd}')
            say(f'[r{rnd}] {fleet.summary()}  covered {len(cov)}, uncovered {len(ALL)-len(cov)}  ({time.time()-tic:.0f} s)')
            if J < best[0]:
                best = (J, rnd)
                for f in (out / f'r{rnd}').glob('route_*.npz'): shutil.copy(f, out / f.name)
                json.dump(dict(round=rnd, J=J, covered=sorted(cov),
                               uncovered=sorted(set(ALL) - set(cov))), open(out / 'result.json', 'w'))
                say(f'[r{rnd}] new best J {J:.4f}')
    say(f'best: round {best[1]} J {best[0]:.4f}')


if __name__ == '__main__':
    main()
