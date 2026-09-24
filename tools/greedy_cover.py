"""Greedy cover: build a fleet one route at a time, each beam run over the targets nobody has visited yet, with a
prize of 1 J on every remaining target (a target nobody flies costs exactly 1 J) and the beam fuel weight taken from
CostModel, so the beam optimises cost(route) - targets(route) directly.  The planner mass --m0 sets the reachable
tank ceiling: CostModel(0.70, 2).tank(m0 - 640, m0), e.g. 805 kg at m0 1000, 971 kg at 1400, 1149 kg at 2000 --
routes deeper than ~36 flybys need a mass ceiling above 1000.
Each accepted route is rebuilt as an impulsive twin (settle_tour) for its true tank, and saved as route_NN.npz.
Usage: greedy_cover.py out_dir [--n 8] [--m0 1400] [--beam 300] [--nproc 8]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, signal, pathlib, argparse, multiprocessing as mp
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


def _timeout(sig, frm):
    raise TimeoutError('settle timeout')


def plan_route(avail, PZ, a, P, grid, nproc, say, tag='', cf=None):
    """Beam over the targets in `avail` (prizes PZ, a len-300 array; only avail get a non-zero prize), then settle the
    diverse top candidates until `a.try_best` of them SUCCEED (a degenerate Lambert chain must not exhaust the budget).
    Returns [(ip, tank, targets, value), ...] sorted by value ascending (value = twin J_i - sum PZ over targets), or []
    if nothing settles.  Shared by greedy_cover.main and tools/fleet_tree.py."""
    avail = sorted(set(avail))
    prize = np.zeros(300)
    for t in avail:
        prize[t - 1] = PZ[t - 1]
    P.prize = prize; P.collect = DeepCollect(a.min_keep)
    val = lambda s: float(CM.cost(s.fuel, s.m0)) - float(sum(PZ[x[0] - 1] for x in s.seq))
    t0 = time.time()
    best, beam = beam_search(RI.eph(), excluded=sorted(set(ALL) - set(avail)), m0=a.m0,
                             t_launch_grid=grid, P=P, n_proc=nproc, verbose=False)
    uniq = {}
    for s in list(P.collect) + list(beam):
        key = frozenset(x[0] for x in s.seq)
        if key not in uniq or s.fuel < uniq[key].fuel:
            uniq[key] = s
    cand = sorted(uniq.values(), key=val)
    if not cand:
        say(f'{tag}: {len(avail)} avail, beam {time.time()-t0:.0f} s, no states'); return []
    if cf:
        for s in cand[:300]:
            r = tour_json(s); r['tag'] = tag; r['targets'] = sorted(x[0] for x in s.seq)
            r['cost_planner'] = float(CM.cost(s.fuel, s.m0)); cf.write(json.dumps(r) + '\n')
        cf.flush()
    say(f'{tag}: {len(avail)} avail, beam {time.time()-t0:.0f} s, best planner {len(cand[0].seq)} flybys '
        f'(fuel {cand[0].fuel:.0f}, tank {CM.tank(cand[0].fuel, cand[0].m0):.0f})')
    # keep a diverse shortlist (<=92% target overlap), then settle until try_best succeed or the budget runs out
    trial = []; seen_sets = []
    for s in cand:
        if len(trial) >= a.scan: break
        ts = set(x[0] for x in s.seq)
        if any(len(ts & o) > 0.92 * max(len(ts), len(o)) for o in seen_sets): continue
        seen_sets.append(ts); trial.append(s)
    out = []; nok = 0; t_settle = time.time()
    for s in trial:
        if nok >= a.try_best or time.time() - t_settle > a.settle_budget: break
        tour = tour_json(s); tg = sorted(x[0] for x in s.seq); tic = time.time()
        try:
            signal.signal(signal.SIGALRM, _timeout); signal.setitimer(signal.ITIMER_REAL, a.settle_timeout)
            ip, miss, lag = settle_tour(RI.eph(), tour, lambda ip: RI.settle(ip, a.iters))
        except Exception as e:                       # singular linearisation / timeout: drop this candidate
            say(f'   {tag} {len(tg)} flybys: settle raised {type(e).__name__}: {repr(e)[:50]} ({time.time()-tic:.0f} s)')
            continue
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
        if ip is None or miss > 150.0:
            say(f'   {tag} {len(tg)} flybys: settle failed (miss {miss:.0f} km, {time.time()-tic:.0f} s)'); continue
        nok += 1
        tank = float(ip.tank()); value = RI.cost(tank) - float(sum(PZ[t - 1] for t in tg))
        say(f'   {tag} {len(tg)} flybys: planner tank {CM.tank(s.fuel, s.m0):.0f} -> twin tank {tank:.0f} '
            f'(J_i {RI.cost(tank):.3f}, dv {ip.dv():.1f} = {ip.dv()/len(tg):.3f}/flyby, value {value:+.3f}, {time.time()-tic:.0f} s)')
        out.append((ip, tank, tg, value))
    out.sort(key=lambda x: x[3])
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--n', type=int, default=8)
    ap.add_argument('--m0', type=float, default=1400.0); ap.add_argument('--beam', type=int, default=300)
    ap.add_argument('--nproc', type=int, default=8); ap.add_argument('--grid', default='0,800,20')
    ap.add_argument('--min-keep', type=int, default=20); ap.add_argument('--max-tank', type=float, default=1300.0)
    ap.add_argument('--iters', type=int, default=100); ap.add_argument('--try-best', type=int, default=5)
    ap.add_argument('--scan', type=int, default=40, help='diverse candidates to consider until try-best settle')
    ap.add_argument('--settle-timeout', type=float, default=180.0)
    ap.add_argument('--settle-budget', type=float, default=600.0, help='wall-clock cap on the settle phase per route')
    ap.add_argument('--wt', type=float, default=1.0); ap.add_argument('--start', default='')
    ap.add_argument('--max-depth', type=int, default=120, help='cap the route length: forces a balanced fleet')
    ap.add_argument('--prizes', default='', help='JSON file with 300 per-target prizes (J units); default all 1')
    ap.add_argument('--vinf', type=float, default=2.0); ap.add_argument('--dvmax', type=float, default=1.2)
    ap.add_argument('--tofmax', type=float, default=400.0); ap.add_argument('--lintofmax', type=float, default=600.0)
    ap.add_argument('--columns', default='', help='also append every kept beam state to this jsonl pool')
    a = ap.parse_args()
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt')
    fleet = RI.IFleet(a.start) if a.start else RI.IFleet()
    RI.eph()
    wf = CM.w_fuel(0.5 * (a.m0 - 640.0), a.m0)
    say(f'greedy cover: n {a.n}, m0 {a.m0:.0f}, w_fuel {wf:.3f}, tank ceiling {CM.tank(a.m0-640.0, a.m0):.0f} kg'
        + (f', start {fleet.summary()}' if a.start else ''))
    g = [float(x) for x in a.grid.split(',')]; grid = np.arange(g[0], g[1] + 1e-9, g[2]) * DAY
    P = BP(beam=a.beam, w_fuel=wf, m_margin=40.0, dv_max=a.dvmax, lin_tofs=np.arange(20, a.lintofmax + 1e-9, 10) * DAY,
           tofs=np.arange(15, a.tofmax + 1e-9, 5) * DAY,
           lin_drmax=0.15 * AU, vinf_cap=a.vinf, tof_refine=True, w_t=a.wt, max_depth=a.max_depth)
    cf = open(a.columns, 'a') if a.columns else None
    PZ = np.ones(300)
    if a.prizes:
        PZ = np.asarray(json.load(open(a.prizes)), float)
        say(f'prizes from {a.prizes}: min {PZ.min():.2f} median {np.median(PZ):.2f} max {PZ.max():.2f}')
    k0 = len(fleet.routes)
    for k in range(k0 + 1, a.n + 1):
        cov = set(fleet.coverage()); R = [t for t in ALL if t not in cov]
        if not R: break
        t0 = time.time()
        got = plan_route(R, PZ, a, P, grid, a.nproc, say, tag=f'route {k}', cf=cf)
        took = next(((ip, tank, tg) for ip, tank, tg, value in got if tank <= a.max_tank), None)
        if took is None:
            say(f'route {k}: no acceptable candidate'); break
        ip, tank, tg = took
        fleet.routes[f'{k:02d}'] = dict(st=RI.ist(ip), tank=tank)
        fleet.save(out, note=f'route {k}')
        say(f'route {k}: ACCEPTED {len(tg)} targets, tank {tank:.0f}; fleet {fleet.summary()}  ({time.time()-t0:.0f} s)')
    say(f'final: {fleet.summary()}')
    left = sorted(set(ALL) - set(fleet.coverage()))
    say(f'uncovered ({len(left)}): {left}')
    json.dump(dict(uncovered=left, covered=sorted(fleet.coverage()), J=fleet.J(),
                   routes={n: sorted(int(x) for x in r['st']['asts']) for n, r in fleet.routes.items()},
                   tanks={n: r['tank'] for n, r in fleet.routes.items()}), open(out / 'result.json', 'w'))


if __name__ == '__main__':
    main()
