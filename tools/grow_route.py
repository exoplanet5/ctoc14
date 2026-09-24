"""Grow one route from scratch with the impulsive whole-trajectory model (cheapest-insertion construction).

Start: launch at --tl days with v_inf --vinf (km/s, 3 comps, ecliptic frame) and no flyby; impulse nodes every 20 d to
the mission end. Loop: close approaches (< --dmax AU) of the current trajectory to every allowed target asteroid, nearest
--trials first; parallel insertion trials (aim-point homotopy + SCP); accept the cheapest if dJ < --lam; stop otherwise.
Output: <out>/route.npz (impulsive state) + grow.log.
Usage: grow_route.py outdir --tl 0 --vinf 0,0,0 [--targets all|file] [--exclude file] [--lam 0.04] [--dmax 0.05] [--trials 16]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse, multiprocessing as mp
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1])); sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import numpy as np
from run_ialns import ist, ipr, w_insert, w_cands, cost, logger, eph
from ctoc14.impulsive import ImpulsiveProblem
from ctoc14.constants import DAY, T_MISSION
from ctoc14.search import UNREACHABLE


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out')
    ap.add_argument('--tl', type=float, default=0.0); ap.add_argument('--vinf', default='0,0,0')
    ap.add_argument('--targets', default='all'); ap.add_argument('--exclude', default=None)
    ap.add_argument('--lam', type=float, default=0.04); ap.add_argument('--dmax', type=float, default=0.05)
    ap.add_argument('--trials', type=int, default=16); ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--max-flybys', type=int, default=60)
    a = ap.parse_args()
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True); say = logger(out / 'grow.log')
    allowed = set(range(1, 301)) - set(UNREACHABLE)
    if a.targets != 'all':
        allowed &= set(int(x) for x in open(a.targets).read().split())
    if a.exclude:
        allowed -= set(int(x) for x in open(a.exclude).read().split())
    E = eph(); tL = a.tl * DAY
    ts = np.arange(tL + 10 * DAY, T_MISSION - 2 * DAY, 20 * DAY)
    ip = ImpulsiveProblem(E, tL, np.array([float(x) for x in a.vinf.split(',')]), ts, np.zeros((len(ts), 3)), np.zeros(0), [])
    st = ist(ip); tank = ip.tank()
    say(f'grow: launch {a.tl} d, v_inf {a.vinf}, {len(allowed)} allowed targets, lambda {a.lam}')
    with mp.get_context('fork').Pool(a.nproc, maxtasksperchild=40) as pool:
        while len(st['asts']) < a.max_flybys:
            todo = sorted(allowed - set(int(x) for x in st['asts']))
            cands = w_cands(('r', st, todo, a.dmax))
            cands = sorted(cands, key=lambda c: c['dist'])[:a.trials]
            if not cands:
                say('no candidate within dmax'); break
            res = pool.map(w_insert, [('r', st, c['ast'], c['t']) for c in cands])
            good = [r for r in res if r['ok']]
            for r in good:
                r['dJ'] = cost(r['tank']) - cost(tank)
            if not good:
                say(f'no successful trial among {len(cands)}'); break
            best = min(good, key=lambda r: r['dJ'])
            if best['dJ'] > a.lam:
                say(f'cheapest insertion {best["ast"]} dJ {best["dJ"]:.4f} > lambda: stop'); break
            st = best['st']; tank = best['tank']
            np.savez(out / 'route.npz', **st)
            say(f'+{best["ast"]} at {best["t"] / DAY:.0f} d (dist {next(c["dist"] for c in cands if c["ast"] == best["ast"]):.3f} AU): '
                f'dJ {best["dJ"]:+.4f}; flybys {len(st["asts"])}, tank {tank:.1f} kg, J_i {cost(tank):.4f} ({len(good)}/{len(cands)} ok)')
    say(f'done: {len(st["asts"])} flybys, tank {tank:.1f} kg, J_i {cost(tank):.4f}, J_i per flyby {cost(tank) / max(1, len(st["asts"])):.4f}')


if __name__ == '__main__':
    main()
