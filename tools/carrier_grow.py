"""Deepen a sparse carrier route by greedy insertion (stage 7 step 4).

Carrier routes (tools/carrier_dp.py) fly ~20 targets for ~5 km/s and are SPARSE in time, the opposite of the time-full
routes on which insertion cost 1.6-9 km/s.  Each round: closest approaches of every unflown target to the current
trajectory (run_ialns.w_cands), insertion trials of the --ntrial nearest in parallel (insert_homotopy + settle), then
accept the cheapest ones that are --sep days apart from each other, one joint re-settle, repeat.

Usage: carrier_grow.py fleet_dir route out_dir [--target 37] [--max-tank 1000]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, time, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from ctoc14.search import UNREACHABLE
from ctoc14.constants import DAY, VE
from ctoc14.globalopt import insert_block


def w_multi(job):
    st, items = job
    ip = RI.ipr(st)
    try:
        d = insert_block(ip, items, log=RI.QUIET)
        if d.max() > 1e4: return None
        miss = RI.settle(ip, 100)
    except Exception:
        return None
    return (RI.ist(ip), float(ip.tank())) if miss <= 150 else None


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('src'); ap.add_argument('route'); ap.add_argument('out')
    ap.add_argument('--target', type=int, default=37); ap.add_argument('--max-tank', type=float, default=1000.0)
    ap.add_argument('--dmax', type=float, default=0.12); ap.add_argument('--ntrial', type=int, default=24)
    ap.add_argument('--naccept', type=int, default=3); ap.add_argument('--sep', type=float, default=250.0)
    ap.add_argument('--gap', type=float, default=12.0, help='min days from an existing flyby')
    ap.add_argument('--max-step', type=float, default=45.0, help='kg: stop when a single insertion costs more (a flyby is worth ~40 kg at N=8)')
    ap.add_argument('--nproc', type=int, default=8); ap.add_argument('--exclude', default='')
    a = ap.parse_args()
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt')
    fl = RI.IFleet(a.src); st = fl.routes[a.route]['st']; tank = fl.routes[a.route]['tank']
    excl = set(UNREACHABLE) | set(int(x) for x in a.exclude.split(',') if x)
    say(f'start {a.route}: {len(st["asts"])} flybys, tank {tank:.0f} kg')
    pool = mp.get_context('fork').Pool(a.nproc); t0 = time.time(); rnd = 0
    while len(st['asts']) < a.target:
        rnd += 1
        have = set(int(x) for x in st['asts']); targets = [t for t in range(1, 301) if t not in have and t not in excl]
        chunks = [targets[i::a.nproc] for i in range(a.nproc)]
        cands = [c for lst in pool.map(RI.w_cands, [('h', st, ch, a.dmax) for ch in chunks]) for c in lst]
        cands = [c for c in cands if np.abs(np.asarray(st['tf']) - c['t']).min() > a.gap * DAY]
        best = {}
        for c in sorted(cands, key=lambda c: c['dist']): best.setdefault(c['ast'], c)
        trial = sorted(best.values(), key=lambda c: c['dist'])[:a.ntrial]
        if not trial: say('no candidates left'); break
        res = [r for r in pool.map(RI.w_insert, [('h', st, c['ast'], c['t']) for c in trial]) if r.get('ok')]
        res.sort(key=lambda r: r['tank'])
        pick = []
        for r in res:
            if r['tank'] > a.max_tank or r['tank'] - tank > a.max_step: break
            if all(abs(r['t'] - q['t']) > a.sep * DAY for q in pick): pick.append(r)
            if len(pick) >= a.naccept: break
        if not pick: say(f'round {rnd}: {len(res)} of {len(trial)} trials ok, none under {a.max_tank:.0f} kg / +{a.max_step:.0f} kg'); break
        new = pick[0]['st'], pick[0]['tank']
        if len(pick) > 1:
            j = w_multi((st, [(r['ast'], r['t']) for r in pick]))
            if j is not None and j[1] <= a.max_tank: new = j
            else: pick = pick[:1]
        st, tk = new
        say(f'round {rnd}: {len(res)}/{len(trial)} ok, singles {" ".join(str(int(r["tank"] - tank)) for r in res[:6])} kg; '
            f'+{[r["ast"] for r in pick]} -> {len(st["asts"])} flybys, tank {tk:.0f} kg (+{tk - tank:.0f}), '
            f'{VE*np.log(tk/601.5)/len(st["asts"]):.3f} km/s/fb  ({time.time()-t0:.0f} s)')
        tank = tk
        f2 = RI.IFleet(); f2.routes[a.route] = dict(st=st, tank=tank); f2.save(out, note='grown carrier route')
    say(f'done: {len(st["asts"])} flybys, tank {tank:.0f} kg, J_i {RI.cost(tank):.3f}')


if __name__ == '__main__':
    main()
