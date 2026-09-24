"""Mass-produce DEEP carrier routes for the set-cover pool (stage 7, section 9).

One route = one carrier from the DP library, its base settled, then greedy insertion growth (accept the first candidate
whose tank cost is within --max-step, nearest first) until nothing is affordable.  Each worker gets a RANDOM SUBSET of
the catalogue (--frac) so that routes built from similar carriers still diverge -- diversity is what the set cover
(tools/route_cover.py) needs, not individual depth.  The measured cost of the best such route is 0.034 J per flyby
against t10d's fleet 0.0415, so a disjoint packing of them is worth J ~ 12.1-12.9.

Usage: carrier_farm.py lib.pkl out_dir [--nroute 80] [--frac 0.6] [--max-tank 1100] [--nproc 8]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pickle, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI, carrier_dp as CD
from ctoc14.globalopt import insert_block
from ctoc14.search import UNREACHABLE
from ctoc14.constants import DAY
ALL = [t for t in range(1, 301) if t not in UNREACHABLE]
_A = None


def base_of(c, allow):
    """settle the DP base restricted to `allow`, with leave-one-out retries."""
    full = [b for b in c['base'] if b[0] in allow]
    if len(full) < 4: return None
    for drop in [None] + list(range(len(full))):
        base = [b for i, b in enumerate(full) if i != drop]
        try:
            ip = CD.seed_ip(RI.eph(), c['tL'], np.array(c['vinf']),
                            np.array([b[1] for b in base]), np.array([b[2] for b in base]))
            if insert_block(ip, [(b[0], b[1]) for b in base], log=RI.QUIET).max() > 1e4: continue
            if RI.settle(ip, 100) <= 150: return ip
        except Exception:
            continue
    return None


KG_OF = lambda d: 4.0 + 600.0 * max(0.0, d - 0.02)      # kg per insertion vs closest approach (measured calibration)


def w_route(job):
    k, c, seed = job; a = _A
    rng = np.random.default_rng(seed)
    PI = _A.pi
    univ = a.allow if a.allow else ALL
    allow = set(int(x) for x in rng.choice(univ, size=int(a.frac * len(univ)), replace=False)) if a.frac < 1 else set(univ)
    allow |= set(b[0] for b in c['base'][:3])
    t0 = time.time()
    ip = base_of(c, allow)
    if ip is None: return k, None
    st = RI.ist(ip); tank = float(ip.tank()); nfail = 0
    while len(st['asts']) < a.depth and time.time() - t0 < a.tmax:
        have = set(int(x) for x in st['asts'])
        left = [t for t in allow if t not in have]
        if not left: break
        cands = RI.w_cands(('h', st, left, a.dmax))
        tf = np.asarray(st['tf'])
        best = {}
        for x in sorted(cands, key=lambda x: x['dist']):
            if np.abs(tf - x['t']).min() > a.gap * DAY: best.setdefault(x['ast'], x)
        if PI is None:
            trial = sorted(best.values(), key=lambda x: x['dist'])[:a.ntrial]
        else:
            # reduced cost: the target's LP dual minus the J the insertion is expected to cost
            dJdkg = (1.0 + 2.0 * (tank - 600.0) / 1400.0) / 1400.0
            sc = {x['ast']: PI.get(x['ast'], 0.0) - dJdkg * KG_OF(x['dist']) for x in best.values()}
            trial = [x for x in sorted(best.values(), key=lambda x: -sc[x['ast']]) if sc[x['ast']] > 0][:a.ntrial]
        if not trial: break
        got = False
        for x in trial:
            r = RI.w_insert(('h', st, x['ast'], x['t']))
            gain = PI.get(x['ast'], 0.0) if PI is not None else 1.0
            lim = a.max_step if PI is None else min(a.max_step_pi, gain / max(1e-9, (1.0 + 2.0 * (tank - 600.0) / 1400.0) / 1400.0))
            if r.get('ok') and r['tank'] <= a.max_tank and r['tank'] - tank <= lim:
                st, tank = r['st'], r['tank']; got = True; break
        if not got:
            nfail += 1
            if nfail >= a.nfail: break
        else:
            nfail = 0
    return k, (st, tank, time.time() - t0)


def main():
    global _A
    ap = argparse.ArgumentParser(); ap.add_argument('lib'); ap.add_argument('out')
    ap.add_argument('--nroute', type=int, default=80); ap.add_argument('--frac', type=float, default=0.6)
    ap.add_argument('--depth', type=int, default=45); ap.add_argument('--max-tank', type=float, default=1100.0)
    ap.add_argument('--max-step', type=float, default=45.0); ap.add_argument('--dmax', type=float, default=0.15)
    ap.add_argument('--ntrial', type=int, default=6); ap.add_argument('--gap', type=float, default=12.0)
    ap.add_argument('--tmax', type=float, default=900.0); ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--seed', type=int, default=0); ap.add_argument('--top', type=int, default=400)
    ap.add_argument('--allow-file', default='', help='json list: grow only on these targets (fix-and-optimise)')
    ap.add_argument('--duals', default=''); ap.add_argument('--pi-floor', type=float, default=0.030,
                    help='J value given to EVERY target on top of its dual (a deep route is worth ~0.034 J per flyby)')
    ap.add_argument('--max-step-pi', type=float, default=120.0)
    ap.add_argument('--nfail', type=int, default=3); ap.add_argument('--by-base', action='store_true',
                    help='rank the library by base depth instead of DP value (deep bases -> deep routes)')
    a = ap.parse_args(); _A = a
    a.allow = json.load(open(a.allow_file)) if a.allow_file else None
    a.pi = {int(k): v + a.pi_floor for k, v in json.load(open(a.duals)).items()} if a.duals else None
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt')
    lib = pickle.load(open(a.lib, 'rb'))
    if a.allow is not None:                    # rank carriers by how much of the REMAINING set they reach
        A = set(a.allow)
        for x in lib:
            x['pval'] = sum(1 for t, _, _ in x['base'] if t in A) + \
                        0.5 * sum(1 for t, (d, _) in x['neigh'].items() if d < 0.06 and t in A)
        lib = sorted(lib, key=lambda x: -x['pval'])[:a.top]
    elif a.pi is not None:                     # rank carriers by the dual value they can reach
        for x in lib:
            x['pval'] = sum(a.pi.get(t, 0.0) for t, _, _ in x['base']) + \
                        sum(a.pi.get(t, 0.0) for t, (d, _) in x['neigh'].items() if d < 0.06)
        lib = sorted(lib, key=lambda x: -x['pval'])[:a.top]
    else:
        lib = sorted(lib, key=(lambda x: -len(x['base'])) if a.by_base else (lambda x: -x['val']))[:a.top]
    rng = np.random.default_rng(a.seed)
    idx = rng.choice(len(lib), size=min(a.nroute, len(lib)), replace=False)
    jobs = [(int(k), lib[int(k)], a.seed * 1000 + i) for i, k in enumerate(idx)]
    say(f'{len(jobs)} routes from {len(lib)} carriers, frac {a.frac}, max-tank {a.max_tank}, step {a.max_step}')
    t0 = time.time(); n = 0
    with mp.get_context('fork').Pool(a.nproc) as pl:
        for k, r in pl.imap_unordered(w_route, jobs):
            if r is None: continue
            st, tank, dt = r; n += 1
            np.savez(out / f'route_c{k:05d}.npz', **st)
            say(f'[{n:3d}/{len(jobs)}] carrier {k}: {len(st["asts"]):2d} fb @ {tank:6.0f} kg, '
                f'{RI.cost(tank)/len(st["asts"]):.4f} J/fb  ({dt:.0f} s, total {time.time()-t0:.0f} s)')
    say(f'done: {n} routes in {time.time()-t0:.0f} s')


if __name__ == '__main__':
    main()
