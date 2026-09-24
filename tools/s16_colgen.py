"""Stage 16: MOVE-COLUMN generation around a fleet, for MILP-coordinated relocation.

s15b_relocate applies single-target moves greedily (one move per route per round, top --top removals, <= 3 hosts per
target) and saves only the removal columns of APPLIED moves.  Here every candidate move is saved as a column so the
exact selector (tools/s16_select.py) can pick the best NON-CONFLICTING set of moves together with any pool column:
  1. removals   every target x of every route k (run_ialns.w_remove: remove + settle) -> OUT/cols/route_<k>m<x>_<n>.npz
  2. insertions every x whose removal saves >= --min-save J, into every OTHER route passing within --dmax AU
                (s15b_close._price_host lin_price screen: lin dJ < --screen x saving, <= --per-target hosts, cheapest
                first), twin insertion s15b_close._insert -> OUT/cols/route_<h>x<x>_<n>.npz
  3. the fleet's own routes are copied to OUT/cols/route_<k>_base_<n>.npz (so the selector always sees them).
Everything is saved as it settles (OUT/moves.jsonl), and a restarted run skips the jobs already in moves.jsonl.

usage: s16_colgen.py OUT FLEET [--nproc 10] [--min-save 0.004] [--dmax 0.2] [--screen 1.5] [--per-target 4]
FLEET = fleet dir | FRONTIER_JSON:KEY | comma list of route files."""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, signal, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)
import numpy as np
import run_ialns as RI
import s15b_fleet as FA
import s15b_close as CL
from greedy_cover import _timeout


def w_rem(job):
    k, st, x, timeout = job
    tic = time.time()
    try:
        signal.signal(signal.SIGALRM, _timeout); signal.setitimer(signal.ITIMER_REAL, timeout)
        r = RI.w_remove((k, st, [x]))
    except Exception as e:
        r = dict(name=k, asts=[x], ok=False, error=type(e).__name__)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    r['sec'] = round(time.time() - tic, 1)
    return r


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('fleet')
    ap.add_argument('--nproc', type=int, default=10); ap.add_argument('--min-save', type=float, default=0.004)
    ap.add_argument('--dmax', type=float, default=0.2); ap.add_argument('--screen', type=float, default=1.5)
    ap.add_argument('--per-target', type=int, default=4); ap.add_argument('--timeout', type=float, default=300.0)
    a = ap.parse_args()
    out = ROOT / a.out; (out / 'cols').mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt'); RI.eph()
    files = FA.fleet_files(a.fleet)
    fl = {k: dict(st=FA.load_st(f), f=f) for k, f in files.items()}
    for k, r in fl.items():
        r['tank'] = float(RI.ipr(r['st']).tank())
        p = out / 'cols' / f'route_{k}_base_{len(r["st"]["asts"])}.npz'
        if not p.exists():
            np.savez(p, **r['st'])
    sj0 = sum(RI.cost(r['tank']) for r in fl.values())
    say(f'colgen {a.fleet}: {len(fl)} routes, sum J_i {sj0:.4f}')
    log = out / 'moves.jsonl'; done = {}
    if log.exists():
        for l in open(log):
            try:
                m = json.loads(l); done[(m['kind'], m['route'], m['ast'])] = m
            except Exception:
                pass
    tic = time.time()
    # ---- 1. removals
    jobs = [(k, r['st'], int(x), a.timeout) for k, r in fl.items() for x in r['st']['asts']
            if ('rem', k, int(x)) not in done]
    say(f'  removals: {len(jobs)} to settle ({sum(1 for q in done if q[0] == "rem")} cached)')
    with mp.get_context('fork').Pool(a.nproc, maxtasksperchild=10) as pool, open(log, 'a') as fo:
        for r in pool.imap_unordered(w_rem, jobs, chunksize=1):
            k = r['name']; x = int(r['asts'][0])
            m = dict(kind='rem', route=k, ast=x, ok=bool(r.get('ok')), sec=r['sec'])
            if m['ok']:
                p = out / 'cols' / f'route_{k}m{x}_{len(r["st"]["asts"])}.npz'; np.savez(p, **r['st'])
                m.update(tank=float(r['tank']), save=float(RI.cost(fl[k]['tank']) - RI.cost(r['tank'])),
                         col=str(p.relative_to(ROOT)))
            fo.write(json.dumps(m) + '\n'); fo.flush(); done[('rem', k, x)] = m
    rem = [m for q, m in done.items() if q[0] == 'rem' and m['ok']]
    sav = sorted([m for m in rem if m['save'] >= a.min_save], key=lambda m: -m['save'])
    say(f'  removals settled: {len(rem)}; saving >= {a.min_save}: {len(sav)}; top ' +
        ', '.join(f'{m["ast"]}@{m["route"]} {m["save"]:.3f}' for m in sav[:15]) + f' ({time.time() - tic:.0f} s)')
    # ---- 2. insertions of the saving targets into the other routes
    by_host = {}
    for m in sav:
        for h in fl:
            if h != m['route']:
                by_host.setdefault(h, set()).add(m['ast'])
    with mp.get_context('fork').Pool(a.nproc, maxtasksperchild=4) as pool:
        P = pool.map(CL._price_host, [(h, fl[h]['st'], sorted(xs), a.dmax) for h, xs in by_host.items()], chunksize=1)
    priced = [c for lst in P for c in lst]
    src = {m['ast']: m for m in sav}
    att = []
    for x, m in src.items():
        cs = sorted([c for c in priced if c['ast'] == x and c['host'] != m['route'] and c['lin_dJ'] < a.screen * m['save']],
                    key=lambda c: c['lin_dJ'])
        seen = set()
        for c in cs:
            if c['host'] in seen:
                continue
            seen.add(c['host'])
            if ('ins', c['host'], x) not in done:
                att.append(c)
            if len(seen) >= a.per_target:
                break
    say(f'  {len(priced)} priced approaches; {len(att)} insertion attempts ({time.time() - tic:.0f} s)')
    jobs = [(c['host'], fl[c['host']]['st'], c['ast'], c['t'], c['dist'], a.timeout, fl[c['host']]['tank']) for c in att]
    nok = 0
    with mp.get_context('fork').Pool(a.nproc, maxtasksperchild=4) as pool, open(log, 'a') as fo:
        for r in pool.imap_unordered(CL._insert, jobs, chunksize=1):
            m = dict(kind='ins', route=r['host'], ast=int(r['ast']), ok=bool(r['ok']), dist=r['dist'], sec=r['sec'])
            if r['ok']:
                p = out / 'cols' / f'route_{r["host"]}x{r["ast"]}_{len(r["st"]["asts"])}.npz'; np.savez(p, **r['st'])
                m.update(tank=float(r['tank']), dJ=float(r['dJ']), col=str(p.relative_to(ROOT)),
                         gain=float(src[int(r['ast'])]['save'] - r['dJ']))
                nok += 1
            fo.write(json.dumps(m) + '\n'); fo.flush(); done[('ins', m['route'], m['ast'])] = m
    ins = sorted([m for q, m in done.items() if q[0] == 'ins' and m['ok']], key=lambda m: -m.get('gain', -9))
    say(f'  insertions settled: {len(ins)} ({nok} new); best single-move gains ' +
        ', '.join(f'{m["ast"]}->{m["route"]} {m["gain"]:+.4f}' for m in ins[:15]) + f' ({time.time() - tic:.0f} s)')
    json.dump(dict(fleet=a.fleet, sumJi=round(sj0, 5), removals=len(rem), insertions=len(ins),
                   positive=[m for m in ins if m.get('gain', -1) > 0]), open(out / 'colgen.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
