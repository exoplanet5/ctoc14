"""s15b: RE-SETTLE every route of a fleet (run_ialns.ipr + run_ialns.settle: restore, L1 SCP, restore, thrust-cap
homotopy) and confirm the tanks: a route is confirmed when the re-settled max miss <= 150 km.  The re-settled fleet
(re-settled state when it settles, else the original state) is saved as OUT/fleet (RI.IFleet) with OUT/resettle.json.

usage: s15b_resettle.py OUT FLEET [--nproc 4] [--iters 100] [--timeout 900]
FLEET = fleet dir | FRONTIER_JSON:K | comma list of route files."""
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
from greedy_cover import ALL, _timeout


def _resettle(job):
    name, f, iters, timeout = job
    st = FA.load_st(f); ip0 = RI.ipr(st); tank0 = float(ip0.tank())
    Yf, _ = ip0.integrate(); dm, _ = ip0.misses(Yf); miss0 = float(np.linalg.norm(dm, axis=1).max())
    tic = time.time(); ip = RI.ipr(st)
    try:
        signal.signal(signal.SIGALRM, _timeout); signal.setitimer(signal.ITIMER_REAL, timeout)
        miss = float(RI.settle(ip, iters))
    except Exception as e:
        return dict(name=name, f=f, ok=False, tank0=tank0, miss0=miss0, err=type(e).__name__, sec=round(time.time() - tic))
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    return dict(name=name, f=f, ok=miss <= 150.0, tank0=tank0, miss0=miss0, tank=float(ip.tank()), miss=miss,
                st=RI.ist(ip), n=len(ip.asts), sec=round(time.time() - tic))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('fleet')
    ap.add_argument('--nproc', type=int, default=4); ap.add_argument('--iters', type=int, default=100)
    ap.add_argument('--timeout', type=float, default=900.0)
    a = ap.parse_args()
    if a.nproc > 8:
        raise SystemExit('CPU budget: --nproc <= 8')
    out = ROOT / a.out; out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt'); RI.eph()
    files = FA.fleet_files(a.fleet)
    with mp.get_context('fork').Pool(min(a.nproc, len(files)), maxtasksperchild=1) as pool:
        R = pool.map(_resettle, [(k, f, a.iters, a.timeout) for k, f in files.items()], chunksize=1)
    F = RI.IFleet(); rows = []
    for r in R:
        if r['ok']:
            F.routes[r['name']] = dict(st=r['st'], tank=r['tank'])
        else:
            st = FA.load_st(r['f']); F.routes[r['name']] = dict(st=st, tank=r['tank0'])
        rows.append({k: v for k, v in r.items() if k != 'st'})
        say(f'  {r["name"]}: {r.get("n", "?")} fb, stored tank {r["tank0"]:.2f} (miss {r["miss0"]:.0f} km) -> re-settled '
            + (f'{r["tank"]:.2f} (miss {r["miss"]:.0f} km) {"OK" if r["ok"] else "NOT settled"}' if 'tank' in r else r.get('err', ''))
            + f' ({r["sec"]} s)')
    cov = F.coverage(); sj = sum(RI.cost(r['tank']) for r in F.routes.values())
    sj0 = sum(RI.cost(r['tank0']) for r in rows)
    F.save(out / 'fleet', note=f're-settled {a.fleet}')
    rep = dict(fleet=a.fleet, craft=len(F.routes), covered=len(cov), misses=sorted(set(ALL) - set(cov)),
               sumJi_stored=round(sj0, 4), sumJi_resettled=round(sj, 4), all_ok=all(r['ok'] for r in rows), routes=rows)
    json.dump(rep, open(out / 'resettle.json', 'w'), indent=1)
    say(f'resettle: {len(F.routes)} craft, covered {len(cov)}, sum J_i stored {sj0:.4f} -> re-settled {sj:.4f}; '
        f'all settled: {rep["all_ok"]}')


if __name__ == '__main__':
    main()
