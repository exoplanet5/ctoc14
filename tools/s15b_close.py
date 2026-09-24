"""s15b: CAP-LADDER twin closer for a 9-craft fleet.

The fleet's misses are offered to every route passing within --dmax AU (run_ialns.w_cands approach minima, not within
10 d of the route's own flybys), priced by the linearised whole-route twin (s14_twinbeam.lin_price, the insertion screen
of s15_isogen.absorb_targets).  Caps (J units per placement) are climbed in order (default 0.03, 0.05, 0.08, 0.12); at
each cap, rounds of:
    attempts  every (host, miss) pair with lin dJ <= --screen x cap + 0.01, at most --per-host per host (cheapest
              first) and --per-miss hosts per miss, all in one fork Pool (<= --nproc); twin insertion = run_ialns.w_insert
              (approach <= 0.08 AU) or s15_isogen.w_insert_long (farther), SIGALRM --timeout per attempt;
    commit    successes sorted by true dJ; a success is committed when dJ <= cap, its host has no placement yet in this
              round and its miss is still open (ONE placement per host per round);
    re-price  only the hosts that changed (SEQUENTIAL pricing per host), then the next round at the same cap.
A cap with no commit hands over to the next cap.  EVERY settled insertion is saved as a column (OUT/cols/route_*.npz)
for the selector; the closed fleet is OUT/fleet (RI.IFleet), the report OUT/close.json, the log OUT/log.txt.

usage: s15b_close.py OUT FLEET [--caps 0.03,0.05,0.08,0.12] [--dmax 0.2] [--nproc 4] [--timeout 300]
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
from ctoc14.constants import DAY, VE


def _price_host(job):
    """Approaches of the misses to one host + lin_price: [dict(host, ast, t, dist, lin_kg, lin_dJ)]."""
    import s14_twinbeam as TB
    host, st, misses, dmax = job
    try:
        apps = RI.w_cands((host, st, misses, dmax))
    except Exception:
        return []
    tf = np.asarray(st['tf'], float)
    apps = [a for a in apps if np.abs(tf - a['t']).min() > 10 * DAY]
    by = {}
    for a in sorted(apps, key=lambda a: a['dist']):
        by.setdefault(a['ast'], [])
        if len(by[a['ast']]) < 3:
            by[a['ast']].append(a)
    ip = RI.ipr(st); tank = float(ip.tank())
    par = dict(st=st, nreg=len(TB.centres(float(st['tL']), float(np.max(st['tf'])))), dv=float(ip.dv()))
    out = []
    for X, lst in by.items():
        for a in lst:
            best = None
            for dt in (0.0, -4.0, 4.0):
                try:
                    dvm, lin, cm = TB.lin_price(par, X, a['t'] + dt * DAY)
                except Exception:
                    continue
                if np.isfinite(dvm) and (best is None or dvm < best[0]):
                    best = (float(dvm), a['t'] + dt * DAY)
            if best is None:
                continue
            kg = float(tank * (np.exp(max(best[0], 0.0) / VE) - 1.0))
            out.append(dict(host=host, ast=int(X), t=best[1], dist=float(a['dist']), lin_kg=kg,
                            lin_dJ=float(RI.cost(tank + kg) - RI.cost(tank))))
    return out


def _insert(job):
    import s15_isogen as IG
    host, st, X, t, dist, timeout, tank0 = job
    tic = time.time()
    try:
        signal.signal(signal.SIGALRM, _timeout); signal.setitimer(signal.ITIMER_REAL, timeout)
        r = IG.w_insert_long((host, st, X, t)) if dist > 0.08 else RI.w_insert((host, st, X, t))
    except Exception as e:
        r = dict(ok=False, error=type(e).__name__)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    out = dict(host=host, ast=int(X), t=float(t), dist=float(dist), ok=bool(r.get('ok')), sec=round(time.time() - tic, 1))
    if out['ok']:
        out.update(st=r['st'], tank=float(r['tank']), miss=float(r['miss']),
                   dJ=float(RI.cost(r['tank']) - RI.cost(tank0)), dkg=float(r['tank'] - tank0))
    else:
        out.update(err=str(r.get('miss', r.get('error', '')))[:40])
    return out


def close(files, out, caps=(0.03, 0.05, 0.08, 0.12), dmax=0.2, nproc=4, timeout=300.0, screen=1.6, per_host=2,
          per_miss=3, say=print, wall=3 * 3600.0):
    out = pathlib.Path(out); (out / 'cols').mkdir(parents=True, exist_ok=True)
    fleet = {k: dict(st=FA.load_st(f), f=f) for k, f in files.items()}
    for k, r in fleet.items():
        r['tank'] = float(RI.ipr(r['st']).tank())
    cov = set(int(a) for r in fleet.values() for a in r['st']['asts'])
    M = sorted(set(ALL) - cov)
    sj0 = sum(RI.cost(r['tank']) for r in fleet.values())
    say(f'close: {len(fleet)} craft, covered {len(cov)}, sum J_i {sj0:.4f}, misses {M}; caps {list(caps)} dmax {dmax}')
    priced = {}; dirty = set(fleet); placed = []; tried = set(); ncol = 0; tic = time.time()
    with mp.get_context('fork').Pool(nproc, maxtasksperchild=4) as pool:
        for cap in caps:
            rnd = 0
            while M and time.time() - tic < wall:
                rnd += 1
                if dirty:
                    R = pool.map(_price_host, [(k, fleet[k]['st'], M, dmax) for k in sorted(dirty)], chunksize=1)
                    for k, lst in zip(sorted(dirty), R):
                        priced[k] = lst
                    dirty = set()
                cands = [c for k in priced for c in priced[k] if c['ast'] in M]
                cands.sort(key=lambda c: c['lin_dJ'])
                att = []; nh = {}; nm = {}
                for c in cands:
                    key = (c['host'], c['ast'], round(c['t'] / DAY), len(fleet[c['host']]['st']['asts']))
                    if key in tried or c['lin_dJ'] > screen * cap + 0.01:
                        continue
                    if nh.get(c['host'], 0) >= per_host or nm.get(c['ast'], 0) >= per_miss:
                        continue
                    if any(a['host'] == c['host'] and a['ast'] == c['ast'] for a in att):
                        continue
                    att.append(c); nh[c['host']] = nh.get(c['host'], 0) + 1; nm[c['ast']] = nm.get(c['ast'], 0) + 1
                    tried.add(key)
                say(f'  cap {cap} round {rnd}: open {M}; {len(cands)} priced approaches, {len(att)} attempts: ' +
                    ', '.join(f'{a["ast"]}->{a["host"]} {a["dist"]:.3f}AU lin {a["lin_dJ"]:.3f}' for a in att))
                if not att:
                    break
                res = pool.map(_insert, [(a['host'], fleet[a['host']]['st'], a['ast'], a['t'], a['dist'], timeout,
                                          fleet[a['host']]['tank']) for a in att], chunksize=1)
                for a, r in zip(att, res):
                    r['lin_dJ'] = a['lin_dJ']
                    if r['ok']:
                        n = len(r['st']['asts'])
                        p = out / 'cols' / f'route_{r["host"]}x{r["ast"]}_{n}.npz'; i = 1
                        while p.exists():
                            p = out / 'cols' / f'route_{r["host"]}x{r["ast"]}_{n}_{i}.npz'; i += 1
                        np.savez(p, **r['st']); r['col'] = str(p.relative_to(ROOT)); ncol += 1
                    say(f'    {r["ast"]} -> {r["host"]} ({r["dist"]:.3f} AU, lin {a["lin_dJ"]:.3f}): ' +
                        (f'OK +{r["dkg"]:.1f} kg dJ {r["dJ"]:.4f} miss {r["miss"]:.0f} km' if r['ok'] else f'fail {r["err"]}')
                        + f' ({r["sec"]:.0f} s)')
                ok = sorted([r for r in res if r['ok']], key=lambda r: r['dJ'])
                used = set(); com = 0
                for r in ok:
                    if r['dJ'] > cap + 1e-9 or r['host'] in used or r['ast'] not in M:
                        continue
                    h = fleet[r['host']]
                    h.update(st=r['st'], tank=r['tank'], f=r['col']); used.add(r['host']); dirty.add(r['host'])
                    M.remove(r['ast']); com += 1
                    placed.append(dict(ast=r['ast'], host=r['host'], dJ=round(r['dJ'], 4), dkg=round(r['dkg'], 1), cap=cap,
                                       col=r['col']))
                    say(f'    COMMIT {r["ast"]} -> {r["host"]} dJ {r["dJ"]:.4f} (cap {cap})')
                if not com:
                    break
    F = RI.IFleet()
    for k, r in fleet.items():
        F.routes[k] = dict(st=r['st'], tank=r['tank'])
    F.save(out / 'fleet', note=f's15b_close caps {list(caps)}')
    sj = sum(RI.cost(r['tank']) for r in fleet.values())
    rep = dict(craft=len(fleet), covered=len(F.coverage()), sumJi=round(sj, 4), sumJi0=round(sj0, 4), misses=M,
               placed=placed, columns=ncol, files={k: r['f'] for k, r in fleet.items()}, wall_s=round(time.time() - tic))
    json.dump(rep, open(out / 'close.json', 'w'), indent=1)
    say(f'close done: covered {rep["covered"]}, sum J_i {sj0:.4f} -> {sj:.4f}, placed {[(p["ast"], p["host"]) for p in placed]}, '
        f'open {M} ({rep["wall_s"]} s)')
    return rep


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('fleet')
    ap.add_argument('--caps', default='0.03,0.05,0.08,0.12'); ap.add_argument('--dmax', type=float, default=0.2)
    ap.add_argument('--nproc', type=int, default=4); ap.add_argument('--timeout', type=float, default=300.0)
    ap.add_argument('--screen', type=float, default=1.6); ap.add_argument('--per-host', type=int, default=2)
    ap.add_argument('--per-miss', type=int, default=3); ap.add_argument('--wall', type=float, default=3 * 3600.0)
    a = ap.parse_args()
    if a.nproc > 8:
        raise SystemExit('CPU budget: --nproc <= 8')
    out = ROOT / a.out; out.mkdir(parents=True, exist_ok=True)
    say = RI.logger(out / 'log.txt'); RI.eph()
    close(FA.fleet_files(a.fleet), out, [float(x) for x in a.caps.split(',')], a.dmax, a.nproc, a.timeout, a.screen,
          a.per_host, a.per_miss, say, a.wall)
