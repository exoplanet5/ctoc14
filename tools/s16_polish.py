"""Stage 16: MULTI-START re-settle of every route of a fleet (same target set, honest regular 20 d grid).

run_ialns.settle (restore -> L1 SCP optimise -> restore -> thrust-cap homotopy) is a LOCAL method: re-settling a
route from its stored state moved tanks by -24..+5 kg in stage 15b.  Here each route is re-settled from several
starts and the lightest settled variant is kept (the selector dedupes by target set and keeps the lightest):
    s0  plain settle (iters 150)             s1  settle iters 300
    pT  impulses scaled by (1 + sigma N(0,1)) per component (sigma 0.1 / 0.25, 2 seeds each), then settle
    pF  flyby epochs jittered by N(0, 1 d) (2 seeds), then settle
    tL  launch epoch shifted by +-3 / +-7 d (>= 0), impulses re-binned onto the grid of the new launch
        (s15b_regrid.regrid), then settle
Every variant that settles (miss <= 150 km) is saved as OUT/cols/route_<k>_<tag>_<n>.npz with its tank, max
impulse/cap and node regularity in OUT/polish.jsonl; OUT/fleet = per route the lightest honest variant (or the
original).  Resumable: variants already in polish.jsonl are skipped.

usage: s16_polish.py OUT FLEET [--nproc 10] [--timeout 600] [--only h01,t1]"""
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
from s15b_regrid import regrid
from greedy_cover import _timeout
from ctoc14.impulsive import ImpulsiveProblem
from ctoc14.constants import DAY, T_MISSION
import zlib

TAGS = ['s0', 's1', 'pT10a', 'pT10b', 'pT25a', 'pT25b', 'pFa', 'pFb', 'tLp3', 'tLp7', 'tLm3', 'tLm7']


def variant(st, tag):
    ip = RI.ipr(st)
    if tag in ('s0', 's1'):
        return ip, (300 if tag == 's1' else 150)
    rng = np.random.default_rng(zlib.crc32(tag.encode()))
    if tag.startswith('pT'):
        sig = int(tag[2:4]) / 100.0
        ip.Ts = ip.Ts * (1.0 + sig * rng.standard_normal(ip.Ts.shape))
        return ip, 150
    if tag.startswith('pF'):
        tf = ip.tf + rng.normal(0.0, 1.0, len(ip.tf)) * DAY
        ip.tf = np.minimum(tf, T_MISSION - 2 * DAY)
        return ip, 150
    if tag.startswith('tL'):
        d = (1 if tag[2] == 'p' else -1) * float(tag[3:]) * DAY
        tL = ip.tL + d
        if tL < 0.0 or tL >= float(np.min(ip.tf)) - 5 * DAY:
            return None, 0
        ip2 = ImpulsiveProblem(ip.eph, tL, ip.vinf, ip.ts, ip.Ts, ip.tf, ip.asts, margin=ip.margin, dtau=ip.h)
        keep = ip2.ts >= tL                        # impulses before the new launch are folded into the first bin
        if not keep.all():
            Ts = ip2.Ts.copy(); Ts[np.argmax(keep)] += Ts[~keep].sum(axis=0)
            ip2 = ImpulsiveProblem(ip.eph, tL, ip.vinf, ip2.ts[keep], Ts[keep], ip.tf, ip.asts, margin=ip.margin, dtau=ip.h)
        return regrid(ip2, ip.h), 150
    raise ValueError(tag)


def walk(st, sign, step_d=1.0, max_steps=15, iters=150, patience=2):
    """Launch-epoch WALK: shift tL by sign x step_d, re-bin onto the new grid, settle; keep walking from the settled
    state while it settles; return the lightest settled state seen (or None) and the number of steps."""
    cur = RI.ipr(st); best = None; best_tank = float(cur.tank()); bad = 0; n = 0
    for n in range(1, max_steps + 1):
        tL = cur.tL + sign * step_d * DAY
        if tL < 0.0 or tL >= float(np.min(cur.tf)) - 5 * DAY:
            break
        ip2 = ImpulsiveProblem(cur.eph, tL, cur.vinf, cur.ts, cur.Ts, cur.tf, cur.asts, margin=cur.margin, dtau=cur.h)
        g = regrid(ip2, cur.h)
        miss = float(RI.settle(g, iters))
        if miss > 150.0:
            break
        t = float(g.tank())
        if t < best_tank - 0.05:
            best_tank = t; best = g; bad = 0
        else:
            bad += 1
            if bad >= patience:
                break
        cur = g
    return best, n


def job(args):
    k, st, tag, timeout = args
    tic = time.time(); out = dict(route=k, tag=tag)
    try:
        signal.signal(signal.SIGALRM, _timeout); signal.setitimer(signal.ITIMER_REAL, timeout)
        if tag.startswith('tLw'):
            ip, nsteps = walk(st, 1 if tag[3] == 'p' else -1, step_d=float(tag[4:]) if len(tag) > 4 else 1.0)
            if ip is None:
                out.update(ok=False, err=f'no gain ({nsteps} steps)'); return out
            Yf, _ = ip.integrate(); dm, _ = ip.misses(Yf); miss = float(np.linalg.norm(dm, axis=1).max())
            out['steps'] = nsteps
        else:
            ip, iters = variant(st, tag)
            if ip is None:
                out.update(ok=False, err='skip'); return out
            miss = float(RI.settle(ip, iters))
        ts = np.asarray(ip.ts); irr = int(np.sum(np.abs(np.diff(ts) - ip.h) > 0.01 * DAY)) if len(ts) > 1 else 0
        cap = float((np.linalg.norm(ip.Ts, axis=1) / ip.tcap).max()) if len(ip.Ts) else 0.0
        out.update(ok=bool(miss <= 150.0), miss=miss, tank=float(ip.tank()), irregular=irr, imp_cap=cap,
                   tL=float(ip.tL) / DAY, st=RI.ist(ip))
    except Exception as e:
        out.update(ok=False, err=type(e).__name__)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        out['sec'] = round(time.time() - tic, 1)
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('fleet')
    ap.add_argument('--nproc', type=int, default=10); ap.add_argument('--timeout', type=float, default=600.0)
    ap.add_argument('--only', default=''); ap.add_argument('--tags', default=','.join(TAGS))
    a = ap.parse_args()
    out = ROOT / a.out; (out / 'cols').mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt'); RI.eph()
    files = FA.fleet_files(a.fleet)
    fl = {k: dict(st=FA.load_st(f), f=f) for k, f in files.items()}
    for r in fl.values():
        r['tank'] = float(RI.ipr(r['st']).tank())
    only = [x for x in a.only.split(',') if x] or sorted(fl)
    log = out / 'polish.jsonl'; done = {}
    if log.exists():
        for l in open(log):
            try:
                m = json.loads(l); done[(m['route'], m['tag'])] = m
            except Exception:
                pass
    tags = [t for t in a.tags.split(',') if t]
    jobs = [(k, fl[k]['st'], t, a.timeout) for k in only for t in tags if (k, t) not in done]
    jobs.sort(key=lambda j: -len(j[1]['asts']))
    say(f'polish {a.fleet}: {len(fl)} routes, sum J_i {sum(RI.cost(r["tank"]) for r in fl.values()):.4f}; {len(jobs)} settles '
        f'({len(done)} cached) on {a.nproc} procs')
    with mp.get_context('fork').Pool(a.nproc, maxtasksperchild=6) as pool, open(log, 'a') as fo:
        for r in pool.imap_unordered(job, jobs, chunksize=1):
            k = r['route']; st = r.pop('st', None)
            if r.get('ok') and st is not None:
                p = out / 'cols' / f'route_{k}_{r["tag"]}_{len(st["asts"])}.npz'; np.savez(p, **st)
                r['col'] = str(p.relative_to(ROOT)); r['dkg'] = round(r['tank'] - fl[k]['tank'], 2)
            fo.write(json.dumps(r) + '\n'); fo.flush(); done[(k, r['tag'])] = r
            say(f'  {k} {r["tag"]}: ' + (f'tank {r["tank"]:.2f} ({r["dkg"]:+.2f} kg) miss {r["miss"]:.0f} km imp/cap '
                                          f'{r["imp_cap"]:.3f} irr {r["irregular"]}' if r.get('ok') else
                                          f'fail {r.get("err", "")} {r.get("miss", "")}') + f' ({r["sec"]} s)')
    F = RI.IFleet(); rep = {}
    for k, r in fl.items():
        best = None
        for (kk, t), m in done.items():
            if kk == k and m.get('ok') and m['irregular'] == 0 and m['imp_cap'] <= 1.02 and m['tank'] < r['tank'] - 1e-6:
                if best is None or m['tank'] < best['tank']:
                    best = m
        if best is not None:
            st = FA.load_st(best['col']); F.routes[k] = dict(st=st, tank=float(RI.ipr(st).tank())); rep[k] = best['tag']
        else:
            F.routes[k] = dict(st=r['st'], tank=r['tank']); rep[k] = 'orig'
    F.save(out / 'fleet', note=f's16_polish {a.fleet}: lightest honest variant per route')
    sj0 = sum(RI.cost(r['tank']) for r in fl.values()); sj = sum(RI.cost(r['tank']) for r in F.routes.values())
    json.dump(dict(fleet=a.fleet, sumJi0=round(sj0, 5), sumJi=round(sj, 5), chosen=rep), open(out / 'polish.json', 'w'), indent=1)
    say(f'polish done: sum J_i {sj0:.4f} -> {sj:.4f}; chosen {rep}')


if __name__ == '__main__':
    main()
