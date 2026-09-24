"""Stage 16 pipeline for ONE candidate fleet over the honest pool:
    close  (s15b_close.close cap ladder; a second ladder at higher caps / dmax if misses stay open)
 -> relocate to convergence (tools/s15b_relocate.py, repeated invocations until a pass makes no move)
 -> re-settle every route (tools/s15b_resettle.py) and keep, per route, the LIGHTER of stored / re-settled state
 -> honest check (irregular nodes, max impulse/cap, max 20 d window/cap, twin miss, coverage) -> OUT/final/fleet.
Every step has its own dir under OUT and is skipped when its report exists, so a killed run resumes where it stopped.

usage: s16_pipe.py OUT FLEET [--nproc 4] [--caps 0.03,0.05,0.08,0.12] [--caps2 0.2,0.3,0.45] [--dmax 0.25]
                   [--kgfb 8] [--top 40] [--rdmax 0.2] [--rounds 4] [--max-reloc 6]
FLEET = fleet dir | FRONTIER_JSON:KEY | comma list of route files (tools/s15b_fleet.fleet_files)."""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse, subprocess, multiprocessing as mp
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)
import numpy as np
import run_ialns as RI
import s15b_fleet as FA
import s15b_close as CL
from greedy_cover import ALL
from ctoc14.constants import DAY

PY = os.path.expanduser('~/.venvs/astro313/bin/python')


def honest_row(job):
    n, st = job
    ip = RI.ipr(st); ts = np.asarray(ip.ts); h = ip.h
    irr = int(np.sum(np.abs(np.diff(ts) - h) > 0.01 * DAY)) if len(ts) > 1 else 0
    nrm = np.linalg.norm(ip.Ts, axis=1) if len(ip.Ts) else np.zeros(0)
    cap = float((nrm / ip.tcap).max()) if len(nrm) else 0.0
    w = 0.0
    for s in ts:
        sel = (ts >= s) & (ts < s + h - 1.0)
        w = max(w, float(nrm[sel].sum() / ip.tcap))
    Yf, _ = ip.integrate(); dm, _ = ip.misses(Yf); miss = float(np.linalg.norm(dm, axis=1).max()) if len(dm) else 0.0
    return dict(name=n, n=len(ip.asts), irregular=irr, max_imp_cap=round(cap, 4), max_window_cap=round(w, 4),
                miss_km=round(miss, 1), tank=round(float(ip.tank()), 3), Ji=round(RI.cost(float(ip.tank())), 5))


def honest_check(fleet_dir, nproc, say):
    F = RI.IFleet(fleet_dir)
    with mp.get_context('fork').Pool(min(nproc, len(F.routes))) as pl:
        R = pl.map(honest_row, [(n, F.routes[n]['st']) for n in sorted(F.routes)])
    cov = set(F.coverage()) - {131, 144}
    ok = all(r['irregular'] == 0 and r['max_window_cap'] <= 1.02 and r['miss_km'] <= 150 for r in R)
    rep = dict(fleet=str(fleet_dir), craft=len(R), covered=len(cov), misses=sorted(set(ALL) - cov),
               sumJi=round(sum(r['Ji'] for r in R), 5), honest=ok, routes=R)
    for r in R:
        say(f'   {r["name"]:>6s} {r["n"]:2d} fb tank {r["tank"]:8.2f} J_i {r["Ji"]:.4f} irregular {r["irregular"]} '
            f'imp/cap {r["max_imp_cap"]:.3f} window/cap {r["max_window_cap"]:.3f} miss {r["miss_km"]:.0f} km')
    say(f'  honest check {fleet_dir}: {len(R)} craft, covered {len(cov)}, sum J_i {rep["sumJi"]:.4f}, honest {ok}')
    return rep


def dedupe(files, out, nproc, say, rounds=8):
    """Targets flown by 2+ routes (the MILP allows overlaps): remove each duplicate from the route where removal
    (run_ialns.w_remove: remove + settle) saves most; one removal per route per round.  -> out/fleet, out/dedupe.json"""
    out = pathlib.Path(out); (out / 'cols').mkdir(parents=True, exist_ok=True)
    fl = {k: dict(st=FA.load_st(f), f=f) for k, f in files.items()}
    for r in fl.values():
        r['tank'] = float(RI.ipr(r['st']).tank())
    sj0 = sum(RI.cost(r['tank']) for r in fl.values()); log = []
    for rnd in range(rounds):
        cov = {}
        for k, r in fl.items():
            for x in r['st']['asts']:
                cov.setdefault(int(x), []).append(k)
        dups = {x: ks for x, ks in cov.items() if len(ks) > 1}
        if not dups:
            break
        jobs = [(k, fl[k]['st'], [x]) for x, ks in sorted(dups.items()) for k in ks]
        with mp.get_context('fork').Pool(min(nproc, len(jobs)), maxtasksperchild=4) as pool:
            R = pool.map(RI.w_remove, jobs, chunksize=1)
        cand = sorted([(RI.cost(fl[r['name']]['tank']) - RI.cost(r['tank']), r) for r in R if r['ok']], key=lambda z: -z[0])
        used = set(); done = set(); n = 0
        for s, r in cand:
            x = int(r['asts'][0]); k = r['name']
            if s <= 0 or k in used or x in done:
                continue
            p = out / 'cols' / f'route_{k}m{x}_{len(r["st"]["asts"])}.npz'; np.savez(p, **r['st'])
            fl[k].update(st=r['st'], tank=r['tank'], f=str(p.relative_to(ROOT))); used.add(k); done.add(x); n += 1
            log.append(dict(ast=x, route=k, save=round(s, 4))); say(f'    dedupe: {x} removed from {k}, saves {s:.4f}')
        if not n:
            break
    F = RI.IFleet()
    for k, r in fl.items():
        F.routes[k] = dict(st=r['st'], tank=r['tank'])
    F.save(out / 'fleet', note='s16_pipe dedupe')
    sj = sum(RI.cost(r['tank']) for r in fl.values())
    json.dump(dict(sumJi0=round(sj0, 5), sumJi=round(sj, 5), removed=log, files={k: r['f'] for k, r in fl.items()}),
              open(out / 'dedupe.json', 'w'), indent=1)
    say(f'  dedupe: sum J_i {sj0:.4f} -> {sj:.4f}, {len(log)} duplicate flybys removed')


def fleet_sumji(spec):
    files = FA.fleet_files(spec)
    return sum(RI.cost(float(RI.ipr(FA.load_st(f)).tank())) for f in files.values()), files


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('fleet')
    ap.add_argument('--nproc', type=int, default=4); ap.add_argument('--caps', default='0.03,0.05,0.08,0.12')
    ap.add_argument('--caps2', default='0.2,0.3,0.45'); ap.add_argument('--dmax', type=float, default=0.25)
    ap.add_argument('--dmax2', type=float, default=0.3)
    ap.add_argument('--kgfb', type=float, default=8.0); ap.add_argument('--top', type=int, default=40)
    ap.add_argument('--rdmax', type=float, default=0.2); ap.add_argument('--rounds', type=int, default=4)
    ap.add_argument('--max-reloc', type=int, default=6); ap.add_argument('--timeout', type=float, default=300.0)
    a = ap.parse_args()
    out = ROOT / a.out; out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt'); RI.eph()
    rel = lambda p: str(pathlib.Path(p).relative_to(ROOT))
    # ---- 0. source fleet
    sj0, files = fleet_sumji(a.fleet)
    cov = set(int(x) for f in files.values() for x in FA.load_st(f)['asts'])
    M = sorted(set(ALL) - cov)
    say(f'pipe {a.fleet}: {len(files)} craft, covered {len(cov)}, sum J_i {sj0:.4f}, misses {M}')
    cur = a.fleet
    ntot = sum(len(FA.load_st(f)['asts']) for f in files.values())
    if ntot > len(cov):
        d = out / 'dedupe'
        if not (d / 'dedupe.json').exists():
            dedupe(files, d, a.nproc, say)
        cur = rel(d / 'fleet')
    # ---- 1. close
    if M:
        for tag, caps, dmax in (('close', a.caps, a.dmax), ('close2', a.caps2, a.dmax2)):
            d = out / tag
            if not (d / 'close.json').exists():
                d.mkdir(parents=True, exist_ok=True); s = RI.logger(d / 'log.txt')
                CL.close(FA.fleet_files(cur), d, [float(x) for x in caps.split(',')], dmax, a.nproc, a.timeout,
                         say=s, wall=2 * 3600.0)
            rep = json.load(open(d / 'close.json'))
            say(f'  {tag}: covered {rep["covered"]}, sum J_i {rep["sumJi0"]:.4f} -> {rep["sumJi"]:.4f}, open {rep["misses"]}')
            cur = rel(d / 'fleet')
            if not rep['misses']:
                break
        rep = json.load(open(ROOT / cur / '..' / 'close.json'))
        if rep['misses']:
            say(f'  NOT CLOSED: open {rep["misses"]}; relocating anyway (columns for the selector)')
    # ---- 2. relocate to convergence
    for k in range(1, a.max_reloc + 1):
        d = out / f'reloc{k}'
        if not (d / 'relocate.json').exists():
            cmd = [PY, 'tools/s15b_relocate.py', rel(d), cur, '--kgfb', str(a.kgfb), '--top', str(a.top), '--dmax',
                   str(a.rdmax), '--nproc', str(min(a.nproc, 8)), '--rounds', str(a.rounds), '--timeout', str(a.timeout)]
            say('  run ' + ' '.join(cmd[1:]))
            with open(out / f'reloc{k}.out', 'a') as fo:
                subprocess.run(cmd, cwd=ROOT, stdout=fo, stderr=subprocess.STDOUT, check=False)
        if not (d / 'relocate.json').exists():
            say(f'  reloc{k} failed (no report)'); break
        rep = json.load(open(d / 'relocate.json'))
        say(f'  reloc{k}: sum J_i {rep["sumJi0"]:.4f} -> {rep["sumJi"]:.4f} ({len(rep["moves"])} moves), covered {rep["covered"]}')
        cur = rel(d / 'fleet')
        if not rep['moves'] or rep['sumJi0'] - rep['sumJi'] < 1e-3:
            break
    # ---- 3. re-settle, keep the lighter settled state per route
    d = out / 'resettle'
    if not (d / 'resettle.json').exists():
        cmd = [PY, 'tools/s15b_resettle.py', rel(d), cur, '--nproc', str(min(a.nproc, 8))]
        say('  run ' + ' '.join(cmd[1:]))
        with open(out / 'resettle.out', 'a') as fo:
            subprocess.run(cmd, cwd=ROOT, stdout=fo, stderr=subprocess.STDOUT, check=False)
    rs = json.load(open(d / 'resettle.json'))
    F = RI.IFleet(); Fr = RI.IFleet(d / 'fleet'); src = FA.fleet_files(cur)
    for r in rs['routes']:
        n = r['name']; st0 = FA.load_st(src[n])
        if r['ok'] and r['tank'] < r['tank0']:
            F.routes[n] = dict(st=Fr.routes[n]['st'], tank=float(RI.ipr(Fr.routes[n]['st']).tank()))
        else:
            F.routes[n] = dict(st=st0, tank=float(RI.ipr(st0).tank()))
    F.save(out / 'final' / 'fleet', note=f's16_pipe {a.fleet}: lighter of stored / re-settled per route')
    say(f'  re-settle: stored {rs["sumJi_stored"]:.4f}, re-settled {rs["sumJi_resettled"]:.4f}; final (lighter per route) '
        f'{sum(RI.cost(r["tank"]) for r in F.routes.values()):.4f}')
    # ---- 4. honest check
    rep = honest_check(out / 'final' / 'fleet', a.nproc, say)
    rep['source'] = a.fleet; rep['sumJi_source'] = round(sj0, 5)
    json.dump(rep, open(out / 'final.json', 'w'), indent=1)
    say(f'PIPE DONE {a.fleet}: {rep["covered"]} covered, sum J_i {sj0:.4f} -> {rep["sumJi"]:.4f}, honest {rep["honest"]}')


if __name__ == '__main__':
    main()
