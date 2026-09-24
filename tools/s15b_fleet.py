"""s15b: fleet analysis -- per-route depth / tank / kg per flyby / J_i, the fleet's misses, and every miss's approaches
(local distance minima <= --dmax AU of the route's Kepler-arc trajectory, run_ialns.w_cands, not within 10 d of one of
the route's own flybys), optionally priced with the linearised whole-route twin (s14_twinbeam.lin_price, as in
s15_isogen.absorb_targets).

FLEET spec: a fleet dir (route_*.npz) | FRONTIER_JSON:K (the picks of frontier K) | comma-separated route files.

usage: s15b_fleet.py FLEET [--dmax 0.3] [--price] [--out JSON] [--nproc 4]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)
import numpy as np
import run_ialns as RI
from greedy_cover import ALL
from ctoc14.constants import DAY, VE


def load_st(f):
    z = np.load(ROOT / f); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL'])
    return st


def fleet_files(spec):
    """{name: route file (repo-relative)} of a fleet spec."""
    if ':' in spec and spec.split(':')[0].endswith('.json'):
        fj, K = spec.rsplit(':', 1)
        fr = json.load(open(ROOT / fj))[K]
        return {f'{i + 1:02d}': p for i, p in enumerate(fr['picks'])}
    p = ROOT / spec
    if p.is_dir():
        return {f.stem[6:]: str(f.relative_to(ROOT)) for f in sorted(p.glob('route_*.npz'))}
    return {f'{i + 1:02d}': x for i, x in enumerate(spec.split(','))}


def _price(job):
    """lin_price of every approach of one host: [(ast, t, dist, lin_kg)]."""
    import s14_twinbeam as TB
    name, st, apps = job
    ip = RI.ipr(st); tank = float(ip.tank())
    par = dict(st=st, nreg=len(TB.centres(float(st['tL']), float(np.max(st['tf'])))), dv=float(ip.dv()))
    out = []
    for a in apps:
        best = None
        for dt in (0.0, -4.0, 4.0):
            try:
                dvm, lin, cm = TB.lin_price(par, a['ast'], a['t'] + dt * DAY)
            except Exception:
                continue
            if np.isfinite(dvm) and (best is None or dvm < best[0]):
                best = (float(dvm), a['t'] + dt * DAY)
        kg = None if best is None else float(tank * (np.exp(max(best[0], 0.0) / VE) - 1.0))
        out.append(dict(host=name, ast=a['ast'], t=a['t'] if best is None else best[1], dist=a['dist'], lin_kg=kg))
    return out


def analyse(files, dmax=0.3, price=False, nproc=4, misses=None):
    sts = {k: load_st(f) for k, f in files.items()}
    tanks = {k: float(RI.ipr(s).tank()) for k, s in sts.items()}
    cov = {}
    for k, s in sts.items():
        for a in s['asts']:
            cov.setdefault(int(a), []).append(k)
    M = sorted(set(ALL) - set(cov)) if misses is None else sorted(misses)
    routes = {}
    for k, s in sts.items():
        n = len(s['asts']); t = tanks[k]
        routes[k] = dict(f=files[k], n=n, tank=round(t, 2), kg_fb=round((t - 600.0) / n, 2), Ji=round(RI.cost(t), 4),
                         t_launch_d=round(float(s['tL']) / DAY, 1), t_last_d=round(float(np.max(s['tf'])) / DAY, 1),
                         targets=sorted(int(a) for a in s['asts']))
    apps = []
    if M:
        with mp.get_context('fork').Pool(min(nproc, len(sts))) as pool:
            R = pool.map(RI.w_cands, [(k, s, M, dmax) for k, s in sts.items()])
        for (k, s), lst in zip(sts.items(), R):
            tf = np.asarray(s['tf'], float)
            apps += [dict(a, host=k) for a in lst if np.abs(tf - a['t']).min() > 10 * DAY]
        if price and apps:
            by = {}
            for a in apps:
                by.setdefault(a['host'], []).append(a)
            with mp.get_context('fork').Pool(min(nproc, len(by))) as pool:
                R = pool.map(_price, [(k, sts[k], v) for k, v in by.items()])
            apps = [x for lst in R for x in lst]
    near = {}
    for X in M:
        lst = sorted([a for a in apps if a['ast'] == X], key=lambda a: a['dist'])
        near[X] = [dict(host=a['host'], dist=round(a['dist'], 4), t_d=round(a['t'] / DAY, 1),
                        lin_kg=None if a.get('lin_kg') is None else round(a['lin_kg'], 1)) for a in lst]
    sj = sum(r['Ji'] for r in routes.values())
    return dict(craft=len(sts), covered=len(cov), sumJi=round(sj, 4), misses=M, routes=routes, near=near,
                dup=sum(len(v) - 1 for v in cov.values()))


def show(A, say=print):
    say(f'{A["craft"]} craft, covered {A["covered"]} (dup {A["dup"]}), sum J_i {A["sumJi"]:.4f}; misses {A["misses"]}')
    for k, r in sorted(A['routes'].items(), key=lambda kv: -kv[1]['n']):
        say(f'   {k:>6s}: {r["n"]:2d} fb @ {r["tank"]:7.1f} kg  {r["kg_fb"]:5.1f} kg/fb  J_i {r["Ji"]:.4f}  '
            f'[{r["t_launch_d"]:.0f}-{r["t_last_d"]:.0f} d]  {r["f"]}')
    for X, lst in A['near'].items():
        say(f'   miss {X:3d}: ' + ('; '.join(f'{a["host"]} {a["dist"]:.3f} AU @{a["t_d"]:.0f}d'
                                             + (f' {a["lin_kg"]:.0f} kg' if a['lin_kg'] is not None else '')
                                             for a in lst[:5]) or 'no approach'))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('fleet'); ap.add_argument('--dmax', type=float, default=0.3)
    ap.add_argument('--price', action='store_true'); ap.add_argument('--out', default=None)
    ap.add_argument('--nproc', type=int, default=4)
    a = ap.parse_args()
    RI.eph()
    A = analyse(fleet_files(a.fleet), a.dmax, a.price, a.nproc)
    show(A)
    if a.out:
        json.dump(A, open(ROOT / a.out, 'w'), indent=1)
