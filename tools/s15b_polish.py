"""s15b: FUEL-POLISH + MISS-ABSORBING re-plan jobs for tools/s15_twintail.py.

For a 9-craft fleet F with misses M, every route r with kg/flyby > --kgfb (or the --routes given) is re-planned ALONE by
the twin-priced beam over pool = targets(r) UNION {m in M passing within --dmax AU of r}, premium p on those misses,
fuel weight w_fuel (one job per value of --wf), every other route fixed (fleet_pre).  The beam re-orders / re-times the
route's own targets for fuel and takes the near misses as DESIGNED-IN targets; its whole front (depth >= n(r) - --slack)
is saved as columns, so the selector can use a lighter route over the same targets or one that adds a miss.

usage: s15b_polish.py QUEUE.json FLEET [FLEET ...] [--kgfb 10.5] [--routes h05,s1] [--dmax 0.15] [--p 0.02]
                      [--wf 1.5,3.0] [--beam 24] [--tries 72] [--slack 3] [--stones] [--append] [--prefix P]
FLEET = fleet dir | FRONTIER_JSON:K | comma list of route files (tools/s15b_fleet.py)."""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, glob, pathlib, argparse, hashlib
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)
import numpy as np
import s15b_fleet as FA
from greedy_cover import ALL


def existing_keys():
    """(replaced route, pool, w_fuel, p) of every twintail job already defined anywhere (s15 + s15b job files)."""
    keys = set()
    for f in glob.glob(str(ROOT / 'results/s15/twin_jobs/*.json')) + glob.glob(str(ROOT / 'results/s15b/jobs/*.json')):
        try:
            J = json.load(open(f))
        except Exception:
            continue
        if J.get('replaced'):
            keys.add((J['replaced'], tuple(sorted(J['pool'])), float(J.get('w_fuel') or 0.52), float(J.get('p', 0))))
    return keys


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('queue'); ap.add_argument('fleets', nargs='+')
    ap.add_argument('--kgfb', type=float, default=10.5); ap.add_argument('--routes', default='')
    ap.add_argument('--dmax', type=float, default=0.15); ap.add_argument('--p', type=float, default=0.02)
    ap.add_argument('--wf', default='1.5,3.0'); ap.add_argument('--beam', type=int, default=24)
    ap.add_argument('--tries', type=int, default=72); ap.add_argument('--slack', type=int, default=3)
    ap.add_argument('--stones', action='store_true'); ap.add_argument('--append', action='store_true')
    ap.add_argument('--prefix', default='P'); ap.add_argument('--max-new', type=int, default=4)
    a = ap.parse_args()
    import run_ialns as RI
    RI.eph()
    qf = ROOT / a.queue
    jobs = json.load(open(qf)) if (a.append and qf.exists()) else []
    have = existing_keys() | {(j['replaced'], tuple(sorted(j['pool'])), float(j['w_fuel']), float(j['p'])) for j in jobs}
    wfs = [float(x) for x in a.wf.split(',') if x]
    want = [x for x in a.routes.split(',') if x]
    for spec in a.fleets:
        files = FA.fleet_files(spec)
        A = FA.analyse(files, dmax=max(a.dmax, 0.3), price=False, nproc=4)
        print(f'fleet {spec}: {A["craft"]} craft, covered {A["covered"]}, sum J_i {A["sumJi"]:.4f}, misses {A["misses"]}')
        cov_all = set().union(*[set(r['targets']) for r in A['routes'].values()])
        for k, r in sorted(A['routes'].items(), key=lambda kv: -kv[1]['kg_fb']):
            if want and k not in want:
                continue
            if not want and r['kg_fb'] <= a.kgfb:
                continue
            near = sorted({X for X, lst in A['near'].items() for x in lst if x['host'] == k and x['dist'] <= a.dmax})
            near = sorted(near, key=lambda X: min(x['dist'] for x in A['near'][X] if x['host'] == k))[:a.max_new]
            pool = sorted(set(r['targets']) | set(near))
            stones = sorted(cov_all - set(pool)) if a.stones else None
            for wf in wfs:
                key = (r['f'], tuple(pool), wf, a.p if near else 0.0)
                if key in have:
                    print(f'   {k}: pool {len(pool)} wf {wf}: already defined'); continue
                have.add(key)
                h = hashlib.md5(f'{r["f"]}|{pool}|{wf}|{a.p}|{bool(stones)}'.encode()).hexdigest()[:5]
                tag = f'{a.prefix}{k}_{h}'
                J = dict(tag=tag, out=f'results/s15b/cols/{tag}', pool=pool, S=near, p=a.p if near else 0.0, steps=1,
                         beam=a.beam, tries=a.tries, nproc=1, min_save=max(3, r['n'] - a.slack), w_fuel=wf,
                         fleet_pre=[f for kk, f in files.items() if kk != k], replaced=r['f'], src_fleet=spec,
                         src_route=dict(n=r['n'], tank=r['tank'], kg_fb=r['kg_fb']))
                if stones:
                    J.update(stones=stones, stone_prize=0.01)
                jobs.append(J)
                print(f'   {tag}: re-plan {k} ({r["n"]} fb @ {r["tank"]:.0f} kg, {r["kg_fb"]:.1f} kg/fb) pool {len(pool)} '
                      f'+misses {near} wf {wf}' + (' +stones' if stones else ''))
    qf.parent.mkdir(parents=True, exist_ok=True)
    json.dump(jobs, open(qf, 'w'), indent=1)
    print(f'{len(jobs)} jobs -> {qf}')


if __name__ == '__main__':
    main()
