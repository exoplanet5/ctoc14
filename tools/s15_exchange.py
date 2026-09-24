"""Stage 15 [E]: single-route EXCHANGE re-plans of a 9-craft fleet with the lintwin twin beam.

The twin beam re-flies a route's own target set at 82-100 % (stage 14, t10d pools), the production beam at 5-84 %.
So for a fleet with misses M: re-plan ONE route r over pool = targets(r) + M with a premium on M (and on r's
isolated targets, which must stay), every other route fixed ('fleet_pre').  A re-plan that keeps r's targets and
adds some of M raises coverage without touching the other craft.

usage: s15_exchange.py OUT_JOBS.json FLEET_DIR [--routes s1,s2,t1] [--p 0.05] [--w-fuel 1.5] [--tag t04] [--name X]"""
import os, sys, json, pathlib, argparse
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
sys.path.insert(0, str(ROOT / 'tools'))
os.chdir(ROOT)
import numpy as np

ISO = [64, 119, 137, 138, 157, 158, 168, 216]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('fleet')
    ap.add_argument('--routes', default=''); ap.add_argument('--p', type=float, default=0.05)
    ap.add_argument('--w-fuel', type=float, default=1.5); ap.add_argument('--tag', default='t04')
    ap.add_argument('--name', default=None); ap.add_argument('--beam', type=int, default=20)
    ap.add_argument('--tries', type=int, default=60); ap.add_argument('--max-pool', type=int, default=40)
    ap.add_argument('--append', action='store_true')
    a = ap.parse_args()
    from greedy_cover import ALL
    d = pathlib.Path(a.fleet); nm = a.name or d.name
    files = {f.stem[6:]: str(f.relative_to(ROOT)) if f.is_absolute() else str(f) for f in sorted((ROOT / d).glob('route_*.npz'))}
    tg = {k: set(int(x) for x in np.load(ROOT / f)['asts']) for k, f in files.items()}
    cov = set().union(*tg.values()); M = sorted(set(ALL) - cov)
    names = [x for x in a.routes.split(',') if x] or sorted(files, key=lambda k: len(tg[k]))
    jobs = json.load(open(a.out)) if (a.append and pathlib.Path(a.out).exists()) else []
    print(f'fleet {nm}: {len(files)} routes, covered {len(cov)}, misses {M}')
    for k in names:
        pool = sorted(tg[k] | set(M))
        if len(pool) > a.max_pool:
            print(f'  {k}: pool {len(pool)} > {a.max_pool}: skipped'); continue
        S = sorted(set(M) | (tg[k] & set(ISO)))
        jt = f'E{nm}_{k}'
        jobs.append(dict(tag=jt, out=f'results/s15/cols/{a.tag}/{jt}', pool=pool, S=S, p=a.p, steps=1, beam=a.beam,
                         tries=a.tries, nproc=1, min_save=max(3, len(tg[k]) - 6), w_fuel=a.w_fuel,
                         fleet_pre=[f for kk, f in files.items() if kk != k], replaced=files[k], src_fleet=str(d)))
        print(f'  {jt}: re-plan {k} ({len(tg[k])} fb) over {len(pool)} with premium on {S}')
    json.dump(jobs, open(a.out, 'w'), indent=1)
    print(f'{len(jobs)} jobs -> {a.out}')


if __name__ == '__main__':
    main()
