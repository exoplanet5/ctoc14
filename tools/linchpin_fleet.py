"""Stage-13 linchpin, step 4: assemble the re-planned fleet and measure it.
For every survivor of a pools file, take the best settled candidate over the given run directories
(<dir>/<route>_*/result.json + cand_k*_a*.npz; best = most pool targets flown, then lightest tank) and write an IFleet
dir (route_<name>.npz + fleet.json).  Prints the per-route table (pool, flown, tank, kg/fb vs the source route), the
fleet totals (covered = union of flown targets, orphans = pool targets nobody flew, sum J_i, J = sum J_i + misses) and
where the orphans come from (own vs extra, hard vs easy, dissolved route of origin, orbit elements, source epochs).
Usage: linchpin_fleet.py pools.json out_fleet_dir rundir [rundir ...] [--src results/s12/cover_clean/fleet]
       [--include-soft]
"""
import sys, json, pathlib, argparse
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from ctoc14.constants import DAY
from ctoc14.search import UNREACHABLE


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('pools'); ap.add_argument('out'); ap.add_argument('dirs', nargs='+')
    ap.add_argument('--src', default=str(ROOT / 'results/s12/cover_clean/fleet'))
    ap.add_argument('--include-soft', action='store_true')
    ap.add_argument('--lam', type=float, default=0.15, help='choose per route the settled candidate minimising J_i - lam * k '
                    '(lam = J value of one more pool target; 0.15 ~ the forced-insertion premium)')
    ap.add_argument('--hard', default=str(ROOT / 'results/n8/prizes_H.json'))
    a = ap.parse_args()
    PJ = json.load(open(a.pools)); pools = {n: set(p) for n, p in PJ['pools'].items()}
    V = PJ.get('dissolved', [])
    src = RI.IFleet(a.src)
    H = np.asarray(json.load(open(a.hard)), float); hard = set(i + 1 for i in range(300) if H[i] >= 1.5)
    ALL = [t for t in range(1, 301) if t not in UNREACHABLE]
    origin = {int(x): n for n, r in src.routes.items() for x in r['st']['asts']}
    src_t = {int(x): float(t) for r in src.routes.values() for x, t in zip(r['st']['asts'], r['st']['tf'])}
    best = {}
    for d in a.dirs:
        for f in sorted(pathlib.Path(d).glob('*/result.json')):
            r = json.load(open(f))
            if r['name'] not in pools or set(r['pool']) != pools[r['name']]:
                continue
            if r.get('soft_prize', -1) >= 0 and not a.include_soft:
                continue
            for s in r['settled']:
                if not s.get('ok'):
                    continue
                key = -(s['J_i'] - a.lam * s['k'])
                if r['name'] not in best or key > best[r['name']][0]:
                    best[r['name']] = (key, s, f.parent / f'cand_k{s["k"]}_a{s["alt"]}.npz', f.parent.name)
    fleet = RI.IFleet()
    print(f'{"route":5s} {"pool":>4s} {"own+ext":>8s} | {"flown":>5s} {"(own/ext)":>9s} {"n":>3s} {"tank":>7s} {"kg/fb":>6s} {"J_i":>6s} '
          f'| {"src fb":>6s} {"src tank":>8s} {"kg/fb":>6s} | run')
    for n in sorted(pools):
        own = set(int(x) for x in src.routes[n]['st']['asts'])
        if n not in best:
            print(f'{n:5s} {len(pools[n]):4d} | nothing settled'); continue
        _, s, npz, run = best[n]
        z = np.load(npz); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL'])
        fleet.routes[n] = dict(st=st, tank=float(RI.ipr(st).tank()))
        t0 = src.routes[n]['tank']; n0 = len(own)
        print(f'{n:5s} {len(pools[n]):4d} {len(own & pools[n]):3d}+{len(pools[n] - own):<3d} | {s["k"]:5d} {s["own_flown"]:4d}/{s["extra_flown"]:<4d} '
              f'{s["n"]:3d} {s["tank"]:7.1f} {s["kg_per_fb"]:6.2f} {s["J_i"]:6.3f} | {n0:6d} {t0:8.1f} {(t0 - 600) / n0:6.2f} | {run}')
    fleet.save(a.out, note=f'linchpin re-plan fleet for {a.pools}')
    cov = set(fleet.coverage())
    target = set().union(*pools.values())
    orph = sorted(target - cov)
    sumJ = sum(RI.cost(r['tank']) for r in fleet.routes.values())
    flown = sum(len(r['st']['asts']) for r in fleet.routes.values())
    prop = sum(r['tank'] - 600 for r in fleet.routes.values())
    print(f'\nFLEET {len(fleet.routes)} craft: covered {len(cov)} of {len(target)} pool targets, orphans {len(orph)}, '
          f'flybys {flown}, propellant {prop:.0f} kg = {prop / max(flown, 1):.2f} kg/flyby, sum J_i {sumJ:.4f}, '
          f'J = sum J_i + {len(orph)} orphans + 2 = {sumJ + len(orph) + 2:.3f}')
    if orph:
        E = RI.eph()
        el = E.elem if hasattr(E, 'elem') else None
        by_src = {}
        for x in orph:
            by_src.setdefault(origin.get(x, '?'), []).append(x)
        home = {x: n for n in pools for x in pools[n]}
        by_home = {}
        for x in orph:
            by_home.setdefault(home[x], []).append(x)
        print(f'orphans: {len(orph)}; hard {len(set(orph) & hard)} ({len(hard & target)} hard in the pools); from a '
              f'DISSOLVED route {sum(1 for x in orph if origin.get(x) in V)} of {sum(1 for x in target if origin.get(x) in V)} '
              f'extra targets; own targets of a survivor {sum(1 for x in orph if origin.get(x) not in V)}')
        print('  by assigned survivor: ' + ', '.join(f'{n}:{len(v)}/{len(pools[n])}' for n, v in sorted(by_home.items())))
        print('  by source route     : ' + ', '.join(f'{n}:{len(v)}' for n, v in sorted(by_src.items())))
        yr = np.array([src_t[x] / DAY / 365.25 for x in orph]); yc = np.array([src_t[x] / DAY / 365.25 for x in cov if x in src_t])
        print(f'  source flyby epoch (yr): orphans median {np.median(yr):.1f} (q25 {np.percentile(yr, 25):.1f}, q75 '
              f'{np.percentile(yr, 75):.1f}); covered median {np.median(yc):.1f}; orphans in last 5 yr '
              f'{(yr > 10).mean():.0%} vs covered {(yc > 10).mean():.0%}')
        if el is not None:
            ai = np.array([[el[x - 1][0], el[x - 1][1], el[x - 1][2]] for x in orph])
            ac = np.array([[el[x - 1][0], el[x - 1][1], el[x - 1][2]] for x in cov])
            print(f'  elements median a/e/i: orphans {np.median(ai[:, 0]):.2f} AU / {np.median(ai[:, 1]):.2f} / '
                  f'{np.median(ai[:, 2]):.1f} deg; covered {np.median(ac[:, 0]):.2f} / {np.median(ac[:, 1]):.2f} / '
                  f'{np.median(ac[:, 2]):.1f}')
        print(f'  orphan ids: {orph}')
    json.dump(dict(pools=a.pools, covered=len(cov), orphans=orph, sum_Ji=sumJ, J=sumJ + len(orph) + 2, flybys=flown,
                   propellant=prop, routes={n: dict(run=best[n][3], k=best[n][1]['k'], tank=best[n][1]['tank'],
                                                   own_flown=best[n][1]['own_flown'], extra_flown=best[n][1]['extra_flown'],
                                                   pool=len(pools[n])) for n in best}),
              open(pathlib.Path(a.out) / 'linchpin_fleet.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
