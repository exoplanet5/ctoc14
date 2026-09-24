"""Stage 14 track A: assemble the best re-flown twin route per t10d route into one route-set directory.

For every route, among the free runs whose variant matches --variants (default lt, lt60: the linearised-twin gate),
take the settled best.npz with the most pool targets re-flown (ties: lighter twin tank).  The t10d pools are disjoint,
so the result is a 10-craft route set; it is written as a run_ialns IFleet (route_<name>.npz + fleet.json) for
inspection only -- NOT an export or a submission.
Usage: s14_bestset.py OUT_DIR [--variants lt,lt60] [--root results/s14/A]
"""
import sys, json, pathlib, argparse
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--variants', default='lt,lt60')
    ap.add_argument('--root', default=str(ROOT / 'results/s14/A'))
    a = ap.parse_args(); R = pathlib.Path(a.root); vs = a.variants.split(',')
    best = {}
    for f in sorted(R.glob('free/*/result.json')):
        v = f.parent.name.rsplit('_', 1)[0]
        if v not in vs:
            continue
        d = json.load(open(f)); r = d['route']
        key = (d['reflown'], -d['tank'])
        if r not in best or key > best[r][0]:
            best[r] = (key, f.parent, d)
    F = RI.IFleet()
    rows = {}
    for r, (key, p, d) in sorted(best.items()):
        z = np.load(p / 'best.npz'); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL'])
        ip = RI.ipr(st); Yf, _ = ip.integrate(); dm, _ = ip.misses(Yf)
        miss = float(np.linalg.norm(dm, axis=1).max())
        F.routes[r] = dict(st=st, tank=float(ip.tank()))
        rows[r] = dict(run=p.name, reflown=d['reflown'], pool=d['pool_size'], tank=round(float(ip.tank()), 2),
                       src_tank=d['src_tank'], miss_km=round(miss, 1), J_i=round(RI.cost(float(ip.tank())), 5))
        print(f'{r}: {p.name:8s} {d["reflown"]:2d}/{d["pool_size"]:2d}  tank {ip.tank():7.1f} (t10d {d["src_tank"]:7.1f})  '
              f'max miss {miss:6.1f} km  J_i {RI.cost(float(ip.tank())):.4f}')
    F.save(a.out, note='s14 track A: best linearised-twin re-fly per t10d route (inspection only, not an export)')
    cov = len(F.coverage()); sj = sum(RI.cost(r['tank']) for r in F.routes.values())
    print(f'route set: {len(F.routes)} routes, covered {cov}, sum J_i {sj:.4f}, raw J (with misses) {sj + 300 - cov:.4f}')
    json.dump(dict(routes=rows, covered=cov, sum_Ji=sj, rawJ=sj + 300 - cov), open(pathlib.Path(a.out) / 'bestset.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
