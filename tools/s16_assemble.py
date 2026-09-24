"""Stage 16: assemble a fleet dir (RI.IFleet, routes r1..r9 by depth) from a fleet spec (FRONTIER_JSON:KEY | dir | list).
usage: s16_assemble.py SPEC OUTDIR [--note TEXT]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, pathlib, argparse
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)
import run_ialns as RI
import s15b_fleet as FA

ap = argparse.ArgumentParser(); ap.add_argument('spec'); ap.add_argument('out'); ap.add_argument('--note', default='')
a = ap.parse_args()
files = FA.fleet_files(a.spec)
sts = sorted((FA.load_st(f) for f in files.values()), key=lambda st: (-len(st['asts']), float(RI.ipr(st).tank())))
F = RI.IFleet()
for i, st in enumerate(sts):
    F.routes[f'r{i + 1}'] = dict(st=st, tank=float(RI.ipr(st).tank()))
F.save(ROOT / a.out, note=a.note or f'assembled from {a.spec}')
json.dump({f'r{i + 1}': f for i, f in enumerate(sorted(files.values(), key=lambda f: (-len(FA.load_st(f)['asts']), float(RI.ipr(FA.load_st(f)).tank()))))},
          open(ROOT / a.out / 'sources.json', 'w'), indent=1)
print(F.summary(), 'sum J_i', round(sum(RI.cost(r['tank']) for r in F.routes.values()), 5))
