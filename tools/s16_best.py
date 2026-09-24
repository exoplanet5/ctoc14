"""Stage 16: keep results/s16/best/fleet = the best CLOSED (298 covered, 9 craft) HONEST fleet seen so far.

Runs the honest check of tools/s16_pipe.honest_check (0 irregular nodes, max 20 d window/cap <= 1.02, twin miss
<= 150 km, coverage) on the candidate fleet; if it passes and its twin sum J_i is below the stored best, the fleet is
copied to results/s16/best/fleet and results/s16/best/best.json records source, sum J_i and the per-route table.

usage: s16_best.py FLEET_DIR [--note TEXT] [--nproc 3]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, shutil, pathlib, argparse
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)
import run_ialns as RI
from s16_pipe import honest_check

BEST = ROOT / 'results/s16/best'


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('fleet'); ap.add_argument('--note', default='')
    ap.add_argument('--nproc', type=int, default=3)
    a = ap.parse_args()
    BEST.mkdir(parents=True, exist_ok=True); say = RI.logger(BEST / 'log.txt'); RI.eph()
    rep = honest_check(ROOT / a.fleet, a.nproc, say)
    old = json.load(open(BEST / 'best.json')) if (BEST / 'best.json').exists() else None
    ok = rep['honest'] and rep['covered'] == 298 and rep['craft'] == 9
    if ok and (old is None or rep['sumJi'] < old['sumJi'] - 1e-6):
        tmp = BEST / 'fleet.tmp'
        if tmp.exists():
            shutil.rmtree(tmp)
        shutil.copytree(ROOT / a.fleet, tmp)
        if (BEST / 'fleet').exists():
            shutil.rmtree(BEST / 'fleet')
        tmp.rename(BEST / 'fleet')
        rep.update(source=a.fleet, note=a.note, when=time.strftime('%m-%d %H:%M'))
        json.dump(rep, open(BEST / 'best.json', 'w'), indent=1)
        say(f'NEW BEST {rep["sumJi"]:.5f} from {a.fleet} ({a.note})')
    else:
        say(f'not stored: honest {rep["honest"]}, covered {rep["covered"]}, sum J_i {rep["sumJi"]:.5f} '
            f'(best {old["sumJi"] if old else None})')


if __name__ == '__main__':
    main()
