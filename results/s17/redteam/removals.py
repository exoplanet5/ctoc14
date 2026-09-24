"""s17 redteam: single-target removal marginals (run_ialns.w_remove = remove + settle) for the s16 best fleet routes.
Resumable: appends to removals.jsonl, skips (route, ast) already done.  Cached r2/r5 rows from results/s16/iter are
merged by q_removal.py, not recomputed.  usage: removals.py r9 r8 ... [--nproc 2]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    sys.path.insert(0, _p)
os.chdir(ROOT)
import numpy as np
import run_ialns as RI
from s16_colgen import w_rem
OUT = ROOT / 'results/s17/redteam/removals.jsonl'
if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('V', nargs='+'); ap.add_argument('--nproc', type=int, default=2)
    a = ap.parse_args(); RI.eph()
    done = set()
    if OUT.exists():
        for l in open(OUT):
            r = json.loads(l); done.add((r['route'], r['ast']))
    jobs = []
    for V in a.V:
        z = np.load(ROOT / f'results/s16/best/fleet/route_{V}.npz'); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL'])
        for x in st['asts']:
            if (V, int(x)) not in done:
                jobs.append((V, st, int(x), 150.0))
    with open(OUT, 'a') as fo, mp.get_context('fork').Pool(a.nproc, maxtasksperchild=8) as pool:
        for (V, st, x, _), r in zip(jobs, pool.imap(w_rem, jobs)):
            tank0 = float(RI.ipr(st).tank())
            row = dict(route=V, ast=x, ok=bool(r.get('ok')), sec=r.get('sec'), tank0=tank0)
            if r.get('ok'):
                row.update(tank=float(r['tank']), dkg=float(tank0 - r['tank']), dJ=float(RI.cost(tank0) - RI.cost(r['tank'])))
            fo.write(json.dumps(row) + '\n'); fo.flush()
            print(row, flush=True)
