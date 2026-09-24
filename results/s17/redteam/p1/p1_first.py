"""P1 frontier restricted to the first NFIRST columns of cols.jsonl (9 fleet routes + pool n>=36; the F2 set), because the
full 1382-column MILP gets no incumbent in 900 s (full2_*.json).  usage: p1_first.py OUT NFIRST job,job,...  job = thr:N:K"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, multiprocessing as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p1_milp as PM
out, nfirst, jobs = sys.argv[1], int(sys.argv[2]), [tuple(j.split(':')) for j in sys.argv[3].split(',')]
keep = set(json.loads(l)['key'] for i, l in zip(range(nfirst), open(PM.P1 / 'cols.jsonl')))


def run(j):
    thr, N, K = float(j[0]), int(j[1]), int(j[2])
    cols = [c for c in PM.load(thr) if c['key'] in keep]
    r = PM.solve(cols, N, K, 150.0, 600.0, maxcov=(K == 0)); r.update(thr=thr, cap=150.0, nfirst=nfirst)
    return r


R = []
with mp.get_context('fork').Pool(len(jobs)) as pool:
    for r in pool.imap_unordered(run, jobs):
        R.append(r); print(json.dumps({k: v for k, v in r.items() if k not in ('columns', 'arcs')}), flush=True)
        json.dump(R, open(PM.P1 / out, 'w'), indent=1)
print('all done', flush=True)
