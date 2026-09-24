"""Why do insertions of route 09's hard targets fail? Trials on the 9 other routes with per-trial diagnostics."""
import os
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'): os.environ.setdefault(v, '1')
import sys, time, numpy as np, multiprocessing as mp
sys.path.insert(0, '.'); sys.path.insert(0, 'tools')
from run_ialns import IFleet, ipr, settle, w_cands, cost
from ctoc14.globalopt import insert_homotopy
from ctoc14.constants import DAY

def trial(job):
    name, st, ast, t, dist = job
    ip = ipr(st); tank0 = ip.tank(); tic = time.time()
    info = dict(host=name, ast=ast, t=round(t / DAY), dist=round(dist, 3)); lines = []
    try:
        d = insert_homotopy(ip, ast, t, stages=10, iters_stage=6, tol_km=100, log=lines.append, verbose=True)
        info['hom_miss'] = round(float(d.max()))
        info['worst_ast'] = ip.asts[int(np.argmax(d))]
        info['n_rej'] = sum('REJ' in l for l in lines)
        info['stages'] = [l.split(':', 1)[1].strip() for l in lines if 'homotopy stage' in l]
        if d.max() <= 1e4:
            info['settle_miss'] = round(settle(ip, 100))
            info['dJ'] = round(cost(ip.tank()) - cost(tank0), 4)
            info['tank'] = f'{tank0:.0f}->{ip.tank():.0f}'
            info['cap_ratio'] = round(float(np.linalg.norm(ip.Ts, axis=1).max() / ip.tcap), 2)
    except Exception as e:
        info['error'] = repr(e)
    info['sec'] = round(time.time() - tic)
    return info

if __name__ == '__main__':
    F = IFleet('results/newgen/ifleet10'); del F.routes['09']
    targets = [8, 39, 168, 251]
    jobs = []
    for n, r in F.routes.items():
        for c in w_cands((n, r['st'], targets, 0.6)):
            jobs.append((n, r['st'], c['ast'], c['t'], c['dist']))
    jobs.sort(key=lambda j: (j[2], j[4]))
    keep = []
    for j in jobs:
        if sum(1 for x in keep if x[2] == j[2]) < 8:
            keep.append(j)
    print(f'{len(keep)} trials of {len(jobs)} candidates', flush=True)
    with mp.get_context('fork').Pool(8) as pool:
        for info in pool.imap_unordered(trial, keep):
            st = info.pop('stages', [])
            print(info, flush=True)
            if st:
                print('   ' + ' | '.join(st[:: max(1, len(st) // 5)]), flush=True)
