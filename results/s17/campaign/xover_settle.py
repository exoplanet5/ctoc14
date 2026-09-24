"""Stage 17 CAMPAIGN test X2: build crossover children as impulsive twins and settle them honestly.

Child twin = A's launch + A's impulses up to t1 + junction impulses (Lambert) at t1 and t2 + B's impulses after t2,
flybys A[< t1] + B[> t2].  Before settle the child is exactly feasible (same states as A before t1 and B after t2),
so estimate == raw twin tank.  Then s15b_regrid.regrid (20 d bins from launch) + RI.settle (tcap enforced) gives the
honest tank.  usage: xover_settle.py [milp_key] [max_children] [nproc]
Outputs results/s17/campaign/settle.json and children/route_<i>.npz
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pickle, pathlib, multiprocessing as mp
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)
import numpy as np
import warnings; warnings.filterwarnings('ignore')
from ctoc14.kepler import propagate_batch
from ctoc14.lambert import lambert
from ctoc14.constants import DAY
from ctoc14.impulsive import ImpulsiveProblem
import run_ialns as RI
import s15b_regrid as RG

OUT = ROOT / 'results/s17/campaign'
GSTEP = 10 * DAY


def load(f):
    z = np.load(f); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL']); return st


def build(D, S, a, b, k, L):
    A = load(D['reps'][a]['f']); B = load(D['reps'][b]['f'])
    t1 = k * GSTEP; t2 = (k + L) * GSTEP
    r1, vA = S[a, k, :3], S[a, k, 3:]; r2, vB = S[b, k + L, :3], S[b, k + L, 3:]
    v1, v2 = lambert(r1[None], r2[None], np.array([t2 - t1]))
    v1 = v1[0]; v2 = v2[0]
    ma = A['ts'] <= t1; mb = B['ts'] > t2
    ts = np.concatenate([A['ts'][ma], [t1 + 1.0, t2 - 1.0], B['ts'][mb]])
    Ts = np.vstack([A['Ts'][ma], (v1 - vA)[None], (vB - v2)[None], B['Ts'][mb]])
    fa = A['tf'] < t1; fb = B['tf'] > t2
    tf = np.concatenate([A['tf'][fa], B['tf'][fb]]); asts = np.concatenate([A['asts'][fa], B['asts'][fb]])
    # drop duplicate targets (keep the first occurrence)
    seen = set(); keep = []
    for i, x in enumerate(asts):
        if int(x) not in seen:
            seen.add(int(x)); keep.append(i)
    tf = tf[keep]; asts = asts[keep]
    # the junction impulse at t1+1 s acts on A's state at t1 propagated 1 s: negligible
    return ImpulsiveProblem(RI.eph(), A['tL'], A['vinf'], ts, Ts, tf, [int(x) for x in asts])


def job(args):
    i, a, b, k, L, est = args
    t0 = time.time()
    D = pickle.load(open(OUT / 'reps.pkl', 'rb')); S = np.load(OUT / 'states.npz')['S']
    ip = build(D, S, a, b, k, L)
    Yf, _ = ip.integrate(); dm, _ = ip.misses(Yf); miss0 = float(np.linalg.norm(dm, axis=1).max())
    tank0 = float(ip.tank())
    g = RG.regrid(ip)
    over0 = float((np.linalg.norm(g.Ts, axis=1) / g.tcap).max())
    try:
        miss = float(RI.settle(g, 100))
    except Exception as e:
        return dict(i=i, err=type(e).__name__, tank_est=est, tank_raw=tank0, miss_raw=miss0)
    over1 = float((np.linalg.norm(g.Ts, axis=1) / g.tcap).max())
    np.savez(OUT / f'children/route_c{i}.npz', **RI.ist(g))
    return dict(i=i, a=a, b=b, k=k, L=L, n=len(g.asts), tank_est=round(est, 2), tank_raw=round(tank0, 2),
                miss_raw_km=round(miss0, 1), over_cap_before=round(over0, 2), tank_settled=round(float(g.tank()), 2),
                miss_km=round(miss, 1), over_cap_after=round(over1, 3), sec=round(time.time() - t0))


if __name__ == '__main__':
    key = sys.argv[1] if len(sys.argv) > 1 else 'with_children_N8'
    mx = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    npr = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    (OUT / 'children').mkdir(exist_ok=True)
    M = json.load(open(OUT / 'milp.json'))[key]
    jobs = [(i, p['a'], p['b'], p['k'], p['L'], p['tank']) for i, p in enumerate(M['picks']) if p['kind'] == 'child'][:mx]
    with mp.get_context('fork').Pool(npr) as pool:
        res = pool.map(job, jobs, chunksize=1)
    json.dump(dict(key=key, results=res), open(OUT / f'settle_{key}.json', 'w'), indent=1)
    for r in res:
        print(r)
