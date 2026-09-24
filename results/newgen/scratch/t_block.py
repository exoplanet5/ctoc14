"""Block insertion prototype: a contiguous time-block of victim V's targets inserted into a host, all aim points moving
together (globalopt.insert_block), with and without ejecting the host's own flybys inside the block window."""
import os
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'): os.environ.setdefault(v, '1')
import sys, time, numpy as np, multiprocessing as mp
sys.path.insert(0, '.'); sys.path.insert(0, 'tools')
from run_ialns import IFleet, ipr, ist, settle, cost, eph, trajectory
from ctoc14.globalopt import insert_block, remove_flybys
from ctoc14.constants import DAY, AU
YR = 365.25 * DAY

def block_job(job):
    host, st, items, eject, tag = job
    ip = ipr(st); tank0 = ip.tank(); tic = time.time()
    ej = []
    if eject:
        t0 = min(t for _, t in items) - 30 * DAY; t1 = max(t for _, t in items) + 30 * DAY
        ej = [int(a) for a, t in zip(ip.asts, ip.tf) if t0 <= t <= t1]
        if ej: remove_flybys(ip, ej)
    try:
        d = insert_block(ip, items, tol_km=100, ds0=0.2, iters_stage=6)
        ok = d.max() <= 1e4
        miss = settle(ip, 100) if ok else float(d.max())
        ok = ok and miss <= 150
    except Exception as e:
        return dict(tag=tag, host=host, ok=False, err=repr(e)[:80], sec=round(time.time() - tic))
    dJ = cost(ip.tank()) - cost(tank0)
    return dict(tag=tag, host=host, n=len(items), eject=len(ej), ok=bool(ok), miss=round(float(miss)), dJ=round(dJ, 4),
                per=round(dJ / len(items), 4), tank=f'{tank0:.0f}->{ip.tank():.0f}', sec=round(time.time() - tic), ejected=ej)

if __name__ == '__main__':
    V = sys.argv[1] if len(sys.argv) > 1 else '09'
    F = IFleet('results/newgen/ifleet10')
    vst = F.routes[V]['st']; o = np.argsort(vst['tf'])
    vt = vst['tf'][o]; va = [int(a) for a in np.array(vst['asts'])[o]]
    print(f'victim {V}: cost {cost(F.routes[V]["tank"]):.4f}, {len(va)} targets from {vt[0]/YR:.1f} to {vt[-1]/YR:.1f} yr')
    E = eph()
    # blocks of 5 consecutive targets
    B = 5; blocks = [list(zip(va[i:i + B], vt[i:i + B])) for i in range(0, len(va), B)]
    jobs = []
    for bi, items in enumerate(blocks):
        if len(items) < 3: continue
        tmid = np.mean([t for _, t in items])
        # rank hosts by mean distance to the block asteroids at their own times
        sc = []
        for n, r in F.routes.items():
            if n == V: continue
            ip = ipr(r['st']); tt, rr = trajectory(ip)
            dd = []
            for a, t in items:
                ra, _ = E.ast_states_at(np.array([a - 1]), np.array([t]))
                i = int(np.clip(np.searchsorted(tt, t), 0, len(tt) - 1))
                dd.append(np.linalg.norm(rr[i] - ra[0]) / AU)
            sc.append((float(np.mean(dd)), n))
        sc.sort()
        print(f'block {bi} ({len(items)} targets {[a for a,_ in items]} at {vt[bi*B]/YR:.1f}-{tmid/YR:.1f} yr): '
              f'closest hosts ' + ', '.join(f'{n} {d:.3f} AU' for d, n in sc[:4]))
        for d, n in sc[:3]:
            jobs.append((n, F.routes[n]['st'], items, False, f'b{bi}'))
            jobs.append((n, F.routes[n]['st'], items, True, f'b{bi}E'))
    print(f'{len(jobs)} block trials', flush=True)
    with mp.get_context('fork').Pool(8, maxtasksperchild=10) as pool:
        for r in pool.imap_unordered(block_job, jobs):
            print(r, flush=True)
