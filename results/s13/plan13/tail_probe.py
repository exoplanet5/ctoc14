"""Measurement (the plan's early-kill test G0.5): keep routes 1..K of a settled pass, re-plan the tail routes K+1..8 with
the TAIL PORTFOLIO (one member per process), take for each route the deepest candidate that settles (twin tank <= 1.15x
planner, miss <= 150 km, tank <= TANKMAX), then the next route.  Reports settled coverage and sum J_i of the new fleet.
usage: tail_probe.py SRC_FLEET OUT_DIR KEEP [TANKMAX]"""
import os, sys, json, time, signal, pathlib, multiprocessing as mp
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
sys.path.insert(0, str(ROOT / 'results/s13/planner_probe'))
os.chdir(ROOT)
import numpy as np
import run_ialns as RI
import ctoc14.impulsive as IM
from greedy_cover import BP, CM, DeepCollect, ALL, _timeout
from ctoc14.search import beam_search, tour_json
from ctoc14.constants import DAY, AU
IM.LIN_FALLBACK = 2.5
MEMBERS = {'T1': dict(dv=1.2, dr=0.15, npt=6, mc=150, g=800, wt=1.0),
           'T2': dict(dv=1.6, dr=0.20, npt=6, mc=150, g=800, wt=1.0),
           'T3': dict(dv=2.0, dr=0.25, npt=2, mc=60, g=800, wt=1.0),
           'T4': dict(dv=1.2, dr=0.15, npt=6, mc=150, g=1500, wt=1.0),
           'T5': dict(dv=1.6, dr=0.20, npt=6, mc=150, g=800, wt=0.5),
           'T6': dict(dv=2.0, dr=0.25, npt=6, mc=150, g=1500, wt=1.0)}

def run_member(job):
    name, avail = job; m = MEMBERS[name]
    prize = np.zeros(300)
    for t in avail: prize[t - 1] = 1.0
    P = BP(beam=100, w_fuel=CM.w_fuel(480, 1600), m_margin=40.0, dv_max=m['dv'], tofs=np.arange(15, 401, 5) * DAY,
           lin_tofs=np.arange(20, 601, 10) * DAY, lin_drmax=m['dr'] * AU, vinf_cap=4.0, tof_refine=True, max_depth=60,
           n_per_target=m['npt'], max_children=m['mc'], w_t=m['wt'])
    P.prize = prize; P.collect = DeepCollect(4); tic = time.time()
    best, beam = beam_search(RI.eph(), excluded=sorted(set(ALL) - set(avail)), m0=1600.0,
                             t_launch_grid=np.arange(0, m['g'] + 1, 20) * DAY, P=P, n_proc=1, verbose=False)
    byn = {}
    for s in list(P.collect) + list(beam):
        k = len(s.seq)
        if k not in byn or s.fuel < byn[k].fuel: byn[k] = s
    return name, [(k, float(CM.tank(byn[k].fuel, 1600.0)), tour_json(byn[k])) for k in sorted(byn, reverse=True)[:8]], time.time() - tic

def settle(job):
    name, k, ptank, tour = job; tic = time.time()
    try:
        signal.signal(signal.SIGALRM, _timeout); signal.setitimer(signal.ITIMER_REAL, 240)
        ip, miss, lag = IM.settle_tour(RI.eph(), tour, lambda ip: RI.settle(ip, 100))
    except Exception:
        ip = None; miss = float('inf')
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    ok = ip is not None and miss <= 150 and float(ip.tank()) <= 1.15 * ptank
    return name, k, ptank, (RI.ist(ip) if ok else None), (float(ip.tank()) if ip is not None else None), miss, time.time() - tic

def main():
    src, out, K = sys.argv[1], pathlib.Path(sys.argv[2]), int(sys.argv[3]); TMAX = float(sys.argv[4]) if len(sys.argv) > 4 else 1050.0
    out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt')
    src_fl = RI.IFleet(src); fl = RI.IFleet()
    names = sorted(src_fl.routes)
    for n in names[:K]: fl.routes[n] = src_fl.routes[n]
    say(f'kept {K}: {fl.summary()}; source tail: ' + ' '.join(f'{n}:{len(src_fl.routes[n]["st"]["asts"])}@{src_fl.routes[n]["tank"]:.0f}' for n in names[K:]))
    with mp.get_context('fork').Pool(6, maxtasksperchild=2) as pool:
        for k in range(K + 1, 9):
            avail = sorted(set(ALL) - set(fl.coverage()))
            res = pool.map(run_member, [(m, avail) for m in MEMBERS])
            cands = []
            for name, lst, sec in res:
                say(f'  route {k:02d} {name}: pool {len(avail)}, front ' + ' '.join(f'{kk}:{t:.0f}' for kk, t, _ in lst[:5]) + f' ({sec:.0f} s)')
                cands += [(kk, t, tour, name) for kk, t, tour in lst if t <= TMAX]
            cands.sort(key=lambda c: (-c[0], c[1]))
            took = None; tried = 0
            while cands and took is None and tried < 18:
                batch = cands[:6]; cands = cands[6:]; tried += len(batch)
                sres = pool.map(settle, [(nm, kk, t, tour) for kk, t, tour, nm in batch])
                for nm, kk, pt, st, tank, miss, sec in sres:
                    say(f'    {nm} {kk} fb planner {pt:.0f} -> twin {tank if tank is None else round(tank)} (miss {miss:.0f}) {"OK" if st else "fail"} ({sec:.0f} s)')
                ok = [r for r in sres if r[3] is not None and r[4] <= TMAX]
                if ok:
                    took = max(ok, key=lambda r: (r[1], -r[4]))
            if took is None:
                say(f'route {k:02d}: nothing settled'); break
            nm, kk, pt, st, tank, miss, sec = took
            fl.routes[f'{k:02d}'] = dict(st=st, tank=tank)
            fl.save(out / 'fleet', note=f'tail probe after route {k}')
            sj = sum(RI.cost(r['tank']) for r in fl.routes.values())
            say(f'route {k:02d}: TOOK {nm} {kk} fb @ {tank:.0f} kg; covered {len(fl.coverage())}, sum J_i {sj:.4f}')
    sj = sum(RI.cost(r['tank']) for r in fl.routes.values())
    say(f'FINAL tail probe: {fl.summary()}; covered {len(fl.coverage())}, sum J_i {sj:.4f}; source was '
        f'{len(src_fl.coverage())} @ {sum(RI.cost(r["tank"]) for r in src_fl.routes.values()):.4f}')

if __name__ == '__main__':
    main()
