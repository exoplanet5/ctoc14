"""Read-only probe: ONE sequential covering pass (8 routes, strict exclusion of covered targets) with the two levers
measured tonight: a small premium delta on a hardness set (gc16's 67 leftovers) and relaxed legs for the late routes.
usage: probe_pass.py OUT_DIR DELTA RELAX_FROM DVMAX_RELAX DRMAX_RELAX [BEAM] [TANKCAP]
Writes OUT_DIR/pass.json and OUT_DIR/route_NN.npz (twin-settled). Planner tank cap per route = TANKCAP (default 990)."""
import os, sys, json, time, signal, pathlib
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
_CWD0 = os.getcwd(); os.chdir(ROOT)
import numpy as np
import run_ialns as RI
import ctoc14.impulsive as IM
from greedy_cover import BP, CM, DeepCollect, ALL, _timeout
from ctoc14.search import beam_search, tour_json
from ctoc14.constants import DAY, AU

IM.LIN_FALLBACK = 2.5
out = pathlib.Path(os.path.join(_CWD0, sys.argv[1])); out.mkdir(parents=True, exist_ok=True)
delta = float(sys.argv[2]); relax_from = int(sys.argv[3]); dv_r = float(sys.argv[4]); dr_r = float(sys.argv[5])
beamw = int(sys.argv[6]) if len(sys.argv) > 6 else 100
tankcap = float(sys.argv[7]) if len(sys.argv) > 7 else 990.0
gc = RI.IFleet('results/newgen/gc16'); cov = set()
for n, r in gc.routes.items(): cov |= set(int(x) for x in r['st']['asts'])
HARD = set(json.load(open(os.environ['HARDFILE']))[os.environ.get('HARDKEY','U1')]) if os.environ.get('HARDFILE') else set(ALL) - cov
NEG = set(json.load(open(os.environ['HARDFILE']))[os.environ['NEGKEY']]) if os.environ.get('NEGKEY') else set()
covered = set(); log = []; routes = []
grid = np.arange(0, 801, 20) * DAY
for k in range(1, 9):
    relaxed = k >= relax_from
    dv = dv_r if relaxed else 1.2; dr = dr_r if relaxed else 0.15
    avail = [t for t in ALL if t not in covered]
    prize = np.zeros(300)
    for t in avail: prize[t - 1] = 1.0 + (delta if t in HARD else 0.0) + (float(os.environ.get('DNEG','0')) if t in NEG else 0.0)
    P = BP(beam=beamw, w_fuel=CM.w_fuel(480, 1600), m_margin=40.0, dv_max=dv,
           tofs=np.arange(15, 401, 5) * DAY, lin_tofs=np.arange(20, 601, 10) * DAY, lin_drmax=dr * AU,
           vinf_cap=4.0, tof_refine=True, max_depth=60, n_per_target=int(os.environ.get('NPT','2')), max_children=int(os.environ.get('MC','60')))
    P.prize = prize; P.collect = DeepCollect(6)
    if os.environ.get('DVT'):
        dvt = np.full(300, dv)
        _key = os.environ.get('DVT_KEY', 'HARD')
        _set = HARD if _key == 'HARD' else (HARD | set(json.load(open(os.environ['HARDFILE']))['U1'])) if _key == 'UNION' else set(json.load(open(os.environ['HARDFILE']))[_key])
        for t in _set: dvt[t - 1] = max(dv, float(os.environ['DVT']))
        P.dv_max_t = dvt
    tic = time.time()
    best, beam = beam_search(RI.eph(), excluded=sorted(covered), m0=1600.0, t_launch_grid=grid, P=P, n_proc=1, verbose=False)
    states = [s for s in list(P.collect) + list(beam) if CM.tank(s.fuel, 1600.0) <= tankcap]
    if not states:
        log.append(dict(k=k, n=0)); break
    s = max(states, key=lambda s: (len(s.seq), -s.fuel))
    tg = [x[0] for x in s.seq]
    ip = None; miss = float('inf'); t1 = time.time()
    try:
        signal.signal(signal.SIGALRM, _timeout); signal.setitimer(signal.ITIMER_REAL, 240)
        ip, miss, lag = IM.settle_tour(RI.eph(), tour_json(s), lambda ip: RI.settle(ip, 100))
    except Exception:
        ip = None
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    ok = ip is not None and miss <= 150
    twin = float(ip.tank()) if ok else None
    if ok: np.savez(out / f'route_{k:02d}.npz', **RI.ist(ip))
    covered |= set(tg)
    rec = dict(k=k, relaxed=relaxed, n=len(tg), hard=len(set(tg) & HARD), planner=round(float(CM.tank(s.fuel, 1600.0))),
               twin=None if twin is None else round(twin), covered=len(covered), beam_s=round(t1 - tic), settle_s=round(time.time() - t1),
               targets=sorted(tg))
    log.append(rec); routes.append(rec)
    print(json.dumps({kk: v for kk, v in rec.items() if kk != 'targets'}), flush=True)
    json.dump(dict(delta=delta, relax_from=relax_from, dv_relax=dv_r, dr_relax=dr_r, beam=beamw, tankcap=tankcap, routes=log,
                   covered=len(covered), hard_left=len(HARD - covered)), open(out / 'pass.json', 'w'), indent=1)
Ji = [RI.cost(r['twin'] if r['twin'] else r['planner'] * 1.02) for r in routes]
print('PASS', json.dumps(dict(covered=len(covered), hard_left=len(HARD - covered), sumJi=round(sum(Ji), 3),
                               depths=[r['n'] for r in routes], tanks=[r['twin'] or r['planner'] for r in routes])), flush=True)
