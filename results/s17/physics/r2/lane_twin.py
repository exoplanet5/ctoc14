"""R2 check (PHYSICS): realise the T1 pin lanes in the twin and compare with the carrier-model price.
Read-only on tools/ and ctoc14/.  Thresholds: THRESHOLDS_R2.txt (written before the run)."""
import os
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'): os.environ.setdefault(v,'1')
import sys, json, time, pathlib, multiprocessing as mp, signal
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT/'tools')); sys.path.insert(0, str(ROOT/'results/s17/physics'))
os.chdir(ROOT)
import numpy as np
import rare_lanes as RL            # sets CD.K_SLOPE = K_PHYS*1.5, loads RARE ids/elem
import carrier_dp as CD
import run_ialns as RI
from s15b_regrid import regrid
from exp_carrier import rv2frame
from ctoc14.constants import DAY, VE
OUT = ROOT/'results/s17/physics/r2'
L = json.load(open(ROOT/'results/s17/physics/rare_lanes.json'))['lanes']

class TO(Exception): pass
def _alarm(*a): raise TO()

def job(k):
    lane = L[k]; tL = lane['tL']*DAY; vinf = np.array(lane['vinf'])
    E = RL.E; rE, vE = E.earth_state(tL); r = np.array(rE, float); v = np.array(vE, float) + vinf
    t, lag, gap, ast, Pc = CD.events_for(dict(tL=tL, r=r, v=v), RL.ids, RL.elem, RL.EPS)
    ac = float(rv2frame(r[None], v[None])[3][0])
    val, seq = CD.dp_route(t, lag, gap, ast, tL, RL.LAM, RL.KGAP, ac=ac)
    seen = set(); seq = [i for i in seq if not (ast[i] in seen or seen.add(ast[i]))]
    tl = np.concatenate([[tL/DAY], t[seq]/DAY]); lg = np.concatenate([[0.0], lag[seq]]); S = np.diff(lg)/np.diff(tl)
    dvp15 = CD.K_SLOPE*np.abs(np.diff(S)).sum(); dvp1 = CD.K_PHYS*np.abs(np.diff(S)).sum()
    dvg = RL.KGAP*np.abs(gap[seq] - (2/3)*ac*S).sum()
    rec = dict(lane=k, n=len(seq), tg=[int(ast[i]) for i in seq], t=[round(float(t[i]/DAY)) for i in seq],
               model_dv=float(dvp15+dvg), model_dv_phase_k1=float(dvp1), model_dv_gap=float(dvg), stored_dv=lane['dv'], stored_n=lane['n_lane'])
    tic = time.time(); signal.signal(signal.SIGALRM, _alarm); signal.alarm(420)
    try:
        from ctoc14.globalopt import insert_block
        ip = CD.seed_ip(E, tL, vinf, t[seq], lag[seq])
        d = insert_block(ip, [(int(ast[i]), float(t[i])) for i in seq], log=RI.QUIET)
        miss = float(np.max(d))
        if miss < 1e4: miss = RI.settle(ip, 100)
        rec.update(miss_raw=miss, tank_raw=float(ip.tank()))
        if miss <= 150:
            g = regrid(ip, 20*DAY); m2 = RI.settle(g, 100)
            rec.update(miss=m2, tank=float(g.tank()), dv_twin=float(VE*np.log(g.tank()/601.5)))
        else:
            rec.update(miss=miss)
    except TO:
        rec.update(miss=float('inf'), note='timeout 420 s')
    except Exception as ex:
        rec.update(miss=float('inf'), note=repr(ex)[:200])
    signal.alarm(0); rec['sec'] = round(time.time()-tic)
    print(json.dumps(rec), flush=True)
    return rec

if __name__ == '__main__':
    ks = sorted(L, key=int)
    with mp.get_context('fork').Pool(2) as pl:
        R = pl.map(job, ks, chunksize=1)
    json.dump(R, open(OUT/'lane_twin.json','w'), indent=1)
    ok = [x for x in R if x.get('miss', 1e9) <= 150 and 'dv_twin' in x]
    print(f'settled {len(ok)}/{len(R)}')
    if ok:
        ratio = np.median([x['dv_twin']/x['model_dv'] for x in ok]); per = np.median([x['dv_twin']/x['n'] for x in ok])
        print(f'median twin/model {ratio:.2f}; median twin dv per pin {per:.3f} km/s; lanes', [(x['lane'], x['n'], round(x['model_dv'],1), round(x['dv_twin'],1), round(x['tank'])) for x in ok])
