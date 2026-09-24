"""Stage 14 track A diagnostic: can the PLANNER's own state chain a source route's legs?

Along the source route's own order and epochs, the planner state is carried exactly as ctoc14.search carries it: at
every step the next source target is looked up among search.expand's children (production limits dv 1.2 / drmax
0.15, then relaxed dv 2.0 / drmax 0.25); the child nearest the source epoch is taken and ITS arrival velocity (Lambert
or LinLeg) becomes the state velocity.  If neither limit admits the leg, the state is re-seeded with the single-rev
Lambert arrival velocity at the source epoch (the chain goes on).  Mass follows the planner (m0 1600).
Compare with s14_twinbeam.py guided (myopic prefix twin state) and linchpin_legcheck (full-route twin state).
Usage: s14_chaincheck.py SRC_FLEET OUT.json [--routes 00,01,...]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from greedy_cover import ALL
from ctoc14.search import State, expand, mask_of, UNREACHABLE
from ctoc14.lambert import lambert
from ctoc14.constants import DAY, VINF_MAX
from s14_twinbeam import params, M0P


def chain(E, st):
    o = np.argsort(st['tf']); asts = [int(st['asts'][i]) for i in o]; tfs = [float(st['tf'][i]) for i in o]
    pool = sorted(asts)
    Pp = params(1.2, 0.15, mc=1000); Pr = params(2.0, 0.25, mc=1000)
    for P in (Pp, Pr):
        P.prize = np.zeros(300)
        for t in pool: P.prize[t - 1] = 1.0
    excl = mask_of((set(ALL) - set(pool)) | UNREACHABLE)
    tL = float(st['tL']); rE, vE = E.earth_state(tL)
    R = E.ast_states_at(np.array(asts) - 1, np.array(tfs))[0]
    v1, v2 = lambert(rE[None], R[0][None], np.array([tfs[0] - tL]))
    vinf = v1[0] - vE; vinf = vinf * min(1.0, VINF_MAX / np.linalg.norm(vinf))
    s = State(t=tfs[0], r=R[0], v=v2[0], m=M0P, visited=excl | mask_of([asts[0]]), seq=((asts[0], tfs[0], 0.0, 0.0),),
              fuel=0.0, t_launch=tL, vinf=vinf, m0=M0P, rare=1.0)
    legs = []
    for k in range(1, len(asts)):
        tgt = asts[k]; rec = dict(k=k, ast=tgt, tof_src_d=round((tfs[k] - tfs[k - 1]) / DAY, 1))
        pick = None
        for tag, P in (('prod', Pp), ('rel', Pr)):
            ch = [c for c in expand(E, s, P) if c.seq[-1][0] == tgt]
            rec[tag] = bool(ch)
            if ch and pick is None:
                pick = min(ch, key=lambda c: abs(c.seq[-1][1] - tfs[k]))
                rec['dv_leg'] = round(pick.seq[-1][2], 3); rec['dt_d'] = round((pick.seq[-1][1] - tfs[k]) / DAY, 1)
        if pick is None:                      # re-seed on the source epoch with the Lambert arrival velocity
            va, vb = lambert(s.r[None], R[k][None], np.array([tfs[k] - s.t]))
            dv = float(np.linalg.norm(va[0] - s.v)) if np.all(np.isfinite(va)) else float('inf')
            rec['reseed_lam1_dv'] = round(dv, 3)
            vnew = vb[0] if np.all(np.isfinite(vb)) else s.v
            s = State(t=tfs[k], r=R[k], v=vnew, m=s.m, visited=s.visited | (1 << (tgt - 1)),
                      seq=s.seq + ((tgt, tfs[k], dv, tfs[k] - s.t),), fuel=s.fuel, t_launch=s.t_launch, vinf=s.vinf,
                      m0=M0P, rare=s.rare + 1)
        else:
            s = pick
        legs.append(rec)
    n = len(legs); runs = []; c = 0
    for l in legs:
        c = c + 1 if l['prod'] else 0; runs.append(c)
    return dict(legs=n, adm_prod=sum(l['prod'] for l in legs), adm_rel=sum(l['rel'] for l in legs),
                longest_prod_run=max(runs) if runs else 0, detail=legs)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('src'); ap.add_argument('out')
    ap.add_argument('--routes', default='')
    a = ap.parse_args()
    E = RI.eph(); F = RI.IFleet(a.src)
    names = a.routes.split(',') if a.routes else sorted(F.routes)
    res = {}; tic = time.time()
    for nme in names:
        r = chain(E, F.routes[nme]['st']); res[nme] = r
        print(f'{nme}: {r["legs"]} legs; planner-state admissible prod {r["adm_prod"]}, relaxed {r["adm_rel"]}; '
              f'longest prod run {r["longest_prod_run"]}', flush=True)
    tot = dict(legs=sum(r['legs'] for r in res.values()), prod=sum(r['adm_prod'] for r in res.values()),
               rel=sum(r['adm_rel'] for r in res.values()))
    print(f'ALL: {tot} ({time.time() - tic:.0f} s)')
    json.dump(dict(total=tot, routes=res), open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()
