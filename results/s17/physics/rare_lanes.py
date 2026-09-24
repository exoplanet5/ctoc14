"""T1: how many of the 95 rare (<= 4 node events) targets can 8 fixed-carrier phase lanes cover disjointly?"""
import os
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'): os.environ.setdefault(v,'1')
import sys, json, time, pathlib, multiprocessing as mp
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT/'tools'))
import numpy as np
import carrier_dp as CD
from exp_carrier import rv2frame, frame, gaps
from ctoc14.kepler import load_mea, Ephemeris
from ctoc14.constants import DAY, VINF_MAX
OUT = ROOT/'results/s17/physics'
EPS, KGAP, KPH, LAM = 0.03, 10.0, 1.5, 0.5
CD.K_SLOPE = CD.K_PHYS * KPH
ev = json.load(open(OUT/'events.json'))
nev = np.zeros(301, int)
for e in ev:
    if e['d'] < 0.05: nev[e['ast']] += 1
ids_all, elem_all = load_mea()
TG = [X for X in range(1,301) if X not in (131,144)]
RARE = [X for X in TG if nev[X] <= 4]
sel = np.array([int(x) in set(RARE) for x in ids_all]); ids = ids_all[sel]; elem = elem_all[sel]
selA = np.array([int(x) in set(TG) for x in ids_all]); idsA = ids_all[selA]; elemA = elem_all[selA]
E = Ephemeris()
G = {}
def lane(job):
    tL, r, v, w, allt = job
    ii, ee = (idsA, elemA) if allt else (ids, elem)
    try:
        t, lag, gap, ast, Pc = CD.events_for(dict(tL=tL, r=r, v=v), ii, ee, EPS)
        if len(t) < 2: return None
        ac = float(rv2frame(r[None], v[None])[3][0])
        val, seq = CD.dp_route(t, lag, gap, ast, tL, LAM, KGAP, ac=ac, w=w)
    except Exception as ex:
        return None
    seen = set(); seq = [i for i in seq if not (ast[i] in seen or seen.add(ast[i]))]
    if not seq: return None
    tl = np.concatenate([[tL/DAY], t[seq]/DAY]); lg = np.concatenate([[0.0], lag[seq]])
    S = np.diff(lg)/np.diff(tl)
    dvp = CD.K_SLOPE*np.abs(np.diff(S)).sum(); dvg = KGAP*np.abs(gap[seq] - (2/3)*ac*S).sum()
    return dict(val=float(val), tg=[int(ast[i]) for i in seq], t=[float(t[i]/DAY) for i in seq], dv=float(dvp+dvg), dvp=float(dvp))
if __name__ == '__main__':
    t0 = time.time(); rng = np.random.default_rng(17); n = 40000
    tLs = rng.uniform(0, 730*DAY, n)
    d = rng.normal(size=(n,3)); d /= np.linalg.norm(d, axis=1)[:,None]
    vinf = d*(VINF_MAX*rng.uniform(0,1,n)**(1/3))[:,None]
    RV = np.array([np.concatenate(E.earth_state(x)) for x in tLs]); r = RV[:,:3]; v = RV[:,3:]+vinf
    pa, na, Ea = frame(elem[:, :5]); cnt = np.zeros(n, int)
    for c0 in range(0, n, 10000):
        s = slice(c0, c0+10000); pc, nc, Ec, ac = rv2frame(r[s], v[s]); g, _ = gaps(pc, nc, Ec, pa, na, Ea)
        cnt[s] = (np.abs(g).min(2) < EPS).sum(1)
    top = np.argsort(-cnt)[:3000]
    print(f'{n} carriers: rare candidates (gap<{EPS}) best {cnt.max()} median {np.median(cnt):.0f} top3000 min {cnt[top].min()}  ({time.time()-t0:.0f}s)', flush=True)
    w1 = np.ones(301)
    with mp.get_context('fork').Pool(2) as pl:
        L = pl.map(lane, [(float(tLs[c]), r[c], v[c], w1, False) for c in top], chunksize=20)
        lib = [(c, x) for c, x in zip(top, L) if x]
        ns = np.array([len(x['tg']) for _, x in lib])
        print(f'library {len(lib)} lanes: rare per single lane best {ns.max()} median {np.median(ns):.0f}; ({time.time()-t0:.0f}s)', flush=True)
        # greedy 8
        held = {}; chosen = []
        for k in range(8):
            got = set(x for s in held.values() for x in s)
            j = max(range(len(lib)), key=lambda j: (len(set(lib[j][1]['tg']) - got), -lib[j][1]['dv']))
            c, x = lib[j]; chosen.append(c); held[k] = [y for y in x['tg'] if y not in got]
        cov = set(x for s in held.values() for x in s)
        print(f'greedy 8 lanes: rare covered {len(cov)} (per lane {[len(s) for s in held.values()]})', flush=True)
        # re-pricing rounds: each slot re-solved over the top-300 library carriers with w=0 on others' targets
        cand = [c for c, _ in sorted(lib, key=lambda cx: -len(cx[1]['tg']))[:300]]
        res = {}
        for rnd in range(3):
            for k in range(8):
                others = set(x for kk, s in held.items() if kk != k for x in s)
                w = np.ones(301); w[list(others)] = -0.5
                R = pl.map(lane, [(float(tLs[c]), r[c], v[c], w, False) for c in cand], chunksize=10)
                best = None
                for c, x in zip(cand, R):
                    if not x: continue
                    new = [y for y in x['tg'] if y not in others]
                    sc = len(new) - 0.05*x['dv']
                    if best is None or sc > best[0]: best = (sc, c, new, x)
                _, c, new, x = best; chosen[k] = c; held[k] = new; res[k] = dict(tL=float(tLs[c]/DAY), vinf=vinf[c].tolist(), tg=new, dv=x['dv'], dvp=x['dvp'], n_lane=len(x['tg']))
            cov = set(x for s in held.values() for x in s)
            print(f'round {rnd+1}: rare covered {len(cov)}/{len(RARE)}  per lane {[len(held[k]) for k in range(8)]}  model dv {[round(res[k]["dv"],1) for k in range(8)]}  ({time.time()-t0:.0f}s)', flush=True)
        # sanity: same 8 carriers on ALL targets
        A = pl.map(lane, [(float(tLs[c]), r[c], v[c], w1, True) for c in chosen])
        print('sanity, same carriers on ALL targets: flybys', [len(a['tg']) if a else 0 for a in A], ' model dv', [round(a['dv'],1) if a else 0 for a in A])
    json.dump(dict(rare=RARE, lanes=res, covered=sorted(cov), missing=sorted(set(RARE)-cov)), open(OUT/'rare_lanes.json','w'), indent=1)
    print('missing rare:', sorted(set(RARE)-cov))
