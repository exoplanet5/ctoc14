"""DRIFTING-CARRIER solver (stage 10) — the architecture the measurements actually point to.

Why not the GTOC-style transfer database (stage 9).  We built one (26 694 encounter nodes, 13 M Lambert edges) and
measured it against a real 37-flyby route that settles at 16.64 km/s:

    same chain, exact epochs, impulsive Lambert junctions : 29.8 km/s   (the MODEL overprices by 1.8x)
    ... with 0.25 d epoch jitter                          : 33.6
    ... with 1.00 d epoch jitter                          : 70.4

The impulsive cost is near-singular in the flyby epoch, so no usable grid exists (GTOC9/GTOC11 databases work because
those are LEO rendezvous problems whose cost is smooth in nodal drift).  What IS smooth here is the craft's ORBIT:
it stays near 1 AU and changes slowly, and the element-space cost of changing it is calibrated to 2 % on our own fleet
(ctoc14/phasemodel).  So the state is the orbit, and the flyby epochs are OUTPUTS -- exactly why carrier routes are the
cheapest we have (0.032-0.034 J per flyby against t10d's fleet 0.0415).

The one thing a FIXED carrier cannot do is serve more than ~20 targets: its cheap set is fixed by geometry.  But
moving to a neighbouring carrier costs

    dv = K_TOTAL [ K_DRIFT |d drift| + K_PLANE |d i_vec| ]  =  0.26-1.5 km/s   (median over the optimised library)

while ONE extra flyby is worth 1.29 km/s at a 900 kg tank.  A drift that unlocks another ~17 cheap targets is therefore
worth roughly ten times what it costs.  This solver searches that: a route is a SEQUENCE OF CARRIERS over the mission,
each serving its own targets, joined by continuous element drifts (no impulsive hop, no epoch grid).

  value(c, s)  = carrier_dp DP over carrier c's events inside time window s  (flybys - lam x phase dv)
  cost(c -> c')= the element metric above
  DP over (segment, carrier), then duplicate targets are blacklisted and the DP re-run until the set is clean.

Usage: carrier_drift.py lib.pkl out_dir [--segs 8] [--lam 1.0] [--settle]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pickle, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import carrier_dp as CD
from ctoc14.constants import DAY, AU, T_MISSION, VE
from ctoc14.kepler import Ephemeris, load_mea
from ctoc14.phasemodel import K_DRIFT, K_PLANE, K_TOTAL, drift as drift_of
from ctoc14.search import UNREACHABLE
from exp_carrier import rv2frame
_S = {}


def w_cell(job):
    """DP value of carrier c inside window [t0, t1): (value, [(target, epoch, lag), ...])"""
    ci, t0, t1, w = job
    g = _S; a = g['a']; c = g['lib'][ci]
    eph = g['eph']
    rE, vE = eph.earth_state(c['tL'])
    r = np.array(rE, float); v = np.array(vE, float) + np.array(c['vinf'], float)
    try:
        t, lag, gap, ast, Pc = CD.events_for(dict(tL=c['tL'], r=r, v=v), g['ids'], g['elem'], a.eps)
    except Exception:
        return ci, t0, 0.0, []
    m = (t >= max(t0, c['tL'] + 20 * DAY)) & (t < t1)
    if m.sum() < 2: return ci, t0, 0.0, []
    t, lag, gap, ast = t[m], lag[m], gap[m], ast[m]
    ac_ = float(rv2frame(r[None], v[None])[3][0])
    try:
        val, seq = CD.dp_route(t, lag, gap, ast, max(t0, c['tL']), a.lam, a.kgap, ac=ac_, w=w)
    except Exception:
        return ci, t0, 0.0, []
    seen = set(); seq = [i for i in seq if not (ast[i] in seen or seen.add(ast[i]))]
    return ci, t0, float(val), [(int(ast[i]), float(t[i]), float(lag[i])) for i in seq]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('lib'); ap.add_argument('out')
    ap.add_argument('--segs', type=int, default=8); ap.add_argument('--lam', type=float, default=1.0)
    ap.add_argument('--kgap', type=float, default=10.0); ap.add_argument('--kphase', type=float, default=1.5)
    ap.add_argument('--eps', type=float, default=0.012); ap.add_argument('--lam-drift', type=float, default=0.8,
                    help='J per km/s charged to a carrier change (1 flyby ~ 1.29 km/s at 900 kg)')
    ap.add_argument('--top', type=int, default=240); ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--rounds', type=int, default=6); ap.add_argument('--nroutes', type=int, default=10)
    ap.add_argument('--settle', action='store_true'); ap.add_argument('--seed', type=int, default=0)
    a = ap.parse_args()
    CD.K_SLOPE = CD.K_PHYS * a.kphase
    eph = Ephemeris(); ids, elem = load_mea()
    keep = np.array([int(x) not in UNREACHABLE for x in ids]); ids = ids[keep]; elem = elem[keep]
    lib = pickle.load(open(a.lib, 'rb'))[:a.top]
    _S.update(a=a, eph=eph, ids=ids, elem=elem, lib=lib)
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    import run_ialns as RI; say = RI.logger(out / 'log.txt')

    # element-space transition cost between carriers (continuous thrust, phasemodel: 2 % on the real fleet)
    R = []; V = []
    for c in lib:
        rE, vE = eph.earth_state(c['tL']); R.append(rE); V.append(np.array(vE) + np.array(c['vinf']))
    p, nv, Ev, sma = rv2frame(np.array(R), np.array(V))
    dr = drift_of(sma * AU)
    D = K_TOTAL * (K_DRIFT * np.abs(dr[:, None] - dr[None, :]) +
                   K_PLANE * np.linalg.norm(nv[:, None, :] - nv[None, :, :], axis=2))
    edges = np.linspace(0.0, T_MISSION, a.segs + 1)
    say(f'{len(lib)} carriers x {a.segs} segments of {edges[1]/DAY:.0f} d; drift cost median {np.median(D):.2f} km/s')

    prize = np.ones(301)
    fleet = []
    with mp.get_context('fork').Pool(a.nproc) as pl:
        for rt in range(a.nroutes):
            taken_by_fleet = set(t for r in fleet for t in r['targets'])
            w = prize.copy()
            for t in taken_by_fleet: w[t] = 0.0
            cells = {}
            jobs = [(ci, edges[s], edges[s + 1], w) for ci in range(len(lib)) for s in range(a.segs)]
            for ci, t0, val, seq in pl.map(w_cell, jobs, chunksize=16):
                cells[(ci, float(t0))] = (val, seq)
            best = None
            for _ in range(a.rounds):
                # DP over (segment, carrier)
                NEG = -1e9
                Vv = np.full((a.segs, len(lib)), NEG); P = np.full((a.segs, len(lib)), -1, int)
                for ci in range(len(lib)):
                    v0, s0 = cells[(ci, float(edges[0]))]
                    Vv[0, ci] = v0 if lib[ci]['tL'] < edges[1] else NEG
                for s in range(1, a.segs):
                    prev = Vv[s - 1]
                    for ci in range(len(lib)):
                        cand = prev - a.lam_drift * D[:, ci]
                        j = int(np.argmax(cand))
                        if cand[j] <= NEG / 2: continue
                        Vv[s, ci] = cand[j] + cells[(ci, float(edges[s]))][0]; P[s, ci] = j
                s = a.segs - 1; ci = int(np.argmax(Vv[s]))
                if Vv[s, ci] <= NEG / 2: break
                chain = []
                while s >= 0:
                    chain.append((s, ci)); ci = P[s, ci]; s -= 1
                    if ci < 0: break
                chain = chain[::-1]
                segs = []; seen = set(); dup = set()
                for s, ci in chain:
                    val, seq = cells[(ci, float(edges[s]))]
                    keepseq = []
                    for tg, te, lg in seq:
                        if tg in seen: dup.add(tg)
                        else: seen.add(tg); keepseq.append((tg, te, lg))
                    segs.append(dict(ci=ci, s=s, base=keepseq))
                dv = sum(D[chain[k - 1][1], chain[k][1]] for k in range(1, len(chain)))
                nfb = len(seen)
                if best is None or nfb > best['nfb']: best = dict(nfb=nfb, segs=segs, dv=float(dv), chain=chain)
                if not dup: break
                for t in dup: w[t] = 0.0           # blacklist duplicates and re-run the DP
                for k in list(cells):
                    pass
                jobs = [(ci, edges[s], edges[s + 1], w) for ci in set(c for _, c in chain) for s in range(a.segs)]
                for ci, t0, val, seq in pl.map(w_cell, jobs, chunksize=8):
                    cells[(ci, float(t0))] = (val, seq)
            if best is None: break
            carr = [s['ci'] for s in best['segs'] if s['base']]
            fleet.append(dict(targets=sorted(set(t for s in best['segs'] for t, _, _ in s['base'])),
                              segs=best['segs'], dv=best['dv'], carriers=carr))
            say(f'route {rt}: {best["nfb"]} flybys over {len(set(carr))} carriers, drift dv {best["dv"]:.2f} km/s, '
                f'segments ' + '+'.join(str(len(s['base'])) for s in best['segs']))
    cov = sorted(set(t for r in fleet for t in r['targets']))
    say(f'fleet: {len(fleet)} routes, {len(cov)} distinct targets')
    json.dump([dict(targets=[int(t) for t in r['targets']], dv=float(r['dv']),
                    segs=[dict(ci=int(s['ci']), base=[[int(x), float(y), float(z)] for x, y, z in s['base']])
                          for s in r['segs']]) for r in fleet],
              open(out / 'routes.json', 'w'), indent=1)
    say(f'-> {out}/routes.json')


if __name__ == '__main__':
    main()
