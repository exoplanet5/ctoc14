"""Choose the 8 carriers JOINTLY so that their cheap neighbourhoods tile the catalogue (stage 7, section 6).

library : top --pool carriers by DP value (no exclusions).  For each: the DP base route and, from the UNSETTLED seed
          trajectory (ballistic carrier + DP phasing impulses), every target's closest approach (distance, epoch).
select  : transportation MILP.  y_c binary, x_tc in [0,1] (integral for fixed y), x_tc <= y_c, sum_c x_tc <= 1,
          sum_t x_tc <= cap, sum_c y_c = N; maximise covered targets, tie-break on distance.  Carriers are pre-filtered
          by randomised greedy max-coverage.
Usage: carrier_tile.py library out.pkl [--n 60000 --pool 4000]
       carrier_tile.py select lib.pkl out_dir [--delta 0.06 --cap 40 --ncraft 8]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pickle, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
from ctoc14.constants import AU, DAY, VINF_MAX, T_MISSION
from ctoc14.kepler import load_mea, Ephemeris
from ctoc14.search import UNREACHABLE
import carrier_dp as CD
from exp_carrier import frame, rv2frame, gaps

_G = {}


def w_lib(c):
    g = _G; a = g['a']; eph = g['eph']
    import run_ialns as RI
    from ctoc14.globalopt import extend_arc
    car = dict(tL=float(g['tLs'][c]), r=g['r'][c], v=g['v'][c])
    t, lag, gap, ast, Pc = CD.events_for(car, g['ids'], g['elem'], a.eps)
    if len(t) < 3: return None
    ac_ = float(rv2frame(car['r'][None], car['v'][None])[3][0])
    val, seq = CD.dp_route(t, lag, gap, ast, car['tL'], a.lam, a.kgap, ac=ac_)
    seen = set(); seq = [i for i in seq if not (ast[i] in seen or seen.add(ast[i]))]
    if len(seq) < a.min_base: return None
    ip = CD.seed_ip(eph, car['tL'], g['vinf'][c], t[seq], lag[seq])
    extend_arc(ip, T_MISSION - 2 * DAY)
    tt, r = RI.trajectory(ip, step=2 * DAY)
    tb = t[seq]; neigh = {}
    free = np.abs(tt[:, None] - tb[None, :]).min(1) > 12 * DAY
    for X in g['ids']:
        ra, _ = eph.ast_states_at(np.full(len(tt), int(X) - 1), tt)
        d = np.linalg.norm(r - ra, axis=1) / AU; d[~free] = 9.0
        k = int(d.argmin())
        if d[k] < a.dmax: neigh[int(X)] = (float(d[k]), float(tt[k]))
    return dict(c=int(c), tL=car['tL'], vinf=g['vinf'][c].tolist(), val=float(val),
                base=[(int(ast[i]), float(t[i]), float(lag[i])) for i in seq], neigh=neigh)


def cmd_library(a):
    rng = np.random.default_rng(a.seed); eph = Ephemeris()
    ids, elem = load_mea(); keep = np.array([int(x) not in UNREACHABLE for x in ids]); ids = ids[keep]; elem = elem[keep]
    CD.K_SLOPE = CD.K_PHYS * a.kphase
    tLs = rng.uniform(0, 365.25 * DAY, a.n)
    d = rng.normal(size=(a.n, 3)); d /= np.linalg.norm(d, axis=1)[:, None]
    vinf = d * (VINF_MAX * rng.uniform(0, 1, a.n) ** (1 / 3))[:, None]
    RV = np.array([np.concatenate(eph.earth_state(x)) for x in tLs]); r = RV[:, :3]; v = RV[:, 3:] + vinf
    CD._G = dict(tLs=tLs, r=r, v=v, ids=ids, elem=elem, a=a)
    _G.update(tLs=tLs, r=r, v=v, vinf=vinf, ids=ids, elem=elem, a=a, eph=eph)
    tic = time.time()
    with mp.get_context('fork').Pool(a.nproc) as pl:
        vals = np.array(pl.map(CD._value, range(a.n), chunksize=64))
        top = np.argsort(-vals[:, 0])[:a.pool]
        print(f'DP over {a.n} carriers: {time.time()-tic:.0f} s; pool value {vals[top,0].min():.1f}..{vals[top,0].max():.1f}', flush=True)
        lib = [x for x in pl.map(w_lib, [int(c) for c in top], chunksize=4) if x]
    pickle.dump(lib, open(a.out, 'wb'))
    nb = [len(x['base']) for x in lib]; nn = [sum(1 for d, _ in x['neigh'].values() if d < 0.06) for x in lib]
    print(f'{len(lib)} carrier routes ({time.time()-tic:.0f} s): base {np.median(nb):.0f} (max {max(nb)}), '
          f'neighbours < 0.06 AU {np.median(nn):.0f} (max {max(nn)})')


def sets_of(lib, delta):
    S = []
    for x in lib:
        s = {t: 0.0 for t, _, _ in x['base']}
        for t, (d, _) in x['neigh'].items():
            if d < delta and t not in s: s[t] = d / delta
        S.append(s)
    return S


def cmd_select(a):
    from scipy.optimize import milp, LinearConstraint, Bounds
    from scipy.sparse import coo_matrix
    lib = pickle.load(open(a.lib, 'rb')); rng = np.random.default_rng(a.seed)
    for delta in [float(x) for x in a.delta.split(',')]:
        S = sets_of(lib, delta)
        # randomised greedy max coverage to pre-filter carriers
        keep = set()
        for rep in range(a.greedy):
            left = set(range(1, 301)) - set(UNREACHABLE); chosen = []
            for k in range(a.ncraft):
                gain = np.array([len(left & s.keys()) for s in S], float) + rng.uniform(0, 3 if rep else 0, len(S))
                j = int(gain.argmax()); chosen.append(j); left -= set(sorted(S[j], key=S[j].get)[:10**6])
            keep |= set(chosen)
            if rep == 0: print(f'delta {delta}: plain greedy union of {a.ncraft} (no cap): {298 - len(left)}')
        K = sorted(keep); T = sorted(set(t for j in K for t in S[j]))
        ti = {t: i for i, t in enumerate(T)}; pairs = [(ci, ti[t], S[j][t]) for ci, j in enumerate(K) for t in S[j]]
        nx = len(pairs); ny = len(K); n = nx + ny
        cvec = np.concatenate([-(1.0 - 0.02 * np.array([p[2] for p in pairs])), np.zeros(ny)])
        R = []; C = []; V = []; lo = []; hi = []; row = 0
        for t in range(len(T)):                                       # each target at most once
            for k, p in enumerate(pairs):
                pass
        bt = [[] for _ in T]; bc = [[] for _ in K]
        for k, p in enumerate(pairs): bt[p[1]].append(k); bc[p[0]].append(k)
        for lst in bt:
            R += [row] * len(lst); C += lst; V += [1.0] * len(lst); lo.append(0); hi.append(1); row += 1
        for ci, lst in enumerate(bc):                                 # capacity, linked to y
            R += [row] * (len(lst) + 1); C += lst + [nx + ci]; V += [1.0] * len(lst) + [-float(a.cap)]; lo.append(-np.inf); hi.append(0); row += 1
        for k, p in enumerate(pairs):                                 # x <= y
            R += [row, row]; C += [k, nx + p[0]]; V += [1.0, -1.0]; lo.append(-np.inf); hi.append(0); row += 1
        R += [row] * ny; C += list(range(nx, n)); V += [1.0] * ny; lo.append(a.ncraft); hi.append(a.ncraft); row += 1
        A = coo_matrix((V, (R, C)), shape=(row, n)).tocsr()
        integ = np.concatenate([np.zeros(nx), np.ones(ny)])
        tic = time.time()
        res = milp(cvec, constraints=LinearConstraint(A, lo, hi), integrality=integ, bounds=Bounds(0, 1),
                   options=dict(time_limit=a.tlim, mip_rel_gap=0.002))
        x = res.x; ys = [K[ci] for ci in range(ny) if x[nx + ci] > 0.5]
        asg = {j: [] for j in ys}
        for k, p in enumerate(pairs):
            if x[k] > 0.5: asg[K[p[0]]].append(T[p[1]])
        cov = sum(len(v) for v in asg.values())
        print(f'delta {delta}: MILP over {ny} carriers ({time.time()-tic:.0f} s, {res.message[:40]}): '
              f'{a.ncraft} carriers tile {cov} of 298 (cap {a.cap}) | ' +
              ' '.join(f'{len(lib[j]["base"])}+{len(asg[j]) - len(set(asg[j]) & set(t for t,_,_ in lib[j]["base"]))}' for j in ys), flush=True)
        out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
        json.dump(dict(delta=delta, cap=a.cap, covered=cov,
                       carriers=[dict(lib=j, tL=lib[j]['tL'], vinf=lib[j]['vinf'], base=lib[j]['base'], pool=sorted(asg[j]),
                                      neigh={str(t): lib[j]['neigh'][t] for t in asg[j] if t in lib[j]['neigh']}) for j in ys]),
                  open(out / f'tile_{delta}.json', 'w'), indent=1)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest='cmd')
    p = sub.add_parser('library'); p.add_argument('out'); p.add_argument('--n', type=int, default=60000)
    p.add_argument('--pool', type=int, default=4000); p.add_argument('--eps', type=float, default=0.012)
    p.add_argument('--lam', type=float, default=1.0); p.add_argument('--kgap', type=float, default=10.0)
    p.add_argument('--kphase', type=float, default=1.5); p.add_argument('--dmax', type=float, default=0.15)
    p.add_argument('--min-base', type=int, default=8); p.add_argument('--seed', type=int, default=0)
    p.add_argument('--nproc', type=int, default=8); p.add_argument('--static-gap', action='store_true')
    p = sub.add_parser('select'); p.add_argument('lib'); p.add_argument('out'); p.add_argument('--delta', default='0.04,0.06,0.1')
    p.add_argument('--cap', type=int, default=40); p.add_argument('--ncraft', type=int, default=8)
    p.add_argument('--greedy', type=int, default=60); p.add_argument('--seed', type=int, default=0); p.add_argument('--tlim', type=float, default=300)
    a = ap.parse_args()
    dict(library=cmd_library, select=cmd_select)[a.cmd](a)
