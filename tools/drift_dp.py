"""PHASE-CONTINUOUS drifting-carrier DP (stage 10, second version).

Why a second version.  `tools/carrier_drift.py` let a route change carrier at fixed segment boundaries and charged the
element-space metric for the change.  Its DP base rose from 17-22 to 28 flybys for ~3.3 km/s of modelled drift, but the
settled twin gave 17-22 flybys at 842-865 kg (0.055-0.072 J/fb) and 3 of 6 routes did not settle at all.  The defect is
structural: each carrier's event lags are measured from ITS OWN launch, so at a segment boundary the DP silently teleported
the craft from phase `lag` on carrier c to phase ~0 on carrier c'.  The real trajectory pays for that teleport (13-14 km/s
actual against ~3 modelled).

The fix here.  Carry the PHASE through the DP, exactly, with a common reference.

  * "the craft is on carrier c with lag l" means: at time t it sits at R_c(t - l), where R_c is c's ballistic trajectory.
    That is precisely how `carrier_dp.events_for` defines the lag of an event, so nothing has to be re-derived.
  * The physical, carrier-independent quantity is the craft's HELIOCENTRIC LONGITUDE.  Write Psi_c(tau) for the unwrapped
    ecliptic longitude of R_c(tau) (monotone: every carrier is prograde with e < 0.25).  Then the craft's longitude at
    time t is Psi_c(t - l).
  * Continuity at a carrier change at time t therefore reads  Psi_c(t - l) = Psi_c'(t - l')  (mod 2 pi), whose solution is

        l'  =  l + D(c, c', tau),        tau = t - l,        D(c, c', tau) = tau - Psi_c'^{-1}( Psi_c(tau) ),

    the branch of the inverse nearest tau.  D is a slowly varying function of tau (it drifts at the two carriers'
    drift-rate difference), so it is tabulated once per ordered pair on a 5-day grid and interpolated.  A carrier change
    is then a RELABELLING of the lag, not a reset of it, and the DP's phase is continuous by construction.

  * The slope state is made ABSOLUTE so that a carrier change also pays for its semi-major-axis change:
        sigma = P_c (1 + s) / P_0 - 1        (P_0 = 365.25 d, s = dl/dt on carrier c)
    is the craft's true relative period offset; changing it costs K_SLOPE |d sigma| = (v/3) |dP/P| km/s, the same
    constant `carrier_dp` uses and the same number as phasemodel.K_DRIFT per deg/yr.  `carrier_drift.py` charged
    K_DRIFT |d drift| at the boundary AND let the DP charge slope changes inside a segment -- here the one sigma state
    does both, with no double count.
  * What a carrier change then still costs is only the part of the orbit a tangential burn cannot buy: the plane and the
    eccentricity vector,  K_TOTAL ( K_PLANE |d i_vec| + kecc K_ECC |d e_vec| ).
  * The residual node gap is charged at the slope the craft ARRIVES with, gap - (2/3) a sigma_local, exactly as in
    carrier_dp (a period offset lifts the whole track).

So the route is ONE smooth trajectory: a phase curve, piecewise linear in (t, longitude), whose breakpoints are flybys
and whose track may morph from carrier to carrier.  There are no segments and no epoch grid; flyby epochs are outputs.

Usage:
    drift_dp.py out_dir [--seeds 24] [--knb 24] [--lam 1.0] [--kgap 10] [--kplane 1.0] [--settle]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pickle, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import carrier_dp as CD
from ctoc14.constants import MU, AU, DAY, T_MISSION, VE, VINF_MAX
from ctoc14.kepler import Ephemeris, load_mea, solve_kepler
from ctoc14.phasemodel import K_DRIFT, K_PLANE, K_ECC, K_TOTAL, drift as drift_of, elements
from ctoc14.search import UNREACHABLE
from exp_carrier import rv2frame

P0 = 365.25 * DAY                     # reference period for the absolute slope coordinate sigma
TAU0 = -900.0 * DAY                   # phase tables cover [TAU0, TAU1] so that tau = t - lag is always inside
TAU1 = T_MISSION + 900.0 * DAY


# ----------------------------------------------------------------------------- carrier phase geometry
def psi_of(r0, v0, taus, t0):
    """Unwrapped heliocentric ecliptic longitude of the ballistic orbit through (r0, v0) at t0, sampled at `taus` [s].

    Analytic Kepler (no universal-variable iteration per sample): elements once, then M -> E -> f -> r."""
    r0 = np.asarray(r0, float); v0 = np.asarray(v0, float)
    rn = np.linalg.norm(r0); v2 = v0 @ v0
    a = 1.0 / (2.0 / rn - v2 / MU)
    h = np.cross(r0, v0); hn = np.linalg.norm(h); hh = h / hn
    ev = np.cross(v0, h) / MU - r0 / rn; e = np.linalg.norm(ev)
    P = ev / e if e > 1e-9 else (r0 / rn)
    Q = np.cross(hh, P)
    f0 = np.arctan2(Q @ r0, P @ r0)
    E0 = 2 * np.arctan2(np.sqrt(1 - e) * np.sin(f0 / 2), np.sqrt(1 + e) * np.cos(f0 / 2))
    M0 = E0 - e * np.sin(E0)
    n = np.sqrt(MU / a ** 3)
    M = M0 + n * (np.asarray(taus, float) - t0)
    E = solve_kepler(M, e)
    f = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
    rr = a * (1 - e * np.cos(E))
    xyz = rr[:, None] * (np.cos(f)[:, None] * P[None, :] + np.sin(f)[:, None] * Q[None, :])
    th = np.arctan2(xyz[:, 1], xyz[:, 0])
    return np.unwrap(th), 2 * np.pi / n, a


class Pool:
    """A set of carriers with their phase tables, pairwise lag-offset tables and switch costs."""

    def __init__(self, lib, idx, eph, dtau=5.0 * DAY):
        self.idx = list(idx)
        self.taus = np.arange(TAU0, TAU1, dtau)
        self.tL = np.array([lib[i]['tL'] for i in self.idx])
        R = []; V = []
        for i in self.idx:
            rE, vE = eph.earth_state(lib[i]['tL'])
            R.append(np.array(rE, float)); V.append(np.array(vE, float) + np.array(lib[i]['vinf'], float))
        self.r = np.array(R); self.v = np.array(V); self.vinf = np.array([lib[i]['vinf'] for i in self.idx], float)
        self.psi = np.empty((len(self.idx), len(self.taus)))
        self.P = np.empty(len(self.idx)); self.a = np.empty(len(self.idx))
        for k in range(len(self.idx)):
            self.psi[k], self.P[k], ak = psi_of(self.r[k], self.v[k], self.taus, self.tL[k])
            self.a[k] = ak / AU
        el = elements(self.r, self.v)
        self.ivec = el['i_vec']; self.evec = el['e_vec']
        # ordered-pair lag offset D(c, c', tau)  [seconds]
        nc = len(self.idx)
        self.D = np.zeros((nc, nc, len(self.taus)))
        for b in range(nc):
            pb = self.psi[b]
            for c in range(nc):
                if b == c: continue
                # unwrapped value of psi_b that is congruent to psi_c(tau) and nearest psi_b(tau)
                d = self.psi[c] - pb
                tgt = pb + (d + np.pi) % (2 * np.pi) - np.pi
                taup = np.interp(tgt, pb, self.taus)          # psi_b is monotone increasing
                self.D[c, b] = self.taus - taup               # l_on_b = l_on_c + D[c, b](tau)

    def delta(self, c, tau):
        """vector over pool carriers b of D(c, b, tau) [s]."""
        j = np.clip(np.searchsorted(self.taus, tau) - 1, 0, len(self.taus) - 2)
        w = (tau - self.taus[j]) / (self.taus[j + 1] - self.taus[j])
        return self.D[c, :, j] * (1 - w) + self.D[c, :, j + 1] * w

    def switch_cost(self, kecc=1.0):
        """km/s of a carrier change: plane + eccentricity vector (the a change is the sigma state, not this)."""
        dp = np.linalg.norm(self.ivec[:, None, :] - self.ivec[None, :, :], axis=2)
        de = np.linalg.norm(self.evec[:, None, :] - self.evec[None, :, :], axis=2)
        return K_TOTAL * (K_PLANE * dp + kecc * K_ECC * de)


# ----------------------------------------------------------------------------- events
def pool_events(pool, ids, elem, eps, images=1):
    """Pooled event list over every carrier of the pool, sorted by time.
    Returns dict of arrays t [s], lag [s], gap [AU], ast, car (pool index), plus per-carrier a [AU]."""
    T = []; L = []; G = []; A = []; C = []
    for k in range(len(pool.idx)):
        car = dict(tL=float(pool.tL[k]), r=pool.r[k], v=pool.v[k])
        try:
            t, lag, gap, ast, Pc = CD.events_for(car, ids, elem, eps, images=images)
        except Exception:
            continue
        if len(t) == 0: continue
        T.append(t); L.append(lag * DAY); G.append(gap); A.append(ast); C.append(np.full(len(t), k))
    if not T: return None
    t = np.concatenate(T); lag = np.concatenate(L); gap = np.concatenate(G)
    ast = np.concatenate(A); car = np.concatenate(C)
    o = np.argsort(t)
    return dict(t=t[o], lag=lag[o], gap=gap[o], ast=ast[o], car=car[o])


# ----------------------------------------------------------------------------- the DP
def dp_drift(ev, pool, SW, lam=1.0, kgap=10.0, smax=0.08, s0max=0.05, sigmax=0.30,
             dtmin=15.0, dtmax=900.0, w=None, maxlab=28, kslope=None):
    """Line-graph DP on the pooled (time, phase) plane.  State = incoming arc (carries sigma exactly).

    Returns (value, [event indices]).  Arc i -> j:
        l_i on carrier(j) = lag[i] + D(car_i, car_j, t_i - lag_i)          (phase continuity, exact)
        s   = (lag[j] - l_i) / (t_j - t_i)                                 (local slope on carrier j)
        sig = P_j (1 + s) / P0 - 1                                         (absolute period offset)
        cost = lam * ( K_SLOPE |sig - sig_in| + kgap |gap_j - (2/3) a_j s| + SW[car_i, car_j] )
    """
    K = CD.K_SLOPE if kslope is None else kslope
    t = ev['t']; lag = ev['lag']; gap = ev['gap']; ast = ev['ast']; car = ev['car']
    n = len(t); td = t / DAY
    gain = np.ones(n) if w is None else np.asarray(w, float)[ast]
    Pc = pool.P; ac = pool.a
    sig_of = lambda j, s: Pc[car[j]] * (1.0 + s) / P0 - 1.0
    # launch labels: carrier's own launch, lag 0 at tL, small initial slope, sigma is free (bought by v_inf)
    lab = [[] for _ in range(n)]
    tl = pool.tL[car]
    s0 = lag / np.maximum(t - tl, 1.0)
    ok0 = (np.abs(s0) <= s0max) & (t - tl < 2 * dtmax * DAY) & (t > tl + 20 * DAY)
    for i in np.nonzero(ok0)[0]:
        sg = sig_of(i, s0[i])          # sigma at launch is bought by v_inf: free, and not bounded by sigmax
        lab[i].append((gain[i] - lam * kgap * abs(gap[i] - (2.0 / 3.0) * ac[car[i]] * s0[i]), sg, -1, -1))
    lo = np.searchsorted(t, t + dtmin * DAY, side='left')
    hi = np.searchsorted(t, t + dtmax * DAY, side='right')
    best = (-1e9, -1, -1)
    for i in range(n):
        L = lab[i]
        if not L: continue
        if len(L) > maxlab:
            L.sort(key=lambda x: -x[0]); L = L[:maxlab]; lab[i] = L
        V = np.array([x[0] for x in L]); S = np.array([x[1] for x in L])
        k = int(V.argmax())
        if V[k] > best[0]: best = (V[k], i, k)
        a0, b0 = lo[i], hi[i]
        if b0 <= a0: continue
        J = np.arange(a0, b0)
        d = pool.delta(int(car[i]), float(t[i] - lag[i]))            # (npool,)
        li = lag[i] + d[car[J]]
        dt = t[J] - t[i]
        s = (lag[J] - li) / dt
        m = (np.abs(s) <= smax) & (ast[J] != ast[i])
        if not m.any(): continue
        J = J[m]; s = s[m]
        sg = Pc[car[J]] * (1.0 + s) / P0 - 1.0
        m2 = np.abs(sg) <= sigmax
        if not m2.any(): continue
        J = J[m2]; s = s[m2]; sg = sg[m2]
        arc = gain[J] - lam * (kgap * np.abs(gap[J] - (2.0 / 3.0) * ac[car[J]] * s) + SW[car[i], car[J]])
        val = V[:, None] - lam * K * np.abs(sg[None, :] - S[:, None]) + arc[None, :]
        kk = val.argmax(0)
        vv = val[kk, np.arange(len(J))]
        for c in range(len(J)):
            lab[J[c]].append((vv[c], sg[c], i, int(kk[c])))
    seq = []; _, i, k = best
    while i >= 0:
        seq.append(i); _, _, pi, pk = lab[i][k]; i, k = pi, pk
    return best[0], seq[::-1]


def route_report(ev, pool, SW, seq, lam, kgap, kslope=None):
    """model dv breakdown of a DP sequence."""
    K = CD.K_SLOPE if kslope is None else kslope
    t = ev['t']; lag = ev['lag']; gap = ev['gap']; car = ev['car']
    sig = []; gp = []; sw = 0.0
    for a, b in zip([None] + list(seq[:-1]), seq):
        if a is None:
            s = lag[b] / max(t[b] - pool.tL[car[b]], 1.0)
        else:
            d = pool.delta(int(car[a]), float(t[a] - lag[a]))
            s = (lag[b] - (lag[a] + d[car[b]])) / (t[b] - t[a])
            sw += SW[car[a], car[b]]
        sig.append(pool.P[car[b]] * (1 + s) / P0 - 1.0)
        gp.append(abs(gap[b] - (2.0 / 3.0) * pool.a[car[b]] * s))
    sig = np.array(sig)
    return dict(dv_phase=float(K * np.abs(np.diff(sig)).sum()), dv_gap=float(kgap * np.sum(gp)),
                dv_switch=float(sw), nswitch=int(sum(1 for a, b in zip(seq[:-1], seq[1:]) if car[a] != car[b])),
                ncar=int(len(set(car[seq].tolist()))), sigma=[float(x) for x in sig])


# ----------------------------------------------------------------------------- driver
_S = {}


def w_route(job):
    seed, w, args = job
    g = _S; pool = g['pools'][seed]; ev = g['evs'][seed]; SW = g['sws'][seed]
    if ev is None: return seed, None
    val, seq = dp_drift(ev, pool, SW, lam=args.lam, kgap=args.kgap, smax=args.smax, s0max=args.s0max,
                        sigmax=args.sigmax, dtmax=args.dtmax, w=w, maxlab=args.maxlab)
    if len(seq) < 4: return seed, None
    seen = set(); seq = [i for i in seq if not (ev['ast'][i] in seen or seen.add(int(ev['ast'][i])))]
    rep = route_report(ev, pool, SW, seq, args.lam, args.kgap)
    c0 = int(ev['car'][seq[0]])
    rep.update(val=float(val), n=len(seq), seed=seed,
               tL=float(pool.tL[c0]), vinf=[float(x) for x in pool.vinf[c0]],
               legs=[[int(ev['ast'][i]), float(ev['t'][i])] for i in seq],
               carriers=[int(pool.idx[int(ev['car'][i])]) for i in seq])
    return seed, rep


def build_pool(job):
    seed, idx, args = job
    eph = Ephemeris(); ids = _S['ids']; elem = _S['elem']; lib = _S['lib']
    pool = Pool(lib, idx, eph)
    ev = pool_events(pool, ids, elem, args.eps, images=args.images)
    return seed, pool, ev, pool.switch_cost(args.kecc)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out'); ap.add_argument('--lib', default='results/s7/libopt.pkl')
    ap.add_argument('--seeds', type=int, default=16, help='number of seed carriers (one pool each)')
    ap.add_argument('--knb', type=int, default=24, help='carriers per pool (seed + nearest by switch cost)')
    ap.add_argument('--top', type=int, default=240)
    ap.add_argument('--eps', type=float, default=0.012); ap.add_argument('--images', type=int, default=1)
    ap.add_argument('--lam', type=float, default=1.0); ap.add_argument('--kgap', type=float, default=10.0)
    ap.add_argument('--kphase', type=float, default=1.5); ap.add_argument('--kecc', type=float, default=1.0)
    ap.add_argument('--smax', type=float, default=0.08); ap.add_argument('--s0max', type=float, default=0.05)
    ap.add_argument('--sigmax', type=float, default=0.30); ap.add_argument('--dtmax', type=float, default=900.0)
    ap.add_argument('--maxlab', type=int, default=28); ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--sdiv', type=float, default=0.30, help='min switch-metric separation between seed carriers')
    ap.add_argument('--rounds', type=int, default=1, help='>1: build a covering pool, re-pricing taken targets to 0')
    ap.add_argument('--nmin', type=int, default=8, help='drop routes shorter than this')
    ap.add_argument('--fixed', action='store_true', help='control: pool = seed carrier only (no drifting)')
    a = ap.parse_args()
    CD.K_SLOPE = CD.K_PHYS * a.kphase
    eph = Ephemeris(); ids, elem = load_mea()
    keep = np.array([int(x) not in UNREACHABLE for x in ids]); ids = ids[keep]; elem = elem[keep]
    lib = pickle.load(open(a.lib, 'rb'))[:a.top]
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    import run_ialns as RI; say = RI.logger(out / 'log.txt')
    _S.update(ids=ids, elem=elem, lib=lib)

    # global switch-cost matrix over the library -> neighbourhoods
    R = []; V = []
    for c in lib:
        rE, vE = eph.earth_state(c['tL']); R.append(rE); V.append(np.array(vE) + np.array(c['vinf']))
    el = elements(np.array(R), np.array(V))
    SWall = K_TOTAL * (K_PLANE * np.linalg.norm(el['i_vec'][:, None, :] - el['i_vec'][None, :, :], axis=2) +
                       a.kecc * K_ECC * np.linalg.norm(el['e_vec'][:, None, :] - el['e_vec'][None, :, :], axis=2))
    vals = np.array([c['val'] for c in lib])
    # seeds: best value first, but spread out in the switch metric (the library holds near-duplicate carriers)
    seeds = []
    for s in np.argsort(-vals):
        if len(seeds) >= a.seeds: break
        if all(SWall[s, q] > a.sdiv for q in seeds): seeds.append(int(s))
    pools = {}
    for s in seeds:
        nb = np.argsort(SWall[s])[:1 if a.fixed else a.knb]
        pools[int(s)] = [int(x) for x in nb]
    say(f'{len(lib)} carriers; switch cost median {np.median(SWall[np.triu_indices(len(lib),1)]):.3f} km/s, '
        f'neighbourhood {a.knb} median {np.median([SWall[s][pools[int(s)]].max() for s in seeds]):.3f} km/s')

    tic = time.time()
    with mp.get_context('fork').Pool(a.nproc) as pl:
        built = pl.map(build_pool, [(int(s), pools[int(s)], a) for s in seeds])
    _S['pools'] = {s: p for s, p, _, _ in built}
    _S['evs'] = {s: e for s, _, e, _ in built}
    _S['sws'] = {s: w for s, _, _, w in built}
    ne = [0 if e is None else len(e['t']) for e in _S['evs'].values()]
    say(f'pools built in {time.time()-tic:.0f} s: {np.median(ne):.0f} events each (median)')

    # --rounds > 1 builds a POOL for the set cover: after each round the targets just taken are priced at zero, so the
    # next round's routes are forced onto the rest of the catalogue.  How slowly the flyby count decays over the rounds
    # is exactly the FREEDOM a partition needs (a fixed carrier: 19.5 10.5 9.0 6.5 ...; a 24-carrier pool: 24.5 20.5
    # 13.0 11.5 ..., i.e. +62% of total supply over 8 rounds).
    prize = np.ones(301)
    routes = []
    tic = time.time()
    with mp.get_context('fork').Pool(a.nproc) as pl:
        for rd in range(a.rounds):
            res = pl.map(w_route, [(int(s), prize, a) for s in seeds])
            new = [r for _, r in res if r and r['n'] >= a.nmin]
            new.sort(key=lambda r: -r['n'])
            for r in new:
                r['round'] = rd
                if a.rounds > 1:                       # greedy: one route per seed per round, then re-price
                    for t, _ in r['legs']: prize[t] = 0.0
            routes += new
            say(f'round {rd}: {len(new)} routes, flybys ' + ' '.join(str(r['n']) for r in new) +
                f'; priced out {int(301 - prize.sum())} targets ({time.time()-tic:.0f} s)')
            if not new: break
    routes.sort(key=lambda r: -r['n'])
    for r in routes[:24]:
        say(f'  seed {r["seed"]:3d} rd {r.get("round",0)}: {r["n"]:2d} flybys over {r["ncar"]} carriers ({r["nswitch"]} switches), '
            f'model dv {r["dv_phase"]:.2f} phase + {r["dv_gap"]:.2f} gap + {r["dv_switch"]:.2f} switch km/s, '
            f'value {r["val"]:.2f}')
    cov = set(t for r in routes for t, _ in r['legs'])
    say(f'{len(routes)} routes, union {len(cov)} of 298 targets')
    json.dump(routes, open(out / 'routes.json', 'w'), indent=1)
    say(f'-> {out}/routes.json')


if __name__ == '__main__':
    main()
