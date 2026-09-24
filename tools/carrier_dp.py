"""Carrier-orbit route generator: exact DP on the (time, lag) event plane, then a real twin settle (stage 7 test).

Model.  A craft stays on one launch-bought carrier orbit and only steers its phase.  For every asteroid whose orbit
passes within --eps of the carrier TRACK (radial gap at the mutual node, tools/exp_carrier.py) each node passage of
the asteroid is an event (t, lag): the craft must reach the node `lag` days later than its ballistic schedule.  A
route is a piecewise-linear lag curve l(t) from (t_launch, 0) through events; the slope s = dl/dt is a period offset
and changing it costs  v/3 |ds| = 9.93 km/s per unit (= phasemodel.K_DRIFT 0.028 km/s per deg/yr).  Each event also
pays K_GAP |gap| for the residual track gap.

Unlike ctoc14/phasedp.py this is NOT a grid: nodes are the events, arcs are event pairs, and the DP runs on the line
graph (state = incoming arc), so the slope is exact and the eccentricity/plane slack is granted to nobody -- the
carrier's e- and i-vectors are fixed by construction.  Value = flybys - lam * dv.

The DP is only a GENERATOR.  The verdict comes from --settle: the sequence goes through impulsive.from_tour and
run_ialns.settle exactly like a planner tour, and the twin tank is the number that counts.

Usage: carrier_dp.py [--n 50000] [--eps 0.01] [--top 5] [--lam 0.5] [--settle]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
from ctoc14.constants import MU, AU, DAY, VINF_MAX, T_MISSION, T_OFF_AST, VE
from ctoc14.kepler import load_mea, Ephemeris
from ctoc14.search import UNREACHABLE
from exp_carrier import frame, rv2frame, gaps

K_SLOPE = 29.78 / 3.0          # km/s per unit of lag slope (day/day); scaled by --kphase for the DP COST only
K_PHYS = 29.78 / 3.0           # the physical value, used for the seed impulses


def f2M(f, e):
    E = 2 * np.arctan2(np.sqrt(1 - e) * np.sin(f / 2), np.sqrt(1 + e) * np.cos(f / 2))
    return E - e * np.sin(E)


def anomaly_of(u, n, Evec):
    """true anomaly of unit direction(s) u in the orbit with normal n and eccentricity vector Evec."""
    e = np.linalg.norm(Evec, axis=-1)
    P = Evec / e[..., None]; Q = np.cross(n, P)
    return np.arctan2(np.einsum('...k,...k->...', Q, u), np.einsum('...k,...k->...', P, u)), e


def events_for(car, ids, elem, eps, images=1):
    """car: dict(tL, r, v).  Returns arrays (t [s], lag [d], gap [AU], ast id)."""
    pc, nc, Ec, ac = rv2frame(car['r'][None], car['v'][None])
    pa, na, Ea = frame(elem[:, :5])
    g, _ = gaps(pc, nc, Ec, pa, na, Ea); g = g[0]                         # (nast, 2)
    Pc = 2 * np.pi * np.sqrt((ac[0] * AU) ** 3 / MU)
    f0, ec = anomaly_of(car['r'] / np.linalg.norm(car['r']), nc[0], Ec[0]); M0c = f2M(f0, ec)
    nm_a = np.sqrt(MU / (elem[:, 0] * AU) ** 3); Pa = 2 * np.pi / nm_a
    ev = []
    for j in range(len(ids)):
        for k, sg in enumerate((1.0, -1.0)):
            if abs(g[j, k]) >= eps: continue
            u = np.cross(nc[0], na[j]); u = sg * u / np.linalg.norm(u)
            fa, ea = anomaly_of(u, na[j], Ea[j]); Mu = f2M(fa, ea)
            t1 = ((Mu - np.radians(elem[j, 5])) / nm_a[j] - T_OFF_AST) % Pa[j]
            ta = np.arange(t1, T_MISSION - 5 * DAY, Pa[j]); ta = ta[ta > car['tL'] + 20 * DAY]
            fc, _ = anomaly_of(u, nc[0], Ec[0])
            tau = car['tL'] + ((f2M(fc, ec) - M0c) / (2 * np.pi) * Pc) % Pc   # first ballistic passage of the node
            for t in ta:
                d = ((t - tau + Pc / 2) % Pc) - Pc / 2
                for im in range(-images, images + 1):
                    ev.append((t, (d + im * Pc) / DAY, g[j, k], int(ids[j])))
    ev.sort()
    E = np.array(ev) if ev else np.zeros((0, 4))
    return E[:, 0], E[:, 1], E[:, 2], E[:, 3].astype(int), Pc


def dp_route(t, lag, gap, ast, tL, lam, kgap, smax=0.08, s0max=0.05, dtmin=15.0, dtmax=900.0, ac=1.0, couple=True, w=None):
    """Line-graph DP.  Returns (value, list of event indices).

    `w`: optional per-target weight (array indexed by asteroid id).  A flyby of target x is worth w[x] instead of 1,
    which is how the joint k-route design prices targets that another route is already taking (tools/carrier_joint.py).
    """
    n = len(t); td = t / DAY
    gain = np.ones(n) if w is None else np.asarray(w, float)[ast]
    # a lag slope s is a period offset, i.e. the whole track sits (2/3) s a further out: the gap an event pays is the
    # one left at the slope the craft ARRIVES with (couple=False: the static carrier gap)
    geff = (lambda j, s: np.abs(gap[j] - (2.0 / 3.0) * ac * s)) if couple else (lambda j, s: np.abs(gap[j]))
    # arcs out of each node
    succ = []
    for i in range(n):
        j = np.nonzero((td > td[i] + dtmin) & (td < td[i] + dtmax) & (ast != ast[i]))[0]
        s = (lag[j] - lag[i]) / (td[j] - td[i]); ok = np.abs(s) <= smax
        succ.append((j[ok], s[ok]))
    # incoming labels per node: lists of (value, slope, pred node, pred label idx)
    lab = [[] for _ in range(n)]
    s0 = lag / np.maximum(td - tL / DAY, 1.0)
    for i in range(n):
        if abs(s0[i]) <= s0max and td[i] - tL / DAY < 2 * dtmax:
            lab[i].append((gain[i] - lam * kgap * geff(i, s0[i]), s0[i], -1, -1))
    best = (-1e9, -1, -1)
    for i in range(n):                                   # events are time sorted
        if not lab[i]: continue
        L = lab[i]
        if len(L) > 60:                                  # keep the best labels per node (slope-diverse)
            L.sort(key=lambda x: -x[0]); L = L[:60]; lab[i] = L
        V = np.array([x[0] for x in L]); S = np.array([x[1] for x in L])
        k = int(V.argmax())
        if V[k] > best[0]: best = (V[k], i, k)
        J, SJ = succ[i]
        if len(J) == 0: continue
        val = V[:, None] - lam * K_SLOPE * np.abs(SJ[None, :] - S[:, None]) + (gain[J] - lam * kgap * geff(J, SJ))[None, :]
        kk = val.argmax(0)
        for c, j in enumerate(J):
            lab[j].append((val[kk[c], c], SJ[c], i, int(kk[c])))
    seq = []; _, i, k = best
    while i >= 0:
        seq.append(i); _, _, pi, pk = lab[i][k]; i, k = pi, pk
    return best[0], seq[::-1]


def seed_ip(eph, tL, vinf, t_ev, lag_ev, dtau=20 * DAY):
    """Impulsive problem WITHOUT flybys: the ballistic carrier plus the DP's phasing impulses (tangential, 1 d after each
    event, v/3 per unit of slope change; the first slope goes into v_inf).  Flybys are added by insert_block."""
    from ctoc14.kepler import propagate_twobody
    from ctoc14.impulsive import ImpulsiveProblem
    td = np.concatenate([[tL], t_ev]) / DAY; lg = np.concatenate([[0.0], lag_ev])
    S = np.diff(lg) / np.diff(td)
    rE, vE = eph.earth_state(tL); v = vE + vinf
    v = v + K_PHYS * S[0] * v / np.linalg.norm(v)
    vi = v - vE; nv = np.linalg.norm(vi)
    if nv > VINF_MAX: vi *= VINF_MAX / nv
    r = np.array(rE, float); v = vE + vi; t = tL; imp_t = []; imp_dv = []
    for k in range(len(t_ev) - 1):
        tk = t_ev[k] + 1.0 * DAY
        r, v = propagate_twobody(r, v, tk - t); t = tk
        dv = K_PHYS * (S[k + 1] - S[k]) * v / np.linalg.norm(v)
        imp_t.append(tk); imp_dv.append(dv); v = v + dv
    edges = np.arange(tL, t_ev.max() + dtau, dtau); ts = 0.5 * (edges[:-1] + edges[1:]); Ts = np.zeros((len(ts), 3))
    ts = np.concatenate([ts, np.array(imp_t)]); Ts = np.concatenate([Ts, np.array(imp_dv).reshape(-1, 3)])
    o = np.argsort(ts); return ImpulsiveProblem(eph, tL, vi, ts[o], Ts[o], np.zeros(0), [], dtau=dtau)


_G = {}


def _value(c):
    g = _G; a = g['a']
    car = dict(tL=float(g['tLs'][c]), r=g['r'][c], v=g['v'][c])
    t, lag, gap, ast, Pc = events_for(car, g['ids'], g['elem'], a.eps)
    if len(t) < 3: return (-1e9, 0)
    ac_ = float(rv2frame(car['r'][None], car['v'][None])[3][0])
    val, seq = dp_route(t, lag, gap, ast, car['tL'], a.lam, a.kgap, ac=ac_, couple=not a.static_gap)
    return (float(val), len(set(int(ast[i]) for i in seq)))


def main():
    global K_SLOPE
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=50000); ap.add_argument('--eps', type=float, default=0.01)
    ap.add_argument('--top', type=int, default=5); ap.add_argument('--lam', type=float, default=0.5)
    ap.add_argument('--kgap', type=float, default=15.0, help='km/s per AU of residual track gap')
    ap.add_argument('--seed', type=int, default=0); ap.add_argument('--settle', action='store_true')
    ap.add_argument('--iters', type=int, default=100); ap.add_argument('--exclude', default='')
    ap.add_argument('--out', default=''); ap.add_argument('--lambert-seed', action='store_true')
    ap.add_argument('--static-gap', action='store_true'); ap.add_argument('--fleet-out', default='')
    ap.add_argument('--pool', type=int, default=0); ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--kphase', type=float, default=1.0, help='multiplier on the phase cost (twin calibration 1.5)')
    a = ap.parse_args(); K_SLOPE = K_SLOPE * a.kphase
    rng = np.random.default_rng(a.seed); eph = Ephemeris()
    ids, elem = load_mea()
    excl = set(UNREACHABLE) | set(int(x) for x in a.exclude.split(',') if x)
    keep = np.array([int(x) not in excl for x in ids]); ids = ids[keep]; elem = elem[keep]
    pa, na, Ea = frame(elem[:, :5])

    # sample launch-reachable carriers (launch in the first year), rank by geometric candidate count
    tLs = rng.uniform(0, 365.25 * DAY, a.n)
    d = rng.normal(size=(a.n, 3)); d /= np.linalg.norm(d, axis=1)[:, None]
    vinf = d * (VINF_MAX * rng.uniform(0, 1, a.n) ** (1 / 3))[:, None]
    RV = np.array([np.concatenate(eph.earth_state(x)) for x in tLs]); r = RV[:, :3]; v = RV[:, 3:] + vinf
    cnt = np.zeros(a.n, int)
    for c0 in range(0, a.n, 10000):
        sl = slice(c0, c0 + 10000)
        pc, nc, Ec, ac = rv2frame(r[sl], v[sl]); g, _ = gaps(pc, nc, Ec, pa, na, Ea)
        cnt[sl] = (np.abs(g).min(2) < a.eps).sum(1)
    order = np.argsort(-cnt)[:a.top]
    print(f'{a.n} carriers, eps {a.eps}: candidate count best {cnt.max()}, median {np.median(cnt):.0f}')

    if a.pool > a.top:
        # rank a pool of carriers by the DP VALUE (timing included), not by the geometric candidate count
        import multiprocessing as mp
        cand = np.argsort(-cnt)[:a.pool]
        global _G; _G = dict(tLs=tLs, r=r, v=v, ids=ids, elem=elem, a=a)
        tic = time.time()
        with mp.get_context('fork').Pool(a.nproc) as pl:
            vals = pl.map(_value, [int(c) for c in cand], chunksize=8)
        vals = np.array(vals); order = cand[np.argsort(-vals[:, 0])[:a.top]]
        print(f'DP value over {a.pool} carriers ({time.time()-tic:.0f} s): flybys best {vals[:,1].max():.0f}, '
              f'median {np.median(vals[:,1]):.0f}; value best {vals[:,0].max():.2f}, median {np.median(vals[:,0]):.2f}')
    if a.settle:
        import run_ialns as RI
        FL = RI.IFleet()
    res = []
    for c in order:
        car = dict(tL=float(tLs[c]), r=r[c], v=v[c])
        t, lag, gap, ast, Pc = events_for(car, ids, elem, a.eps)
        tic = time.time()
        ac_ = float(rv2frame(car['r'][None], car['v'][None])[3][0])
        val, seq = dp_route(t, lag, gap, ast, car['tL'], a.lam, a.kgap, ac=ac_, couple=not a.static_gap)
        seen = set(); seq = [i for i in seq if not (ast[i] in seen or seen.add(ast[i]))]
        S = np.diff(np.concatenate([[0.0], lag[seq]])) / np.diff(np.concatenate([[car['tL'] / DAY], t[seq] / DAY]))
        dv_s = K_SLOPE * np.abs(np.diff(S)).sum()
        dv_g = a.kgap * np.abs(gap[seq] - (0.0 if a.static_gap else (2.0 / 3.0) * ac_ * S)).sum()
        print(f'\ncarrier tL {car["tL"]/DAY:5.0f} d |vinf| {np.linalg.norm(vinf[c]):.2f}  P {Pc/DAY:.1f} d: {cnt[c]} candidates,'
              f' {len(t)} events -> DP {len(seq)} flybys, model dv {dv_s:.2f} (phase) + {dv_g:.2f} (gap) km/s'
              f'  ({time.time()-tic:.1f} s)')
        print('   lag [d]: ' + ' '.join(f'{x:.0f}' for x in lag[seq]))
        tour = dict(t_launch=car['tL'], vinf=[float(x) for x in vinf[c]],
                    legs=[dict(ast=int(ast[i]), t_flyby=float(t[i])) for i in seq])
        rec = dict(tour=tour, n=len(seq), dv_model=float(dv_s + dv_g), dv_phase=float(dv_s), gap_sum=float(dv_g / a.kgap), eps=a.eps, lam=a.lam)
        if a.settle and len(seq) >= 3:
            from ctoc14.impulsive import settle_tour
            tic = time.time()
            ip, miss, lg = settle_tour(eph, tour, lambda ip: RI.settle(ip, a.iters)) if a.lambert_seed else (None, np.inf, None)
            if ip is None:
                # carrier seed: ballistic carrier + DP phasing impulses, all flybys brought in by one joint homotopy
                from ctoc14.globalopt import insert_block
                ip = seed_ip(eph, car['tL'], vinf[c], t[seq], lag[seq])
                d = insert_block(ip, [(int(ast[i]), float(t[i])) for i in seq], log=RI.QUIET)
                miss = float(np.max(d))
                if miss < 1e4: miss = RI.settle(ip, a.iters)
                if miss > 150.0: ip = None
            if ip is None:
                print(f'   twin: NOT settled (best miss {miss:.0f} km)  ({time.time()-tic:.0f} s)')
            else:
                tank = float(ip.tank()); dv = VE * np.log(tank / 601.5)
                print(f'   twin: {len(seq)} flybys, tank {tank:.0f} kg, dv {dv:.2f} km/s = {dv/len(seq):.3f} km/s per flyby,'
                      f' J_i {1 + (tank-600)/1400 + ((tank-600)/1400)**2:.3f}  ({time.time()-tic:.0f} s)')
                rec.update(tank=tank)
                if a.fleet_out: FL.routes[f'c{len(FL.routes):02d}'] = dict(st=RI.ist(ip), tank=tank)
        res.append(rec)
    if a.out: json.dump(res, open(a.out, 'w'), indent=1)
    if a.settle and a.fleet_out and FL.routes: FL.save(a.fleet_out, note='carrier routes (may overlap)')


if __name__ == '__main__':
    main()
