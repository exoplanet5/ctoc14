"""Beam search on the precomputed transfer graph (stage 9).

State = the last EDGE (o->p), because a flyby fixes no velocity: the junction impulse at p is
|v_dep(p->q) - v_arr(o->p)|, and both vectors are already in the table (tools/transfer_db.py).  So expanding a
partial route is a gather and a subtraction -- no Lambert, no propagation, no ephemeris.  Everything the old beam
(ctoc14/search.py) spent its time on is precomputed, which is what turns route generation from minutes into
milliseconds and lets the compute go into the fleet combinatorics instead.

Launch: source labels are Earth at a grid of epochs; the first leg is free up to |v_inf| <= 4 km/s.
Covered targets are kept as zero-gain STEPPING STONES (the 2026-09-16 waypoint lesson), never excluded.

Usage: graph_beam.py tdb.npz [--beam 3000] [--lam 1.2] [--depth 60] [--exclude 1,2,...]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
from ctoc14.constants import DAY, VE, VINF_MAX, T_MISSION, M_DRY
from ctoc14.kepler import Ephemeris

NW = 5                                        # 298 target bits in 5 uint64 words


class Graph:
    def __init__(self, path):
        z = np.load(path)
        self.src = z['src']; self.dst = z['dst']; self.v1 = z['v1'].astype(np.float64)
        self.v2 = z['v2'].astype(np.float64); self.ptr = z['ptr']; self.A = z['A']
        self.ids = z['ids']
        self.t = self.A[:, 0]; self.tgt = self.A[:, 1].astype(np.int32); self.R = self.A[:, 4:7]
        self.n = len(self.A)

    def sources(self, eph, t0=0.0, t1=900.0, step=5.0, dtmax=500.0, vinf=VINF_MAX):
        """Earth -> first encounter: returns (edge-like) arrays of arrival node, arrival velocity, cost."""
        from ctoc14.lambert import lambert_all
        tl = np.arange(t0, t1, step) * DAY
        RE = np.array([eph.earth_state(x)[0] for x in tl]); VE_ = np.array([eph.earth_state(x)[1] for x in tl])
        I = []; J = []
        for k, t in enumerate(tl):
            lo = np.searchsorted(self.t, t + 15 * DAY); hi = np.searchsorted(self.t, t + dtmax * DAY)
            I.append(np.full(hi - lo, k)); J.append(np.arange(lo, hi))
        I = np.concatenate(I); J = np.concatenate(J)
        v1, v2, _, _ = lambert_all(RE[I], self.R[J], self.t[J] - tl[I], nrev_max=1)
        dv = np.linalg.norm(v1 - VE_[I][None], axis=-1)
        dv = np.where(np.isfinite(dv), dv, np.inf)
        k = np.argmin(dv, axis=0); g = np.arange(len(I))
        cost = np.maximum(0.0, dv[k, g] - vinf)
        ok = np.isfinite(cost) & (cost < 0.5)
        return dict(tl=tl[I[ok]], node=J[ok], varr=v2[k, g][ok], cost=cost[ok], vinf=v1[k, g][ok] - VE_[I[ok]])


def beam(G, S, a, prize=None):
    """returns the best routes found: list of dict(nodes, targets, dv, tL)"""
    tgt = G.tgt; ptr = G.ptr; dst = G.dst; v1 = G.v1; v2 = G.v2
    pz = np.ones(301) if prize is None else prize
    # level 0 labels: one per source arc
    m = np.argsort(S['cost'])[:a.beam * 4]
    node = S['node'][m].astype(np.int64); varr = S['varr'][m]; dv = S['cost'][m].copy()
    cnt = pz[tgt[node]].copy()
    mask = np.zeros((len(node), NW), np.uint64)
    for i, nd in enumerate(node):
        t = int(tgt[nd]); mask[i, t >> 6] |= np.uint64(1) << np.uint64(t & 63)
    par = [(-1, np.arange(len(node)), node.copy(), S['tl'][m].copy())]
    best = None; hist = []
    for lev in range(1, a.depth):
        deg = (ptr[node + 1] - ptr[node]).astype(np.int64)
        if deg.sum() == 0: break
        lab = np.repeat(np.arange(len(node)), deg)
        off = np.concatenate([np.arange(ptr[p], ptr[p + 1]) for p in node]) if len(node) else np.zeros(0, int)
        q = dst[off]
        step = np.linalg.norm(v1[off] - varr[lab], axis=1)
        keep = step < a.dvmax
        lab, off, q, step = lab[keep], off[keep], q[keep], step[keep]
        if len(lab) == 0: break
        ndv = dv[lab] + step
        tq = tgt[q].astype(np.uint64)
        bit = np.uint64(1) << (tq & np.uint64(63))
        seen = (mask[lab, (tq >> np.uint64(6)).astype(np.intp)] & bit) != 0
        gain = np.where(seen, 0.0, pz[tq.astype(np.intp)])
        ncnt = cnt[lab] + gain
        tank = M_DRY * 1.0025 * np.exp(ndv / VE)
        alive = tank < a.max_tank
        lab, off, q, ndv, ncnt, seen, tq, bit = (x[alive] for x in (lab, off, q, ndv, ncnt, seen, tq, bit))
        if len(lab) == 0: break
        score = ncnt - a.lam * ndv
        # keep the best a.beam labels, at most --per-node per arrival node
        o = np.argsort(-score)
        o = o[:a.beam * 12]
        seenc = {}; sel = []
        capby = tgt if getattr(a, 'cap_by', 'target') == 'target' else np.arange(len(tgt))
        for i in o:
            k = int(capby[q[i]])
            c = seenc.get(k, 0)
            if c >= a.per_node: continue
            seenc[k] = c + 1; sel.append(i)
            if len(sel) >= a.beam: break
        sel = np.array(sel, int)
        lab, off, q, ndv, ncnt, tq = (x[sel] for x in (lab, off, q, ndv, ncnt, tq))
        nmask = mask[lab].copy()
        nmask[np.arange(len(tq)), (tq >> np.uint64(6)).astype(np.intp)] |= np.uint64(1) << (tq & np.uint64(63))
        par.append((lev, lab, q.copy(), None))
        node, varr, dv, cnt, mask = q, v2[off], ndv, ncnt, nmask
        j = int(np.argmax(cnt - a.lam2 * dv))
        hist.append((lev + 1, float(cnt[j]), float(dv[j]), float(M_DRY * 1.0025 * np.exp(dv[j] / VE))))
        if best is None or cnt[j] - a.lam2 * dv[j] > best[0]:
            best = (cnt[j] - a.lam2 * dv[j], lev, j)
    return par, hist, best, (node, dv, cnt)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('tdb')
    ap.add_argument('--beam', type=int, default=3000); ap.add_argument('--lam', type=float, default=1.2)
    ap.add_argument('--lam2', type=float, default=0.35); ap.add_argument('--depth', type=int, default=60)
    ap.add_argument('--dvmax', type=float, default=1.2); ap.add_argument('--max-tank', type=float, default=1500.0)
    ap.add_argument('--per-node', type=int, default=3); ap.add_argument('--tl0', type=float, default=0.0)
    ap.add_argument('--tl1', type=float, default=900.0)
    ap.add_argument('--cap-by', default='target', choices=['target', 'node'])
    a = ap.parse_args()
    G = Graph(a.tdb); eph = Ephemeris()
    tic = time.time()
    S = G.sources(eph, a.tl0, a.tl1)
    print(f'{len(S["node"])} launch arcs (cost < 0.5 km/s beyond v_inf)  ({time.time()-tic:.0f} s)')
    tic = time.time()
    par, hist, best, fin = beam(G, S, a)
    print(f'beam {a.beam}, depth {a.depth}: {time.time()-tic:.1f} s')
    print('  level  flybys   dv     tank')
    for k, c, d, t in hist[::4]: print(f'  {k:5d}  {c:5.0f}  {d:6.2f}  {t:6.0f} kg')
    n, dv, cnt = fin
    J = lambda tk: 1 + (tk - 600) / 1400 + ((tk - 600) / 1400) ** 2
    tank = M_DRY * 1.0025 * np.exp(dv / VE)
    o = np.argsort([J(t) / max(c, 1) for t, c in zip(tank, cnt)])[:6]
    print('\n  best final labels by J per flyby:')
    for i in o:
        print(f'    {cnt[i]:.0f} flybys, dv {dv[i]:.2f} km/s, tank {tank[i]:.0f} kg, J_i {J(tank[i]):.3f}, '
              f'{J(tank[i])/max(cnt[i],1):.4f} J/fb')
    print('  reference: t10d best route 0.0346, our best carrier route 0.0339, HIT fleet ~0.0353')


if __name__ == '__main__':
    main()
