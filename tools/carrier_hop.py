"""Multi-carrier routes: design the HOP into the route (stage 8, the leaders' missing ingredient).

Arithmetic that forces this.  A fleet covering 298 targets needs every target at ~12 kg of tank (t10d achieves 11,
HIT 11.8).  A carrier's OWN targets cost ~13 kg -- competitive -- but one carrier only ever has ~20 of them, and every
target past those costs 45-90 kg (stage 7 section 10.4).  The leaders fly 37-43 per craft.  37-43 is two or three
carriers' worth, so their craft must change carrier mid-mission.

Inserting a second carrier's targets into a finished route does not work (0 of 24 block homotopies settled: the route
is pinned by its own flybys).  But a hop needs no insertion machinery at all -- it is structurally identical to the
launch.  At any time t_h the craft is at r(t_h) with velocity v; adding Dv there puts it on a NEW carrier defined by
(t_h, r, v + Dv), whose events and lag structure follow from that state exactly as the launch carrier's follow from
(tL, r_E, v_E + vinf).  So a two-segment route is 8 free numbers (tL, vinf, t_h, Dv) and three segments is 12, scored
by the same exact DP per segment and searched by the same pattern search.

Usage: carrier_hop.py out_dir [--segs 2] [--starts 64] [--settle]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pickle, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import carrier_dp as CD
from ctoc14.constants import DAY, VINF_MAX, T_MISSION, VE
from ctoc14.kepler import load_mea, Ephemeris, propagate_twobody
from ctoc14.search import UNREACHABLE
from exp_carrier import rv2frame
_S = {}


def segments(x, nseg):
    """x = [tL_d, vinf(3), (t_hop_d, dv(3)) x (nseg-1)] -> list of (t_start, r, v, dv_cost)."""
    eph = _S['eph']
    tL = float(np.clip(x[0], 0.0, 700.0)) * DAY
    vi = np.array(x[1:4], float); nv = np.linalg.norm(vi)
    if nv > VINF_MAX: vi = vi * (VINF_MAX / nv)
    rE, vE = eph.earth_state(tL)
    segs = [[tL, np.array(rE, float), np.array(vE, float) + vi, 0.0]]
    t = tL; r = segs[0][1]; v = segs[0][2]
    for k in range(nseg - 1):
        th = float(np.clip(x[4 + 4 * k], (t + 200 * DAY) / DAY, (T_MISSION - 400 * DAY) / DAY)) * DAY
        if th <= t + 100 * DAY: return None
        r, v = propagate_twobody(r, v, th - t); t = th
        dv = np.array(x[5 + 4 * k:8 + 4 * k], float)
        n = np.linalg.norm(dv)
        if n > _S['a'].dvmax: dv = dv * (_S['a'].dvmax / n); n = _S['a'].dvmax
        v = v + dv
        segs.append([t, np.array(r, float), np.array(v, float), float(n)])
    return segs


def route_of(x, nseg, w=None):
    """DP each segment on its own carrier over its own time window; value = flybys - lam (phase dv + hop dv)."""
    a = _S['a']; ids = _S['ids']; elem = _S['elem']
    segs = segments(x, nseg)
    if segs is None: return None
    out = []; val = 0.0; taken = set()
    for k, (t0, r, v, dvc) in enumerate(segs):
        t1 = segs[k + 1][0] if k + 1 < len(segs) else T_MISSION
        keep = np.array([int(t) not in taken for t in ids])
        if keep.sum() < 3: continue
        try:
            t, lag, gap, ast, Pc = CD.events_for(dict(tL=t0, r=r, v=v), ids[keep], elem[keep], a.eps)
        except Exception:
            return None
        m = (t > t0 + 20 * DAY) & (t < t1 - 5 * DAY)
        if m.sum() < 3: continue
        t, lag, gap, ast = t[m], lag[m], gap[m], ast[m]
        ac_ = float(rv2frame(r[None], v[None])[3][0])
        try:
            f, seq = CD.dp_route(t, lag, gap, ast, t0, a.lam, a.kgap, ac=ac_, w=w)
        except Exception:
            return None
        seen = set(); seq = [i for i in seq if not (ast[i] in seen or seen.add(ast[i]))]
        seq = [i for i in seq if int(ast[i]) not in taken]
        if not seq: continue
        taken |= set(int(ast[i]) for i in seq)
        val += f - a.lam * a.kphase_hop * dvc
        out.append(dict(t0=float(t0), dv=float(dvc), r=r.tolist(), v=v.tolist(),
                        base=[(int(ast[i]), float(t[i]), float(lag[i])) for i in seq]))
    if not out: return None
    return dict(val=float(val), nfb=len(taken), segs=out, x=[float(q) for q in x])


def w_opt(job):
    x0, nseg, seed = job; a = _S['a']
    rng = np.random.default_rng(seed)
    x = np.array(x0, float); r0 = route_of(x, nseg); f = r0['val'] if r0 else -1e9
    step = np.array([30.0, 0.6, 0.6, 0.6] + [200.0, 0.4, 0.4, 0.4] * (nseg - 1))
    for it in range(a.iters):
        best = None
        for y in x + step * rng.normal(size=(a.batch, len(x))):
            r = route_of(y, nseg)
            if r and (best is None or r['val'] > best[0]): best = (r['val'], y, r)
        if best and best[0] > f: f, x, r0 = best[0], best[1], best[2]; step *= 1.15
        else: step *= 0.7
        if step[0] < 0.5: break
    return r0


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--segs', type=int, default=2)
    ap.add_argument('--starts', type=int, default=64); ap.add_argument('--iters', type=int, default=150)
    ap.add_argument('--batch', type=int, default=8); ap.add_argument('--n', type=int, default=20000)
    ap.add_argument('--eps', type=float, default=0.012); ap.add_argument('--lam', type=float, default=1.0)
    ap.add_argument('--kgap', type=float, default=10.0); ap.add_argument('--kphase', type=float, default=1.5)
    ap.add_argument('--kphase-hop', type=float, default=1.0, help='J per km/s of hop dv in the DP objective')
    ap.add_argument('--dvmax', type=float, default=3.0); ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--seed', type=int, default=0); ap.add_argument('--settle', action='store_true')
    a = ap.parse_args()
    CD.K_SLOPE = CD.K_PHYS * a.kphase
    eph = Ephemeris(); ids, elem = load_mea()
    keep = np.array([int(x) not in UNREACHABLE for x in ids]); ids = ids[keep]; elem = elem[keep]
    _S.update(a=a, eph=eph, ids=ids, elem=elem)
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    import run_ialns as RI; say = RI.logger(out / 'log.txt')
    rng = np.random.default_rng(a.seed)
    nseg = a.segs
    X0 = []
    for _ in range(a.n):
        x = [rng.uniform(0, 700)] + list(rng.normal(size=3) * 1.2)
        t = rng.uniform(700, 4200)
        for k in range(nseg - 1):
            x += [t] + list(rng.normal(size=3) * 0.8); t = min(t + rng.uniform(600, 2000), 4800)
        X0.append(x)
    tic = time.time()
    with mp.get_context('fork').Pool(a.nproc) as pl:
        vals = np.array(pl.map(_score, [(X0[i::a.nproc], nseg) for i in range(a.nproc)], chunksize=1), dtype=object)
        sc = np.concatenate([np.asarray(v) for v in vals]); idx = np.concatenate([np.arange(a.n)[i::a.nproc] for i in range(a.nproc)])
        order = idx[np.argsort(-sc)][:a.starts]
        say(f'{a.n} random {nseg}-segment starts scored in {time.time()-tic:.0f} s: best {sc.max():.2f}')
        res = [r for r in pl.map(w_opt, [(X0[j], nseg, a.seed * 131 + i) for i, j in enumerate(order)]) if r]
    res.sort(key=lambda r: -r['val'])
    say(f'optimised {len(res)} routes in {time.time()-tic:.0f} s')
    for r in res[:8]:
        say('  value %.2f: %d flybys, segments ' % (r['val'], r['nfb'])
            + ' + '.join(f'{len(s["base"])}@{s["t0"]/DAY:.0f}d' + (f' (hop {s["dv"]:.2f} km/s)' if s['dv'] else '')
                         for s in r['segs']))
    json.dump(res[:64], open(out / 'routes.json', 'w'), indent=1)
    say(f'-> {out}/routes.json')


def _score(job):
    X, nseg = job
    o = []
    for x in X:
        r = route_of(x, nseg)
        o.append(r['val'] if r else -1e9)
    return o


if __name__ == '__main__':
    main()
