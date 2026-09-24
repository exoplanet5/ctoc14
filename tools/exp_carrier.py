"""Carrier-orbit geometry screen (stage 7 idea test).

Hypothesis: the leaders' regime (N=10 at ~70 kg of propellant per craft, 0.14 km/s per flyby against our 0.49) is not
a better search over the same routes but a different kind of route: the craft sits on ONE carrier orbit that launch
v_inf buys for free, never changes its plane or eccentricity vector, and only steers its PHASE (0.028 km/s per deg/yr
of drift change, phasemodel.K_DRIFT).  An asteroid is then cheap for that carrier iff the two ORBITS nearly intersect,
which is pure geometry: at the mutual node the radial gap  r_ast(u) - r_car(u)  must be ~0.

This script measures how many asteroids a launch-reachable carrier can have within eps, and how many carriers cover
the catalogue (greedy set cover).  No dynamics, no time: if this count is small the idea is dead in a minute.

Usage: exp_carrier.py [--n 200000] [--eps 0.003] [--k 10] [--seed 0] [--out carriers.json]
"""
import sys, json, pathlib, argparse
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import numpy as np
from ctoc14.constants import MU, AU, VINF_MAX, EARTH_ELEMENTS
from ctoc14.kepler import load_mea
from ctoc14.search import UNREACHABLE


def frame(el):
    """elements (a AU, e, i, Om, om [deg]) rows -> p [AU], unit normal n, eccentricity vector E (3-D)."""
    a, e = el[:, 0], el[:, 1]
    i, Om, om = np.radians(el[:, 2]), np.radians(el[:, 3]), np.radians(el[:, 4])
    cO, sO, ci, si, co, so = np.cos(Om), np.sin(Om), np.cos(i), np.sin(i), np.cos(om), np.sin(om)
    P = np.stack([cO * co - sO * so * ci, sO * co + cO * so * ci, so * si], 1)
    n = np.stack([sO * si, -cO * si, ci], 1)
    return a * (1 - e * e), n, P * e[:, None]


def rv2frame(r, v):
    """states [km, km/s] rows -> p [AU], n, E, a [AU]."""
    h = np.cross(r, v); hn = np.linalg.norm(h, axis=1)
    rn = np.linalg.norm(r, axis=1)
    E = np.cross(v, h) / MU - r / rn[:, None]
    a = 1.0 / (2.0 / rn - np.einsum('ij,ij->i', v, v) / MU)
    return hn * hn / MU / AU, h / hn[:, None], E, a / AU


def gaps(pc, nc, Ec, pa, na, Ea):
    """radial gap [AU] at both mutual nodes: (ncar, nast, 2)."""
    u = np.cross(nc[:, None, :], na[None, :, :])
    s = np.linalg.norm(u, axis=2); u = u / np.maximum(s, 1e-12)[:, :, None]
    out = []
    for sg in (1.0, -1.0):
        ra = pa[None, :] / (1 + sg * np.einsum('ajk,jk->aj', u, Ea))
        rc = pc[:, None] / (1 + sg * np.einsum('ajk,ak->aj', u, Ec))
        g = ra - rc
        g[(ra <= 0) | (ra > 3.0)] = 9.0
        out.append(g)
    return np.stack(out, 2), s


def earth_state(L, el=EARTH_ELEMENTS):
    """Earth heliocentric state at true longitude-like angle L (rad, measured from Earth's perihelion)."""
    a, e = el[0] * AU, el[1]
    i, Om, om = np.radians(el[2:5])
    p = a * (1 - e * e); r = p / (1 + e * np.cos(L))
    rp = np.stack([r * np.cos(L), r * np.sin(L), 0 * L], 1)
    vp = np.sqrt(MU / p) * np.stack([-np.sin(L), e + np.cos(L), 0 * L], 1)
    cO, sO, ci, si, co, so = np.cos(Om), np.sin(Om), np.cos(i), np.sin(i), np.cos(om), np.sin(om)
    R = np.array([[cO * co - sO * so * ci, -cO * so - sO * co * ci, sO * si],
                  [sO * co + cO * so * ci, -sO * so + cO * co * ci, -cO * si],
                  [so * si, co * si, ci]])
    return rp @ R.T, vp @ R.T


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=200000); ap.add_argument('--eps', type=float, default=0.003)
    ap.add_argument('--k', type=int, default=10); ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out', default=''); ap.add_argument('--chunk', type=int, default=20000)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    ids, elem = load_mea()
    keep = np.array([int(t) not in UNREACHABLE for t in ids])
    el = elem[keep, :5]; ids = ids[keep]
    pa, na, Ea = frame(el)

    # reference: Earth's own orbit as the carrier
    pe, ne, Ee = frame(EARTH_ELEMENTS[None, :5])
    g, s = gaps(pe, ne, Ee, pa, na, Ea); gm = np.abs(g).min(2)[0]
    print('carrier = Earth orbit: asteroids with node gap < eps')
    for eps in (0.001, 0.002, 0.003, 0.005, 0.01, 0.02, 0.05):
        print(f'   eps {eps:6.3f} AU: {(gm < eps).sum():4d}')

    # launch-reachable carriers: Earth state at a random longitude + v_inf anywhere in the 4 km/s ball
    hit = []; par = []
    for c0 in range(0, a.n, a.chunk):
        m = min(a.chunk, a.n - c0)
        L = rng.uniform(0, 2 * np.pi, m)
        d = rng.normal(size=(m, 3)); d /= np.linalg.norm(d, axis=1)[:, None]
        vinf = d * (VINF_MAX * rng.uniform(0, 1, m) ** (1 / 3))[:, None]
        r, v = earth_state(L); v = v + vinf
        pc, nc, Ec, ac = rv2frame(r, v)
        g, s = gaps(pc, nc, Ec, pa, na, Ea)
        ok = np.abs(g).min(2) < a.eps
        hit.append(np.packbits(ok, axis=1)); par.append(np.column_stack([L, vinf, ac, np.linalg.norm(Ec, axis=1),
                                                                        np.degrees(np.arccos(np.clip(nc[:, 2], -1, 1)))]))
    H = np.unpackbits(np.concatenate(hit), axis=1)[:, :len(ids)].astype(bool); par = np.concatenate(par)
    cnt = H.sum(1)
    print(f'\n{a.n} launch-reachable carriers, eps {a.eps} AU: asteroids within eps per carrier')
    print('   percentiles  50/90/99/max: ' + ' '.join(f'{np.percentile(cnt, q):.0f}' for q in (50, 90, 99, 100)))
    print(f'   asteroids reachable by at least one carrier: {H.any(0).sum()} of {len(ids)}')

    # greedy set cover
    left = np.ones(len(ids), bool); chosen = []
    for k in range(a.k):
        gain = (H & left).sum(1); j = int(gain.argmax())
        chosen.append(j); left &= ~H[j]
        L, vx, vy, vz, ac, ec, ic = par[j]
        print(f'   carrier {k+1:2d}: +{gain[j]:3d} (own {cnt[j]:3d})  a {ac:.3f} e {ec:.3f} i {ic:5.2f} deg'
              f'  |vinf| {np.linalg.norm([vx, vy, vz]):.2f}  -> uncovered {left.sum()}')
    print(f'   uncovered after {a.k}: {sorted(int(x) for x in ids[left])}')
    if a.out:
        json.dump(dict(eps=a.eps, ids=[int(x) for x in ids],
                       carriers=[dict(L=float(par[j][0]), vinf=[float(x) for x in par[j][1:4]],
                                      targets=[int(x) for x in ids[H[j]]]) for j in chosen]), open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()
