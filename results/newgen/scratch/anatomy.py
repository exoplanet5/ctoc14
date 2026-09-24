"""Where does the Delta-v of the 10-craft fleet go, how many crossing events do asteroids offer, how often do routes cross."""
import sys, pathlib, numpy as np
sys.path.insert(0, '.'); sys.path.insert(0, 'tools')
from run_ialns import IFleet, ipr, trajectory, eph
from ctoc14.constants import AU, DAY, T_MISSION, VE
from ctoc14.search import UNREACHABLE
E = eph()
F = IFleet(sys.argv[1] if len(sys.argv) > 1 else 'results/newgen/ifleet10')
print('== A. fleet Delta-v anatomy')
tot_dv = 0; nfb = 0; rows = []
allfb = []
for n, r in sorted(F.routes.items()):
    ip = ipr(r['st']); Yf, Ys = ip.integrate()
    dvm = np.linalg.norm(ip.Ts, axis=1); dv = dvm.sum(); tot_dv += dv; nfb += len(ip.tf)
    # RTN decomposition of impulses at the pre-impulse states
    rr = Ys[:, :3]; vv = Ys[:, 3:6]
    R = rr / np.linalg.norm(rr, axis=1)[:, None]; N = np.cross(rr, vv); N /= np.linalg.norm(N, axis=1)[:, None]; T = np.cross(N, R)
    cr = np.abs((ip.Ts * R).sum(1)).sum(); ct = np.abs((ip.Ts * T).sum(1)).sum(); cn = np.abs((ip.Ts * N).sum(1)).sum()
    # orbital elements along the route (post-impulse states)
    h = np.cross(rr, vv); hn = np.linalg.norm(h, axis=1); inc = np.degrees(np.arccos(h[:, 2] / hn))
    v2 = (vv ** 2).sum(1); rn = np.linalg.norm(rr, axis=1); mu = 1.32712440018e11
    a = 1 / (2 / rn - v2 / mu); ev = np.linalg.norm(np.cross(vv, h) / mu - rr / rn[:, None], axis=1)
    # flyby geometry: asteroid offset from the 1 AU circle, relative speed
    ra, va = E.ast_states_at(np.array(ip.asts) - 1, ip.tf)
    rxy = np.hypot(ra[:, 0], ra[:, 1]) / AU - 1; z = ra[:, 2] / AU; vrel = np.linalg.norm(Yf[:, 3:6] - va, axis=1)
    # per-leg dv: impulses between consecutive flybys
    edges = np.concatenate([[ip.tL], np.sort(ip.tf)]); leg = np.searchsorted(edges, ip.ts, side='right') - 1
    legdv = np.bincount(np.clip(leg, 0, len(ip.tf) - 1), dvm, minlength=len(ip.tf))
    o = np.argsort(ip.tf)
    for j, k in enumerate(o):
        allfb.append((n, ip.asts[k], ip.tf[k] / DAY / 365.25, legdv[j], rxy[k], z[k], vrel[k]))
    rows.append((n, len(ip.tf), dv, cr, ct, cn, a.min() / AU, a.max() / AU, ev.max(), inc.max(), np.linalg.norm(ip.vinf), ip.tank()))
    print(f'route {n}: {len(ip.tf):2d} flybys, dv {dv:6.2f} km/s ({dv / len(ip.tf):.3f}/flyby) RTN split R {cr / dv:.2f} T {ct / dv:.2f} N {cn / dv:.2f} '
          f'a {a.min() / AU:.3f}-{a.max() / AU:.3f} e<={ev.max():.3f} i<={inc.max():.2f} vinf {np.linalg.norm(ip.vinf):.2f} tank {ip.tank():.0f}')
print(f'fleet: {nfb} flybys, {tot_dv:.1f} km/s, {tot_dv / nfb:.3f} km/s per flyby')
allfb = np.array(allfb, dtype=object)
legdv = allfb[:, 3].astype(float); rxy = allfb[:, 4].astype(float); z = allfb[:, 5].astype(float); vrel = allfb[:, 6].astype(float)
print('per-leg dv quantiles (km/s): ' + ' '.join(f'p{q}={np.percentile(legdv, q):.2f}' for q in (10, 25, 50, 75, 90, 100)))
print(f'share of fleet dv in the costliest 20% of legs: {np.sort(legdv)[::-1][:len(legdv) // 5].sum() / legdv.sum():.2f}')
print('flyby geometry: |r_xy-1AU| quantiles ' + ' '.join(f'p{q}={np.percentile(np.abs(rxy), q):.3f}' for q in (50, 75, 90, 100)) +
      ' | |z| quantiles ' + ' '.join(f'p{q}={np.percentile(np.abs(z), q):.3f}' for q in (50, 75, 90, 100)) +
      ' | vrel quantiles ' + ' '.join(f'p{q}={np.percentile(vrel, q):.1f}' for q in (50, 75, 90, 100)))
print('costliest legs: ' + ', '.join(f'{r[0]}:{r[1]}@{r[2]:.1f}y {r[3]:.2f}km/s (dr {r[4]:+.3f}, z {r[5]:+.3f})' for r in sorted(allfb.tolist(), key=lambda r: -r[3])[:15]))
print('== B. crossing events per asteroid (local minima of the distance to the 1 AU circle)')
tt = np.arange(0, T_MISSION, 0.5 * DAY)
ids = np.array([k for k in range(1, 301) if k not in UNREACHABLE])
cnt = {0.03: [], 0.05: [], 0.1: []}; ev_per = {}
elem = {}
for k in ids:
    ra, va = E.ast_states_at(np.full(len(tt), k - 1), tt)
    d = np.hypot(np.hypot(ra[:, 0], ra[:, 1]) / AU - 1, ra[:, 2] / AU)
    i = np.where((d[1:-1] <= d[:-2]) & (d[1:-1] < d[2:]))[0] + 1
    ev_per[k] = d[i]
    for th in cnt: cnt[th].append(int((d[i] < th).sum()))
    b = E.ast[k - 1]; elem[k] = (b.a / AU, b.e, np.degrees(np.arccos(b.R[2, 2])), d[i].min())
for th in cnt:
    c = np.array(cnt[th]); print(f'events with d < {th} AU: total {c.sum()}, per asteroid min {c.min()} median {np.median(c):.0f} max {c.max()}, asteroids with 0: {int((c == 0).sum())}')
hard = [8, 39, 168, 251, 7, 171, 208, 34]
print('unplaceable targets: ' + '; '.join(f'{k}: a {elem[k][0]:.2f} e {elem[k][1]:.2f} i {elem[k][2]:.1f} dmin {elem[k][3]:.3f} n<0.05 {cnt[0.05][list(ids).index(k)]}' for k in hard))
inc = np.array([elem[k][2] for k in ids]); print('inclination quantiles (deg): ' + ' '.join(f'p{q}={np.percentile(inc, q):.1f}' for q in (25, 50, 75, 90, 100)))
print('== C. route pair crossings (distance < 0.05 AU and relative speed < 0.5 km/s)')
tr = {}
for n, r in sorted(F.routes.items()):
    ip = ipr(r['st']); t, p = trajectory(ip, step=2 * DAY); tr[n] = (t, p)
names = sorted(tr); npair = 0; nx = 0
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        ta, pa = tr[names[i]]; tb, pb = tr[names[j]]
        t0 = max(ta[0], tb[0]); t1 = min(ta[-1], tb[-1])
        if t1 <= t0: continue
        ia = (ta >= t0) & (ta <= t1); ib = (tb >= t0) & (tb <= t1); m = min(ia.sum(), ib.sum())
        d = np.linalg.norm(pa[ia][:m] - pb[ib][:m], axis=1) / AU
        vv = np.linalg.norm(np.gradient(pa[ia][:m], axis=0) - np.gradient(pb[ib][:m], axis=0), axis=1) / (2 * DAY)
        k = np.where((d[1:-1] <= d[:-2]) & (d[1:-1] < d[2:]) & (d[1:-1] < 0.05) & (vv[1:-1] < 0.5))[0] + 1
        npair += 1; nx += len(k)
        if len(k): print(f'  {names[i]}x{names[j]}: ' + ', '.join(f'{ta[ia][kk] / DAY / 365.25:.1f}y d{d[kk]:.3f} v{vv[kk]:.2f}' for kk in k))
print(f'{nx} crossings over {npair} pairs')
