"""Rigid-orbit coverage test: how many asteroids can K in-ecliptic craft on fixed (a, e, omega) orbits cover if every
flyby happens at an asteroid's ecliptic-node crossing (z = 0)? Events = node crossings with |r - 1 AU| < 0.08 AU.
Coverage: |phase error| < PHI deg and |r_craft - r_event| < DR AU, t_event >= launch time (craft phase = Earth's at launch).
Compared with the same test for craft with inclination <= 3 deg matching min-distance events (current fleet style)."""
import sys, numpy as np
sys.path.insert(0, '.')
from ctoc14.kepler import Ephemeris
from ctoc14.constants import AU, DAY, T_MISSION, MU
from ctoc14.search import UNREACHABLE
E = Ephemeris(); rng = np.random.default_rng(1)
YR = 365.25 * DAY; ids = np.array([k for k in range(1, 301) if k not in UNREACHABLE])
tt = np.arange(0, T_MISSION, 0.25 * DAY)
# Earth mean longitude reference (linear): lambda_E(t) ~ L0 + nE t
rE, vE = E.earth_state(np.array([0.0, YR / 4]))
nE = 2 * np.pi / (365.256 * DAY); L0 = np.arctan2(rE[0, 1], rE[0, 0])
def wrap(x): return (x + np.pi) % (2 * np.pi) - np.pi
ev_t, ev_lam, ev_r, ev_z, ev_id, ev_kind = [], [], [], [], [], []
for k in ids:
    ra, va = E.ast_states_at(np.full(len(tt), k - 1), tt)
    z = ra[:, 2] / AU; rxy = np.hypot(ra[:, 0], ra[:, 1]) / AU; lam = np.arctan2(ra[:, 1], ra[:, 0])
    i = np.where(np.sign(z[1:]) != np.sign(z[:-1]))[0]
    for j in i:
        f = z[j] / (z[j] - z[j + 1]); r = rxy[j] + f * (rxy[j + 1] - rxy[j])
        if abs(r - 1) < 0.08:
            ev_t.append(tt[j] + f * (tt[1] - tt[0])); ev_lam.append(lam[j] + f * wrap(lam[j + 1] - lam[j])); ev_r.append(r); ev_z.append(0.0); ev_id.append(k); ev_kind.append(0)
    d = np.hypot(rxy - 1, z)
    i = np.where((d[1:-1] <= d[:-2]) & (d[1:-1] < d[2:]) & (d[1:-1] < 0.08))[0] + 1
    for j in i:
        ev_t.append(tt[j]); ev_lam.append(lam[j]); ev_r.append(rxy[j]); ev_z.append(z[j]); ev_id.append(k); ev_kind.append(1)
ev = dict(t=np.array(ev_t), lam=np.array(ev_lam), r=np.array(ev_r), z=np.array(ev_z), id=np.array(ev_id), kind=np.array(ev_kind))
np.savez('results/newgen/scratch/events.npz', **ev)
kn = ev['kind'] == 0
print(f'node events |r-1|<0.08: {kn.sum()} ({len(set(ev["id"][kn]))} asteroids), |r-1|<0.05: {(kn & (np.abs(ev["r"] - 1) < 0.05)).sum()} '
      f'({len(set(ev["id"][kn & (np.abs(ev["r"] - 1) < 0.05)]))} asteroids); min-distance events d<0.08: {(~kn).sum()}')
miss_node = sorted(set(ids) - set(ev['id'][kn & (np.abs(ev['r'] - 1) < 0.05)]))
print(f'asteroids without a node crossing within 0.05 AU of 1 AU: {miss_node}')

def cover_curve(sel, PHI, DR, DZ, incl, nsamp=400000, K=14, label=''):
    t = ev['t'][sel]; lam = ev['lam'][sel]; r = ev['r'][sel]; z = ev['z'][sel]; aid = ev['id'][sel]
    aidx = np.searchsorted(ids, aid)
    # sample rigid craft: launch time tL, delta-a, e, omega (+ i, Omega)
    tL = rng.uniform(0, 12 * YR, nsamp); da = rng.uniform(-0.12, 0.12, nsamp); e = rng.uniform(0, 0.12, nsamp)
    om = rng.uniform(-np.pi, np.pi, nsamp)
    inc = rng.uniform(0, np.radians(incl), nsamp) if incl > 0 else np.zeros(nsamp); Om = rng.uniform(-np.pi, np.pi, nsamp)
    a = 1 + da; n = nE * a ** -1.5
    lamL = L0 + nE * tL                                          # craft longitude at launch = Earth's
    # launch consistency: r = 1 at launch -> e cos(lamL - om) = 1 - 1/a  (adjust e sign via omega choice: keep if |mismatch| < 0.03)
    ok = np.abs(a * (1 - e * np.cos(lamL - om)) - 1) < 0.03
    tL, da, e, om, inc, Om, a, n, lamL = (x[ok] for x in (tL, da, e, om, inc, Om, a, n, lamL))
    nl = len(tL); covered = np.zeros(len(ids), bool); chosen = []
    cov_sets = None
    for kk in range(K):
        best = (-1, None, None)
        for c0 in range(0, nl, 4000):
            s = slice(c0, c0 + 4000)
            lc = lamL[s][:, None] + n[s][:, None] * (t[None, :] - tL[s][:, None])
            ph = np.abs(wrap(lc - lam[None, :])) < np.radians(PHI)
            rc = a[s][:, None] * (1 - e[s][:, None] * np.cos(lc - om[s][:, None]))
            zc = a[s][:, None] * np.sin(inc[s][:, None]) * np.sin(lc - Om[s][:, None])
            m = ph & (np.abs(rc - r[None, :]) < DR) & (np.abs(zc - z[None, :]) < DZ) & (t[None, :] >= tL[s][:, None])
            # new asteroids per craft
            newm = m & ~covered[aidx][None, :]
            # count distinct asteroids: use bincount over (row, aidx) pairs
            rows, cols = np.nonzero(newm)
            if len(rows) == 0: continue
            pairs = np.unique(rows * 1000 + aidx[cols])
            cnt = np.bincount(pairs // 1000, minlength=newm.shape[0])
            j = int(np.argmax(cnt))
            if cnt[j] > best[0]:
                best = (int(cnt[j]), c0 + j, np.unique(aidx[np.nonzero(m[j])[0]]))
        if best[0] <= 0: break
        cnt, j, cov = best; covered[cov] = True; chosen.append(j)
        print(f'  {label} craft {kk + 1}: +{cnt} new (covers {len(cov)} total on its own) tL {tL[j] / YR:.1f}y da {da[j]:+.3f} e {e[j]:.3f} i {np.degrees(inc[j]):.1f} -> covered {covered.sum()}/{len(ids)}')
    return covered

sel_node = (ev['kind'] == 0) & (np.abs(ev['r'] - 1) < 0.06)
print('== in-ecliptic rigid craft, node events, phase +-2.5 deg, dr 0.02 AU')
c1 = cover_curve(sel_node, 2.5, 0.02, 1e9, 0.0, label='node')
print('== in-ecliptic rigid craft, node events, phase +-4 deg, dr 0.03 AU (cheap corrections allowed)')
c2 = cover_curve(sel_node, 4.0, 0.03, 1e9, 0.0, label='node4')
print('== inclined (<=3 deg) rigid craft, min-distance events, phase +-2.5 deg, dr 0.02, dz 0.015')
c3 = cover_curve(ev['kind'] == 1, 2.5, 0.02, 0.015, 3.0, label='incl')
