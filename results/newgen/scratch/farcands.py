"""For the costliest legs of the 10-craft fleet: does the target have a cheaper event (closer to the 1 AU circle) that
some route could reach by slow phasing (phase gap at the event time, lead time since the route's previous flyby)?"""
import sys, numpy as np
sys.path.insert(0, '.'); sys.path.insert(0, 'tools')
from run_ialns import IFleet, ipr, trajectory, eph
from ctoc14.constants import AU, DAY, T_MISSION
E = eph(); YR = 365.25 * DAY
F = IFleet('results/newgen/ifleet10'); ev = np.load('results/newgen/scratch/events.npz')
d_ev = np.hypot(ev['r'] - 1, ev['z'])
# route longitude vs time (2 d sampling) and flyby times
R = {}
for n, r in sorted(F.routes.items()):
    ip = ipr(r['st']); t, p = trajectory(ip); R[n] = dict(t=t, lam=np.arctan2(p[:, 1], p[:, 0]), rxy=np.hypot(p[:, 0], p[:, 1]) / AU, z=p[:, 2] / AU, tf=np.sort(ip.tf), ip=ip)
def wrap(x): return (x + np.pi) % (2 * np.pi) - np.pi
legs = []
for n, r in R.items():
    ip = r['ip']; dvm = np.linalg.norm(ip.Ts, axis=1)
    edges = np.concatenate([[ip.tL], np.sort(ip.tf)]); leg = np.clip(np.searchsorted(edges, ip.ts, side='right') - 1, 0, len(ip.tf) - 1)
    legdv = np.bincount(leg, dvm, minlength=len(ip.tf)); o = np.argsort(ip.tf)
    for j, k in enumerate(o):
        legs.append((legdv[j], n, ip.asts[k], ip.tf[k]))
legs.sort(reverse=True); top = legs[:60]
print(f'top-60 legs carry {sum(l[0] for l in top):.1f} of {sum(l[0] for l in legs):.1f} km/s')
n_cheap = 0; n_reach = 0; dv_reach = 0.0
for dv, n, a, tf in top:
    m = ev['id'] == a
    cur = np.hypot(np.hypot(*E.ast_states_at(np.array([a - 1]), np.array([tf]))[0][0, :2]) / AU - 1, E.ast_states_at(np.array([a - 1]), np.array([tf]))[0][0, 2] / AU)
    good = np.where(m & (d_ev < 0.03))[0]
    if len(good) == 0:
        print(f'{n}:{a} {dv:.2f} km/s at d {cur:.3f}: no event with d<0.03 (best {d_ev[m].min():.3f})'); continue
    n_cheap += 1; best = None
    for g in good:
        te = ev['t'][g]
        for rn, rr in R.items():
            if te < rr['t'][0] + 60 * DAY or te > rr['t'][-1]: continue
            i = np.searchsorted(rr['t'], te); gap = np.degrees(abs(wrap(rr['lam'][min(i, len(rr['t']) - 1)] - ev['lam'][g])))
            prev = rr['tf'][rr['tf'] < te - 30 * DAY]; lead = (te - (prev[-1] if len(prev) else rr['ip'].tL)) / YR
            # rough phasing cost: 2 * dphi * v / (3 n T) with T = lead (>= 0.25 yr) + offset costs (e, i) ~ v/2 dr + v dz
            est = 2 * np.radians(gap) * 29.8 / (3 * 2 * np.pi * max(lead, 0.25)) + 0.5 * 29.8 * abs(ev['r'][g] - rr['rxy'][min(i, len(rr['t']) - 1)]) + 29.8 * abs(ev['z'][g] - rr['z'][min(i, len(rr['t']) - 1)])
            if best is None or est < best[0]: best = (est, rn, te / YR, gap, lead, d_ev[g])
    if best and best[0] < dv:
        n_reach += 1; dv_reach += dv - best[0]
    print(f'{n}:{a} {dv:.2f} km/s at d {cur:.3f}: {len(good)} events d<0.03; best alt est {best[0]:.2f} km/s via route {best[1]} at {best[2]:.1f}y (gap {best[3]:.0f} deg, lead {best[4]:.1f} yr, d {best[5]:.3f})')
print(f'{n_cheap}/60 targets have an event within 0.03 AU; {n_reach} have a rough-estimate cheaper alternative (sum of estimated savings {dv_reach:.1f} km/s)')
