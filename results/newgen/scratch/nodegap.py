"""Are the fleet's flybys already at the asteroids' node crossings? For each flyby: |z|, and the nearest event of the
same asteroid with |z| <= 0.01 AU (time gap and how much plane Delta-v it would save)."""
import sys, numpy as np
sys.path.insert(0, '.'); sys.path.insert(0, 'tools')
from run_ialns import IFleet, ipr, eph
from ctoc14.constants import AU, DAY
from ctoc14.phasemodel import plane_cost
E = eph(); F = IFleet('results/newgen/ifleet10')
ev = dict(np.load('results/newgen/scratch/events2.npz')); ev['ast_id'] = np.asarray(ev['ids'])[ev['ast']]
gaps = []; zs = []; saves = []
for n, r in sorted(F.routes.items()):
    ip = ipr(r['st'])
    ra, _ = E.ast_states_at(np.array(ip.asts) - 1, ip.tf)
    z = np.abs(ra[:, 2]) / AU
    for X, t, zz in zip(ip.asts, ip.tf, z):
        m = np.nonzero((ev['ast_id'] == X) & (np.abs(ev['z']) <= 0.01))[0]
        zs.append(zz)
        if len(m) == 0:
            gaps.append(np.nan); continue
        dt = np.abs(ev['t'][m] - t) / DAY
        gaps.append(dt.min())
        saves.append(plane_cost(zz) - plane_cost(abs(ev['z'][m][np.argmin(dt)])))
g = np.array(gaps, float); zs = np.array(zs); sv = np.array(saves)
print(f'flybys: {len(zs)}; |z| at flyby p50 {np.percentile(zs,50):.4f} p90 {np.percentile(zs,90):.4f} AU')
print(f'targets with a node-crossing event (|z|<=0.01): {np.isfinite(g).sum()} of {len(g)}')
print(f'time gap to the nearest node crossing: p10 {np.nanpercentile(g,10):.0f} p50 {np.nanpercentile(g,50):.0f} '
      f'p90 {np.nanpercentile(g,90):.0f} days')
print(f'flybys already within 30 d of a node crossing: {(g < 30).sum()}; within 100 d: {(g < 100).sum()}')
print(f'plane Delta-v that node crossings would save (round-trip model): total {sv.sum():.1f} km/s '
      f'over {len(sv)} flybys, median {np.median(sv):.3f} km/s')
print(f'fleet actual plane Delta-v: 62.6 km/s')
