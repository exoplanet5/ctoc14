"""Stage 17 CAMPAIGN diagnostic: phase structure and internal tail-exchange (2-opt*) junctions of the s16a twin fleet.

For the 9 routes of results/s16/best/fleet (read-only): heliocentric longitude of every craft on a 10 d grid,
pairwise longitude separation statistics, and every Lambert junction A(t) -> B(t+T) <= 1.0 km/s (T 30..150 d).
Output results/s17/campaign/fleet_phase.json
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, pathlib, itertools
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    sys.path.insert(0, _p)
os.chdir(ROOT)
import numpy as np
import warnings; warnings.filterwarnings('ignore')
from ctoc14.kepler import propagate_batch
from ctoc14.lambert import lambert
from ctoc14.constants import DAY, AU, VE, T_MISSION
import run_ialns as RI

GSTEP = 10 * DAY; GRID = np.arange(0, T_MISSION, GSTEP); K = len(GRID)
names = [f'r{i}' for i in range(1, 10)]
S = np.full((9, K, 6), np.nan); DVC = np.full((9, K), np.nan); fl = []; DVT = []
for i, n in enumerate(names):
    z = np.load(f'results/s16/best/fleet/route_{n}.npz'); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL'])
    ip = RI.ipr(st); Yf, Ys = ip.integrate(); nrm = np.linalg.norm(ip.Ts, axis=1)
    post_r = np.vstack([ip.rE[None], Ys[:, :3]]); post_v = np.vstack([(ip.vE + ip.vinf)[None], Ys[:, 3:6] + ip.Ts])
    node_t = np.concatenate([[ip.tL], ip.ts]); cum = np.concatenate([[0.0], np.cumsum(nrm)])
    ks = np.nonzero(GRID >= ip.tL)[0]; m = np.searchsorted(node_t, GRID[ks], side='right') - 1
    r, v = propagate_batch(post_r[m], post_v[m], GRID[ks] - node_t[m])
    S[i, ks, :3] = r; S[i, ks, 3:] = v; DVC[i, ks] = cum[m]
    o = np.argsort(ip.tf); fl.append((np.asarray(ip.tf)[o], np.asarray(ip.asts)[o])); DVT.append(nrm.sum())
lon = np.degrees(np.arctan2(S[:, :, 1], S[:, :, 0]))
out = dict(pairs={})
seps = []
for a, b in itertools.combinations(range(9), 2):
    d = np.abs((lon[a] - lon[b] + 180) % 360 - 180); d = d[np.isfinite(d)]
    dist = np.linalg.norm(S[a, :, :3] - S[b, :, :3], axis=1) / AU; dist = dist[np.isfinite(dist)]
    seps.append(np.median(d))
    out['pairs'][f'{names[a]}-{names[b]}'] = dict(med_lon_sep=round(float(np.median(d)), 1),
                                                  frac_within_30deg=round(float((d < 30).mean()), 3),
                                                  frac_within_0p05AU=round(float((dist < 0.05).mean()), 4))
out['median_pair_lon_sep_deg'] = round(float(np.median(seps)), 1)
# all-pairs junctions A(k) -> B(k+L)
J = []
for L in (3, 5, 8, 11, 15):
    T = L * GSTEP
    for a, b in itertools.permutations(range(9), 2):
        ks = np.nonzero(np.isfinite(S[a, :K - L, 0]) & np.isfinite(S[b, L:, 0]))[0]
        if len(ks) == 0:
            continue
        v1, v2 = lambert(S[a, ks, :3], S[b, ks + L, :3], np.full(len(ks), T))
        d = np.linalg.norm(v1 - S[a, ks, 3:], axis=1) + np.linalg.norm(S[b, ks + L, 3:] - v2, axis=1)
        for k, dv in zip(ks[d <= 1.0], d[d <= 1.0]):
            t1 = GRID[k]; t2 = GRID[k + L]
            npre = int((fl[a][0] < t1).sum()); nsuf = int((fl[b][0] > t2).sum())
            tank = 601.5 * np.exp((DVC[a, k] + dv + DVT[b] - DVC[b, k + L]) / VE)
            J.append(dict(a=names[a], b=names[b], day=round(t1 / DAY), T=L * 10, dvj=round(float(dv), 3),
                          depth=npre + nsuf, tank=round(float(tank), 1)))
out['n_junctions_le_1kms'] = len(J)
out['n_junctions_le_0p5kms'] = sum(j['dvj'] <= 0.5 for j in J)
out['pairs_with_junction'] = sorted(set(f"{j['a']}->{j['b']}" for j in J))
out['junctions_best'] = sorted(J, key=lambda j: j['dvj'])[:40]
json.dump(out, open('results/s17/campaign/fleet_phase.json', 'w'), indent=1)
print('median pair lon sep', out['median_pair_lon_sep_deg'], 'junctions<=1', len(J), '<=0.5', out['n_junctions_le_0p5kms'])
print('pairs', out['pairs_with_junction'])
for j in out['junctions_best'][:15]:
    print(j)
