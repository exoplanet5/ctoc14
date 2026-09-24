"""Per-route timeline of the s16 best twin fleet: flybys and delta-v per 2-yr bin, per-leg dv, suffix budgets."""
import sys, json, pathlib
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from ctoc14.constants import DAY, VE
RI.eph()
F = {}
for f in sorted((ROOT / 'results/s16/best/fleet').glob('route_*.npz')):
    z = np.load(f); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL']); F[f.stem[6:]] = st
bins = np.arange(0, 5479 + 730, 730)
out = {}
print('route   n  tank  dv   | flybys per 2yr bin                | dv (km/s) per 2yr bin              | top-3 legs dv share')
for k, st in F.items():
    ip = RI.ipr(st); tank = float(ip.tank())
    ts = np.asarray(st['ts']) / DAY; T = np.linalg.norm(np.asarray(st['Ts']).reshape(-1, 3), axis=1)
    tf = np.sort(np.asarray(st['tf']) / DAY)
    nb = np.histogram(tf, bins)[0]; dvb = np.histogram(ts, bins, weights=T)[0]
    # per leg dv (impulses in (tf[k-1], tf[k]]; first leg from launch)
    edges = np.concatenate([[st['tL'] / DAY - 1], tf])
    leg = np.array([T[(ts > edges[j]) & (ts <= edges[j + 1])].sum() for j in range(len(tf))])
    after = T[ts > tf[-1]].sum()
    vinf = float(np.linalg.norm(st['vinf']))
    s = np.sort(leg)[::-1]
    out[k] = dict(n=len(tf), tank=round(tank, 1), dv=round(float(T.sum()), 3), vinf=round(vinf, 3), fb_per_2yr=nb.tolist(),
                  dv_per_2yr=[round(float(x), 2) for x in dvb], leg_dv=[round(float(x), 3) for x in leg],
                  top3_share=round(float(s[:3].sum() / max(T.sum(), 1e-9)), 3), leg_median=round(float(np.median(leg)), 3),
                  legs_gt_1=int((leg > 1.0).sum()), gap_median=round(float(np.median(np.diff(edges))), 1))
    print(f'{k:4s} {len(tf):3d} {tank:6.1f} {T.sum():5.2f} | {" ".join(f"{x:2d}" for x in nb)} | {" ".join(f"{x:4.2f}" for x in dvb)} | {out[k]["top3_share"]:.2f} med {out[k]["leg_median"]:.2f} >1:{out[k]["legs_gt_1"]} gap {out[k]["gap_median"]}')
tot_n = sum(np.array(v['fb_per_2yr']) for v in out.values()); tot_dv = sum(np.array(v['dv_per_2yr']) for v in out.values())
print('fleet', tot_n.tolist(), [round(float(x), 1) for x in tot_dv])
json.dump(out, open(ROOT / 'results/s17/redteam/timeline.json', 'w'), indent=1)
