"""Settle drifting-carrier routes (tools/carrier_drift.py) in the impulsive twin and report J per flyby.
A drifting route is physically one smooth trajectory, so it is settled like any tour: Lambert chain seed with the
linearised-leg fallback for the long/degenerate legs (ctoc14.impulsive.settle_tour), then the SCP.
Usage: drift_settle.py routes.json out_dir"""
import os
for _v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'): os.environ.setdefault(_v,'1')
import sys, json, pathlib, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np, run_ialns as RI, pickle
from ctoc14.impulsive import settle_tour
from ctoc14.constants import VE, DAY

lib = pickle.load(open('results/s7/libopt.pkl', 'rb'))


def w(job):
    k, r = job
    ev = sorted([(float(t), int(x)) for s in r['segs'] for x, t, _ in s['base']])
    if len(ev) < 4: return k, None
    ci = r['segs'][0]['ci']
    tour = dict(t_launch=float(lib[ci]['tL']), vinf=[float(q) for q in lib[ci]['vinf']],
                legs=[dict(ast=x, t_flyby=t) for t, x in ev])
    try:
        ip, miss, lag = settle_tour(RI.eph(), tour, lambda ip: RI.settle(ip, 120))
    except Exception as e:
        return k, ('error', repr(e)[:60])
    if ip is None: return k, ('miss', float(miss))
    return k, (RI.ist(ip), float(ip.tank()), len(ip.asts))


R = json.load(open(sys.argv[1])); out = pathlib.Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
with mp.get_context('fork').Pool(8) as pl: res = dict(pl.map(w, list(enumerate(R))))
fl = RI.IFleet()
for k, v in sorted(res.items()):
    n = len(R[k]['targets'])
    if v is None or v[0] in ('miss', 'error'):
        print(f'{k:02d}: {n} flybys -> NOT settled ({v[0] if v else "short"} {v[1] if v else ""})'); continue
    st, tk, nf = v
    fl.routes[f'{k:02d}'] = dict(st=st, tank=tk)
    np.savez(out / f'route_d{k:02d}.npz', **st)
    print(f'{k:02d}: {nf} flybys, tank {tk:.0f} kg, dv {VE*np.log(tk/601.5):.2f} km/s, '
          f'J_i {RI.cost(tk):.3f}, {RI.cost(tk)/nf:.4f} J/fb')
if fl.routes:
    fl.save(out / 'fleet', note='drifting carrier routes')
    print(f'\nfleet: {len(fl.routes)} routes, {len(fl.coverage())} covered, sum J_i '
          f'{sum(RI.cost(r["tank"]) for r in fl.routes.values()):.3f}')
print('reference: best carrier route 0.0339 J/fb, t10d fleet 0.0415, HIT fleet ~0.0353')
