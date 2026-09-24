import os
os.environ['OMP_NUM_THREADS']='1'
import sys, numpy as np
sys.path.insert(0,'.'); sys.path.insert(0,'tools')
from run_ialns import IFleet, ipr, trajectory, eph, cost
from ctoc14.constants import AU
fleet = IFleet(sys.argv[1]); E = eph()
traj = {n: trajectory(ipr(r['st'])) for n, r in fleet.routes.items()}
cov = fleet.coverage()
for V, r in sorted(fleet.routes.items(), key=lambda kv: len(kv[1]['st']['asts'])):
    best = []
    for X in [int(x) for x in r['st']['asts']]:
        dmin = 9; host = None
        for n, (tt, rr) in traj.items():
            if n == V: continue
            ra, _ = E.ast_states_at(np.full(len(tt), X - 1), tt)
            d = np.linalg.norm(rr - ra, axis=1).min() / AU
            if d < dmin: dmin, host = d, n
        best.append((dmin, X, host))
    best.sort()
    ds = np.array([b[0] for b in best])
    print(f'{V}: {len(best)} targets, cost {cost(r["tank"]):.3f}; nearest-host dist median {np.median(ds):.3f}, #>0.1 AU {int((ds>0.1).sum())}, #>0.2 {int((ds>0.2).sum())}; worst ' +
          ', '.join(f'{X}({d:.2f}@{h})' for d, X, h in best[-4:]))
