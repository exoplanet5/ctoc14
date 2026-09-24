"""Full element-space model vs the real Delta-v of the 10-craft fleet."""
import sys, numpy as np
sys.path.insert(0, '.'); sys.path.insert(0, 'tools')
from run_ialns import IFleet, ipr
from ctoc14.phasemodel import elements, path_cost, drift
F = IFleet('results/newgen/ifleet10')
print('route | actual dv | predicted | drift part | plane part | ratio')
rows = []
for n, r in sorted(F.routes.items()):
    ip = ipr(r['st']); Yf, Ys = ip.integrate()
    el = elements(Ys[:, :3], Ys[:, 3:6])
    pred, parts = path_cost(el)
    rows.append((ip.dv(), pred))
    print(f'{n} | {ip.dv():6.2f} | {pred:6.2f} | {parts["drift"]:6.2f} | {parts["plane"]:6.2f} | {pred / ip.dv():.3f}')
rows = np.array(rows)
print(f'corr {np.corrcoef(rows[:,0], rows[:,1])[0,1]:.3f}, mean ratio {np.mean(rows[:,1]/rows[:,0]):.3f}, '
      f'spread {np.std(rows[:,1]/rows[:,0]):.3f}')
