import sys, pathlib, json, numpy as np
sys.path.insert(0,'/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0,'/Users/mickey/solarsystem/ctoc14/tools')
import run_ialns as RI
from ctoc14.constants import DAY
F = RI.IFleet('results/s16/best/fleet')
out = {}
for n, r in sorted(F.routes.items()):
    st = r['st']; o = np.argsort(st['tf']); tf = st['tf'][o]/DAY; a = st['asts'][o]
    gaps = np.diff(np.r_[st['tL']/DAY, tf])
    print(f"{n}: n={len(a)} tank={r['tank']:.1f} tL={st['tL']/DAY:.0f} first={tf[0]:.0f} last={tf[-1]:.0f} medgap={np.median(gaps):.0f} maxgap={gaps.max():.0f}")
    print('   ', ' '.join(f'{int(x)}@{t:.0f}' for x, t in zip(a, tf)))
    out[n] = dict(tank=r['tank'], tL=st['tL']/DAY, asts=[int(x) for x in a], tf=[float(t) for t in tf])
json.dump(out, open('results/s17/lns/anatomy.json','w'), indent=1)
