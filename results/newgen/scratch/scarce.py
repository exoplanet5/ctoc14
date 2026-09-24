import sys, numpy as np
sys.path.insert(0, '.'); sys.path.insert(0, 'tools')
from run_ialns import IFleet, ipr
F = IFleet('results/newgen/ifleet10'); ev = np.load('results/newgen/scratch/events.npz')
d = np.hypot(ev['r'] - 1, ev['z']); m = ev['kind'] == 1
n03 = {a: int(((ev['id'] == a) & m & (d < 0.03)).sum()) for a in range(1, 301)}
n05 = {a: int(((ev['id'] == a) & m & (d < 0.05)).sum()) for a in range(1, 301)}
print('route: flybys | targets with 0 events <0.05 AU (hard) | with <3 events <0.03 AU (scarce) | tank')
for n, r in sorted(F.routes.items()):
    ip = ipr(r['st']); A = [int(a) for a in ip.asts]
    hard = [a for a in A if n05[a] == 0]; scarce = [a for a in A if n03[a] < 3]
    print(f'{n}: {len(A):2d} | hard {len(hard):2d} {hard} | scarce {len(scarce):2d} | tank {ip.tank():.0f} dv/flyby {ip.dv() / len(A):.2f}')
