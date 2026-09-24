import os
os.environ['OMP_NUM_THREADS']='1'
import sys, time, numpy as np
sys.path.insert(0,'.'); sys.path.insert(0,'tools')
from run_ialns import IFleet, ipr, settle, cost, ist
from ctoc14.impulsive import ImpulsiveProblem
from ctoc14.constants import DAY
fleet = IFleet(sys.argv[1]); tot = 0.0
for n, r in sorted(fleet.routes.items()):
    st = r['st']; base = cost(r['tank']); res = []
    for dl in [-30, -15, 15, 30]:
        tL = st['tL'] + dl * DAY
        if tL < 0 or tL >= min(st['tf']) - 5 * DAY: continue
        keep = st['ts'] > tL + 1 * DAY
        ip = ImpulsiveProblem(ipr(st).eph, tL, st['vinf'], st['ts'][keep], st['Ts'][keep], st['tf'], [int(a) for a in st['asts']])
        try:
            miss = settle(ip, 150)
        except Exception as e:
            continue
        if miss <= 150: res.append((cost(ip.tank()) - base, dl))
    best = min(res) if res else (0, 0)
    tot += min(0, best[0])
    print(f'{n}: tL {st["tL"]/DAY:.0f} d, J_i {base:.4f}; shifts ' + ', '.join(f'{d:+d}d {g:+.4f}' for g, d in sorted(res, key=lambda x: x[1])), flush=True)
print('sum of best gains', tot)
