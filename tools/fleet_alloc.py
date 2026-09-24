"""Optimal flyby allocation over per-route Pareto curves (stage 6).

Each skeleton has a measured curve k -> tank(k) (deepest tours of a full-pool fill, settled into twins).  The fleet
problem "cover 298 targets with N routes" is then: choose k_r on each curve with sum k_r >= 298 minimising
sum J_i(tank_r(k_r)).  That is a separable resource allocation: sort the marginal costs J_i(k+1) - J_i(k) of every
route and buy the cheapest flybys until the target is met.  It gives the BEST J this skeleton set can reach if the
packing were perfect, i.e. an optimistic bound to compare against the G1/G2 lines before spending hours on packing.

Usage: fleet_alloc.py 'results/n8s6/G/*.out' [--need 298] [--twin]
"""
import sys, re, glob, argparse
import numpy as np

ap = argparse.ArgumentParser(); ap.add_argument('logs', nargs='+'); ap.add_argument('--need', type=int, default=298)
ap.add_argument('--twin', action='store_true', help='use the settled twin tanks (default: planner Pareto)')
ap.add_argument('--tank-cal', type=float, default=1.0)
a = ap.parse_args()
J = lambda tk: 1.0 + (tk - 600.0) / 1400.0 + ((tk - 600.0) / 1400.0) ** 2

curves = {}
for pat in a.logs:
    for f in sorted(glob.glob(pat)):
        for ln in open(f):
            m = re.search(r'(\d\d): Pareto (.*)$', ln)
            if m:
                d = {}
                for tok in m.group(2).split():
                    k, tk = tok.split(':'); d[int(k)] = float(tk) * a.tank_cal
                curves[m.group(1)] = d
            m2 = re.search(r'(\d+) fb: planner\s+(\d+) -> twin\s+(\d+) kg', ln)
            if m2 and a.twin:
                name = None
                for nm in curves:
                    pass
        if a.twin:                                   # twin tanks are logged without the route name; take them in order
            name = None
            for ln in open(f):
                m = re.search(r'(\d\d): Pareto', ln)
                if m: name = m.group(1)
                m2 = re.search(r'\s+(\d+) fb: planner\s+\d+ -> twin\s+(\d+) kg', ln)
                if m2 and name: curves.setdefault(name, {})[int(m2.group(1))] = float(m2.group(2))
if not curves: print('no Pareto lines found'); sys.exit(1)

print(f'{len(curves)} routes, target {a.need} flybys, tanks = {"twin" if a.twin else "planner"} x {a.tank_cal}')
for n in sorted(curves):
    d = curves[n]
    print(f'  {n}: ' + ' '.join(f'{k}:{d[k]:.0f}' for k in sorted(d)) + f'   (cap {max(d)})')
cap = sum(max(d) for d in curves.values())
print(f'\n  capacity {cap} vs need {a.need}: slack {cap - a.need:+d} ({100.0*(cap-a.need)/a.need:+.1f} %)')
if cap < a.need:
    print('  -> INFEASIBLE with this skeleton set (upper bound below the requirement)')

# greedy over marginal costs, starting from each route's cheapest listed depth
cur = {n: min(d) for n, d in curves.items()}
tot = sum(cur.values()); steps = []
while tot < a.need:
    best = None
    for n, d in curves.items():
        k = cur[n]
        if k + 1 in d:
            dc = J(d[k + 1]) - J(d[k])
            if best is None or dc < best[0]: best = (dc, n)
    if best is None:
        print(f'  ran out of curve at {tot} flybys'); break
    dc, n = best; cur[n] += 1; tot += 1; steps.append((n, cur[n], dc))
sumJ = sum(J(curves[n][cur[n]]) for n in curves)
print(f'\n  optimal allocation ({tot} flybys): ' + ' '.join(f'{n}:{cur[n]}@{curves[n][cur[n]]:.0f}' for n in sorted(cur)))
print(f'  sum J_i {sumJ:.3f} -> raw J {sumJ + 2:.3f}   (shown {0.99643*(sumJ+2):.3f} if submitted 09-27)')
print(f'  lines: G2 13.00, G1 13.86 (banked shown 13.81)')
if steps:
    print(f'  last flybys bought cost {steps[-1][2]:.4f} J each; total marginal spend {sum(s[2] for s in steps):.3f} J')
