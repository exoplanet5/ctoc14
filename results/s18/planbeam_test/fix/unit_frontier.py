"""Unit check of the frontier_select cut (reviewer finding: the deepest columns were dropped) and of candidates()'s
launch / continuation fields on a synthetic plan."""
import sys, os
sys.path.insert(0, '/Users/mickey/solarsystem/ctoc14/tools'); sys.path.insert(0, '/Users/mickey/solarsystem/ctoc14')
import numpy as np
import s18_planbeam as PB
from ctoc14.constants import DAY

# 6 depths x 5 candidates, distinct target sets
cands = []
for d in range(1, 7):
    for j in range(5):
        legs = [(100 * d + 10 * j + i, (1650 + 60 * i) * DAY) for i in range(d)]
        cands.append(dict(legs=legs, pot=j % 3, pot_cont=j % 3 + 1, fuel_c=10.0 * d + j, first=legs[0][0], n_new=d,
                          plan_n=d + 2, plan_fuel=100.0, plan_t_end=3000 * DAY, score=0.0, t_launch=0.0, vinf=np.zeros(3)))
ch = PB.frontier_select(cands, dict(PB.PB_DEF))
cnt = {}
for c in ch:
    cnt[c['n_new']] = cnt.get(c['n_new'], 0) + 1
print('kept per depth (max_cols 16, per_depth 3, 6 depths):', dict(sorted(cnt.items())), 'total', len(ch))
assert len(ch) == 16
assert cnt[6] == 3 and cnt[5] == 3 and cnt[1] >= 2, cnt        # every depth keeps ranks 0 and 1, deepest get rank 2
# the same with 4 depths (the measured case): unchanged from before (3 per depth, 12 + rest)
c4 = [c for c in cands if c['n_new'] <= 4]
ch4 = PB.frontier_select(c4, dict(PB.PB_DEF))
cnt4 = {}
for c in ch4:
    cnt4[c['n_new']] = cnt4.get(c['n_new'], 0) + 1
print('kept per depth with 4 depths:', dict(sorted(cnt4.items())), 'total', len(ch4))
assert len(ch4) == 16 and all(cnt4[d] >= 3 for d in range(1, 5))

# candidates(): synthetic states with two launches sharing a first leg
from ctoc14.search import State
T = 1620 * DAY; H = 540 * DAY; L = 360 * DAY
def st(tl, seq, fuel):
    return State(t=seq[-1][1], r=np.zeros(3), v=np.zeros(3), m=1600 - fuel, visited=0, seq=tuple(seq), fuel=fuel,
                 t_launch=tl, vinf=np.array([tl / DAY, 0, 0]), m0=1600.0, rare=0.0)
seq_a = [(5, 1700 * DAY, 0.5, 0), (7, 1900 * DAY, 0.4, 0), (9, 2100 * DAY, 0.3, 0), (11, 2300 * DAY, 0.3, 0)]
states = [st(1630 * DAY, seq_a[:1], 20.0), st(1650 * DAY, seq_a[:1], 30.0), st(1650 * DAY, seq_a, 90.0)]
res = dict(states=states, n_pre=0)
cs = PB.candidates(res, T, H, L, None, {}, dict(PB.PB_DEF), print)
by = {frozenset(a for a, _ in c['legs']): c for c in cs}
c5 = by[frozenset([5])]; c57 = by[frozenset([5, 7])]; c579 = by[frozenset([5, 7, 9])]
print('{5}: launch', c5['t_launch'] / DAY, 'pot', c5['pot'], 'cont', c5['pot_cont'], 'fuel_c', c5['fuel_c'])
print('{5,7}: pot', c57['pot'], 'cont', c57['pot_cont'], '{5,7,9}: pot', c579['pot'], 'cont', c579['pot_cont'])
assert c5['t_launch'] == 1650 * DAY and c5['pot_cont'] == 3 and c5['pot'] == 1     # the deep plan wins the dedupe by continuation
assert c57['pot_cont'] == 2 and c579['pot_cont'] == 1 and c579['pot'] == 1
assert abs(c5['fuel_c'] - 30.0) < 1e-9                                            # fuel of the same-launch ancestor, not the 1630 d one
print('unit_frontier OK')
