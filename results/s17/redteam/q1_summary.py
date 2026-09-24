"""Q1: dissolve-by-insertion price of each route V (lin_price screen, results/s17/redteam/dissolve_*.jsonl).
Per target: cheapest host by lin_dJ among approaches with a TRUSTED linear step (res_km <= RES) and unfiltered.
Greedy host assignment with convex host cost (host dkg summed, J recomputed), bins by price."""
import json, glob, collections, sys
import numpy as np
M_DRY = 600.0
def cost(t): x = (t - M_DRY) / 1400.0; return 1 + x + x * x
R = []
for f in glob.glob('/Users/mickey/solarsystem/ctoc14/results/s17/redteam/dissolve_*.jsonl'):
    R += [json.loads(l) for l in open(f)]
fleet = json.load(open('/Users/mickey/solarsystem/ctoc14/results/s16/best/fleet/fleet.json'))['routes']
import numpy as np
out = {}
for V in sorted(fleet):
    z = np.load(f'/Users/mickey/solarsystem/ctoc14/results/s16/best/fleet/route_{V}.npz'); asts = [int(a) for a in z['asts']]
    rows = [r for r in R if r['V'] == V]
    if not rows: continue
    by = collections.defaultdict(list)
    for r in rows: by[r['ast']].append(r)
    res = {}
    for RES in (3e7, 1e12):
        best = {}
        for a in asts:
            L = [r for r in by.get(a, []) if r['res_km'] <= RES]
            if L: best[a] = min(L, key=lambda r: r['lin_dJ'])
        # greedy convex: order targets by their min price, assign to host minimising incremental J with host's accumulated kg
        add = collections.defaultdict(float); tot = 0.0
        for a in sorted(best, key=lambda a: best[a]['lin_dJ']):
            opts = [r for r in by[a] if r['res_km'] <= RES]
            def inc(r):
                h = r['host']; t0 = fleet[h]['tank'] + add[h]
                return cost(t0 + r['lin_kg']) - cost(t0)
            r = min(opts, key=inc); tot += inc(r); add[r['host']] += r['lin_kg']
        p = np.array([best[a]['lin_dJ'] for a in best])
        bins = dict(le03=int((p <= 0.03).sum()), le06=int(((p > 0.03) & (p <= 0.06)).sum()), le15=int(((p > 0.06) & (p <= 0.15)).sum()), gt15=int((p > 0.15).sum()))
        res['trusted' if RES < 1e11 else 'all'] = dict(placeable=len(best), unplaceable=len(asts) - len(best), sum_min_dJ=round(float(p.sum()), 3),
                                                       greedy_convex_dJ=round(tot, 3), bins=bins,
                                                       cheap_dJ_sum=round(float(p[p <= 0.06].sum()), 3),
                                                       host_kg={h: round(v, 1) for h, v in add.items()})
    Ji = fleet[V]['J_i']
    out[V] = dict(n=len(asts), J_i=round(Ji, 4), **res)
    t = res['all']
    print(f"{V}: n {len(asts)} J_i {Ji:.3f} | all: placeable {t['placeable']} sum {t['sum_min_dJ']:.2f} convex {t['greedy_convex_dJ']:.2f} bins {t['bins']} cheap-sum {t['cheap_dJ_sum']}"
          f" | trusted(res<=3e7km): placeable {res['trusted']['placeable']} sum {res['trusted']['sum_min_dJ']:.2f} bins {res['trusted']['bins']}")
json.dump(out, open('/Users/mickey/solarsystem/ctoc14/results/s17/redteam/q1_dissolve.json', 'w'), indent=1)
