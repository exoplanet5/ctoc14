"""Stage 17 CAMPAIGN: who can cross over with whom?  Junction statistics by parent target-set overlap and by epoch.
Reads reps.pkl, junctions.npy (a, b, k, L, dvj, npre, nsuf); writes junction_stats.json."""
import json, pickle, pathlib, collections
import numpy as np
OUT = pathlib.Path('/Users/mickey/solarsystem/ctoc14/results/s17/campaign')
D = pickle.load(open(OUT / 'reps.pkl', 'rb')); J = np.load(OUT / 'junctions.npy')
sets = [frozenset(r['asts']) for r in D['reps']]; R = len(sets)
a = J[:, 0].astype(int); b = J[:, 1].astype(int); k = J[:, 2].astype(int); dv = J[:, 4]
# distinct ordered route pairs with at least one junction <= x
res = {'n_junctions': int(len(J)), 'n_reps': R}
for thr in (1.0, 0.5, 0.3):
    m = dv <= thr
    pairs = set(zip(a[m].tolist(), b[m].tolist()))
    jac = np.array([len(sets[p] & sets[q]) / len(sets[p] | sets[q]) for p, q in pairs]) if pairs else np.zeros(0)
    ov = np.array([len(sets[p] & sets[q]) for p, q in pairs]) if pairs else np.zeros(0)
    part = collections.Counter(p for p, q in pairs)
    res[f'dv<={thr}'] = dict(junctions=int(m.sum()), route_pairs=len(pairs),
                             frac_all_ordered_pairs=round(len(pairs) / (R * (R - 1)), 5),
                             pair_jaccard_pct=[round(float(x), 3) for x in np.percentile(jac, [10, 50, 90])] if len(jac) else None,
                             pair_overlap_targets_pct=[int(x) for x in np.percentile(ov, [10, 50, 90])] if len(ov) else None,
                             frac_pairs_disjoint=round(float((ov == 0).mean()), 4) if len(ov) else None,
                             partners_per_route_pct=[int(x) for x in np.percentile([part.get(i, 0) for i in range(R)], [10, 50, 90])],
                             junction_day_pct=[int(x) for x in np.percentile(k[m] * 10, [10, 50, 90])])
# control: overlap of random rep pairs
rng = np.random.default_rng(0); p = rng.integers(0, R, 20000); q = rng.integers(0, R, 20000); ok = p != q
ovr = np.array([len(sets[i] & sets[j]) for i, j in zip(p[ok], q[ok])])
res['random_pairs'] = dict(frac_disjoint=round(float((ovr == 0).mean()), 4), overlap_pct=[int(x) for x in np.percentile(ovr, [10, 50, 90])])
json.dump(res, open(OUT / 'junction_stats.json', 'w'), indent=1)
print(json.dumps(res, indent=1))
