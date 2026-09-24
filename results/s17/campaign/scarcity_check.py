"""Stage 17 CAMPAIGN: does pool scarcity (families flying a target, encounter windows) predict which targets the best
N<=8 selection misses?  Output scarcity_check.json"""
import json, pickle, pathlib, collections
import numpy as np
OUT = pathlib.Path('/Users/mickey/solarsystem/ctoc14/results/s17/campaign')
D = pickle.load(open(OUT / 'reps.pkl', 'rb')); M = json.load(open(OUT / 'milp.json'))['parents_N8']
fam = collections.Counter(int(x) for a in D['fl_a'] for x in set(a.tolist()))
allv = np.array([fam.get(t, 0) for t in range(1, 301) if t not in (131, 144)])
mv = np.array([fam.get(t, 0) for t in M['misses']])
rank = [float((allv < v).mean()) for v in mv]
res = dict(all_families_pct=[int(x) for x in np.percentile(allv, [10, 25, 50, 75, 90])], misses_families=mv.tolist(),
           misses_percentile_in_all=[round(r, 2) for r in rank], median_percentile=round(float(np.median(rank)), 2),
           n_targets_below_min_miss=int((allv < mv.min()).sum()))
json.dump(res, open(OUT / 'scarcity_check.json', 'w'), indent=1); print(res)
