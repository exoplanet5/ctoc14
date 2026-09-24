"""Stage 17 CAMPAIGN: encounter windows (pool flyby epochs, clustered at 60 d) of the 12 misses of the N<=8 286 selection,
and the number of pool families flying each.  Output miss_windows.json"""
import json, pickle, pathlib, collections
import numpy as np
OUT = pathlib.Path('/Users/mickey/solarsystem/ctoc14/results/s17/campaign')
D = pickle.load(open(OUT / 'reps.pkl', 'rb')); M = json.load(open(OUT / 'milp.json'))['parents_N8']
res = {}
for m in M['misses']:
    e = sorted(t / 86400 for tf, a in zip(D['fl_t'], D['fl_a']) for t, x in zip(tf, a) if int(x) == m)
    wins = []
    for t in e:
        if wins and t - wins[-1][-1] <= 60:
            wins[-1].append(t)
        else:
            wins.append([t])
    res[m] = dict(families=len(e), windows=[(int(w[0]), len(w)) for w in wins])
json.dump(res, open(OUT / 'miss_windows.json', 'w'), indent=1)
for m, v in res.items():
    print(m, v['families'], v['windows'])
