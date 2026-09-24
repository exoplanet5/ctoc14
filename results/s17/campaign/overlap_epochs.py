"""Stage 17 CAMPAIGN: are route overlaps SAME-encounter (same epoch) or DIFFERENT-encounter conflicts?
For every target: flyby epochs across the 1220 pool families; pairwise |dt| distribution; distinct encounter windows
(epochs clustered at 60 d).  Output overlap_epochs.json"""
import json, pickle, pathlib, collections
import numpy as np
OUT = pathlib.Path('/Users/mickey/solarsystem/ctoc14/results/s17/campaign')
D = pickle.load(open(OUT / 'reps.pkl', 'rb'))
ep = collections.defaultdict(list)
for tf, asts in zip(D['fl_t'], D['fl_a']):
    for t, a in zip(tf, asts):
        ep[int(a)].append(t / 86400.0)
dts = []; nwin = {}; same = 0; tot = 0
rng = np.random.default_rng(0)
for a, e in ep.items():
    e = np.sort(np.array(e))
    # distinct encounter windows: gaps > 60 d separate windows
    nwin[a] = int(1 + (np.diff(e) > 60).sum()) if len(e) else 0
    if len(e) > 1:
        i, j = np.triu_indices(len(e), 1)
        d = np.abs(e[i] - e[j]); tot += len(d); same += int((d < 60).sum())
        dts.append(rng.choice(d, min(len(d), 2000), replace=False))
dts = np.concatenate(dts)
w = np.array(list(nwin.values()))
res = dict(pair_overlaps=tot, frac_same_encounter_lt60d=round(same / tot, 4),
           dt_days_pct=[int(x) for x in np.percentile(dts, [10, 25, 50, 75, 90])],
           encounter_windows_per_target_pct=[int(x) for x in np.percentile(w, [10, 25, 50, 75, 90])],
           targets_with_1_window=int((w == 1).sum()), targets_with_le2_windows=int((w <= 2).sum()),
           n_targets=len(w))
json.dump(res, open(OUT / 'overlap_epochs.json', 'w'), indent=1); print(res)
