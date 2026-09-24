"""Reviewer check: contract / honesty / horizon invariants of every planbeam RESULT under a test run.
usage: check_contract.py RUN_DIR [RUN_DIR ...]   (each with slices/sKK/jobs/*.json and <job out>/result.json)
Prints one line per column with every violated invariant, then a summary.  Optionally --regrid N: re-settle N columns
(cheapest per n_new of the first craft job) with C.regrid_settle and report the tank delta."""
import os, sys, json, pathlib, time
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_v, '1')
sys.path.insert(0, '/Users/mickey/solarsystem/ctoc14/tools'); sys.path.insert(0, '/Users/mickey/solarsystem/ctoc14')
import numpy as np
import s18_common as C
from ctoc14.constants import DAY

args = [a for a in sys.argv[1:] if not a.startswith('--')]
n_regrid = 0
if '--regrid' in sys.argv:
    n_regrid = int(sys.argv[sys.argv.index('--regrid') + 1])
W = C.load_windows()
bad = 0; n_cols = 0; regrid_todo = []
for run in args:
    rd = C.run_dir(run)
    for jf in sorted(rd.glob('slices/s*/jobs/s*_*.json')):
        job = C.jload(jf)
        rf = rd / job['out'] / 'result.json'
        if not rf.exists():
            print(f'{jf.name}: no result.json'); continue
        R = C.jload(rf)
        T, H, L = float(job['T']), float(job['H']), float(job['L'])
        rem = set(int(t) for t in job['remaining']); quota = job.get('quota') or {}
        print(f'== {run} {R["job_id"]} kind {R["kind"]} ok {R["ok"]} columns {len(R["columns"])} retried {R.get("retried")} '
              f'wall {R.get("stats", {}).get("wall_s")} end {R.get("stats", {}).get("end")} error {(R.get("error") or "")[-120:]!r}')
        if not R['ok']:
            continue
        for k in ('job_id', 'kind', 'craft', 'ok', 'error', 'columns', 'stats'):
            assert k in R, f'RESULT missing {k}'
        st = R['stats']
        for k in ('levels', 'n_states', 'n_columns', 'settles', 'settle_ok', 'price_s', 'settle_s', 'wall_s', 'end', 'hist'):
            if k not in st:
                print(f'  stats missing {k}'); bad += 1
        prefix = []
        if job['kind'] == 'craft':
            pst = C.load_twin(rd / job['route']); o = C.order_of(pst)
            prefix = [int(pst['asts'][i]) for i in o]
            tank_prev_twin = float(C.twin_of(pst).tank())
            if abs(tank_prev_twin - float(job['tank_prev'])) > 0.5:
                print(f'  prefix twin tank {tank_prev_twin:.1f} != job tank_prev {job["tank_prev"]:.1f}')
        ids = set(); sets = set()
        for c in R['columns']:
            n_cols += 1; errs = []
            try:
                C.check_column(c, T, H)
            except AssertionError as e:
                errs.append(f'check_column: {e}')
            for k in ('id', 'craft', 'kind', 'targets', 'epochs', 'all_targets', 'n_new', 'npz', 'tank', 'tank_prev', 'dtank',
                      'miss', 't_launch', 't_end', 'potential', 'counts', 'depth', 'score', 'job_id'):
                if k not in c:
                    errs.append(f'missing {k}')
            if c['id'] in ids: errs.append('duplicate id')
            ids.add(c['id'])
            key = frozenset(c['targets'])
            if key in sets: errs.append('duplicate target set')
            sets.add(key)
            if not set(c['targets']) <= rem: errs.append(f'targets not in remaining: {sorted(set(c["targets"]) - rem)}')
            if c['kind'] != job['kind'] or c['craft'] != job.get('craft'): errs.append('kind/craft mismatch')
            if c['job_id'] != job['job_id']: errs.append('job_id mismatch')
            if abs(float(c['tank_prev']) - float(job['tank_prev'])) > 1e-6: errs.append('tank_prev != job tank_prev')
            if abs(float(c['dtank']) - (float(c['tank']) - float(c['tank_prev']))) > 1e-6: errs.append('dtank mismatch')
            for cl, n in (quota or {}).items():
                if c['counts'].get(cl, 0) > n: errs.append(f'quota {cl}: {c["counts"][cl]} > head-room {n}')
            if c['counts'] != C.count_classes(W, c['targets']): errs.append('counts mismatch')
            if c['all_targets'] != prefix + list(c['targets']): errs.append(f'all_targets != prefix + targets')
            if c['epochs'] and abs(float(c['t_end']) - max(c['epochs'])) > 1.0: errs.append(f't_end {c["t_end"]} != last epoch')
            if job['kind'] == 'root':
                if not (T - 1e-6 <= float(c['t_launch']) <= T + H + 1e-6): errs.append('root launch outside [T, T+H]')
                if job.get('grid') and not (float(job['grid'][0]) - 1e-6 <= float(c['t_launch']) <= float(job['grid'][1]) + 1e-6):
                    errs.append('root launch outside the job grid')
            # the twin itself
            try:
                tw = C.load_twin(rd / c['npz']); ip = C.twin_of(tw); o = C.order_of(tw)
                asts = [int(tw['asts'][i]) for i in o]; tfs = [float(tw['tf'][i]) for i in o]
                if asts != c['all_targets']: errs.append(f'twin asts {asts} != all_targets')
                new_e = [e for a, e in zip(asts, tfs) if a not in set(prefix)]
                if any(abs(x - y) > 1.0 for x, y in zip(new_e, c['epochs'])): errs.append('twin epochs != column epochs')
                if any(not (T < e <= T + H + 1e-6) for e in new_e): errs.append('twin new epoch outside (T, T+H]')
                if len(new_e) != c['n_new']: errs.append('lookahead flyby in twin')
                miss = C.miss_of(ip)
                if miss > C.TOL_MISS: errs.append(f'twin miss {miss:.0f} km > 150')
                if abs(miss - float(c['miss'])) > 5: errs.append(f'twin miss {miss:.1f} != column miss {c["miss"]:.1f}')
                if abs(float(ip.tank()) - float(c['tank'])) > 0.05: errs.append(f'twin tank {ip.tank():.2f} != column tank {c["tank"]:.2f}')
                hc = C.honest_check(ip)
                if not hc['ok']: errs.append(f'honest_check fails: irregular {hc["irregular"]} cap {hc["max_window_cap"]} miss {hc["miss_km"]}')
                if abs(float(tw['tL']) - float(c['t_launch'])) > 1.0: errs.append('twin tL != t_launch')
                pl = st.get('plan', {}).get(c['id'], {})
                if pl and c['potential'] > pl.get('plan_n', 99) - c['n_new']: errs.append('potential > planned flybys beyond commit')
            except Exception as e:
                errs.append(f'twin load: {e!r}')
            tag = 'OK ' if not errs else 'BAD'
            bad += bool(errs)
            print(f'  {tag} {c["id"]:12s} n {c["n_new"]} dtank {c["dtank"]:6.1f} pot {c["potential"]} miss {c["miss"]:5.1f} '
                  f'epochs {[round(e / DAY) for e in c["epochs"]]} tL {c["t_launch"] / DAY:.0f} d' + ('  ' + '; '.join(errs) if errs else ''))
            if n_regrid and job['kind'] == 'craft':
                regrid_todo.append((rd, c))
print(f'SUMMARY: {n_cols} columns, {bad} with violations')
if n_regrid and regrid_todo:
    # cheapest column per n_new across the collected craft columns, up to n_regrid
    best = {}
    for rd, c in regrid_todo:
        if c['n_new'] not in best or c['dtank'] < best[c['n_new']][1]['dtank']:
            best[c['n_new']] = (rd, c)
    picks = [best[n] for n in sorted(best, reverse=True)][:n_regrid]
    for rd, c in picks:
        tw = C.load_twin(rd / c['npz']); ip = C.twin_of(tw); tic = time.time()
        r = C.regrid_settle(ip)
        print(f'REGRID {c["id"]} n {c["n_new"]}: column tank {c["tank"]:.2f} -> honest {r["tank"]:.2f} (delta {r["tank"] - c["tank"]:+.2f} kg), '
              f'ok {r["ok"]} miss {r["miss"]:.1f} km irregular {r["honest"].get("irregular")} cap {r["honest"].get("max_window_cap")} '
              f'({time.time() - tic:.1f} s)')
