"""H test: grow a hard-rich skeleton route (twin state from run_colgen twin_states) by cheapest insertion of EASY targets.
Per round: candidates = closest approaches (< dmax AU) of every uncovered easy target to the route (w_cands), the best
`--per-round` by distance are tried with insert_homotopy + settle (w_insert) in parallel, the cheapest tank increase is
accepted. Reports the marginal Delta-v per inserted target.
Usage: grow_skel.py twin_state.npz [--rounds 8] [--per-round 24] [--dmax 0.3] [--nproc 4]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from ctoc14.constants import VE

ap = argparse.ArgumentParser(); ap.add_argument('state'); ap.add_argument('--rounds', type=int, default=8)
ap.add_argument('--per-round', type=int, default=24); ap.add_argument('--dmax', type=float, default=0.3); ap.add_argument('--nproc', type=int, default=4)
ap.add_argument('--easy-only', action='store_true', default=True)
a = ap.parse_args()
z = np.load(a.state); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL'])
# hardness: number of dense (>= 35) recosted pool columns containing the target
cnt = np.zeros(301, int)
for line in open(ROOT / 'results/newgen/recost35.jsonl'):
    r = json.loads(line)
    if r['ok']:
        for t in r['targets']: cnt[t] += 1
hard = set(t for t in range(1, 301) if cnt[t] <= 5) | {131, 144}
easy = [t for t in range(1, 301) if t not in hard]
ip = RI.ipr(st); tank0 = float(ip.tank()); dv0 = float(ip.dv())
print(f'skeleton: {len(ip.asts)} targets ({len(set(ip.asts) & hard)} hard), tank {tank0:.1f}, dv {dv0:.2f}; {len(easy)} easy targets available', flush=True)
RI.eph()
with mp.get_context('fork').Pool(a.nproc) as pool:
    for rnd in range(a.rounds):
        tic = time.time()
        have = set(int(x) for x in st['asts'])
        targets = [t for t in easy if t not in have]
        cands = RI.w_cands(('skel', st, targets, a.dmax))
        cands = sorted(cands, key=lambda c: c['dist'])
        # at most 2 events per target, best per-round by distance
        per = {}; sel = []
        for c in cands:
            if per.get(c['ast'], 0) < 2:
                per[c['ast']] = per.get(c['ast'], 0) + 1; sel.append(c)
            if len(sel) >= a.per_round: break
        jobs = [('skel', st, c['ast'], c['t']) for c in sel]
        res = [r for r in pool.map(RI.w_insert, jobs) if r['ok']]
        if not res:
            print(f'round {rnd}: no feasible insertion among {len(sel)} trials'); break
        res.sort(key=lambda r: r['tank'])
        best = res[0]; dtank = best['tank'] - float(RI.ipr(st).tank())
        ddv = VE * np.log(best['tank'] / float(RI.ipr(st).tank()))
        st = best['st']
        print(f'round {rnd}: {len(sel)} trials, {len(res)} feasible; inserted {best["ast"]} at {best["t"]/86400/365.25:.2f} yr: '
              f'tank +{dtank:.1f} kg (+{ddv:.3f} km/s) -> {best["tank"]:.1f}; next best +{VE*np.log(res[1]["tank"]/(best["tank"]-dtank)) if len(res) > 1 else float("nan"):.3f}; '
              f'{len(st["asts"])} targets  ({time.time()-tic:.0f} s)', flush=True)
ipf = RI.ipr(st)
print(f'final: {len(ipf.asts)} targets, tank {ipf.tank():.1f} (J_i {RI.cost(ipf.tank()):.3f}), dv {ipf.dv():.2f} = {ipf.dv()/len(ipf.asts):.3f} per flyby; '
      f'marginal {(ipf.dv()-dv0)/max(1,len(ipf.asts)-len(ip.asts)):.3f} km/s per inserted target')
np.savez(pathlib.Path(a.state).with_name(pathlib.Path(a.state).stem + '_grown.npz'), **st)
