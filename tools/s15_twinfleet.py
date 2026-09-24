"""Stage 15 [T']: assemble and close the 9-craft fleets of finished twin-tail chains (tools/s15_twintail.py).

Fleet of a twin job = the fixed heads (--heads DIR, routes --head-names) + the job's pre-taken tail-1 (job 'tail1', if
any) + every route the chain TOOK.  Reports coverage / sum J_i / misses, then (--close) offers every miss, nearest
first, to every fleet route passing within --close-d AU (s15_isogen.run_close: lin_price screen, twin insertion,
sequential greedy absorption per host, every intermediate route saved as a column under results/s15/cols/<rr>/tclose/).

usage: s15_twinfleet.py OUT JOBDIR [JOBDIR ...] [--heads results/s14/probe_cover/N9/fleet] [--close] [--nproc 8]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)
import numpy as np
import run_ialns as RI
import s15_isogen as IG
from greedy_cover import ALL


def load_st(f):
    z = np.load(ROOT / f); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL'])
    return st


def fleet_files(jobdir, heads, head_names):
    meta = json.load(open(ROOT / jobdir / 'meta.json'))
    if meta['job'].get('fleet_pre'):                 # hybrid chains: the job carries its whole non-twin fleet
        files = {f'f{i + 1}': f for i, f in enumerate(meta['job']['fleet_pre'])}
        for t in meta['chain'].get('taken', []):
            files[f's{t["k"]}'] = t['f']
        return files, meta
    files = {f'h{n}': str(pathlib.Path(heads) / f'route_{n}.npz') for n in head_names}
    t1 = meta['job'].get('tail1')
    if t1:
        files['t1'] = t1
    for i, f in enumerate(meta['job'].get('pre') or []):
        files[f'p{i + 1}'] = f
    for t in meta['chain'].get('taken', []):
        files[f's{t["k"]}'] = t['f']
    return files, meta


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('jobs', nargs='+')
    ap.add_argument('--heads', default='results/s14/probe_cover/N9/fleet'); ap.add_argument('--head-names', default='01,02,03,04,05,06')
    ap.add_argument('--close', action='store_true'); ap.add_argument('--close-d', type=float, default=0.25)
    ap.add_argument('--max-abs', type=int, default=4); ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--tag', default='t02')
    a = ap.parse_args()
    if a.nproc > 8:
        raise SystemExit('CPU budget: --nproc <= 8')
    out = ROOT / a.out; out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt')
    RI.eph()
    hn = [x for x in a.head_names.split(',') if x]
    fleets = {}
    for jd in a.jobs:
        if not (ROOT / jd / 'meta.json').exists():
            say(f'{jd}: not finished'); continue
        files, meta = fleet_files(jd, a.heads, hn)
        F = RI.IFleet()
        for k, f in files.items():
            st = load_st(f); F.routes[k] = dict(st=st, tank=float(RI.ipr(st).tank()))
        cov = set(F.coverage()); M = sorted(set(ALL) - cov)
        sj = sum(RI.cost(r['tank']) for r in F.routes.values())
        dup = sum(len(h) - 1 for h in F.coverage().values())
        nm = pathlib.Path(jd).name
        say(f'fleet {nm}: {len(F.routes)} craft, covered {len(cov)} (dup {dup}), sum J_i {sj:.4f}; misses {len(M)} {M}; '
            + ' '.join(f'{k}:{len(r["st"]["asts"])}@{r["tank"]:.0f}' for k, r in F.routes.items()))
        F.save(out / 'fleets' / nm, note=f'twin fleet {nm}')
        fleets[nm] = dict(files=files, covered=len(cov), sumJi=round(sj, 4), misses=M, craft=len(F.routes))
    rep = dict(fleets=fleets)
    if a.close:
        jobs = []; seen = {}
        for nm, fl in fleets.items():
            if not fl['misses']:
                continue
            scan = [(k, load_st(f), fl['misses'], a.close_d) for k, f in fl['files'].items()]
            with mp.get_context('fork').Pool(min(a.nproc, len(scan))) as pool:
                R = pool.map(RI.w_cands, scan)
            fl['approach'] = {}
            for (k, st, _, _), c in zip(scan, R):
                tf = np.asarray(st['tf']); d = {}
                for x in c:
                    if np.abs(tf - x['t']).min() > 10 * 86400.0 and (x['ast'] not in d or x['dist'] < d[x['ast']]):
                        d[x['ast']] = x['dist']
                fl['approach'][k] = {int(X): round(v, 4) for X, v in d.items()}
                Xs = tuple(X for X, _ in sorted(d.items(), key=lambda kv: kv[1])[:a.max_abs + 1])
                if not Xs:
                    continue
                key = (fl['files'][k], Xs)
                if key in seen:
                    seen[key]['fleets'].append(nm); continue
                tag = f'{nm}_{k}'
                j = dict(kind='close', tag=tag, out=f'results/s15/cols/{a.tag}/tclose/{tag}', items=[(fl['files'][k], list(Xs))],
                         absorb_dmax=a.close_d, absorb_kg=150.0, item_trials=3, max_abs=a.max_abs, wall_budget=2400.0, fleets=[nm])
                seen[key] = j; jobs.append(j)
            near = sorted({X for k in fl['approach'] for X in fl['approach'][k]})
            say(f'fleet {nm}: misses with an approach <= {a.close_d} AU: {len(near)} of {len(fl["misses"])}')
        say(f'{len(jobs)} closer host jobs')
        if jobs:
            with mp.get_context('fork').Pool(min(a.nproc, len(jobs)), maxtasksperchild=1) as pool:
                res = pool.map(IG.run_job, [{k: v for k, v in j.items() if k != 'fleets'} for j in jobs], chunksize=1)
            for j, r in zip(jobs, res):
                say(f'  {j["tag"]}: ' + '; '.join(f'tried {h["tried"]} got {h["got"]} (+{h["dkg"]} kg)' for h in r.get('hosts', []))
                    + (f' ERROR {r.get("error", "")[-120:]}' if r.get('error') else ''))
                for h in r.get('hosts', []):
                    for nm in j['fleets']:
                        fleets[nm].setdefault('closed', []).append(dict(host=j['tag'].split('_')[-1], got=h['got'], dkg=h['dkg']))
        for nm, fl in fleets.items():
            got = sorted({X for h in fl.get('closed', []) for X in h['got']})
            fl['placed'] = got; fl['unplaced'] = sorted(set(fl['misses']) - set(got))
            say(f'fleet {nm}: {fl["covered"]} covered @ {fl["sumJi"]:.4f}; closer placed {got}; never placed {fl["unplaced"]}')
    IG.jdump(rep, out / 'twinfleet.json')


if __name__ == '__main__':
    main()
