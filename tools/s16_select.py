"""Stage 16: exact FRONTIER selection over the HONEST pool (results/s16/honest), parallel over (N, K), with
alternative solutions per K and an ESTIMATED CLOSED cost per solution.

Why.  Every fleet selection of stages 13-15 ranked routes on optimistic planner-twin tanks (per-impulse thrust cap,
irregular nodes).  tools/s16_honest_pool.py re-priced every route (regrid + settle); this selector re-runs the
s15_select frontier MILP on those honest prices:
    min sum J_i  s.t.  sum y <= N,  coverage >= K   (s15_select.frontier / s15b_select.solve formulation)
Alternatives: K < 298 -> no-good cut on the miss set (each alternative leaves a DIFFERENT miss set for the closer);
K = 298 -> no-good cut on the pick set.  Each solution's misses are priced into its routes (s15b_close._price_host,
lin_price screen, one miss per host, greedy) -> est = sum J_i + sum lin dJ ('None' when a miss has no approach).

Columns: s13_master.load (miss <= 150 km, tank <= 1400), dedupe by target set (lightest), dominance pruning -- exactly
s15_select.load_columns, with a (path, mtime, size) cache in results/s16/colcache.pkl.

Outputs OUT/frontier.json  {"N9K298": {covered, sumJi, misses, picks, depths, tanks, est, plan}, "N9K298a1": ...}
        (usable as FLEET spec OUT/frontier.json:N9K298 for tools/s15b_fleet / s15b_close / s15b_relocate),
        OUT/log.txt, and with --swap FLEET a 1-swap table of that fleet over the pool (OUT/swap.json).

usage: s16_select.py OUT [--pats P ...] --jobs 9:292,9:293,...,8:282 [--alt 3] [--ftlim 900] [--nproc 10]
                     [--dmax 0.2] [--swap results/s15b/reloc_h2/fleet]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, glob, time, pickle, fnmatch, pathlib, argparse, collections, multiprocessing as mp
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import csr_matrix, hstack, identity
import run_ialns as RI
import s13_master as SM
import s15_select as SS
import s15b_fleet as FA
import s15b_close as CL

ALL = SS.ALL; TI = SS.TI; NT = SS.NT
PATS = ['results/s16/honest/route_*.npz', 'results/s16/**/cols/route_*.npz']
CACHE = ROOT / 'results/s16/colcache.pkl'
COLS = None


def load_cached(pats, nproc, say, max_tank=1400.0, exclude=('*/relocated/*',)):
    files = sorted(set(f for p in pats for f in glob.glob(p, recursive=True)
                       if not any(fnmatch.fnmatch(f, e) for e in exclude)))
    try:
        cache = pickle.load(open(CACHE, 'rb'))
    except Exception:
        cache = {}
    keys = {}
    for f in files:
        try:
            s = os.stat(f); keys[f] = (f, s.st_mtime, s.st_size)
        except OSError:
            pass
    files = [f for f in files if f in keys]
    todo = [f for f in files if keys[f] not in cache]
    if todo:
        with mp.get_context('fork').Pool(nproc) as pool:
            res = pool.map(SM.load, [(f, max_tank) for f in todo], chunksize=4)
        for f, r in zip(todo, res):
            cache[keys[f]] = r
        tmp = CACHE.with_suffix(f'.tmp{os.getpid()}'); pickle.dump(cache, open(tmp, 'wb')); tmp.replace(CACHE)
    res = [cache[keys[f]] for f in files]
    why = collections.Counter(r for _, r in res)
    best = {}
    for c, _ in res:
        if c is None:
            continue
        k = frozenset(c['asts'])
        if k not in best or c['tank'] < best[k]['tank']:
            best[k] = dict(c)
    cols = sorted(best.values(), key=lambda x: (-len(x['asts']), x['tank']))
    for c in cols:
        c['cost'] = RI.cost(c['tank']); c['set'] = frozenset(c['asts'])
    by_t = collections.defaultdict(list)
    for i, c in enumerate(cols):
        for t in c['asts']:
            by_t[t].append(i)
    keep = []; dom = 0
    for j, c2 in enumerate(cols):
        t0 = min(c2['asts'], key=lambda t: len(by_t[t]))
        if any(i != j and cols[i]['cost'] <= c2['cost'] + 1e-12 and len(cols[i]['set']) >= len(c2['set'])
               and c2['set'] <= cols[i]['set'] and (cols[i]['cost'] < c2['cost'] - 1e-12 or len(cols[i]['set']) > len(c2['set']) or i < j)
               for i in by_t[t0]):
            dom += 1; continue
        keep.append(c2)
    say(f'load: {len(files)} files ({len(todo)} new; {dict(why)}), {len(cols)} distinct, {dom} dominated -> {len(keep)} columns')
    return keep


def solve(cols, N, K, cuts_miss, cuts_pick, tlim):
    A, c = SS.matrices(cols); nr = len(cols)
    Mcov = hstack([A, -identity(NT, format='csr')]).tocsr()
    cons = [LinearConstraint(Mcov, np.zeros(NT), np.full(NT, np.inf)),
            LinearConstraint(csr_matrix(np.concatenate([np.ones(nr), np.zeros(NT)])[None, :]), -np.inf, N),
            LinearConstraint(csr_matrix(np.concatenate([np.zeros(nr), np.ones(NT)])[None, :]), K, np.inf)]
    for ms in cuts_miss:                                  # at least one of this miss set covered
        row = np.zeros(nr + NT)
        for t in ms:
            row[nr + TI[t]] = 1.0
        cons.append(LinearConstraint(csr_matrix(row[None, :]), 1, np.inf))
    for pk in cuts_pick:                                  # not exactly this pick set again
        row = np.zeros(nr + NT)
        for j in pk:
            row[j] = 1.0
        cons.append(LinearConstraint(csr_matrix(row[None, :]), -np.inf, len(pk) - 1))
    r = milp(np.concatenate([c, np.zeros(NT)]), constraints=cons, integrality=np.ones(nr + NT), bounds=Bounds(0, 1),
             options=dict(time_limit=tlim, mip_rel_gap=1e-5))
    if r.x is None:
        return None, r.message
    return [j for j in range(nr) if r.x[j] > 0.5], r.message


def w_solve(job):
    N, K, nalt, tlim = job
    cols = COLS; out = []; cm = []; cp = []; t0 = time.time()
    for i in range(nalt + 1):
        pick, msg = solve(cols, N, K, cm, cp, tlim)
        if pick is None:
            out.append(dict(N=N, K=K, i=i, fail=msg[:80], sec=round(time.time() - t0))); break
        J, cv, sj = SS.evaluate(cols, pick)
        miss = sorted(set(ALL) - set().union(*[cols[j]['set'] for j in pick]))
        srt = sorted(pick, key=lambda j: (-len(cols[j]['asts']), cols[j]['tank']))
        out.append(dict(N=N, K=K, i=i, covered=cv, sumJi=round(sj, 5), misses=miss, msg=msg[:40],
                        picks=[cols[j]['f'] for j in srt], depths=[len(cols[j]['asts']) for j in srt],
                        tanks=[round(cols[j]['tank'], 2) for j in srt], sec=round(time.time() - t0)))
        if miss and K < 298:
            cm.append(miss)
        else:
            cp.append(pick)
    return out


def estimate(rows, nproc, dmax, say):
    jobs = []; idx = []
    for ri, r in enumerate(rows):
        if r.get('fail') or not r['misses']:
            continue
        for h, f in enumerate(r['picks']):
            jobs.append((f'{h + 1:02d}', FA.load_st(f), r['misses'], dmax)); idx.append(ri)
    if jobs:
        with mp.get_context('fork').Pool(nproc, maxtasksperchild=8) as pool:
            P = pool.map(CL._price_host, jobs, chunksize=1)
    else:
        P = []
    per = collections.defaultdict(list)
    for ri, lst in zip(idx, P):
        per[ri] += lst
    for ri, r in enumerate(rows):
        if r.get('fail'):
            continue
        if not r['misses']:
            r['est'] = r['sumJi']; r['plan'] = []; continue
        pr = sorted(per[ri], key=lambda c: c['lin_dJ'])
        used = set(); est = r['sumJi']; plan = []
        for X in r['misses']:
            cs = [c for c in pr if c['ast'] == X and c['host'] not in used]
            if not cs:
                est = float('inf'); plan.append((X, None, None, None)); continue
            c = cs[0]; used.add(c['host']); est += c['lin_dJ']
            plan.append((X, c['host'], round(c['dist'], 3), round(c['lin_dJ'], 4)))
        r['est'] = round(est, 5) if np.isfinite(est) else None; r['plan'] = plan


def swap_table(cols, spec, say, top=5):
    files = FA.fleet_files(spec)
    sts = {k: FA.load_st(f) for k, f in files.items()}
    sets = {k: frozenset(int(a) for a in s['asts']) - {131, 144} for k, s in sts.items()}
    tanks = {k: float(RI.ipr(s).tank()) for k, s in sts.items()}
    base = sum(RI.cost(t) for t in tanks.values())
    say(f'swap table of {spec}: sum J_i {base:.4f}, covered {len(set().union(*sets.values()))}')
    rep = {}
    for k in sorted(sets):
        rest = set().union(*[sets[q] for q in sets if q != k])
        U = set(ALL) - rest
        cand = []
        for c in cols:
            if c['set'] == sets[k]:
                continue
            left = len(U - c['set'])
            if left > 2:
                continue
            d = c['cost'] - RI.cost(tanks[k])
            cand.append((left, d, c))
        cand.sort(key=lambda z: (z[0], z[1]))
        best0 = [z for z in cand if z[0] == 0][:top]; best1 = [z for z in cand if z[0] == 1][:top]
        rep[k] = dict(n=len(sets[k]), tank=round(tanks[k], 2), unique=len(U),
                      swaps0=[dict(dJ=round(d, 4), n=len(c['asts']), tank=round(c['tank'], 1), f=c['f']) for _, d, c in best0],
                      swaps1=[dict(dJ=round(d, 4), n=len(c['asts']), tank=round(c['tank'], 1), f=c['f'],
                                   left=sorted(U - c['set'])) for _, d, c in best1])
        say(f'  {k}: {len(sets[k])} fb @ {tanks[k]:.1f}, {len(U)} unique; full swaps ' +
            (', '.join(f'{d:+.4f} ({len(c["asts"])}fb)' for _, d, c in best0) or 'none') + '; 1-left swaps ' +
            (', '.join(f'{d:+.4f} ({len(c["asts"])}fb, miss {sorted(U - c["set"])})' for _, d, c in best1) or 'none'))
    return rep


def main():
    global COLS
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--pats', nargs='*', default=PATS)
    ap.add_argument('--jobs', default='9:292,9:293,9:294,9:295,9:296,9:297,9:298')
    ap.add_argument('--alt', type=int, default=0); ap.add_argument('--ftlim', type=float, default=900.0)
    ap.add_argument('--nproc', type=int, default=10); ap.add_argument('--dmax', type=float, default=0.2)
    ap.add_argument('--swap', default=''); ap.add_argument('--no-est', action='store_true')
    ap.add_argument('--exclude-src', default='', help='comma fnmatch patterns on the ORIGINAL path (honest index src) of '
                    'results/s16/honest columns to leave out, e.g. "results/s15b/reloc*,results/s15b/honest/*"')
    a = ap.parse_args()
    out = ROOT / a.out; out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt'); RI.eph()
    tic = time.time()
    COLS = load_cached(a.pats, a.nproc, say)
    if a.exclude_src:
        pats = [x for x in a.exclude_src.split(',') if x]; src = {}
        for l in open(ROOT / 'results/s16/honest/index.jsonl'):
            try:
                r = json.loads(l); src[f'results/s16/honest/route_{r["key"]}.npz'] = r['src']
            except Exception:
                pass
        n0 = len(COLS)
        COLS = [c for c in COLS if not any(fnmatch.fnmatch(src.get(c['f'], ''), q) for q in pats)]
        say(f'exclude-src {pats}: {n0} -> {len(COLS)} columns')
    if a.swap:
        json.dump(swap_table(COLS, a.swap, say), open(out / 'swap.json', 'w'), indent=1)
    jobs = [(int(x.split(':')[0]), int(x.split(':')[1]), a.alt, a.ftlim) for x in a.jobs.split(',') if x]
    rows = []
    if jobs:
        with mp.get_context('fork').Pool(min(a.nproc, len(jobs))) as pool:
            for lst in pool.imap_unordered(w_solve, jobs):
                for r in lst:
                    if r.get('fail'):
                        say(f'  N<={r["N"]} K>={r["K"]} alt {r["i"]}: {r["fail"]} ({r["sec"]} s)')
                    else:
                        say(f'  N<={r["N"]} K>={r["K"]} alt {r["i"]}: {len(r["picks"])} craft {r["covered"]} covered '
                            f'sum J_i {r["sumJi"]:.4f} ({r["msg"][:24]}, {r["sec"]} s); misses {r["misses"]}; depths {r["depths"]}')
                rows += lst
    if not a.no_est:
        estimate(rows, a.nproc, a.dmax, say)
    fr = {}
    for r in rows:
        if r.get('fail'):
            continue
        fr[f'N{r["N"]}K{r["K"]}' + (f'a{r["i"]}' if r['i'] else '')] = r
    old = {}
    if (out / 'frontier.json').exists():
        try:
            old = json.load(open(out / 'frontier.json'))
        except Exception:
            old = {}
    old.update(fr)
    json.dump(old, open(out / 'frontier.json', 'w'), indent=1)
    ok = sorted([r for r in rows if not r.get('fail')], key=lambda r: (r.get('est') is None, r.get('est') or 99))
    say('ranking by est closed: ' + ' | '.join(f'N{r["N"]}K{r["K"]}a{r["i"]} {r["sumJi"]:.4f}->{r.get("est")}' for r in ok[:20])
        + f' ({time.time() - tic:.0f} s)')


if __name__ == '__main__':
    main()
