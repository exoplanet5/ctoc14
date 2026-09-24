"""Stage 15 [L]: the miss-driven column loop (resumable).  Each round:
  (a) build a batch of column jobs (tools/s15_isogen.py) and run them on a fork Pool(--nproc <= 8), longest first:
        F  full-pool (298) heads steered to S (subsets of the CURRENT miss set of the best N<=9 fleet, grouped by
           archive co-occurrence = epoch compatibility), premiums {0.005, 0.01, 0.015, 0.02}, member H, seeds;
        E  full-pool pure depth (S empty);
        C  chains = tail columns on RESIDUAL pools: all 298 minus a head subset of the best fleet (and, from round 2 on,
           minus the best new steered head + the best-fleet routes disjoint from it), S = misses in the pool;
        A  archive absorption: twin insertion (lin_price screened) of each miss into existing beam-born columns that
           pass within --absorb-d AU of its window;
  (b) tools/s13_master.py over the archive + every s15 column with --N 8,9,10 (exact selection, LP bounds);
  (c) the new best N<=9 fleet's misses become the next round's targets;
  (d) log per round: best N<=9 covered / sum J_i / misses / LP bound (and N8, N10), columns generated, wall time.
Stops when the best N<=9 covers all 298 at sum J_i <= --target (11.90), after --rounds rounds, or on results/s15/STOP.
State: results/s15/loop.json (+ every job's meta.json / chain.json); a restart continues where it stopped.

usage: s15_loop.py [--rounds 1] [--nproc 8] [--start results/s14/probe_cover/N9/fleet] [--target 11.90]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, glob, time, pathlib, argparse, subprocess, itertools, collections, multiprocessing as mp
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)
import numpy as np
import s15_isogen as IG
import run_ialns as RI
from greedy_cover import ALL

OUT = ROOT / 'results/s15'
ARCHIVE = ['results/s12/cover_clean/fleet/route_*.npz', 'results/newgen/**/route_*.npz', 'results/s13/**/route_*.npz',
           'results/s14/**/route_*.npz']
S15COLS = 'results/s15/cols/**/route_*.npz'
PREMS = [0.005, 0.01, 0.015, 0.02]
BAR = dict(covered=290, sumJi=11.7045, J=21.7045, lp=14.5935)          # results/s14/probe_cover N<=9


def say(s):
    line = f'[{time.strftime("%m-%d %H:%M:%S")}] {s}'
    print(line, flush=True)
    with open(OUT / 'loop_log.txt', 'a') as f:
        f.write(line + '\n')


def jload(p, default=None):
    p = pathlib.Path(p)
    return json.load(open(p)) if p.exists() else default


# ------------------------------------------------------------------------------------------------ archive helpers
def distinct_columns(pats):
    files = sorted(set(f for p in pats for f in glob.glob(str(p), recursive=True) if '/relocated/' not in f))
    seen = {}
    for f in files:
        try:
            z = np.load(f); k = frozenset(int(a) for a in z['asts'])
        except Exception:
            continue
        if k and k not in seen:
            seen[k] = f
    return seen


def s_groups(M, cols, k_max=8):
    """Premium sets from the miss set M: triples/quads/pairs of M that archive routes already fly together
    (co-occurrence = compatible encounter epochs for one craft), ranked by support, plus singletons."""
    M = sorted(M); Ms = set(M)
    sup = collections.Counter()
    for k in cols:
        s = tuple(sorted(Ms & set(k)))
        for r in (2, 3, 4):
            for c in itertools.combinations(s, r):
                sup[c] += 1
    groups = []
    for r, need in ((3, 2), (4, 2), (2, 4)):
        for c, n in sorted(((c, n) for c, n in sup.items() if len(c) == r and n >= need), key=lambda x: -x[1]):
            groups.append((c, n))
    out = []; pair_used = collections.Counter(); tgt_used = collections.Counter()
    for c, n in groups:
        if any(pair_used[p] >= 2 for p in itertools.combinations(c, 2)):
            continue
        out.append(list(c))
        for p in itertools.combinations(c, 2):
            pair_used[p] += 1
        for t in c:
            tgt_used[t] += 1
        if len(out) >= k_max - 1:
            break
    for t in M:                                  # every miss in >= 1 set: singletons for the uncovered ones
        if tgt_used[t] == 0 and len(out) < k_max + 4:
            out.append([t])
    if not out:
        out = [[t] for t in M]
    return out[:max(k_max, len([t for t in M if tgt_used[t] == 0]))]


def fleet_info(d):
    F = RI.IFleet(d)
    names = sorted(F.routes, key=lambda n: (-len(F.routes[n]['st']['asts']), F.routes[n]['tank']))
    cov = set(F.coverage())
    return F, names, cov, sorted(set(ALL) - cov)


# ------------------------------------------------------------------------------------------------------ round jobs
def build_round(r, best, args, M_hard=None, grown=None):
    """Job list of round r from the current best N<=9 fleet (dir), the premium set M_hard (default: the fleet's
    misses) and, from round 2 on, the grown heads (best-fleet route name -> grown column file)."""
    F, names, cov, M = fleet_info(best['fleet'])
    MH = sorted(set(M_hard) if M_hard else set(M))
    grown = grown or {}
    cols = distinct_columns(ARCHIVE + [S15COLS])
    base = f'results/s15/cols/r{r:02d}'
    jobs = []
    def head_file(n):
        return grown.get(n) or str(pathlib.Path(best['fleet']) / f'route_{n}.npz')
    def file_targets(f):
        z = np.load(ROOT / f if not str(f).startswith('/') else f); return set(int(a) for a in z['asts'])
    # ---- F: steered heads over the full pool
    groups = args.round1_groups if (r == 1 and args.round1_groups) else s_groups(MH, cols.keys(), k_max=args.groups)
    prems = PREMS if r == 1 else [0.01, 0.02]
    for gi, S in enumerate(groups):
        for pi, p in enumerate(prems):
            tag = f'F{gi:02d}p{int(round(p * 1000)):02d}'
            jobs.append(dict(kind='column', tag=tag, out=f'{base}/{tag}', S=S, p=p, seed=100 * r + 10 * gi + pi,
                             member='H', beam=args.beam, absorb=MH, wall_budget=2400.0, _est=9.0))
    # ---- E: pure depth
    for e in range(args.pure if r == 1 else 0):
        tag = f'E{e:02d}'
        jobs.append(dict(kind='column', tag=tag, out=f'{base}/{tag}', S=[], p=0.0, seed=1000 * r + e, member='H',
                         beam=args.beam, absorb=MH, wall_budget=2400.0, _est=8.0))
    # ---- C: chains on residual pools of head subsets of the best fleet (grown heads from round 2 on)
    heads = names[:6]
    subsets = [(heads[:6], 0.02), (heads[:6], 0.01), (heads[:5], 0.015), (heads[:4] + heads[5:6], 0.015),
               (heads[:4], 0.01), (heads[:3], 0.02), (heads[:3] + heads[4:5], 0.015), (heads[:2] + heads[3:5], 0.015)]
    for ci, (H, p) in enumerate(subsets[:args.chains]):
        hf = [head_file(n) for n in H]
        hc = set()
        for f in hf:
            hc |= file_targets(f)
        pool = sorted(set(ALL) - hc)
        steps = 9 - len(H)
        tag = f'C{ci:02d}h{len(H)}'
        jobs.append(dict(kind='chain', tag=tag, out=f'{base}/{tag}', pool=pool, S=sorted(set(MH) & set(pool)), p=p,
                         seed=500 * r + ci, steps=steps, members=['H', 'T2'], members_small=['T2', 'T3', 'T6'],
                         beam=args.beam, absorb=None, joint_last=True, chain_cap=1150.0,
                         step_grow=args.step_grow if r > 1 else 0, step_grow_d=0.05, step_grow_kg=25.0,
                         heads=hf, wall_budget=steps * 1500.0, _est=steps * (12.0 if len(pool) < 150 else 14.0)))
    # ---- K: complement chains of the best NEW steered heads of the previous round (round >= 2)
    if r > 1:
        cand = []
        for m in sorted(glob.glob(f'results/s15/cols/r{r - 1:02d}/F*/meta.json')):
            for c in (jload(m) or {}).get('cols', []):
                if c['n'] >= 35 and c.get('tank', 9e9) <= 1100 and set(c['targets']) & set(MH):
                    cand.append((len(set(c['targets']) & set(MH)), c['n'], -c['tank'], c))
        cand.sort(key=lambda x: x[:3], reverse=True)
        used = set()
        for q in cand:
            c = q[3]; key = frozenset(c['targets'])
            if key in used:
                continue
            used.add(key)
            H = [c['f']]; hc = set(c['targets'])
            for n in names:                      # best-fleet routes nearly disjoint from the new head
                f = head_file(n); t = file_targets(f)
                if len(t & hc) <= 2 and len(H) < 6:
                    H.append(f); hc |= t
            pool = sorted(set(ALL) - hc)
            steps = 9 - len(H)
            if steps < 1 or len(pool) < 10:
                continue
            tag = f'K{len(used) - 1:02d}h{len(H)}'
            jobs.append(dict(kind='chain', tag=tag, out=f'{base}/{tag}', pool=pool, S=sorted(set(MH) & set(pool)), p=0.015,
                             seed=700 * r + len(used), steps=steps, members=['H', 'T2'], members_small=['T2', 'T3', 'T6'],
                             joint_last=True, chain_cap=1150.0, step_grow=args.step_grow, step_grow_d=0.05, step_grow_kg=25.0,
                             beam=args.beam, heads=H, wall_budget=steps * 1500.0, _est=steps * 12.0))
            if len(used) >= args.complements:
                break
    # ---- A: archive absorption of each hard target into columns passing close to its window
    ap = scan_approaches(cols, MH, args)
    for X in MH:
        items = [(f, X) for d, n, f in ap.get(X, [])[:args.absorb_items]]
        if items:
            tag = f'A{X:03d}'
            jobs.append(dict(kind='absorb', tag=tag, out=f'{base}/{tag}', items=items, item_trials=1,
                             absorb_dmax=args.absorb_d, absorb_kg=80.0, wall_budget=1800.0, _est=len(items) * 0.8))
    return jobs


def grow_jobs(r, best, args):
    """Round >= 2 pre-phase: grow each of the best fleet's 6 heads by cheap nearby targets of the residual of the
    6 heads (twin insertion into beam-born hosts, <= 0.05 AU, <= 25 kg each, <= 4)."""
    F, names, cov, M = fleet_info(best['fleet'])
    heads = names[:6]; hc = set()
    for n in heads:
        hc |= set(int(a) for a in F.routes[n]['st']['asts'])
    R = sorted(set(ALL) - hc)
    return [dict(kind='grow', tag=f'G{n}', out=f'results/s15/cols/r{r:02d}/grow/G{n}',
                 items=[(str(pathlib.Path(best['fleet']) / f'route_{n}.npz'), R)], grow_d=0.05, grow_kg=25.0, grow_n=4,
                 wall_budget=1500.0, _est=8.0, head=n) for n in heads]


def hard_set_next(r, L, args):
    """Premium set for round r+1: misses of the round's N<=9 pick + misses of every candidate chain fleet
    (>= args.min_cov covered, from the round's fleet closer report)."""
    MH = set()
    R = L['rounds'][f'{r:02d}']
    n9 = (R.get('N9') or {}).get('fleet')
    if n9:
        MH |= set(fleet_info(n9)[3])
    fc = jload(OUT / f'fclose/r{r:02d}/fleetclose.json', {}) or {}
    for nm, f in (fc.get('summary') or {}).items():
        if f.get('covered', 0) >= args.min_cov:
            MH |= set(f.get('unplaced', f.get('misses', [])))
    return sorted(MH)


def _w_scan(job):
    f, M, dmax = job
    try:
        z = np.load(f); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL'])
        c = RI.w_cands(('x', st, M, dmax))
    except Exception:
        return f, [], 0
    tf = np.asarray(st['tf'])
    return f, [x for x in c if np.abs(tf - x['t']).min() > 10 * 86400.0], len(st['asts'])


def scan_approaches(cols, M, args):
    """Per miss X: archive columns (not containing X) passing within args.absorb_d AU of X, nearest first,
    preferring deep hosts: [(dist, n, file)]."""
    files = [f for k, f in cols.items() if len(k) >= 12]
    with mp.get_context('fork').Pool(args.nproc) as pool:
        R = pool.map(_w_scan, [(f, M, args.absorb_d) for f in files], chunksize=8)
    by = collections.defaultdict(list)
    for f, c, n in R:
        d = {}
        for x in c:
            if x['ast'] not in d or x['dist'] < d[x['ast']]:
                d[x['ast']] = x['dist']
        for X, dd in d.items():
            by[X].append((dd, n, f))
    return {X: sorted(v, key=lambda x: (x[0] > 0.03, -x[1] if x[0] <= 0.03 else x[0])) for X, v in by.items()}


def run_jobs(jobs, nproc):
    todo = [j for j in jobs if not (ROOT / j['out'] / 'meta.json').exists()]
    say(f'  {len(jobs)} jobs ({len(jobs) - len(todo)} already done); est {sum(j["_est"] for j in todo):.0f} proc-min '
        f'on {nproc} procs')
    todo.sort(key=lambda j: -j['_est'])
    if not todo:
        return
    tic = time.time(); k = 0
    with mp.get_context('fork').Pool(nproc, maxtasksperchild=1) as pool:
        for res in pool.imap_unordered(IG.run_job, [{kk: v for kk, v in j.items() if kk != '_est'} for j in todo],
                                       chunksize=1):
            k += 1
            if res.get('error'):
                say(f'  job {res.get("tag")} ERROR: {res["error"][-200:]}')
                continue
            jb = res.get('job', {}); cols = res.get('cols', [])
            iso = max([len(c.get('iso', [])) for c in cols] + [0])
            deep = max([c['n'] for c in cols] + [0])
            say(f'  [{k}/{len(todo)}] {jb.get("tag")} ({res.get("kind", "column")}): {len(cols)} columns, deepest {deep}, '
                f'max iso {iso}, {res.get("wall_s")} s (elapsed {(time.time() - tic) / 60:.0f} min)')
            if (OUT / 'STOP').exists():
                say('  STOP file: terminating the pool'); pool.terminate(); break


def master(r, args, suffix=''):
    d = OUT / f'master/r{r:02d}{suffix}'
    if (d / 'master.json').exists():
        return jload(d / 'master.json')
    cmd = [sys.executable, 'tools/s13_master.py', str(d)] + ARCHIVE + [S15COLS, '--N', '8,9,10', '--nproc', '4', '--tlim', str(args.tlim)]
    say(f'  master: {" ".join(cmd[1:4])} ... --N 8,9,10')
    with open(OUT / f'master_r{r:02d}{suffix}.out', 'w') as f:
        subprocess.run(['nice', '-n', '5'] + cmd, cwd=ROOT, stdout=f, stderr=subprocess.STDOUT, check=False)
    return jload(d / 'master.json')


def close_jobs(r, fleet_dir, args):
    """Closer jobs for a picked fleet: every miss is offered, nearest first, to each fleet route passing within
    args.close_d AU of it (sequential greedy absorption per host, run_close)."""
    F, names, cov, M = fleet_info(fleet_dir)
    if not M:
        return []
    with mp.get_context('fork').Pool(min(args.nproc, len(names))) as pool:
        R = pool.map(RI.w_cands, [(n, F.routes[n]['st'], M, args.close_d) for n in names])
    jobs = []
    for n, c in zip(names, R):
        tf = np.asarray(F.routes[n]['st']['tf'])
        d = {}
        for x in c:
            if np.abs(tf - x['t']).min() > 10 * 86400.0 and (x['ast'] not in d or x['dist'] < d[x['ast']]):
                d[x['ast']] = x['dist']
        Xs = [X for X, _ in sorted(d.items(), key=lambda kv: kv[1])][:4]
        if Xs:
            tag = f'Z{n}'
            jobs.append(dict(kind='close', tag=tag, out=f'results/s15/cols/r{r:02d}/close/{tag}',
                             items=[(str(pathlib.Path(fleet_dir).relative_to(ROOT) / f'route_{n}.npz'), Xs)],
                             absorb_dmax=args.close_d, absorb_kg=150.0, item_trials=3, max_abs=3, wall_budget=2400.0,
                             _est=6.0 * len(Xs)))
    return jobs


def summarize(r, jobs, m):
    """Per-round record: master results + generation statistics (gate (i) numbers)."""
    cols = []; chains = []
    for j in jobs:
        mm = jload(ROOT / j['out'] / 'meta.json', {})
        for c in mm.get('cols', []):
            c = dict(c); c['job'] = j['tag']; c['jkind'] = j['kind']; c['S'] = j.get('S', []); c['p'] = j.get('p')
            cols.append(c)
        if j['kind'] == 'chain' and mm:
            chains.append(dict(tag=j['tag'], pool=mm.get('pool_n'), covered=len(mm.get('chain', {}).get('covered', [])),
                               taken=[(t['n'], t['tank'], sorted(set(t['targets']) & set(IG.ISO))) for t in mm.get('chain', {}).get('taken', [])]))
    steer = [c for c in cols if c['jkind'] == 'column' and c['S']]
    g1 = [c for c in steer if c['n'] >= 35 and c['tank'] <= 1100 and set(c['iso']) & set(c['S'])]
    g1_any = [c for c in cols if c['n'] >= 35 and c['tank'] <= 1100 and c['iso']]
    rec = dict(round=r, columns=len(cols), jobs=len(jobs),
               steered_cols=len(steer), gate_i=len(g1), gate_i_any_iso=len(g1_any),
               gate_i_examples=sorted([(c['n'], round(c['tank']), c['iso'], c['job']) for c in g1], key=lambda x: (-len(x[2]), -x[0]))[:12],
               absorbed=sum(1 for c in cols if 'absorbed' in c['kind']), chains=chains)
    if m:
        for N in ('8', '9', '10'):
            x = m['N'].get(N) or m['N'].get(int(N)) if isinstance(m.get('N'), dict) else None
            if x:
                rec[f'N{N}'] = dict(covered=x['covered'], sumJi=x['sumJi'], J=x['J'], lp=x.get('lp'), fleet=x.get('fleet'),
                                    depths=[p['n'] for p in x.get('picks', [])])
        rec['master_columns'] = m.get('columns'); rec['master_files'] = m.get('files')
    return rec


def report(r, R, L):
    say(f'=== round {r} done ({R["wall_s"] / 60:.0f} min): {R["columns"]} columns '
        f'({R["absorbed"]} absorbed); gate(i) steered iso columns at depth>=35 & tank<=1100: {R["gate_i"]}')
    for N in ('8', '9', '10'):
        x = R.get(f'N{N}')
        if x:
            say(f'    N<={N}: {x["covered"]} covered, sum J_i {x["sumJi"]:.4f}, J {x["J"]:.4f}, LP {x["lp"]}; depths {x["depths"]}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--rounds', type=int, default=1); ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--start', default='results/s14/probe_cover/N9/fleet'); ap.add_argument('--target', type=float, default=11.90)
    ap.add_argument('--beam', type=int, default=100); ap.add_argument('--groups', type=int, default=8)
    ap.add_argument('--pure', type=int, default=2); ap.add_argument('--chains', type=int, default=8)
    ap.add_argument('--complements', type=int, default=4)
    ap.add_argument('--absorb-d', type=float, default=0.05); ap.add_argument('--absorb-items', type=int, default=5)
    ap.add_argument('--close', type=int, default=1, help='fleet closer on the round\'s chain fleets before the master')
    ap.add_argument('--close-d', type=float, default=0.25); ap.add_argument('--min-cov', type=int, default=280)
    ap.add_argument('--grow', type=int, default=1); ap.add_argument('--step-grow', type=int, default=3)
    ap.add_argument('--tlim', type=float, default=1200.0)
    a = ap.parse_args()
    if a.nproc > 8:
        raise SystemExit('CPU budget: --nproc <= 8')
    # round-1 premium sets: archive co-occurrence groups of the 8 isolated targets (t10d's own groupings first)
    a.round1_groups = [[137, 158, 216], [157, 168], [119, 138], [64], [64, 137, 157], [64, 119, 137], [157, 158, 216],
                       [137, 138, 157, 158]]
    OUT.mkdir(parents=True, exist_ok=True)
    L = jload(OUT / 'loop.json', dict(rounds={}, best=dict(fleet=a.start, covered=BAR['covered'], sumJi=BAR['sumJi'],
                                                           J=BAR['J'], lp=BAR['lp'], src='s14 probe_cover N9')))
    say(f's15_loop start: {vars(a)}; best {L["best"]}')
    RI.eph()
    for r in range(1, a.rounds + 1):
        key = f'{r:02d}'
        if (OUT / 'STOP').exists():
            say('STOP file present: exit'); break
        R = L['rounds'].setdefault(key, {})
        tic = time.time()
        # ---- premium set: misses of the previous round's N<=9 pick and of its candidate chain fleets
        if r > 1 and 'M_hard' not in R:
            R['M_hard'] = hard_set_next(r - 1, L, a); IG.jdump(L, OUT / 'loop.json')
        # ---- phase G: grow the best fleet's heads by cheap nearby residual targets
        if r > 1 and a.grow and 'grown' not in R:
            gj = R.get('grow_jobs') or grow_jobs(r, L['best'], a)
            R['grow_jobs'] = gj; IG.jdump(L, OUT / 'loop.json')
            say(f'=== round {r} phase G: grow {len(gj)} heads of {L["best"]["fleet"]}')
            run_jobs(gj, a.nproc)
            grown = {}
            for g in gj:
                mm = jload(ROOT / g['out'] / 'meta.json', {}) or {}
                for src, fin in (mm.get('final') or {}).items():
                    if fin != src:
                        grown[g['head']] = fin
                for h in mm.get('hosts', []):
                    say(f'    head {g["head"]}: +{h["got"]} for +{h["dkg"]} kg')
            R['grown'] = grown; IG.jdump(L, OUT / 'loop.json')
        # ---- phase J: column jobs
        if 'jobs' not in R:
            R['jobs'] = build_round(r, L['best'], a, R.get('M_hard'), R.get('grown'))
            R['started'] = time.strftime('%m-%d %H:%M:%S'); R['best_in'] = L['best']; IG.jdump(L, OUT / 'loop.json')
        jobs = R['jobs']
        if not R.get('jobs_done'):
            say(f'=== round {r}: best in {L["best"]["covered"]} covered @ {L["best"]["sumJi"]:.4f} (J {L["best"]["J"]:.4f}); '
                f'premium set {len(R.get("M_hard") or [])}; {len(jobs)} jobs: ' +
                ', '.join(f'{k} {v}' for k, v in collections.Counter(j['kind'] for j in jobs).items()))
            run_jobs(jobs, a.nproc)
            if (OUT / 'STOP').exists():
                say('STOP file present after jobs: exit'); break
            R['jobs_done'] = time.strftime('%m-%d %H:%M:%S'); IG.jdump(L, OUT / 'loop.json')
        # ---- phase Z: fleet closer on the round's candidate chain fleets
        if a.close and not R.get('fclose_done'):
            say(f'  fleet closer: tools/s15_fleetclose.py --round {r} --close-d {a.close_d}')
            with open(OUT / f'fclose_r{r:02d}.out', 'a') as f:
                subprocess.run(['nice', '-n', '5', sys.executable, 'tools/s15_fleetclose.py', f'results/s15/fclose/r{r:02d}',
                                '--round', str(r), '--close-d', str(a.close_d), '--nproc', str(a.nproc),
                                '--min-cov', str(a.min_cov)], cwd=ROOT, stdout=f, stderr=subprocess.STDOUT, check=False)
            fc = jload(OUT / f'fclose/r{r:02d}/fleetclose.json', {}) or {}
            for nm, f in (fc.get('summary') or {}).items():
                say(f'    chain fleet {nm}: {f["covered"]} covered @ {f["sumJi"]:.4f}; closer placed {f.get("placed")}; '
                    f'unplaced {f.get("unplaced")}')
            R['fclose_done'] = time.strftime('%m-%d %H:%M:%S'); IG.jdump(L, OUT / 'loop.json')
        # ---- phase M: exact selection
        if not R.get('done'):
            m = master(r, a)
            rec = summarize(r, jobs, m)
            rec['wall_s'] = round(time.time() - tic) + R.get('wall_prev', 0)
            R.update(rec); R['done'] = time.strftime('%m-%d %H:%M:%S'); IG.jdump(L, OUT / 'loop.json')
            report(r, R, L)
        n9 = R.get('N9')
        if n9 and n9['J'] < L['best']['J'] - 1e-6:
            L['best'] = dict(fleet=n9['fleet'], covered=n9['covered'], sumJi=n9['sumJi'], J=n9['J'], lp=n9['lp'],
                             src=f'round {r}')
        IG.jdump(L, OUT / 'loop.json')
        F, names, cov, M = fleet_info(L['best']['fleet'])
        say(f'=== round {r} complete: best N<=9 now {L["best"]["covered"]} @ {L["best"]["sumJi"]:.4f} (J {L["best"]["J"]:.4f}, '
            f'{L["best"]["src"]}); misses {M}')
        if L['best']['covered'] >= 298 and L['best']['sumJi'] <= a.target:
            say(f'TARGET reached: 298 covered at sum J_i {L["best"]["sumJi"]:.4f} <= {a.target}'); break
    say('s15_loop exit')


if __name__ == '__main__':
    main()
