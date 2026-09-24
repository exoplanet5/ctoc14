"""Stage 15 [T]: TWIN TAILS -- the lintwin twin beam of tools/s14_twinbeam.py (state = a SETTLED impulsive twin;
children = every remaining pool target at its single-leg encounter epochs, priced by the linearised whole-prefix twin;
kept children settled) over a residual pool WITH PRIZES (1 on the pool, 1 + p on the premium set S).

Why: the production beam starves in small sparse pools (round 1, fixed-head 65-pool chain: 27 / 17 / 8 flybys;
stage 14 re-fly of t10d's tail pools: production-leg gate 5 of 24-28, lintwin gate 21-24 of 24-27).  The fixed-head
9-craft problem is exactly that: six deep heads leave 65 sparse targets for three tails.

Score  t_end/T + w_fuel F/1400 - sum(prize over flown targets)   (s14 convention, prizes added; F = 1600 (1-exp(-dv/ve)))
Front  per depth: the prize-score best, the lightest and the most-S settled state (every level); at the end every
       front state of depth >= min_save is saved as a column route_<tag>_<n>[_i].npz (run_ialns ist format).
Chain  steps k = 1..steps: pool_k = pool - (targets of the routes TAKEN at steps < k), S_k = S & pool_k; the step takes
       the front state maximising (n + w_iso |ISO & route|, n - (tank - 600)/kappa) with tank <= cap.
State  OUT/chain.json + OUT/s<k>/ckpt.pkl (every level); a restarted job resumes at the level it stopped.

usage: s15_twintail.py run JOB.json            (JOB: tag, out, pool, S, p, steps, beam, tries, nproc, ...; see DEF)
       s15_twintail.py batch JOBS.json --procs 8   (runs jobs as subprocesses, sum of job nproc <= --procs)"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pickle, pathlib, argparse, warnings, subprocess
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)
import numpy as np
import run_ialns as RI
import s14_twinbeam as TB
from greedy_cover import ALL
from ctoc14.search import mask_of, UNREACHABLE
from ctoc14.search import roots as search_roots
from ctoc14.constants import DAY, VE, T_MISSION, FUEL_MAX
warnings.filterwarnings('ignore')

ISO = [64, 119, 137, 138, 157, 158, 168, 216]
DEF = dict(beam=20, tries=60, roots=300, dv=1.2, drmax=0.15, mc=60, k_epochs=3, lin_cap=3.0, res_max=1.5e8,
           grid='0,800,20', iters=100, timeout=240.0, p=0.02, S=[], w_fuel=None, w_t=None, cap=1150.0, kappa=90.0,
           w_iso=1.0, steps=1, nproc=2, wall_budget=8 * 3600.0, min_save=6, save_cap=1300.0, level_budget=None)


def jdump(obj, path):
    path = pathlib.Path(path); tmp = path.with_suffix(path.suffix + f'.tmp{os.getpid()}')
    with open(tmp, 'w') as f:
        json.dump(obj, f, indent=1, default=lambda o: o.tolist() if isinstance(o, np.ndarray) else
                  (float(o) if isinstance(o, np.floating) else (int(o) if isinstance(o, np.integer) else sorted(o))))
    os.replace(tmp, path)


def prize_vec(pool, S, p):
    pz = np.zeros(300)
    for t in pool:
        pz[t - 1] = 1.0
    for t in S:
        if pz[t - 1] > 0:
            pz[t - 1] = 1.0 + p
    return pz


def pscore(ts_, P, prize):
    F = TB.M0P * (1.0 - np.exp(-ts_['dv'] / VE))
    return P.w_t * ts_['t_end'] / T_MISSION + P.w_fuel * F / FUEL_MAX - float(sum(prize[a - 1] for a in ts_['asts']))


def front_add(front, ts_, P, prize, Sset, real=None):
    """front[n] = dict(score=state, light=state, s=state); n = number of REAL pool targets flown (stepping stones,
    J['stones'], are flybys of targets other routes already cover: prize eps, not counted)."""
    n = len(ts_['asts']) if real is None else len(real & set(ts_['asts'])); f = front.setdefault(n, {})
    sc = pscore(ts_, P, prize); ns = len(Sset & set(ts_['asts']))
    if 'score' not in f or sc < f['score'][0]:
        f['score'] = (sc, ts_)
    if 'light' not in f or ts_['tank'] < f['light'][0]:
        f['light'] = (ts_['tank'], ts_)
    if ns > 0 and ('s' not in f or (ns, -sc) > f['s'][0]):
        f['s'] = ((ns, -sc), ts_)


def twin_beam(J, pool, S, outdir, say):
    """Prize-weighted lintwin twin beam over `pool` (all other targets excluded).  Returns the checkpoint dict
    (front, hist, end, ...).  Resumes from outdir/ckpt.pkl."""
    outdir.mkdir(parents=True, exist_ok=True)
    E = RI.eph(); P = TB.params(J['dv'], J['drmax'], J['mc'])
    if J.get('w_fuel') is not None:
        P.w_fuel = float(J['w_fuel'])
    if J.get('w_t') is not None:
        P.w_t = float(J['w_t'])
    prize = prize_vec(pool, S, J['p']); Sset = set(S)
    stones = sorted(set(J.get('stones') or []) - set(pool)); real = set(pool) if stones else None
    for t in stones:
        prize[t - 1] = float(J.get('stone_prize', 0.01))
    P.prize = prize
    bpool = sorted(set(pool) | set(stones))           # the beam's candidate set (real targets + stepping stones)
    excl = sorted(set(ALL) - set(bpool)); excl_mask = mask_of(set(excl) | UNREACHABLE)
    g = [float(x) for x in J['grid'].split(',')]; grid = np.arange(g[0], g[1] + 1e-9, g[2]) * DAY
    ck = outdir / 'ckpt.pkl'; nproc = int(J['nproc']); tic0 = time.time()
    if ck.exists():
        C = pickle.load(open(ck, 'rb'))
        say(f'  RESUME level {C["level"]}: beam {len(C["beam"])}, deepest {max(C["front"]) if C["front"] else 0}')
    else:
        rts = search_roots(E, grid, TB.M0P, excl, P)
        rts.sort(key=lambda s: s.score(P)); seen = set(); b2 = []
        for s in rts:
            key = (s.seq[0][0], int(s.t_launch // (20 * DAY)))
            if key in seen:
                continue
            seen.add(key); b2.append(s)
        b2 = b2[:J['roots']]
        beam = [TB.root_state(s.t_launch, s.vinf, s.seq[0][0], s.seq[0][1]) for s in b2]
        C = dict(level=1, beam=beam, front={}, hist=[], end=None, settles=0, settle_ok=0, wall=0.0)
        for b in beam:
            front_add(C['front'], b, P, prize, Sset, real)
        say(f'  roots: {len(rts)} launch legs -> {len(beam)} kept; pool {len(pool)}, S {sorted(Sset)} p {J["p"]}')
        TB.save_ckpt(outdir, C)
    t_end_wall = time.time() + J['wall_budget'] - C.get('wall', 0.0)
    while C['end'] is None:
        if time.time() > t_end_wall:
            C['end'] = dict(kind='wall_budget', depth=C['level']); break
        tic = time.time(); beam = C['beam']
        res = TB._pool_map(TB._price_parent, [(ts_, bpool, J['k_epochs'], J['lin_cap'], P.w_t, P.w_fuel, J['res_max'])
                                              for ts_ in beam], nproc)
        kids = []
        for i, (ts_, lst) in enumerate(zip(beam, res)):
            npar = len(ts_['asts']); ppar = float(sum(prize[a - 1] for a in ts_['asts']))
            vm = excl_mask | mask_of(ts_['asts'])
            for sc, X, t, dvm, lin, cm, c1 in lst:
                sc2 = sc + (npar + 1) - (ppar + prize[int(X) - 1])
                kids.append((sc2, i, int(X), float(t), vm | (1 << (int(X) - 1)), lin, dvm))
        t_gen = time.time() - tic
        kids.sort(key=lambda x: x[0])
        seen = set(); cand = []
        for kd in kids:
            key = (kd[4], int(kd[3] // (5 * DAY)))
            if key in seen:
                continue
            seen.add(key); cand.append(kd)
        new = []; tried = 0; nk = len(kids); nu = len(cand)
        while cand and len(new) < J['beam'] and tried < J['tries']:
            take = cand[:max(1, min(nproc, J['tries'] - tried))]; cand = cand[len(take):]
            jobs = [(beam[i], X, t, J['iters'], J['timeout'], lin) for sc, i, X, t, vm, lin, dvm in take]
            out = TB._pool_map(TB.settle_child, jobs, nproc)
            tried += len(take)
            for tw, info in out:
                C['settles'] += 1
                if tw is not None and tw['tank'] <= J['save_cap']:
                    C['settle_ok'] += 1; new.append(tw)
        new.sort(key=lambda t: pscore(t, P, prize))
        seen = set(); nb = []
        for t in new:
            key = (mask_of(t['asts']), int(t['t_end'] // (5 * DAY)))
            if key in seen:
                continue
            seen.add(key); nb.append(t)
        nb = nb[:J['beam']]
        depth = C['level'] + 1
        for t in nb:
            front_add(C['front'], t, P, prize, Sset, real)
        rec = dict(depth=depth, parents=len(beam), children=nk, unique=nu, tried=tried, settled=len(new), kept=len(nb),
                   wall_s=round(time.time() - tic, 1), gen_s=round(t_gen, 1))
        if nb:
            rec.update(t_end_d=[round(min(b['t_end'] for b in nb) / DAY), round(max(b['t_end'] for b in nb) / DAY)],
                       tank=[round(min(b['tank'] for b in nb)), round(max(b['tank'] for b in nb))],
                       maxS=max(len(Sset & set(b['asts'])) for b in nb))
        C['hist'].append(rec)
        say(f'  depth {depth:2d}: {len(beam)} parents -> {nk} kids ({nu} uniq) -> tried {tried}, settled {len(new)}, kept '
            f'{len(nb)}' + (f'; t_end {rec["t_end_d"][0]}-{rec["t_end_d"][1]} d, tank {rec["tank"][0]}-{rec["tank"][1]}, '
                            f'max S {rec["maxS"]}' if nb else '') + f' ({rec["wall_s"]:.0f} s, gen {t_gen:.0f} s)')
        if not nb:
            C['end'] = dict(kind='no_children' if nk == 0 else 'settle_wall', depth=C['level'])
        elif max(len(set(pool) & set(b['asts'])) for b in nb) >= len(pool):
            C['end'] = dict(kind='pool_exhausted', depth=depth)
        C['beam'] = nb if nb else beam; C['level'] = depth if nb else C['level']
        C['wall'] = C.get('wall', 0.0) + (time.time() - tic)
        TB.save_ckpt(outdir, C)
    C['wall_this'] = round(time.time() - tic0)
    TB.save_ckpt(outdir, C)
    return C


def save_front(C, outdir, tag, min_save, Sset, real=None):
    """Every distinct front state of depth >= min_save -> column.  Returns meta list."""
    cols = []; seen = set()
    for n in sorted(C['front'], reverse=True):
        if n < min_save:
            continue
        for kind in ('score', 'light', 's'):
            if kind not in C['front'][n]:
                continue
            ts_ = C['front'][n][kind][1]
            key = (frozenset(ts_['asts']), round(ts_['tank'], 1))
            if key in seen:
                continue
            seen.add(key)
            p = outdir / f'route_{tag}_{n}.npz'; i = 1
            while p.exists():
                p = outdir / f'route_{tag}_{n}_{i}.npz'; i += 1
            np.savez(p, **ts_['st'])
            tg = sorted(int(a) for a in ts_['asts'])
            cols.append(dict(f=str(p.relative_to(ROOT)), n=len(tg), n_real=len(tg) if real is None else len(set(tg) & real),
                             tank=round(ts_['tank'], 2), miss=round(ts_['miss'], 1),
                             kind=kind, targets=tg, iso=sorted(set(tg) & set(ISO)), s_taken=sorted(set(tg) & Sset),
                             t_end_d=round(ts_['t_end'] / DAY, 1), t_launch_d=round(float(ts_['st']['tL']) / DAY, 1)))
    return cols


def run(J):
    J = dict(DEF, **J)
    out = ROOT / J['out']; out.mkdir(parents=True, exist_ok=True)
    if (out / 'meta.json').exists():
        return json.load(open(out / 'meta.json'))
    say = RI.logger(out / 'log.txt'); tic = time.time()
    RI.eph()
    ckj = out / 'chain.json'
    CH = json.load(open(ckj)) if ckj.exists() else dict(steps=[], covered=[], cols=[], taken=[])
    pool0 = sorted(J['pool'])
    say(f'twintail {J["tag"]}: pool {len(pool0)}, S {sorted(set(J["S"]) & set(pool0))}, p {J["p"]}, steps {J["steps"]}, '
        f'beam {J["beam"]}/{J["tries"]}, nproc {J["nproc"]}; resume at step {len(CH["steps"]) + 1}')
    while len(CH['steps']) < J['steps']:
        k = len(CH['steps']) + 1
        pool = sorted(set(pool0) - set(CH['covered']))
        if len(pool) < 3:
            break
        S = sorted(set(J['S']) & set(pool))
        sd = out / f's{k}'
        C = twin_beam(J, pool, S, sd, say)
        cols = save_front(C, out, f'{J["tag"]}s{k}', J['min_save'] if k == 1 else 3, set(S),
                          set(pool) if J.get('stones') else None)
        for c in cols:
            c['step'] = k
        CH['cols'] += cols
        ok = [c for c in cols if c['tank'] <= J['cap']]
        step = dict(k=k, pool=len(pool), S=S, end=C['end'], deepest=max(C['front']) if C['front'] else 0,
                    settles=C['settles'], settle_ok=C['settle_ok'], wall_s=round(C['wall']), ncols=len(cols))
        if not ok:
            step['took'] = None; CH['steps'].append(step); jdump(CH, ckj)
            say(f'  step {k}: nothing under the cap -> chain ends'); break
        took = max(ok, key=lambda c: (c['n_real'] + J['w_iso'] * len(c['iso']), c['n_real'] - (c['tank'] - 600.0) / J['kappa']))
        step['took'] = dict(f=took['f'], n=took['n'], tank=took['tank'], iso=took['iso'], s_taken=took['s_taken'])
        CH['taken'].append(dict(k=k, f=took['f'], n=took['n'], tank=took['tank'], targets=took['targets']))
        CH['covered'] = sorted(set(CH['covered']) | set(took['targets']))
        step['covered'] = len(CH['covered'])
        CH['steps'].append(step); jdump(CH, ckj)
        say(f'  step {k}: pool {len(pool)} deepest {step["deepest"]} ({C["end"]}); TOOK {took["n"]} fb @ {took["tank"]:.0f} kg '
            f'iso {took["iso"]}; chain covered {len(CH["covered"])}/{len(pool0)} ({step["wall_s"]} s)')
    meta = dict(job={kk: v for kk, v in J.items() if kk != 'pool'}, kind='twintail', pool_n=len(pool0), chain=CH,
                cols=CH['cols'], covered=len(CH['covered']), wall_s=round(time.time() - tic))
    jdump(meta, out / 'meta.json')
    say(f'twintail {J["tag"]} done: {len(CH["cols"])} columns, covered {len(CH["covered"])}/{len(pool0)} '
        f'({time.time() - tic:.0f} s)')
    return meta


def busy_python():
    """System-wide number of CPU-bound Python processes (macOS shows the framework 'Python.app' binary)."""
    try:
        out = subprocess.run(['ps', '-Ao', 'pcpu,command'], capture_output=True, text=True).stdout.splitlines()
        return sum(1 for l in out if 'Python.app' in l and float(l.split()[0]) > 50)
    except Exception:
        return 0


def batch(jobs_file, procs, global_cap=None):
    """Run twintail jobs as subprocesses; the sum of running jobs' nproc stays <= procs, and (global_cap) a job starts
    only while the system-wide count of CPU-bound Python processes + its nproc stays <= global_cap."""
    jobs = json.load(open(jobs_file))
    todo = [j for j in jobs if not (ROOT / j['out'] / 'meta.json').exists()]
    run_ = []; jdir = ROOT / 'results/s15/twin_jobs'; jdir.mkdir(parents=True, exist_ok=True)
    while todo or run_:
        run_ = [(j, p) for j, p in run_ if p.poll() is None]
        used = sum(int(j.get('nproc', DEF['nproc'])) for j, _ in run_)
        while (todo and used + int(todo[0].get('nproc', DEF['nproc'])) <= procs and
               (global_cap is None or busy_python() + int(todo[0].get('nproc', DEF['nproc'])) <= global_cap)):
            j = todo.pop(0); jf = jdir / f'{j["tag"]}.json'; json.dump(j, open(jf, 'w'))
            (ROOT / j['out']).mkdir(parents=True, exist_ok=True)
            p = subprocess.Popen(['nice', '-n', '5', sys.executable, str(ROOT / 'tools/s15_twintail.py'), 'run', str(jf)],
                                 cwd=ROOT, stdout=open(ROOT / j['out'] / 'stdout.txt', 'a'), stderr=subprocess.STDOUT)
            run_.append((j, p)); used += int(j.get('nproc', DEF['nproc']))
            print(f'[{time.strftime("%H:%M:%S")}] started {j["tag"]} (nproc {j.get("nproc", DEF["nproc"])}; {used}/{procs})',
                  flush=True)
            if global_cap is not None:
                time.sleep(40)                      # let the new process show up as busy before the next check
        time.sleep(10)
    print(f'[{time.strftime("%H:%M:%S")}] batch done', flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('mode', choices=['run', 'batch']); ap.add_argument('file')
    ap.add_argument('--procs', type=int, default=8)
    ap.add_argument('--global-cap', type=int, default=None, help='start a job only while system-wide busy Python + nproc <= this')
    a = ap.parse_args()
    if a.mode == 'run':
        m = run(json.load(open(a.file)))
        print(json.dumps(dict(covered=m.get('covered'), cols=len(m.get('cols', [])), wall_s=m.get('wall_s'))))
    else:
        if a.procs > 8:
            raise SystemExit('CPU budget: --procs <= 8')
        if a.global_cap is not None and a.global_cap > 8:
            raise SystemExit('CPU budget: --global-cap <= 8')
        batch(a.file, a.procs, a.global_cap)
