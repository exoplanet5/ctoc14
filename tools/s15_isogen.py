"""Stage 15 [G]: ONE column = one deep FREE-beam route steered to absorb a chosen subset S of (isolated) targets.

Pipeline of one column job (runs in ONE process; the loop runs 8 of them at a time):
  1. PLAN   ctoc14.search.beam_search (m0 1600, production leg model) over `pool` (default all 298, strict exclusion of
            everything else); prize 1 on the pool, 1 + p on S (p in (0, 0.025]), a tiny seeded all-target jitter
            (+-jit, default 0.001) and a seeded launch-grid offset for diversity; legs INTO S get the leg cap
            max(member dv, dvS) (the planner cannot even generate 2 of the 8 isolated legs at 2.0 km/s: stage-15
            diagnostic, results/s15/diag/mispricing.json).
            Front = per depth the min-fuel state AND the min-fuel state among those taking the most of S (the plain
            min-fuel front is biased against S: prizes are not fuel), deepest `front` depths.
  2. SETTLE deepest first (ties: more of S, lower planner tank); a candidate is kept when miss <= 150 km, twin tank <=
            1.15 x planner and <= save_cap; EVERY kept candidate is saved as a column (the selector gets the whole
            Pareto front); strided walk: the first `dense` successes consecutive, then every `stride`-th depth.
  3. REPAIR if the deepest candidate FAILED to settle and nothing within 2 flybys of it settled: s13_bisect.failing_leg
            (the longest settling prefix is saved as a column), ban that leg's target, re-plan once, settle.
  4. ABSORB (twin pricing of the S legs -- the diagnostic's action): for the best settled hosts, every target of
            `absorb` (default S + the loop's current misses) not on the host is screened at the host's approach minima
            (run_ialns.w_cands, <= dmax AU) with s14_twinbeam.lin_price (linearised whole-route twin, all impulses /
            epochs / v_inf free; lintwin/removal price ratio median 1.01 on t10d), and the cheapest are inserted in the
            twin (run_ialns.w_insert); settled insertions are saved as extra columns (a second absorption is tried on
            the first absorbed column).
Outputs (OUT = the job dir): route_<tag>_<n>[_<i>].npz columns (run_ialns ist format), meta.json (job, front, every
settle, which targets of S / of the ISO set each column takes, timings).  A job whose meta.json exists is done.

CHAIN jobs (tail columns on residual pools): steps k = 1..L, step k plans over pool_k = pool - (targets of the chosen
routes of steps < k), with S_k = the chain's premium set restricted to pool_k; each step runs 1-4 (tail portfolio
members) as above and TAKES the deepest settled candidate (value tie-break n - (tank-600)/kappa); state chain.json.

usage (test): s15_isogen.py JOB.json            (JOB = dict(kind='column'|'chain', tag, out, ...); see run_job)"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, signal, pathlib, traceback, warnings
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools'), str(ROOT / 'results/s13/planner_probe')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)
import numpy as np
import run_ialns as RI
import s13_bisect as SB
import s14_twinbeam as TB
import ctoc14.impulsive as IM
from greedy_cover import BP, CM, DeepCollect, ALL, _timeout
from ctoc14.search import beam_search, tour_json
from ctoc14.constants import DAY, AU, VE
warnings.filterwarnings('ignore')
IM.LIN_FALLBACK = 2.5

ISO = [64, 119, 137, 138, 157, 158, 168, 216]
M0 = 1600.0
GUARD = 1.15
MISS_OK = 150.0
# leg-model members (tools/s13_pass.py MEMBERS; plan section 4.1 / M-D)
MEMBERS = {'H':  dict(dv=1.2, dr=0.15, npt=6, mc=150, g=800, wt=1.0),
           'T1': dict(dv=1.2, dr=0.15, npt=6, mc=150, g=800, wt=1.0),
           'T2': dict(dv=1.6, dr=0.20, npt=6, mc=150, g=800, wt=1.0),
           'T3': dict(dv=2.0, dr=0.25, npt=2, mc=60, g=800, wt=1.0),
           'T5': dict(dv=1.6, dr=0.20, npt=6, mc=150, g=800, wt=0.5),
           'T6': dict(dv=2.0, dr=0.25, npt=6, mc=150, g=1500, wt=1.0)}

DEFAULTS = dict(pool=None, S=[], p=0.01, seed=0, beam=100, member='H', dvS=2.0, jit=0.001, front=24,
                cap=1100.0, save_cap=1300.0, max_settles=10, ok_target=6, dense=3, stride=3,
                repair=True, repair_settles=5, absorb=None, absorb_hosts=2, absorb_trials=4, absorb_dmax=0.20,
                absorb_kg=80.0, settle_timeout=240.0, kappa=90.0, wall_budget=2400.0)


def now():
    return time.strftime('%H:%M:%S')


def _jd(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (set, frozenset)):
        return sorted(o)
    raise TypeError(type(o))


def jdump(obj, path):
    path = pathlib.Path(path); tmp = path.with_suffix(path.suffix + f'.tmp{os.getpid()}')
    with open(tmp, 'w') as f:
        json.dump(obj, f, indent=1, default=_jd)
    os.replace(tmp, path)


class Log:
    def __init__(self, path):
        self.f = open(path, 'a')

    def __call__(self, s):
        self.f.write(f'[{now()}] {s}\n'); self.f.flush()


# -------------------------------------------------------------------------------------------------------- 1. plan
def params(m, beam, prize, S, dvS):
    P = BP(beam=beam, w_fuel=CM.w_fuel(480, 1600), m_margin=40.0, dv_max=m['dv'], tofs=np.arange(15, 401, 5) * DAY,
           lin_tofs=np.arange(20, 601, 10) * DAY, lin_drmax=m['dr'] * AU, vinf_cap=4.0, tof_refine=True, max_depth=60,
           n_per_target=m['npt'], max_children=m['mc'], w_t=m['wt'])
    P.prize = np.asarray(prize, float)
    if S and dvS and dvS > m['dv']:
        d = np.full(300, m['dv'])
        for t in S:
            d[t - 1] = max(m['dv'], dvS)
        P.dv_max_t = d
    return P


def prize_vector(pool, S, p, seed, jit):
    rng = np.random.default_rng(10007 + int(seed))
    pz = np.zeros(300)
    for t in pool:
        pz[t - 1] = 1.0 + (rng.uniform(-jit, jit) if jit > 0 else 0.0)
    for t in S:
        if t in set(pool):
            pz[t - 1] = 1.0 + p
    return pz


def plan(pool, S, p, member, beam, seed, jit, dvS, front, excl_extra=()):
    """Member beam over `pool` (targets in excl_extra banned).  Returns dict(front=[cand], sec, nstates, deepest_S)."""
    m = MEMBERS[member] if isinstance(member, str) else member
    tic = time.time()
    pool = sorted(set(pool) - set(excl_extra))
    S = [t for t in S if t in set(pool)]
    P = params(m, beam, prize_vector(pool, S, p, seed, jit), S, dvS)
    P.collect = DeepCollect(4)
    off = (int(seed) % 4) * 5.0
    grid = (np.arange(0, m['g'] + 1e-9, 20) + off) * DAY
    best, fin = beam_search(RI.eph(), excluded=sorted(set(ALL) - set(pool)), m0=M0, t_launch_grid=grid, P=P,
                            n_proc=1, verbose=False)
    states = list(P.collect) + list(fin)
    Sset = set(S)
    byn = {}; bys = {}
    for s in states:
        k = len(s.seq); ns = sum(1 for x in s.seq if x[0] in Sset)
        if k not in byn or s.fuel < byn[k].fuel:
            byn[k] = s
        if ns > 0:
            o = bys.get(k)
            if o is None or (ns, -s.fuel) > (sum(1 for x in o.seq if x[0] in Sset), -o.fuel):
                bys[k] = s
    depths = sorted(byn, reverse=True)[:front]
    cands = []
    for k in depths:
        for s in ([byn[k]] + ([bys[k]] if k in bys and bys[k] is not byn[k] else [])):
            tr = tour_json(s)
            tg = sorted(int(l['ast']) for l in tr['legs'])
            cands.append(dict(n=k, pt=float(CM.tank(s.fuel, M0)), tour=tr, targets=tg,
                              s_taken=sorted(Sset & set(tg)), member=member if isinstance(member, str) else 'X'))
    dS = max((k for k in bys), default=0)
    return dict(front=cands, sec=time.time() - tic, nstates=len(states), deepest_S=dS,
                deepest=max(byn) if byn else 0, S=S, pool_n=len(pool))


# ------------------------------------------------------------------------------------------------------ 2. settle
def settle(tour, timeout):
    return SB.settle_one(tour, 'stock', timeout=timeout)


def ok_res(r, pt, save_cap):
    return (r.get('st') is not None and r.get('tank') is not None and r['miss'] <= MISS_OK and
            r['tank'] <= GUARD * pt + 1e-9 and r['tank'] <= save_cap + 1e-9)


def walk(cands, J, log, t_end):
    """Strided deepest-first settle walk.  Returns the list of settle records (cand + res + ok)."""
    order = sorted(cands, key=lambda c: (-c['n'], -len(c['s_taken']), c['pt']))
    order = [c for c in order if c['pt'] <= J['save_cap']]
    recs = []; nok = 0; last_ok_n = None; skip_to = None; cap_ok = False
    need_cap = J.get('need_cap_ok', False)
    for c in order:
        if len(recs) >= J['max_settles'] or time.time() > t_end:
            break
        if nok >= J['ok_target'] and (cap_ok or not need_cap):
            break
        if skip_to is not None and c['n'] > skip_to:
            continue
        r = settle(c['tour'], J['settle_timeout'])
        ok = ok_res(r, c['pt'], J['save_cap'])
        rec = dict(c, res=r, ok=ok, failed=not (r.get('st') is not None and r['miss'] <= MISS_OK and
                                                    r['tank'] is not None and r['tank'] <= GUARD * c['pt'] + 1e-9))
        recs.append(rec)
        log(f'    settle {c["member"]} {c["n"]} fb (S {c["s_taken"]}) planner {c["pt"]:.0f} -> twin '
            f'{"-" if r["tank"] is None else round(r["tank"])} miss {r["miss"]:.0f} {"OK" if ok else "no"} ({r["sec"]:.0f} s)')
        if ok:
            nok += 1; last_ok_n = c['n']; cap_ok = cap_ok or r['tank'] <= J['cap'] + 1e-9
            if nok >= J['dense']:
                skip_to = c['n'] - J['stride']
        elif skip_to is not None:
            skip_to = min(skip_to, c['n'])          # a failure: go on with the next candidate (same depth allowed)
    return recs


# ------------------------------------------------------------------------------------------------------ 4. absorb
def nreg_of(st):
    return len(TB.centres(float(st['tL']), float(np.max(st['tf']))))


def absorb_targets(st, tank, X_list, J, log, t_end, tag=''):
    """Twin insertion of targets X_list into the settled route st, screened by the linearised twin (lin_price).
    Returns list of dict(st, tank, miss, ast, t, lin_kg, dist)."""
    have = set(int(a) for a in st['asts'])
    X_list = [x for x in X_list if x not in have]
    if not X_list:
        return []
    try:
        apps = RI.w_cands(('h', st, X_list, J['absorb_dmax']))
    except Exception as e:
        log(f'    absorb{tag}: w_cands failed {type(e).__name__}'); return []
    tf = np.asarray(st['tf'], float)
    apps = [a for a in apps if np.abs(tf - a['t']).min() > 10 * DAY]
    if not apps:
        log(f'    absorb{tag}: no approach of {X_list} within {J["absorb_dmax"]} AU'); return []
    by = {}
    for a in sorted(apps, key=lambda a: a['dist']):
        by.setdefault(a['ast'], [])
        if len(by[a['ast']]) < 3:
            by[a['ast']].append(a)
    ip = RI.ipr(st); par = dict(st=st, nreg=nreg_of(st), dv=float(ip.dv()))
    priced = []
    for X, lst in by.items():
        for a in lst:
            best = None
            for dt in (0.0, -4.0, 4.0):
                try:
                    dvm, lin, cm = TB.lin_price(par, X, a['t'] + dt * DAY)
                except Exception:
                    continue
                if not np.isfinite(dvm):
                    continue
                if best is None or dvm < best[0]:
                    best = (float(dvm), a['t'] + dt * DAY)
            if best is not None:
                kgp = float(tank * (np.exp(max(best[0], 0.0) / VE) - 1.0))
                priced.append(dict(ast=X, t=best[1], lin_kg=kgp, dist=a['dist']))
    priced.sort(key=lambda d: d['lin_kg'])
    log(f'    absorb{tag}: {len(apps)} approaches of {sorted(by)}; lintwin ' +
        ', '.join(f'{d["ast"]}@{d["t"] / DAY:.0f}d {d["dist"]:.3f}AU {d["lin_kg"]:.0f}kg' for d in priced[:6]))
    out = []; tried = 0; done_X = set()
    for d in priced:
        if tried >= J['absorb_trials'] or time.time() > t_end:
            break
        if d['lin_kg'] > J['absorb_kg'] or d['ast'] in done_X:
            continue
        tried += 1; tic = time.time()
        try:
            signal.signal(signal.SIGALRM, _timeout); signal.setitimer(signal.ITIMER_REAL, J['settle_timeout'])
            r = (w_insert_long(('h', st, d['ast'], d['t'])) if d['dist'] > J.get('long_d', 0.08)
                 else RI.w_insert(('h', st, d['ast'], d['t'])))
        except Exception as e:
            r = dict(ok=False, error=type(e).__name__)
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
        ok = bool(r.get('ok')) and r.get('tank', 1e9) <= J['save_cap']
        log(f'    absorb{tag}: insert {d["ast"]} at {d["t"] / DAY:.0f} d (lintwin {d["lin_kg"]:.0f} kg) -> '
            + (f'OK tank {r["tank"]:.0f} (+{r["tank"] - tank:.1f} kg)' if ok else f'fail {r.get("miss", r.get("error", ""))}')
            + f' ({time.time() - tic:.0f} s)')
        if ok:
            done_X.add(d['ast'])
            out.append(dict(st=r['st'], tank=float(r['tank']), miss=float(r['miss']), ast=d['ast'], t=d['t'],
                            lin_kg=d['lin_kg'], dist=d['dist'], dkg=float(r['tank'] - tank)))
    return out


def w_insert_long(job, stages=20, iters=8):
    """run_ialns.w_insert with a 2x longer aim-point homotopy (far targets: 0.1-0.3 AU = 15-45e6 km; plain 10-stage
    homotopies stall, stage-13 M-C closer)."""
    from ctoc14.globalopt import insert_homotopy
    name, st, ast, t = job
    ip = RI.ipr(st)
    try:
        d = insert_homotopy(ip, ast, t, stages=stages, iters_stage=iters, tol_km=100, log=RI.QUIET)
        if d.max() > 1e4:
            return dict(name=name, ast=ast, t=t, ok=False, miss=float(d.max()))
        miss = RI.settle(ip, 100)
    except Exception as e:
        return dict(name=name, ast=ast, t=t, ok=False, miss=np.inf, error=repr(e)[:80])
    return dict(name=name, ast=ast, t=t, ok=miss <= 150, st=RI.ist(ip), tank=float(ip.tank()), miss=miss)


def grow(st, tank, X_list, J, log, t_end, dmax, kg_cap, n_max, tag=''):
    """Sequential greedy twin absorption: absorb the cheapest target of X_list (lin_price screened, within dmax AU)
    if it costs <= kg_cap, re-screen the rest on the grown route, at most n_max.  Returns (st, tank, got, steps)."""
    got = []; steps = []; left = [int(x) for x in X_list]
    while left and len(got) < n_max and time.time() < t_end - 60:
        res = absorb_targets(st, tank, left, dict(J, absorb_dmax=dmax, absorb_kg=kg_cap, absorb_trials=3), log, t_end,
                             tag=f'{tag}+{got}')
        res = [a for a in res if a['dkg'] <= kg_cap]
        if not res:
            break
        b = min(res, key=lambda a: a['dkg'])
        for a in res:                       # the other settled insertions of this iteration: columns too
            if a is not b:
                steps.append(dict(st=a['st'], tank=a['tank'], ast=a['ast'], dkg=a['dkg'], extra=True, base=list(got)))
        st, tank = b['st'], b['tank']; got.append(b['ast']); left.remove(b['ast'])
        steps.append(dict(st=st, tank=tank, ast=b['ast'], dkg=b['dkg'], extra=False, base=list(got[:-1])))
    return st, tank, got, steps


# ------------------------------------------------------------------------------------------------------ column
def save_col(outdir, tag, st, n, meta_cols, info):
    p = outdir / f'route_{tag}_{n}.npz'; i = 1
    while p.exists():
        p = outdir / f'route_{tag}_{n}_{i}.npz'; i += 1
    np.savez(p, **st)
    tg = sorted(int(a) for a in st['asts'])
    info = dict(info, f=str(p.relative_to(ROOT)), n=len(tg), targets=tg, iso=sorted(set(tg) & set(ISO)))
    meta_cols.append(info)
    return p


def column_step(J, pool, S, outdir, tag, log, t_end, member=None, excl_extra=()):
    """Plan (+ repair) + settle walk for ONE member.  Returns (recs, plan_info)."""
    member = member or J['member']
    pl = plan(pool, S, J['p'], member, J['beam'], J['seed'], J['jit'], J['dvS'], J['front'], excl_extra)
    log(f'  plan {member}: pool {pl["pool_n"]}, S {pl["S"]}, deepest {pl["deepest"]}, deepest with S {pl["deepest_S"]}; front '
        + ' '.join(f'{c["n"]}:{c["pt"]:.0f}{"*" * len(c["s_taken"])}' for c in pl['front'][:10]) + f' ({pl["sec"]:.0f} s)')
    recs = walk(pl['front'], J, log, t_end)
    info = dict(member=member, sec=round(pl['sec']), deepest=pl['deepest'], deepest_S=pl['deepest_S'],
                front=[(c['n'], round(c['pt']), c['s_taken']) for c in pl['front']])
    # repair: the deepest candidate FAILED and nothing within 2 of it settled
    if J['repair'] and recs and recs[0]['failed'] and time.time() < t_end - 400:
        top = recs[0]['n']
        if not any(r['ok'] and r['n'] >= top - 2 for r in recs):
            lo, hi, leg, bi = SB.failing_leg(recs[0]['tour'], max_settles=J['repair_settles'], timeout=J['settle_timeout'])
            info['repair'] = dict(k_ok=lo, k_fail=hi, leg=None if leg is None else int(leg['ast']), settles=bi['settles'])
            log(f'  repair: deepest {top} failed; last settling prefix {lo}, failing leg {hi} '
                f'({None if leg is None else leg["ast"]})')
            if bi['st_ok'] is not None and lo >= 8:
                recs.append(dict(n=lo, pt=float('nan'), tour=None, targets=sorted(int(a) for a in bi['st_ok']['asts']),
                                 s_taken=sorted(set(S) & set(int(a) for a in bi['st_ok']['asts'])), member=member + 'pre',
                                 res=dict(st=bi['st_ok'], tank=float(RI.ipr(bi['st_ok']).tank()), miss=0.0, sec=0.0),
                                 ok=True, failed=False))
            if leg is not None and time.time() < t_end - 300:
                ban = int(leg['ast'])
                pl2 = plan(pool, S, J['p'], member, J['beam'], J['seed'], J['jit'], J['dvS'], J['front'],
                           tuple(excl_extra) + (ban,))
                log(f'  repair plan (ban {ban}): deepest {pl2["deepest"]}; front ' +
                    ' '.join(f'{c["n"]}:{c["pt"]:.0f}{"*" * len(c["s_taken"])}' for c in pl2['front'][:8]) + f' ({pl2["sec"]:.0f} s)')
                nok = max([r['n'] for r in recs if r['ok']] + [0])
                J2 = dict(J, max_settles=4, ok_target=2)
                r2 = walk([c for c in pl2['front'] if c['n'] > nok], J2, log, t_end)
                for r in r2:
                    r['member'] = member + 'r'
                recs += r2
                info['repair']['ban'] = ban; info['repair']['replan_deepest'] = pl2['deepest']
    return recs, info


def run_column(J):
    """One column job.  J = DEFAULTS updated with the job dict (tag, out, pool, S, p, seed, ...)."""
    J = dict(DEFAULTS, **J)
    outdir = ROOT / J['out']; outdir.mkdir(parents=True, exist_ok=True)
    if (outdir / 'meta.json').exists():
        return json.load(open(outdir / 'meta.json'))
    log = Log(outdir / 'log.txt'); tic = time.time(); t_end = tic + J['wall_budget']
    pool = sorted(J['pool']) if J['pool'] else list(ALL)
    S = sorted(set(J['S']) & set(pool))
    ab = sorted(set(J['absorb'] if J['absorb'] is not None else S) & set(pool))
    log(f'column {J["tag"]}: pool {len(pool)}, S {S}, p {J["p"]}, member {J["member"]}, seed {J["seed"]}, absorb {ab}')
    meta_cols = []
    recs, info = column_step(J, pool, S, outdir, J['tag'], log, t_end)
    kept = []
    for r in recs:
        if r['ok']:
            st = r['res']['st']
            save_col(outdir, J['tag'], st, r['n'], meta_cols,
                     dict(kind='settled', member=r['member'], planner=None if not np.isfinite(r['pt']) else round(r['pt'], 1),
                          tank=round(r['res']['tank'], 2), miss=round(r['res']['miss'], 1), s_taken=r['s_taken']))
            kept.append(dict(st=st, tank=r['res']['tank'], n=r['n']))
    # absorb: best hosts = deepest settled with tank <= cap, more isolated targets first, distinct target sets
    isoc = lambda k: len(set(int(a) for a in k['st']['asts']) & set(ISO))
    hosts = []; hsets = []
    for k in sorted([k for k in kept if k['tank'] <= J['cap']], key=lambda k: (-k['n'], -isoc(k), k['tank'])):
        ts_ = frozenset(int(a) for a in k['st']['asts'])
        if ts_ in hsets:
            continue
        hosts.append(k); hsets.append(ts_)
        if len(hosts) >= J['absorb_hosts']:
            break
    if not hosts:
        hosts = sorted(kept, key=lambda k: (-k['n'], k['tank']))[:1]
    ab_log = []
    for hi_, h in enumerate(hosts):
        if not ab or time.time() > t_end - 120:
            break
        res = absorb_targets(h['st'], h['tank'], ab, J, log, t_end, tag=f'[h{h["n"]}]')
        for a in res:
            save_col(outdir, J['tag'] + f'a{a["ast"]}', a['st'], len(a['st']['asts']), meta_cols,
                     dict(kind='absorbed', host_n=h['n'], host_tank=round(h['tank'], 2), ast=a['ast'],
                          t_d=round(a['t'] / DAY, 1), lin_kg=round(a['lin_kg'], 1), dkg=round(a['dkg'], 1),
                          tank=round(a['tank'], 2), miss=round(a['miss'], 1)))
            ab_log.append(dict(host_n=h['n'], ast=a['ast'], dkg=round(a['dkg'], 1), lin_kg=round(a['lin_kg'], 1)))
        if res and time.time() < t_end - 120:             # a second absorption on the cheapest absorbed column
            b = min(res, key=lambda a: a['dkg'])
            rest = [x for x in ab if x not in set(int(q) for q in b['st']['asts'])]
            res2 = absorb_targets(b['st'], b['tank'], rest, dict(J, absorb_trials=2), log, t_end, tag=f'[h{h["n"]}+{b["ast"]}]')
            for a in res2:
                save_col(outdir, J['tag'] + f'a{b["ast"]}a{a["ast"]}', a['st'], len(a['st']['asts']), meta_cols,
                         dict(kind='absorbed2', host_n=h['n'], ast=[b['ast'], a['ast']], dkg=round(a['tank'] - h['tank'], 1),
                              tank=round(a['tank'], 2), miss=round(a['miss'], 1)))
                ab_log.append(dict(host_n=h['n'], ast=[b['ast'], a['ast']], dkg=round(a['tank'] - h['tank'], 1)))
    meta = dict(job={k: v for k, v in J.items() if k != 'pool'}, pool_n=len(pool), S=S, plan=info,
                settles=[dict(member=r['member'], n=r['n'], planner=None if not np.isfinite(r['pt']) else round(r['pt'], 1),
                              tank=None if r['res']['tank'] is None else round(r['res']['tank'], 1),
                              miss=round(r['res']['miss'], 1) if np.isfinite(r['res']['miss']) else None,
                              ok=bool(r['ok']), s_taken=r['s_taken'], sec=round(r['res'].get('sec', 0), 1)) for r in recs],
                absorbed=ab_log, cols=meta_cols, wall_s=round(time.time() - tic), done=now())
    jdump(meta, outdir / 'meta.json')
    log(f'column {J["tag"]} done: {len(meta_cols)} columns ({sum(1 for c in meta_cols if c["kind"] != "settled")} absorbed); '
        f'{time.time() - tic:.0f} s')
    return meta


# ------------------------------------------------------------------------------------------------------ chain
def value(n, tank, kappa):
    return n - (tank - 600.0) / kappa


def joint_pair(pool, S, p, seed, beam, dvS, log, t_end, timeout=240.0, save_cap=1300.0):
    """Joint K=2 beam (results/s13/planner_probe/jointsym, order-free key, wait 60 d, member T7 of s13_pass) over a
    small pool with premium S.  The deepest 3 coverage levels are settled (wait-aware seed, stock fallback); returns
    the settled pairs [(cov, [dict(st, tank, n)])] (both craft settled)."""
    import jointsym
    m = dict(dv=1.6, dr=0.20, npt=2, mc=60, g=800, wt=1.0)
    S = [t for t in S if t in set(pool)]
    P = params(m, beam, prize_vector(pool, S, p, seed, 0.001), S, dvS)
    P.max_depth = 120; P.w_rare = 1.0; P.wait = 60.0 * DAY; P.collect_joint = []
    tic = time.time()
    best = jointsym.joint_beam_search(RI.eph(), [M0, M0], np.arange(0, 800 + 1e-9, 20) * DAY, P=P,
                                      excluded=sorted(set(ALL) - set(pool)), n_proc=1, verbose=False, canon=True)
    front = {}
    for j in list(P.collect_joint) + [best]:
        rs = [s for s in j.craft if s is not None and s.n() > 0]
        if len(rs) < 2:
            continue
        cov = len(set().union(*[{x[0] for x in s.seq} for s in rs]))
        tanks = [float(CM.tank(s.fuel, s.m0)) for s in rs]
        sj = sum(RI.cost(t) for t in tanks)
        if cov not in front or sj < front[cov][0]:
            front[cov] = (sj, [tour_json(s) for s in rs], tanks)
    levels = sorted(front, reverse=True)[:3]
    log(f'  joint K=2: pool {len(pool)}, levels ' + ' '.join(f'{c}={"+".join(str(len(t["legs"])) for t in front[c][1])}'
                                                              f'@{"/".join(f"{x:.0f}" for x in front[c][2])}' for c in levels)
        + f' ({time.time() - tic:.0f} s)')
    out = []
    for c in levels:
        if time.time() > t_end:
            break
        craft = []
        for tr, pt in zip(front[c][1], front[c][2]):
            r = SB.settle_one(tr, 'wait', timeout=timeout)
            if not ok_res(r, pt, save_cap):
                r2 = SB.settle_one(tr, 'stock', timeout=timeout)
                if ok_res(r2, pt, save_cap):
                    r = r2
            ok = ok_res(r, pt, save_cap)
            log(f'    joint level {c}: {len(tr["legs"])} fb planner {pt:.0f} -> twin '
                f'{"-" if r["tank"] is None else round(r["tank"])} miss {r["miss"]:.0f} {"OK" if ok else "no"}')
            if ok:
                craft.append(dict(st=r['st'], tank=r['tank'], n=len(tr['legs'])))
            else:
                break
        if len(craft) == 2:
            out.append((c, craft))
    return out


def run_chain(J):
    """Tail columns over a residual pool, chained.  J: tag, out, pool (residual), S (premium set), p, steps,
    members (list), plus column settings.  Each step k: all members plan over pool_k, one merged strided walk, the
    deepest settled (value tie-break) is TAKEN, every settled candidate is saved; absorb on the taken route."""
    J = dict(DEFAULTS, **J)
    outdir = ROOT / J['out']; outdir.mkdir(parents=True, exist_ok=True)
    if (outdir / 'meta.json').exists():
        return json.load(open(outdir / 'meta.json'))
    log = Log(outdir / 'log.txt'); tic = time.time(); t_end = tic + J['wall_budget']
    ck = outdir / 'chain.json'
    C = json.load(open(ck)) if ck.exists() else dict(steps=[], covered=[], cols=[])
    pool0 = sorted(J['pool'])
    log(f'chain {J["tag"]}: pool {len(pool0)}, S {sorted(set(J["S"]) & set(pool0))}, p {J["p"]}, steps {J["steps"]}, '
        f'members {J["members"]}; resume at step {len(C["steps"]) + 1}')
    while len(C['steps']) < J['steps'] and time.time() < t_end:
        k = len(C['steps']) + 1
        pool = sorted(set(pool0) - set(C['covered']))
        if len(pool) < 3:
            break
        S = sorted(set(J['S']) & set(pool))
        if (J.get('joint_last') and J['steps'] - k == 1 and len(pool) <= J.get('joint_max_pool', 55)
                and not C.get(f'joint{k}') and time.time() < t_end - 900):
            try:
                prs = joint_pair(pool, S, J['p'], J['seed'] + 31 * k, J['beam'], J['dvS'], log, t_end,
                                 J['settle_timeout'], J['save_cap'])
            except Exception as e:
                prs = []; log(f'  joint K=2 failed: {type(e).__name__}: {str(e)[:80]}')
            for c, craft in prs:
                for ci, cr in enumerate(craft):
                    save_col(outdir, f'{J["tag"]}j{k}c{c}{"ab"[ci]}', cr['st'], cr['n'], C['cols'],
                             dict(kind='joint', step=k, level=c, tank=round(cr['tank'], 2)))
            C[f'joint{k}'] = [(c, [cr['n'] for cr in craft], [round(cr['tank'], 1) for cr in craft]) for c, craft in prs] or 'none'
            jdump(C, ck)
        recs = []; infos = []
        mems = J['members'] if len(pool) >= 150 else J.get('members_small', J['members'])
        for m in mems:
            if time.time() > t_end:
                break
            pl = plan(pool, S, J['p'], m, J['beam'], J['seed'] + 17 * k, J['jit'], J['dvS'], J['front'])
            log(f'  step {k} plan {m}: pool {pl["pool_n"]}, S {pl["S"]}, deepest {pl["deepest"]} (with S {pl["deepest_S"]}); '
                'front ' + ' '.join(f'{c["n"]}:{c["pt"]:.0f}{"*" * len(c["s_taken"])}' for c in pl['front'][:8]) +
                f' ({pl["sec"]:.0f} s)')
            infos.append(dict(member=m, sec=round(pl['sec']), deepest=pl['deepest'], deepest_S=pl['deepest_S']))
            recs.append(pl['front'])
        allc = [c for f in recs for c in f]
        seen = set(); cands = []
        for c in sorted(allc, key=lambda c: (-c['n'], -len(c['s_taken']), c['pt'])):
            key = tuple(c['targets'])
            if key in seen:
                continue
            seen.add(key); cands.append(c)
        Jw = dict(J, max_settles=J.get('chain_settles', 12), ok_target=J.get('chain_ok', 5), need_cap_ok=True,
                  cap=J.get('chain_cap', J['cap']))
        rs = walk(cands, Jw, log, t_end)
        oks = [r for r in rs if r['ok']]
        step = dict(k=k, pool=len(pool), S=S, plans=infos, settles=len(rs), ok=len(oks))
        tagk = f'{J["tag"]}s{k}'
        for r in oks:
            save_col(outdir, tagk, r['res']['st'], r['n'], C['cols'],
                     dict(kind='chain', step=k, member=r['member'], planner=round(r['pt'], 1), tank=round(r['res']['tank'], 2),
                          miss=round(r['res']['miss'], 1), s_taken=r['s_taken']))
        ok_cap = [r for r in oks if r['res']['tank'] <= Jw['cap']]
        # REPAIR (stage-15 round 2: C02h5 and K03h4 died here -- every front candidate shared one unsettleable leg):
        # bisect the deepest failing candidate, ban its failing leg's target, re-plan the step once and walk again
        if not ok_cap and J.get('chain_repair', True) and time.time() < t_end - 600:
            fails = [r for r in rs if r.get('failed') and r.get('tour') is not None]
            if fails:
                top = max(fails, key=lambda r: (r['n'], -r['pt']))
                try:
                    lo, hi, leg, bi = SB.failing_leg(top['tour'], max_settles=J['repair_settles'], timeout=J['settle_timeout'])
                except Exception as e:
                    lo, hi, leg, bi = 0, None, None, dict(st_ok=None, settles=0); log(f'  step {k} repair bisect failed: {e}')
                step['repair'] = dict(n=top['n'], k_ok=lo, k_fail=hi, leg=None if leg is None else int(leg['ast']))
                log(f'  step {k} repair: deepest failing {top["n"]} fb; last settling prefix {lo}, failing leg {hi} '
                    f'({None if leg is None else leg["ast"]})')
                if bi.get('st_ok') is not None and lo >= 6:
                    stp = bi['st_ok']; tk = float(RI.ipr(stp).tank())
                    save_col(outdir, tagk + 'pre', stp, len(stp['asts']), C['cols'],
                             dict(kind='chain_prefix', step=k, tank=round(tk, 2)))
                if leg is not None and time.time() < t_end - 400:
                    ban = int(leg['ast']); rec2 = []
                    for m in mems[:2]:
                        pl = plan(pool, S, J['p'], m, J['beam'], J['seed'] + 17 * k + 1, J['jit'], J['dvS'], J['front'], (ban,))
                        log(f'  step {k} repair plan {m} (ban {ban}): deepest {pl["deepest"]}; front ' +
                            ' '.join(f'{c["n"]}:{c["pt"]:.0f}{"*" * len(c["s_taken"])}' for c in pl['front'][:6]) + f' ({pl["sec"]:.0f} s)')
                        rec2 += pl['front']
                    seen2 = set(); c2 = []
                    for c in sorted(rec2, key=lambda c: (-c['n'], -len(c['s_taken']), c['pt'])):
                        if tuple(c['targets']) not in seen2:
                            seen2.add(tuple(c['targets'])); c2.append(c)
                    rs2 = walk(c2, dict(Jw, max_settles=8), log, t_end)
                    for r in rs2:
                        if r['ok']:
                            r['member'] = r['member'] + 'r'
                            save_col(outdir, tagk + 'r', r['res']['st'], r['n'], C['cols'],
                                     dict(kind='chain', step=k, member=r['member'], planner=round(r['pt'], 1),
                                          tank=round(r['res']['tank'], 2), miss=round(r['res']['miss'], 1), s_taken=r['s_taken']))
                    rs = rs + rs2; oks = [r for r in rs if r['ok']]
                    ok_cap = [r for r in oks if r['res']['tank'] <= Jw['cap']]
                    step['repair'].update(ban=ban, ok=len([r for r in rs2 if r['ok']]))
        if not ok_cap:
            step['took'] = None; C['steps'].append(step); jdump(C, ck)
            log(f'  step {k}: nothing settled under the cap -> chain ends'); break
        # an isolated target taken now saves the sparser later steps a hard target: it counts w_iso extra
        wiso = J.get('w_iso', 1.0)
        took = max(ok_cap, key=lambda r: (r['n'] + wiso * len(set(r['targets']) & set(ISO)),
                                          value(r['n'], r['res']['tank'], J['kappa'])))
        st = took['res']['st']; tank = took['res']['tank']
        step['took'] = dict(member=took['member'], n=took['n'], tank=round(tank, 1), s_taken=took['s_taken'])
        # absorb the rest of S (twin-priced) into the taken route, before the next step
        rest = [x for x in S if x not in set(took['targets'])]
        if rest and time.time() < t_end - 300:
            res = absorb_targets(st, tank, rest, J, log, t_end, tag=f'[s{k}]')
            for a in res:
                save_col(outdir, tagk + f'a{a["ast"]}', a['st'], len(a['st']['asts']), C['cols'],
                         dict(kind='chain_absorbed', step=k, host_n=took['n'], ast=a['ast'], dkg=round(a['dkg'], 1),
                              lin_kg=round(a['lin_kg'], 1), tank=round(a['tank'], 2), miss=round(a['miss'], 1)))
            if res:
                b = min(res, key=lambda a: a['dkg'])
                if b['dkg'] <= J.get('chain_absorb_take_kg', 40.0):
                    st = b['st']; tank = b['tank']
                    step['took'].update(absorbed=b['ast'], tank=round(tank, 1), dkg=round(b['dkg'], 1))
        # grow the taken route by CHEAP nearby residual targets (insertion into a beam-born host at d < 0.05 AU:
        # 80 % settle at a median 17.5 kg, stage 13), so the sparser later steps inherit fewer targets
        if J.get('step_grow', 0) > 0 and time.time() < t_end - 300:
            rest_pool = [x for x in pool if x not in set(int(a) for a in st['asts'])]
            st2, tank2, got, gsteps = grow(st, tank, rest_pool, J, log, t_end, J.get('step_grow_d', 0.05),
                                           J.get('step_grow_kg', 25.0), J['step_grow'], tag=f'[s{k} grow]')
            for g in gsteps:
                save_col(outdir, tagk + 'g' + ''.join(f'x{x}' for x in g['base'] + [g['ast']]), g['st'],
                         len(g['st']['asts']), C['cols'], dict(kind='chain_grown', step=k, ast=g['base'] + [g['ast']],
                                                               dkg=round(g['dkg'], 1), tank=round(g['tank'], 2)))
            if got:
                step['took'].update(grown=got, tank=round(tank2, 1)); st, tank = st2, tank2
        C['covered'] = sorted(set(C['covered']) | set(int(a) for a in st['asts']))
        step['covered'] = len(C['covered'])
        C.setdefault('taken', []).append(dict(k=k, n=len(st['asts']), tank=round(tank, 2),
                                              targets=sorted(int(a) for a in st['asts'])))
        C['steps'].append(step); jdump(C, ck)
        log(f'  step {k}: TOOK {step["took"]}; chain covered {len(C["covered"])}/{len(pool0)}')
    meta = dict(job={kk: v for kk, v in J.items() if kk != 'pool'}, kind='chain', pool_n=len(pool0), chain=C,
                cols=C['cols'], wall_s=round(time.time() - tic), done=now())
    jdump(meta, outdir / 'meta.json')
    log(f'chain {J["tag"]} done: {len(C["cols"])} columns, covered {len(C["covered"])}/{len(pool0)} ({time.time() - tic:.0f} s)')
    return meta


def run_absorb(J):
    """Archive absorption: J['items'] = [(route file, target X)]: twin insertion of X into existing (beam-born) columns
    that pass close to X's window, screened by lin_price (absorb_targets).  Every settled insertion is a new column."""
    J = dict(DEFAULTS, **J)
    outdir = ROOT / J['out']; outdir.mkdir(parents=True, exist_ok=True)
    if (outdir / 'meta.json').exists():
        return json.load(open(outdir / 'meta.json'))
    log = Log(outdir / 'log.txt'); tic = time.time(); t_end = tic + J['wall_budget']
    cols = []; tried = []
    for f, X in J['items']:
        if time.time() > t_end:
            break
        z = np.load(ROOT / f); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL'])
        tank = float(RI.ipr(st).tank())
        res = absorb_targets(st, tank, [int(X)], dict(J, absorb_trials=J.get('item_trials', 2)), log, t_end,
                             tag=f'[{pathlib.Path(f).parent.name}/{pathlib.Path(f).stem}]')
        tried.append(dict(f=f, X=int(X), n=len(st['asts']), tank=round(tank, 1), ok=len(res)))
        for a in res:
            save_col(outdir, f'{J["tag"]}x{a["ast"]}', a['st'], len(a['st']['asts']), cols,
                     dict(kind='archive_absorbed', src=f, ast=a['ast'], host_n=len(st['asts']), host_tank=round(tank, 2),
                          dkg=round(a['dkg'], 1), lin_kg=round(a['lin_kg'], 1), tank=round(a['tank'], 2), miss=round(a['miss'], 1)))
    meta = dict(job={k: v for k, v in J.items() if k != 'items'}, kind='absorb', items=tried, cols=cols,
                wall_s=round(time.time() - tic), done=now())
    jdump(meta, outdir / 'meta.json')
    log(f'absorb {J["tag"]} done: {len(cols)} columns from {len(tried)} items ({time.time() - tic:.0f} s)')
    return meta


def run_close(J):
    """Closer on a picked fleet: J['items'] = [(route file, [misses near it])].  Per host, greedy sequential twin
    absorption (lin_price screen, w_insert), re-screening the remaining misses on the grown route after each success,
    at most J['max_abs'] absorptions; every intermediate route is saved as a column (the selector picks)."""
    J = dict(DEFAULTS, **J)
    outdir = ROOT / J['out']; outdir.mkdir(parents=True, exist_ok=True)
    if (outdir / 'meta.json').exists():
        return json.load(open(outdir / 'meta.json'))
    log = Log(outdir / 'log.txt'); tic = time.time(); t_end = tic + J['wall_budget']
    cols = []; hosts = []
    for f, Xs in J['items']:
        if time.time() > t_end:
            break
        z = np.load(ROOT / f); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL'])
        tank0 = tank = float(RI.ipr(st).tank()); got = []
        left = [int(x) for x in Xs]
        while left and len(got) < J.get('max_abs', 3) and time.time() < t_end:
            res = absorb_targets(st, tank, left, dict(J, absorb_trials=J.get('item_trials', 3)), log, t_end,
                                 tag=f'[{pathlib.Path(f).stem}+{got}]')
            if not res:
                break
            b = min(res, key=lambda a: a['dkg'])
            st, tank = b['st'], b['tank']; got.append(b['ast']); left.remove(b['ast'])
            save_col(outdir, f'{J["tag"]}{pathlib.Path(f).stem[-2:]}' + ''.join(f'x{g}' for g in got), st,
                     len(st['asts']), cols, dict(kind='closed', src=f, ast=list(got), host_tank=round(tank0, 2),
                                                 dkg=round(tank - tank0, 1), tank=round(tank, 2), miss=round(b['miss'], 1)))
        hosts.append(dict(f=f, tried=[int(x) for x in Xs], got=got, dkg=round(tank - tank0, 1)))
    meta = dict(job={k: v for k, v in J.items() if k != 'items'}, kind='close', hosts=hosts, cols=cols,
                wall_s=round(time.time() - tic), done=now())
    jdump(meta, outdir / 'meta.json')
    log(f'close {J["tag"]} done: {len(cols)} columns ({time.time() - tic:.0f} s)')
    return meta


def run_grow(J):
    """Grow existing routes by cheap nearby targets: J['items'] = [(route file, [candidate targets])]; per item
    grow() (lin_price screen + twin insertion, within J['grow_d'] AU, <= J['grow_kg'] kg each, <= J['grow_n']); every
    intermediate route is a column; meta 'final' maps each source file to its final grown column file."""
    J = dict(DEFAULTS, **J)
    outdir = ROOT / J['out']; outdir.mkdir(parents=True, exist_ok=True)
    if (outdir / 'meta.json').exists():
        return json.load(open(outdir / 'meta.json'))
    log = Log(outdir / 'log.txt'); tic = time.time(); t_end = tic + J['wall_budget']
    cols = []; final = {}; hosts = []
    for f, Xs in J['items']:
        z = np.load(ROOT / f); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL'])
        tank0 = float(RI.ipr(st).tank())
        st2, tank2, got, steps = grow(st, tank0, Xs, J, log, t_end, J.get('grow_d', 0.05), J.get('grow_kg', 25.0),
                                      J.get('grow_n', 4), tag=f'[{pathlib.Path(f).stem}]')
        last = None
        for g in steps:
            pth = save_col(outdir, f'{J["tag"]}' + ''.join(f'x{x}' for x in g['base'] + [g['ast']]), g['st'],
                           len(g['st']['asts']), cols, dict(kind='grown', src=f, ast=g['base'] + [g['ast']],
                                                            dkg=round(g['tank'] - tank0, 1), tank=round(g['tank'], 2)))
            if not g['extra']:
                last = pth
        final[f] = str(last.relative_to(ROOT)) if last is not None else f
        hosts.append(dict(f=f, got=got, dkg=round(tank2 - tank0, 1), tank0=round(tank0, 1)))
    meta = dict(job={k: v for k, v in J.items() if k != 'items'}, kind='grow', hosts=hosts, final=final, cols=cols,
                wall_s=round(time.time() - tic), done=now())
    jdump(meta, outdir / 'meta.json')
    log(f'grow {J["tag"]} done: {hosts} ({time.time() - tic:.0f} s)')
    return meta


def run_job(J):
    """Pool worker entry: never raises (an error is written to OUT/error.txt and returned)."""
    try:
        RI.eph()
        k = J.get('kind')
        fn = dict(chain=run_chain, absorb=run_absorb, close=run_close, grow=run_grow).get(k, run_column)
        return fn(J)
    except Exception:
        outdir = ROOT / J['out']; outdir.mkdir(parents=True, exist_ok=True)
        (outdir / 'error.txt').write_text(traceback.format_exc())
        return dict(error=traceback.format_exc()[-400:], tag=J.get('tag'))


if __name__ == '__main__':
    J = json.load(open(sys.argv[1]))
    r = run_job(J)
    print(json.dumps({k: v for k, v in r.items() if k in ('wall_s', 'error', 'absorbed', 'S')}, default=_jd)[:2000])
