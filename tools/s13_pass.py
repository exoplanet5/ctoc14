"""Stage 13 (PCC-8G) [P]: ONE covering pass = N sequential FREE beams over the uncovered pool (never pinned), with the
tail portfolio for routes >= --tail-from.  docs/stage13_solver_plan.md sections 3.1, 3.2 and 4.1.
Built from results/s13/seed_passes/scripts/probe_pass5.py (head routes) and results/s13/plan13/tail_probe.py (tail
portfolio); joint member T7 from results/s13/planner_probe/{jointsym,twin_wait}.py.  Shared tools are imported only.

usage: s13_pass.py OUT [--start FLEET --keep 01,02,03,04,05] [--prem PREM.json] [--delta 0.01] [--dvt 1.8]
                   [--beam 100] [--npt 6] [--mc 150] [--tail-from 6] [--members T1,T2,T3,T5,T6,T7] [--kappa 90]
                   [--cap-head 1050] [--cap-tail 1080] [--grid 0,800,20] [--pilot 0|M] [--seed S] [--jitter 0]
                   [--nproc 4] [--no-repair] [--n 8] [--cache results/s13/gen/cache]

Per route k:
  1. plan: every member beam (search.beam_search, m0 1600, strict exclusion of covered targets) returns its front =
     min-fuel state per depth, deepest 8.  Head routes (k < tail-from) use member H (= T1 legs, --npt/--mc/--grid);
     tail routes use the portfolio T1 T2 T3 T5 T6 (T4 = T1 was dropped, M-D); T7 = joint K=2 beam over the pool of
     route N-1 (the last pair), settled with the wait-aware seed, taken only if its pair beats sequential N-1 + N.
  2. walk: candidates with planner tank <= cap, deepest first (ties: lower planner tank), batches of 6, at most 18
     settles; stop at the first batch with a SETTLED candidate (miss <= 150 km, twin <= 1.15 x planner, twin <= cap).
  3. choice: the deepest settled candidate of that batch; ties by value = n - (tank - 600)/kappa (= lower tank).
  4. repair (bisect-ban): if the deepest candidate FAILED to settle (not merely over the cap) and nothing as deep
     settled, s13_bisect.failing_leg finds the breaking leg; its target is banned for this route and that member is
     re-planned once; its candidates deeper than the choice are settled (one batch) and replace it if one settles.
  5. pilot (--pilot M, routes 3-6): the M best diverse settled candidates (<= 85% overlap) are rolled out to route N
     with beam-30 T1 npt-2 planner-only beams; the one with the most rollout coverage (tie: lower sum J_i) is taken.
  6. library: every other settled candidate -> OUT/lib/route_<k>_<member>_<n>.npz, all listed in OUT/cols.jsonl.
Premiums: prize_t = clip(1 + lambda_t, 1.00, 1.02); lambda = PREM.json (or --delta on the gc16 hard set) + seeded
uniform jitter, lambda clipped to [-0.01, 0.02]; legs to the premium set (base lambda > 0) get dv cap max(dv, --dvt).
Restart-safe: pass.json + route files are the state; every beam / settle / rollout result is cached on disk."""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse, hashlib, warnings, multiprocessing as mp
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools'), str(ROOT / 'results/s13/planner_probe')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
_CWD0 = os.getcwd()
os.chdir(ROOT)
import numpy as np
import s13_bisect as SB
import run_ialns as RI
import ctoc14.impulsive as IM
from greedy_cover import BP, CM, DeepCollect, ALL
from ctoc14.search import beam_search, tour_json
from ctoc14.constants import DAY, AU
warnings.filterwarnings('ignore')
IM.LIN_FALLBACK = 2.5

MEMBERS = {'H':  dict(dv=1.2, dr=0.15, npt=6, mc=150, g=800, wt=1.0),
           'T1': dict(dv=1.2, dr=0.15, npt=6, mc=150, g=800, wt=1.0),
           'T2': dict(dv=1.6, dr=0.20, npt=6, mc=150, g=800, wt=1.0),
           'T3': dict(dv=2.0, dr=0.25, npt=2, mc=60, g=800, wt=1.0),
           'T4': dict(dv=1.2, dr=0.15, npt=6, mc=150, g=1500, wt=1.0),
           'T5': dict(dv=1.6, dr=0.20, npt=6, mc=150, g=800, wt=0.5),
           'T6': dict(dv=2.0, dr=0.25, npt=6, mc=150, g=1500, wt=1.0),
           'T7': dict(dv=1.6, dr=0.20, npt=2, mc=60, g=800, wt=1.0, joint=2, wait=60.0, wrare=1.0),
           'R':  dict(dv=1.2, dr=0.15, npt=2, mc=60, g=800, wt=1.0)}
ORDER = ['T7', 'T6', 'T3', 'T5', 'T2', 'T4', 'T1', 'H']        # expected cost, longest first (pool packing)
M0 = 1600.0
FRONT = 8
BATCH = 6
MAX_SETTLES = 18
OVERLAP = 0.85


def value(n, tank, kappa):
    return n - (tank - 600.0) / kappa


def relp(p):
    p = pathlib.Path(p)
    return p.relative_to(ROOT) if str(p).startswith(str(ROOT) + '/') else p


def overlap(a, b):
    a, b = set(a), set(b)
    return len(a & b) / max(len(a), len(b), 1)


def hard_set():
    gc = RI.IFleet(ROOT / 'results/newgen/gc16'); cov = set()
    for n, r in gc.routes.items():
        cov |= set(int(x) for x in r['st']['asts'])
    return set(ALL) - cov


# -------------------------------------------------------------------------------------------------- planner jobs
def _bp(m, beam, prize, prem, dvt):
    P = BP(beam=beam, w_fuel=CM.w_fuel(480, 1600), m_margin=40.0, dv_max=m['dv'], tofs=np.arange(15, 401, 5) * DAY,
           lin_tofs=np.arange(20, 601, 10) * DAY, lin_drmax=m['dr'] * AU, vinf_cap=4.0, tof_refine=True, max_depth=60,
           n_per_target=m['npt'], max_children=m['mc'], w_t=m['wt'])
    P.prize = np.asarray(prize, float)
    if prem and dvt and dvt > m['dv']:
        d = np.full(300, m['dv'])
        for t in prem:
            d[t - 1] = max(m['dv'], dvt)
        P.dv_max_t = d
    return P


def _prize_vec(avail, prize_full):
    pz = np.zeros(300)
    for t in avail:
        pz[t - 1] = prize_full[t - 1]
    return pz


def beam_key(job):
    avail = job['avail']
    return SB.Cache.key(dict(v=SB.CODE_TAG, kind=job['kind'], m=job['m'], beam=job['beam'], avail=avail,
                             prize=[round(float(job['prize'][t - 1]), 12) for t in avail],
                             prem=sorted(set(job['prem']) & set(avail)), dvt=job['dvt'],
                             **({'covered': job['covered'], 'k0': job['k0'], 'N': job['N'], 'caps': job['caps']}
                                if job['kind'] == 'rollout' else {})))


def run_member(job):
    """One single-route member beam over `avail`.  Returns dict(member, front=[(n, planner tank, tour)], extra, sec)."""
    m = job['m']; avail = job['avail']; tic = time.time()
    P = _bp(m, job['beam'], _prize_vec(avail, job['prize']), job['prem'], job['dvt'])
    P.collect = DeepCollect(4)
    best, beam = beam_search(RI.eph(), excluded=sorted(set(ALL) - set(avail)), m0=M0,
                             t_launch_grid=np.arange(0, m['g'] + 1e-9, 20) * DAY, P=P, n_proc=1, verbose=False)
    states = list(P.collect) + list(beam)
    byn = {}
    for s in states:
        k = len(s.seq)
        if k not in byn or s.fuel < byn[k].fuel:
            byn[k] = s
    front = [(k, float(CM.tank(byn[k].fuel, M0)), tour_json(byn[k])) for k in sorted(byn, reverse=True)[:FRONT]]
    # diverse extras for the pilot: deepest levels, <= 85% target overlap with every kept state
    extra = []; sets = []
    if byn:
        nmax = max(byn)
        for s in sorted((s for s in states if len(s.seq) >= nmax - 3), key=lambda s: (-len(s.seq), s.fuel)):
            ts = [x[0] for x in s.seq]
            if any(overlap(ts, o) > OVERLAP for o in sets):
                continue
            sets.append(ts); extra.append((len(s.seq), float(CM.tank(s.fuel, M0)), tour_json(s)))
            if len(extra) >= 12:
                break
    return dict(member=job['member'], front=front, extra=extra, sec=time.time() - tic)


def run_joint(job):
    """Member T7: joint K=2 beam (jointsym, order-free key, wait 60 d, w_rare 1) over `avail`.  Returns the deepest
    3 coverage levels of the harvested joint states, each as its min planner-sum-J_i pair of tours."""
    import jointsym
    m = job['m']; avail = job['avail']; tic = time.time()
    P = _bp(m, job['beam'], _prize_vec(avail, job['prize']), job['prem'], job['dvt'])
    P.max_depth = 120
    P.w_rare = m.get('wrare', 1.0); P.wait = m.get('wait', 60.0) * DAY; P.collect_joint = []
    K = int(m.get('joint', 2))
    best = jointsym.joint_beam_search(RI.eph(), [M0] * K, np.arange(0, m['g'] + 1e-9, 20) * DAY, P=P,
                                      excluded=sorted(set(ALL) - set(avail)), n_proc=1, verbose=False, canon=True)
    front = {}
    for j in list(P.collect_joint) + [best]:
        rs = [s for s in j.craft if s is not None and s.n() > 0]
        if len(rs) < K:
            continue
        cov = len(set().union(*[{x[0] for x in s.seq} for s in rs]))
        tanks = [float(CM.tank(s.fuel, s.m0)) for s in rs]
        sj = sum(RI.cost(t) for t in tanks)
        if cov not in front or sj < front[cov][0]:
            front[cov] = (sj, [tour_json(s) for s in rs], tanks)
    levels = [dict(cov=c, sumJ=front[c][0], tours=front[c][1], ptanks=front[c][2]) for c in sorted(front, reverse=True)[:3]]
    return dict(member='T7', levels=levels, sec=time.time() - tic)


def run_rollout(job):
    """Pilot rollout: planner-only sequential beams (member R, beam 30) for routes k0+1..N from `covered`."""
    m = job['m']; covered = set(job['covered']); tic = time.time(); depths = []; sj = 0.0
    for k in range(job['k0'] + 1, job['N'] + 1):
        avail = sorted(set(ALL) - covered)
        if not avail:
            break
        cap = job['caps'][0] if k < job['caps'][2] else job['caps'][1]
        P = _bp(m, job['beam'], _prize_vec(avail, job['prize']), job['prem'], job['dvt'])
        P.collect = DeepCollect(4)
        best, beam = beam_search(RI.eph(), excluded=sorted(covered), m0=M0,
                                 t_launch_grid=np.arange(0, m['g'] + 1e-9, 20) * DAY, P=P, n_proc=1, verbose=False)
        st = [s for s in list(P.collect) + list(beam) if CM.tank(s.fuel, M0) <= cap]
        if not st:
            break
        s = max(st, key=lambda s: (len(s.seq), -s.fuel))
        covered |= {x[0] for x in s.seq}; depths.append(len(s.seq)); sj += RI.cost(float(CM.tank(s.fuel, M0)))
    return dict(cov=len(covered), sumJ=sj, depths=depths, sec=time.time() - tic)


def w_job(job):
    k = job['kind']
    if k == 'beam':
        return run_member(job)
    if k == 'joint':
        return run_joint(job)
    if k == 'rollout':
        return run_rollout(job)
    if k == 'settle':
        return SB.w_settle(job)
    raise ValueError(k)


def run_jobs(jobs, pmap, cache):
    """Planner jobs through the cache (main process reads/writes it), misses in parallel, longest first."""
    out = [None] * len(jobs); todo = []
    for i, j in enumerate(jobs):
        j['_key'] = beam_key(j)
        v = cache.get(j['kind'], j['_key'])
        if v is not None:
            v = dict(v); v['cached'] = True; out[i] = v
        else:
            todo.append(i)
    rank = {m: r for r, m in enumerate(ORDER)}
    todo.sort(key=lambda i: rank.get(jobs[i]['member'], 99))
    if todo:
        res = pmap(w_job, [jobs[i] for i in todo])
        for i, r in zip(todo, res):
            cache.put(jobs[i]['kind'], jobs[i]['_key'], r)
            r = dict(r); r['cached'] = False; out[i] = r
    return out


# ---------------------------------------------------------------------------------------------------------- pass
class Pass:
    def __init__(self, a):
        self.a = a
        self.out = pathlib.Path(a.out) if os.path.isabs(a.out) else pathlib.Path(_CWD0) / a.out
        self.out.mkdir(parents=True, exist_ok=True); (self.out / 'lib').mkdir(exist_ok=True)
        self.say = RI.logger(self.out / 'log.txt')
        self.cache = SB.Cache(a.cache if os.path.isabs(a.cache) else ROOT / a.cache)
        self.members = [m.strip() for m in a.members.split(',') if m.strip()]
        self.tail_members = [m for m in self.members if m != 'T7']
        g = [float(x) for x in a.grid.split(',')]
        self.head_m = dict(MEMBERS['H'], npt=a.npt, mc=a.mc, g=g[1])
        if g[0] != 0 or g[2] != 20:
            raise SystemExit('--grid must be 0,G,20 (the tail probe grid)')
        self.pool = None; self.t0 = time.time()

    # -------------------------------------------------------------------- premiums
    def premiums(self):
        a = self.a
        base = np.zeros(300)
        if a.prem:
            d = json.load(open(a.prem if os.path.isabs(a.prem) else pathlib.Path(_CWD0) / a.prem))
            d = d.get('lambda', d) if isinstance(d, dict) else d
            if isinstance(d, dict):
                for t, v in d.items():
                    base[int(t) - 1] = float(v)
            else:
                base[:] = np.asarray(d, float)[:300]
        elif a.delta:
            for t in hard_set():
                base[t - 1] = a.delta
        base = np.clip(base, -0.01, 0.02)
        lam = base.copy()
        if a.jitter > 0:
            rng = np.random.default_rng(a.seed)
            lam = np.clip(lam + rng.uniform(-a.jitter, a.jitter, 300), -0.01, 0.02)
        self.lam = lam
        self.prize = np.clip(1.0 + lam, 1.00, 1.02)
        self.prem = sorted(t for t in ALL if base[t - 1] > 1e-9)
        self.lam_digest = dict(n_prem=len(self.prem), n_pos=int((lam > 1e-9).sum()), n_neg=int((lam < -1e-9).sum()),
                               mean=round(float(lam.mean()), 5), max=round(float(lam.max()), 4),
                               prize_min=round(float(self.prize[np.array(ALL) - 1].min()), 4),
                               prize_max=round(float(self.prize[np.array(ALL) - 1].max()), 4),
                               sha1=hashlib.sha1(np.round(lam, 9).tobytes()).hexdigest()[:12],
                               src=a.prem or f'delta {a.delta} on gc16 hard set', jitter=a.jitter, seed=a.seed)

    def pmap(self, f, jobs):
        if not jobs:
            return []
        if self.pool is None:
            return [f(j) for j in jobs]
        return self.pool.map(f, jobs, chunksize=1)

    def job(self, kind, member, avail, m=None, beam=None, **kw):
        a = self.a
        m = m or (self.head_m if member == 'H' else MEMBERS[member])
        j = dict(kind=kind, member=member, m=m, beam=beam or a.beam, avail=list(avail), prize=self.prize.tolist(),
                 prem=self.prem, dvt=a.dvt)
        j.update(kw)
        return j

    def save_state(self):
        json.dump(self.rec, open(self.out / 'pass.json.tmp', 'w'), indent=1, default=SB._jdefault)
        os.replace(self.out / 'pass.json.tmp', self.out / 'pass.json')

    def add_col(self, path, k, member, n, tank, planner, miss, targets, chosen):
        with open(self.out / 'cols.jsonl', 'a') as f:
            f.write(json.dumps(dict(f=str(path), k=k, member=member, n=n, tank=round(tank, 3), planner=round(planner, 2),
                                    miss=round(miss, 2), chosen=chosen, targets=sorted(targets))) + '\n')

    def lib_save(self, k, c, tag=''):
        p = self.out / 'lib' / f'route_{k:02d}_{c["member"]}{tag}_{c["n"]}.npz'
        i = 1
        while p.exists():
            p = self.out / 'lib' / f'route_{k:02d}_{c["member"]}{tag}_{c["n"]}_{i}.npz'; i += 1
        np.savez(p, **c['res']['st'])
        self.add_col(relp(p), k, c['member'], c['n'], c['res']['tank'],
                     c['pt'], c['res']['miss'], [int(x) for x in c['res']['st']['asts']], False)
        return p

    # -------------------------------------------------------------------- settle walk
    def settle_batch(self, cands, cap):
        jobs = [dict(kind='wait' if c.get('wait') else 'stock', tour=c['tour'], id=i, timeout=self.a.settle_timeout)
                for i, c in enumerate(cands)]
        res = SB.settle_many(jobs, self.pmap, self.cache)
        for c, r in zip(cands, res):
            c['res'] = r; c['ok'] = SB.is_ok(r, c['pt'], cap); c['failed'] = SB.failed_settle(r, c['pt'])
            self.say(f'    {c["member"]} {c["n"]} fb planner {c["pt"]:.0f} -> twin '
                     f'{r["tank"] if r["tank"] is None else round(r["tank"])} (miss {r["miss"]:.0f}) '
                     f'{"OK" if c["ok"] else ("over cap" if not c["failed"] else "fail")} ({r["sec"]:.0f} s'
                     f'{", cached" if r.get("cached") else ""})')
        return cands

    def walk(self, cands, cap):
        tries = []; took = None; i = 0
        while i < len(cands) and took is None and len(tries) < MAX_SETTLES:
            b = cands[i:i + BATCH]; i += len(b)
            self.settle_batch(b, cap); tries += b
            ok = [c for c in b if c['ok']]
            if ok:
                took = max(ok, key=self.choice_key)
        return took, tries

    def choice_key(self, c):
        """'lex' (default, plan 4.1): the deepest settled candidate, ties by value = n - (tank - 600)/kappa (= lower tank);
        'value': the settled Pareto point of maximum value (an extra flyby must cost less than kappa kg)."""
        v = value(c['n'], c['res']['tank'], self.a.kappa)
        return (v,) if self.a.choice == 'value' else (c['n'], v)

    def cand_list(self, results, cap, extra=False):
        cands = []; seen = set()
        for r in results:
            for n, pt, tour in (r['front'] + (r['extra'] if extra else [])):
                if pt > cap:
                    continue
                key = (round(tour['t_launch'] / DAY, 4),) + tuple((l['ast'], round(l['t_flyby'] / DAY, 4)) for l in tour['legs'])
                if key in seen:
                    continue
                seen.add(key)
                cands.append(dict(n=n, pt=pt, tour=tour, member=r['member'],
                                  targets=sorted(int(l['ast']) for l in tour['legs'])))
        cands.sort(key=lambda c: (-c['n'], c['pt']))
        return cands

    # -------------------------------------------------------------------- one route
    def route(self, k, fl):
        a = self.a; tic = time.time()
        cap = a.cap_head if k < a.tail_from else a.cap_tail
        members = ['H'] if k < a.tail_from else self.tail_members
        joint = ('T7' in self.members and k == a.n - 1 and k >= a.tail_from and self.rec.get('t7') is None)
        avail = sorted(set(ALL) - set(fl.coverage()))
        jobs = [self.job('beam', m, avail) for m in members]
        if joint:
            jobs.append(self.job('joint', 'T7', avail))
        res = run_jobs(jobs, self.pmap, self.cache)
        single = [r for r in res if r['member'] != 'T7']
        fronts = {}
        for r in single:
            fronts[r['member']] = [(n, round(pt)) for n, pt, _ in r['front'][:5]]
            self.say(f'  route {k:02d} {r["member"]}: pool {len(avail)}, front ' +
                     ' '.join(f'{n}:{pt:.0f}' for n, pt, _ in r['front'][:5]) +
                     f' ({r["sec"]:.0f} s{", cached" if r.get("cached") else ""})')
        t7res = next((r for r in res if r['member'] == 'T7'), None)
        if t7res is not None:
            self.say(f'  route {k:02d} T7 (joint K=2): pool {len(avail)}, levels ' +
                     ' '.join(f'{L["cov"]}={"+".join(str(len(t["legs"])) for t in L["tours"])}'
                              f'@{"/".join(f"{p:.0f}" for p in L["ptanks"])}' for L in t7res['levels']) +
                     f' ({t7res["sec"]:.0f} s{", cached" if t7res.get("cached") else ""})')
        cands = self.cand_list(single, cap)
        rec = dict(k=k, cap=cap, pool=len(avail), members=members, fronts=fronts,
                   beam_s={r['member']: round(r['sec']) for r in res})
        took, tries = self.walk(cands, cap)
        rec['walk'] = dict(took=None if took is None else f'{took["member"]}:{took["n"]}', tries=len(tries))
        # ---- repair: the deepest candidate failed and nothing as deep settled
        if a.repair and cands and cands[0].get('failed') and (took is None or took['n'] < cands[0]['n']):
            took = self.repair(k, cands[0], took, avail, cap, rec, tries)
        # ---- pilot
        if a.pilot > 0 and 3 <= k <= 6 and k < a.n and took is not None:
            took = self.pilot(k, took, tries, cands, single, fl, cap, rec)
        oks = [c for c in tries if c.get('ok')]
        if took is None:
            self.say(f'route {k:02d}: nothing settled ({len(tries)} settles)')
            rec.update(n=0, sec=round(time.time() - tic)); return None, rec, t7res, avail
        for c in oks:
            if c is not took:
                self.lib_save(k, c)
        rec.update(member=took['member'], n=took['n'], planner=round(took['pt'], 1), twin=round(took['res']['tank'], 2),
                   miss=round(took['res']['miss'], 1), targets=sorted(int(x) for x in took['res']['st']['asts']),
                   tries_detail=[(c['member'], c['n'], round(c['pt']), None if c['res']['tank'] is None else round(c['res']['tank']),
                                  round(c['res']['miss']) if np.isfinite(c['res']['miss']) else None, bool(c['ok']))
                                 for c in tries],
                   settles=len(tries), sec=round(time.time() - tic))
        return took, rec, t7res, avail

    def repair(self, k, deep, took, avail, cap, rec, tries):
        a = self.a; tic = time.time()
        self.say(f'  route {k:02d}: REPAIR -- deepest {deep["member"]} {deep["n"]} failed; bisecting')
        width = max(1, a.nproc)
        lo, hi, leg, info = SB.failing_leg(deep['tour'], max_settles=3 * width if width > 1 else 6, pmap=self.pmap,
                                           width=width, cache=self.cache, say=self.say, timeout=a.settle_timeout)
        rep = dict(member=deep['member'], n=deep['n'], k_ok=lo, k_fail=hi, leg=leg, settles=info['settles'])
        if info['st_ok'] is not None and lo >= 4:
            p = self.out / 'lib' / f'route_{k:02d}_{deep["member"]}pre_{lo}.npz'
            np.savez(p, **info['st_ok'])
            ip = RI.ipr(info['st_ok']); tank = float(ip.tank())
            self.add_col(relp(p), k, deep['member'] + 'pre', lo, tank, float('nan'), float('nan'),
                         [int(x) for x in info['st_ok']['asts']], False)
        if leg is None:
            self.say(f'  route {k:02d}: repair found no failing leg (probes {info["probes"]})')
            rec['repair'] = rep; return took
        ban = int(leg['ast'])
        self.say(f'  route {k:02d}: failing leg {hi} (target {ban}); last settling prefix {lo}; banning {ban} and '
                 f're-planning {deep["member"]} ({info["settles"]} settles, {time.time() - tic:.0f} s)')
        av2 = [t for t in avail if t != ban]
        r2 = run_jobs([self.job('beam', deep['member'], av2)], self.pmap, self.cache)[0]
        self.say(f'  route {k:02d} {deep["member"]}-ban{ban}: front ' + ' '.join(f'{n}:{pt:.0f}' for n, pt, _ in r2['front'][:5]) +
                 f' ({r2["sec"]:.0f} s{", cached" if r2.get("cached") else ""})')
        nmin = took['n'] if took is not None else 0
        c2 = [c for c in self.cand_list([r2], cap) if c['n'] > nmin][:BATCH]
        for c in c2:
            c['member'] = deep['member'] + 'r'
        if c2:
            self.settle_batch(c2, cap); tries += c2
            ok = [c for c in c2 if c['ok']]
            if ok:
                new = max(ok, key=self.choice_key)
                was = 'none' if took is None else f'{took["n"]} @ {took["res"]["tank"]:.0f}'
                self.say(f'  route {k:02d}: repair TOOK {new["member"]} {new["n"]} @ {new["res"]["tank"]:.0f} (was {was})')
                rep.update(took=f'{new["member"]}:{new["n"]}', tank=round(new['res']['tank'], 1)); took = new
        rep['sec'] = round(time.time() - tic); rec['repair'] = rep
        return took

    def pilot(self, k, took, tries, cands, single, fl, cap, rec):
        a = self.a; tic = time.time(); M = a.pilot
        oks = sorted([c for c in tries if c.get('ok')], key=self.choice_key, reverse=True)
        sel = []
        for c in oks:
            if all(overlap(c['targets'], s['targets']) <= OVERLAP for s in sel):
                sel.append(c)
            if len(sel) >= M:
                break
        # settle more diverse candidates (fronts + extras) until M are settled or the settle budget is spent
        ext = self.cand_list(single, cap, extra=True)
        tried_keys = {tuple(c['targets']) for c in tries}
        while len(sel) < M and len(tries) < MAX_SETTLES + BATCH:
            more = [c for c in ext if tuple(c['targets']) not in tried_keys and c['n'] >= took['n'] - 3 and
                    all(overlap(c['targets'], s['targets']) <= OVERLAP for s in sel)][:BATCH]
            if not more:
                break
            self.settle_batch(more, cap); tries += more
            for c in more:
                tried_keys.add(tuple(c['targets']))
                if c['ok'] and all(overlap(c['targets'], s['targets']) <= OVERLAP for s in sel) and len(sel) < M:
                    sel.append(c)
        if len(sel) <= 1:
            rec['pilot'] = dict(n_sel=len(sel), note='fewer than 2 diverse settled candidates')
            return took
        cov0 = set(fl.coverage())
        jobs = []
        for c in sel:
            cv = sorted(cov0 | set(c['targets']))
            jobs.append(self.job('rollout', 'R', sorted(set(ALL) - set(cv)), m=MEMBERS['R'], beam=30, covered=cv, k0=k,
                                 N=a.n, caps=[a.cap_head, a.cap_tail, a.tail_from]))
        res = run_jobs(jobs, self.pmap, self.cache)
        best = None
        for c, r in zip(sel, res):
            key = (r['cov'], -(RI.cost(c['res']['tank']) + r['sumJ']))
            self.say(f'  route {k:02d} pilot: {c["member"]} {c["n"]} @ {c["res"]["tank"]:.0f} -> rollout depths {r["depths"]}, '
                     f'covered {r["cov"]}, planner sum J_i (tail) {r["sumJ"]:.3f} ({r["sec"]:.0f} s)')
            if best is None or key > best[0]:
                best = (key, c)
        rec['pilot'] = dict(n_sel=len(sel), chose=f'{best[1]["member"]}:{best[1]["n"]}', rollout_cov=best[0][0],
                            was=f'{took["member"]}:{took["n"]}', sec=round(time.time() - tic))
        return best[1]

    # -------------------------------------------------------------------- T7
    def t7_settle(self, t7res, cap):
        a = self.a; lv = []
        jobs = []
        for li, L in enumerate(t7res['levels']):
            for ci, (t, pt) in enumerate(zip(L['tours'], L['ptanks'])):
                if pt > cap:
                    continue
                jobs.append(dict(li=li, ci=ci, pt=pt, tour=t, n=len(t['legs'])))
        if not jobs:
            return None
        sj = [dict(kind='wait', tour=j['tour'], id=i, timeout=a.settle_timeout) for i, j in enumerate(jobs)]
        res = SB.settle_many(sj, self.pmap, self.cache)
        for j, r in zip(jobs, res):
            j['res'] = r; j['ok'] = SB.is_ok(r, j['pt'], cap)
        redo = [j for j in jobs if not j['ok']]
        if redo:
            res2 = SB.settle_many([dict(kind='stock', tour=j['tour'], id=i, timeout=a.settle_timeout) for i, j in enumerate(redo)],
                                  self.pmap, self.cache)
            for j, r in zip(redo, res2):
                j['seed'] = 'wait+stock'
                if SB.is_ok(r, j['pt'], cap):
                    j['res'] = r; j['ok'] = True; j['seed'] = 'stock'
        for j in jobs:
            r = j['res']
            self.say(f'    T7 level {t7res["levels"][j["li"]]["cov"]} craft {j["ci"]}: {j["n"]} fb planner {j["pt"]:.0f} -> twin '
                     f'{r["tank"] if r["tank"] is None else round(r["tank"])} (miss {r["miss"]:.0f}) {"OK" if j["ok"] else "fail"}'
                     f' [{j.get("seed", "wait")}]')
        for li, L in enumerate(t7res['levels']):
            js = [j for j in jobs if j['li'] == li]
            if len(js) == len(L['tours']) and all(j['ok'] for j in js):
                return dict(cov=L['cov'], craft=[dict(n=j['n'], pt=j['pt'], tank=j['res']['tank'], miss=j['res']['miss'],
                                                      st=j['res']['st'], seed=j.get('seed', 'wait')) for j in js])
        return None

    # -------------------------------------------------------------------- main loop
    def run(self):
        a = self.a; say = self.say
        self.premiums()
        pj = self.out / 'pass.json'
        self.rec = json.load(open(pj)) if pj.exists() else None
        if self.rec is not None and self.rec.get('done'):
            say(f'pass already done: {self.rec.get("summary")}'); return self.rec
        fl = RI.IFleet()
        if self.rec is None:
            self.rec = dict(args=vars(a), lam=self.lam_digest, kept=[], routes={}, t7=None, done=False,
                            started=time.strftime('%Y-%m-%d %H:%M:%S'))
            if a.start:
                src = RI.IFleet(a.start if os.path.isabs(a.start) else pathlib.Path(_CWD0) / a.start)
                keep = [x.strip() for x in a.keep.split(',') if x.strip()] if a.keep else sorted(src.routes)
                for i, name in enumerate(keep):
                    fl.routes[f'{i + 1:02d}'] = src.routes[name]
                    st = src.routes[name]['st']
                    self.rec['kept'].append(dict(k=i + 1, src=f'{a.start}/route_{name}.npz', n=len(st['asts']),
                                                 twin=round(src.routes[name]['tank'], 2)))
                fl.save(self.out, note=f'pass start: kept {keep} of {a.start}')
                for i, name in enumerate(keep):
                    self.add_col(relp(self.out / f'route_{i + 1:02d}.npz'), i + 1, 'kept',
                                 len(src.routes[name]['st']['asts']), src.routes[name]['tank'], float('nan'), float('nan'),
                                 [int(x) for x in src.routes[name]['st']['asts']], True)
            self.save_state()
            say(f'pass start: {vars(a)}; lambda {self.lam_digest}; kept {fl.summary() if fl.routes else "none"}')
        else:
            fl = RI.IFleet(self.out)
            K0 = len(self.rec['kept'])
            done = {int(x) for x, r in self.rec['routes'].items() if r.get('n')}
            for name in list(fl.routes):          # a route file written after the last pass.json update is re-planned
                if int(name) > K0 and int(name) not in done and not (self.rec.get('t7') or {}).get('taken'):
                    del fl.routes[name]
            say(f'pass RESUME: {fl.summary()}; routes done {sorted(self.rec["routes"])}')
        K0 = len(self.rec['kept'])
        if a.nproc > 1:
            self.pool = mp.get_context('fork').Pool(a.nproc, maxtasksperchild=2)
        try:
            k = K0 + 1
            while k <= a.n:
                if str(k) in self.rec['routes']:
                    k += 1; continue
                took, rec, t7res, avail = self.route(k, fl)
                if t7res is not None:
                    tr = self.t7_settle(t7res, a.cap_tail)
                    self.rec['t7'] = dict(k=k, pool=len(avail), levels=[(L['cov'], [len(t['legs']) for t in L['tours']],
                                                                         [round(p) for p in L['ptanks']]) for L in t7res['levels']],
                                          settled=None)
                    if tr is not None:
                        paths = []
                        for ci, c in enumerate(tr['craft']):
                            p = self.out / 'lib' / f'route_{k:02d}_T7{"ab"[ci]}_{c["n"]}.npz'
                            np.savez(p, **c['st']); paths.append(str(relp(p)))
                            self.add_col(relp(p), k, 'T7', c['n'], c['tank'], c['pt'], c['miss'],
                                         [int(x) for x in c['st']['asts']], False)
                        self.rec['t7']['settled'] = dict(cov=tr['cov'], n=[c['n'] for c in tr['craft']],
                                                         tank=[round(c['tank'], 1) for c in tr['craft']], files=paths)
                        depths = '+'.join(str(c['n']) for c in tr['craft'])
                        tanks = '/'.join(f'{c["tank"]:.0f}' for c in tr['craft'])
                        say(f'  T7 pair settled: {tr["cov"]} targets = {depths} @ {tanks} kg')
                    else:
                        say('  T7: no level settled as a pair')
                self.rec['routes'][str(k)] = rec
                if took is None:
                    self.save_state(); break
                name = f'{k:02d}'
                fl.routes[name] = dict(st=took['res']['st'], tank=float(took['res']['tank']))
                fl.save(self.out, note=f's13_pass after route {k}')
                self.add_col(relp(self.out / f'route_{name}.npz'), k, took['member'], took['n'],
                             took['res']['tank'], took['pt'], took['res']['miss'],
                             [int(x) for x in took['res']['st']['asts']], True)
                sj = sum(RI.cost(r['tank']) for r in fl.routes.values())
                say(f'route {k:02d}: TOOK {took["member"]} {took["n"]} fb @ {took["res"]["tank"]:.0f} kg; covered '
                    f'{len(fl.coverage())}, sum J_i {sj:.4f} ({rec["sec"]} s)')
                self.save_state()
                k += 1
            self.finish_t7(fl)
        finally:
            if self.pool is not None:
                self.pool.close(); self.pool.join()
        sj = sum(RI.cost(r['tank']) for r in fl.routes.values())
        names = sorted(fl.routes)
        summ = dict(covered=len(fl.coverage()), sumJi=round(sj, 4), J=round(sj + 300 - len(fl.coverage()), 4),
                    depths=[len(fl.routes[n]['st']['asts']) for n in names], tanks=[round(fl.routes[n]['tank']) for n in names],
                    n_routes=len(names), wall_s=round(time.time() - self.t0))
        self.rec.update(done=True, summary=summ, finished=time.strftime('%Y-%m-%d %H:%M:%S'))
        self.save_state()
        say('PASS ' + json.dumps(summ))
        return self.rec

    def finish_t7(self, fl):
        a = self.a; t7 = self.rec.get('t7')
        if not t7 or not t7.get('settled') or t7.get('decided'):
            return
        kA, kB = f'{a.n - 1:02d}', f'{a.n:02d}'
        seq = [fl.routes[x] for x in (kA, kB) if x in fl.routes]
        sn = sum(len(r['st']['asts']) for r in seq); sv = sum(value(len(r['st']['asts']), r['tank'], a.kappa) for r in seq)
        s = t7['settled']
        tn = sum(s['n']); tv = sum(value(n, t, a.kappa) for n, t in zip(s['n'], s['tank']))
        take = (tn, tv) > (sn, sv)
        t7.update(decided=True, taken=take, seq=dict(n=sn, value=round(sv, 3)), pair=dict(n=tn, value=round(tv, 3)))
        self.say(f'T7 vs sequential routes {kA}+{kB}: joint {tn} (value {tv:.2f}) vs sequential {sn} (value {sv:.2f}) -> '
                 f'{"TAKE T7" if take else "keep sequential"}')
        if take:
            for x in (kA, kB):
                if x in fl.routes:           # keep the sequential routes as library columns
                    r = fl.routes[x]; n = len(r['st']['asts'])
                    p = self.out / 'lib' / f'route_{x}_seq_{n}.npz'
                    np.savez(p, **r['st'])
                    self.add_col(relp(p), int(x), 'seq', n, r['tank'], float('nan'), float('nan'),
                                 [int(q) for q in r['st']['asts']], False)
            for x, f in zip((kA, kB), s['files']):
                z = np.load(ROOT / f); st = {q: z[q] for q in z.files}; st['tL'] = float(st['tL'])
                fl.routes[x] = dict(st=st, tank=float(RI.ipr(st).tank()))
            fl.save(self.out, note='s13_pass: T7 pair replaces the last two routes')
        self.save_state()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out'); ap.add_argument('--start', default=''); ap.add_argument('--keep', default='')
    ap.add_argument('--prem', default=''); ap.add_argument('--delta', type=float, default=0.01)
    ap.add_argument('--dvt', type=float, default=1.8); ap.add_argument('--beam', type=int, default=100)
    ap.add_argument('--npt', type=int, default=6); ap.add_argument('--mc', type=int, default=150)
    ap.add_argument('--tail-from', type=int, default=6); ap.add_argument('--members', default='T1,T2,T3,T5,T6,T7')
    ap.add_argument('--kappa', type=float, default=90.0); ap.add_argument('--cap-head', type=float, default=1050.0)
    ap.add_argument('--cap-tail', type=float, default=1080.0); ap.add_argument('--grid', default='0,800,20')
    ap.add_argument('--pilot', type=int, default=0); ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--jitter', type=float, default=0.0); ap.add_argument('--nproc', type=int, default=4)
    ap.add_argument('--no-repair', dest='repair', action='store_false'); ap.add_argument('--n', type=int, default=8)
    ap.add_argument('--cache', default='results/s13/gen/cache')
    ap.add_argument('--choice', default='lex', choices=['lex', 'value'])
    ap.add_argument('--settle-timeout', type=float, default=SB.SETTLE_TIMEOUT,
                    help='SIGALRM per settle [s]; the plan value is 240 (G0 reproduction of the contended tail probe: 40)')
    a = ap.parse_args()
    if a.nproc > 4:
        raise SystemExit('CPU budget: --nproc must be <= 4')
    RI.eph()
    Pass(a).run()


if __name__ == '__main__':
    main()
