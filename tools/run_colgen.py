"""Column generation + dive for a 10-11 craft fleet (ctoc14/colgen.py, docs/design_colgen.md).

Usage:
  run_colgen.py outdir --pools results/phasing/camp*/pool.jsonl results/phasing/pool/*.jsonl
                [--homotopy 12,11,10] [--rounds 3,3,6] [--beams 2] [--beam 300] [--nproc 5] [--no-refine]
                [--cost-s 1.15] [--reserve 20] [--wt 0.9] [--dive] [--dive-beams 1] [--final-mip 300]
Writes outdir/cg.log, outdir/columns.jsonl (every generated column, inline tour + end state), outdir/duals_N*_r*.json,
outdir/mip_N*/ (integer incumbents over the working set), outdir/dive/ (dive fleet), outdir/final/ (best fleet).
Re-running with the same outdir reloads columns.jsonl (resume).
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse, glob, multiprocessing as mp
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from ctoc14.kepler import Ephemeris
from ctoc14.colgen import CostModel, ColumnStore, RMP, PricingSpec, price, fleet_summary, write_fleet, TARGETS
from ctoc14.constants import cost_sc, M0_MAX
from ctoc14.impulsive import from_tour, settle_tour
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import run_ialns as RI

ap = argparse.ArgumentParser()
ap.add_argument('outdir'); ap.add_argument('--pools', nargs='*', default=[])
ap.add_argument('--homotopy', default='12,11,10'); ap.add_argument('--rounds', default='3,3,6')
ap.add_argument('--beams', type=int, default=2, help='pricing beams per round')
ap.add_argument('--beam', type=int, default=300); ap.add_argument('--nproc', type=int, default=5)
ap.add_argument('--no-refine', action='store_true'); ap.add_argument('--cost-s', type=float, default=1.15); ap.add_argument('--reserve', type=float, default=20.0)
ap.add_argument('--wt', type=float, default=0.9); ap.add_argument('--wfuel', type=float, default=None)
ap.add_argument('--dive', action='store_true'); ap.add_argument('--dive-beams', type=int, default=1)
ap.add_argument('--final-mip', type=float, default=300.0); ap.add_argument('--mip-time', type=float, default=120.0)
ap.add_argument('--alpha', type=float, default=0.5, help='dual smoothing weight on the previous prize vector')
ap.add_argument('--min-targets', type=int, default=5)
ap.add_argument('--m0s', default='1000', help='search masses rotated over the pricing beams, e.g. 1000,1100')
ap.add_argument('--dive-launch-max', type=float, default=720.0)
ap.add_argument('--hard-cap', type=float, default=0.0, help='leg cap (km/s) for legs into targets priced >= 0.5 (0 = off)')
ap.add_argument('--dvmax', type=float, default=1.2)
ap.add_argument('--dive-exclude', action='store_true', help='old behaviour: covered targets excluded in the dive pricing (no waypoints)')
ap.add_argument('--skip-root', action='store_true', help='no root pricing rounds (dive straight from the loaded columns)')
ap.add_argument('--craft-cost', type=float, default=0.0, help='fixed cost per craft in the master (use with a loose --homotopy cap): pricing threshold mu + craft cost')
ap.add_argument('--prize-scale', type=float, default=1.0, help='beam prizes = scale x LP duals (smaller: the beam keeps building long routes)')
ap.add_argument('--twin-recost', type=int, default=0, help='workers: recost the LP support with the impulsive twin after every master solve (0 = off)')
ap.add_argument('--twin-files', nargs='*', default=[], help='recost_pool.py outputs to preload (matched by target set)')
a = ap.parse_args()
out = pathlib.Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
logf = open(out / 'cg.log', 'a')
def say(s):
    line = f'[{time.strftime("%H:%M:%S")}] {s}'; print(line, flush=True); logf.write(line + '\n'); logf.flush()

tic0 = time.time()
cm = CostModel(s=a.cost_s, reserve=a.reserve)
store = ColumnStore(cm)
for pat in a.pools:
    for f in sorted(glob.glob(pat)):
        n = store.load_jsonl(f, tag=pathlib.Path(f).parent.name + '/' + pathlib.Path(f).stem, min_targets=a.min_targets)
        say(f'pool {f}: +{n} columns')
colfile = out / 'columns.jsonl'
if colfile.exists():
    say(f'resume: +{store.load_columns_file(colfile)} generated columns from {colfile}')
say(f'store: {len(store)} columns ({time.time()-tic0:.0f} s); cost model s={cm.s} reserve={cm.reserve}; w_fuel(360 kg @1000) = {cm.w_fuel():.2f}')

eph = Ephemeris()
rmp = RMP(store, craft_cost=a.craft_cost)
all_rows = TARGETS - 1

# ----------------------------------------------------------------------------------------------- twin recost hook
twin_done = {}                     # column id -> twin tank (None = failed)
twin_dir = out / 'twin_states'; twin_dir.mkdir(exist_ok=True)
twin_log = open(out / 'twin.jsonl', 'a')
def _key_of(targets):
    mk = np.zeros(304, np.uint8); mk[np.array(sorted(set(int(t) for t in targets))) - 1] = 1
    return np.packbits(mk).tobytes()
def _apply_twin(j, tank, st=None):
    twin_done[j] = tank
    if tank is None:
        store.cost[j] = float(store.cost[j]) * 1.2          # failed restoration: distrust the planner estimate
    else:
        store.cost[j] = float(cost_sc(tank)) if tank <= M0_MAX else 10.0
        if st is not None:
            np.savez(twin_dir / f'col_{j}.npz', **st)
for f in a.twin_files:
    k0 = 0
    for line in open(f):
        r = json.loads(line)
        j = store.key.get(_key_of(r['targets']))
        if j is None or j in twin_done: continue
        ok = r.get('ok', r.get('tank') is not None)          # recost_pool.py output or a previous run's twin.jsonl
        st = {k: np.asarray(v) for k, v in r['st'].items()} if ok and 'st' in r else None
        if ok and st is None and (twin_dir / f'col_{r["col"]}.npz').exists() and f.startswith(str(out)):
            z = np.load(twin_dir / f'col_{r["col"]}.npz'); st = {k: z[k] for k in z.files}
        _apply_twin(j, r['tank'] if ok else None, st); k0 += 1
    say(f'twin preload {f}: {k0} columns')
ap_twin_timeout = 150.0        # s per column (a pathological restoration can otherwise stall a whole round)
def _twin_alarm(signum, frame):
    raise TimeoutError('twin recost time limit')
def _twin_work(job):
    j, tour = job
    import signal
    signal.signal(signal.SIGALRM, _twin_alarm); signal.alarm(int(ap_twin_timeout))
    try:
        ip, miss, lag = settle_tour(eph, tour, lambda ip: RI.settle(ip, 100))
        if ip is None: return j, None, None, f'miss {miss:.0f} km'
        signal.alarm(0)
        return j, float(ip.tank()), RI.ist(ip), f'ok lag {lag / 86400:.2f} d'
    except Exception as e:  # noqa
        signal.alarm(0)
        return j, None, None, repr(e)[:120]
def twin_recost(cols):
    """Recost the given columns with the twin (those not done yet). Returns the number recosted."""
    for j, t in twin_done.items():                      # price() may have reset a twin cost to a cheaper planner estimate
        if t is not None: store.cost[j] = float(cost_sc(t)) if t <= M0_MAX else 10.0
    todo = [int(j) for j in cols if int(j) not in twin_done]
    if not todo or not a.twin_recost: return 0
    tic = time.time(); nf = 0
    with mp.get_context('fork').Pool(a.twin_recost) as pool:
        for j, tank, st, why in pool.imap_unordered(_twin_work, [(j, store.tour(j)) for j in todo]):
            c0 = float(store.cost[j]); _apply_twin(j, tank, st); nf += tank is None
            twin_log.write(json.dumps(dict(col=j, n=int(store.n[j]), tank=tank, cost_planner=c0, cost=float(store.cost[j]), why=why,
                                           targets=store.targets(j))) + '\n'); twin_log.flush()
    say(f'  twin recost: {len(todo)} columns ({nf} failed) in {time.time()-tic:.0f} s; '
        f'ratio twin/planner median {np.median([store.cost[j] / max(1e-9, c) for j, c in [(j, store.cm.cost(store.fuel[j], store.m0s[j])) for j in todo if twin_done[j] is not None]] or [np.nan]):.3f}')
    return len(todo)
def solve_twin(rows, N, fixed_cost=0.0, verbose=None):
    res = rmp.solve(rows, N, fixed_cost=fixed_cost, verbose=verbose)
    if twin_recost(res.cols[res.y > 0.02]):                     # one pass: support with real weight only
        res = rmp.solve(rows, N, fixed_cost=fixed_cost)
    return res

def spec_for(k, tag_round, excluded_n=0):
    """Rotating diversification of the pricing beams."""
    m0list = [float(x) for x in a.m0s.split(',')]
    base = dict(beam=a.beam, nproc=a.nproc, tof_refine=not a.no_refine, w_t=a.wt, w_fuel=a.wfuel, seed=k,
                m0=m0list[(k // 4) % len(m0list)], hard_cap=a.hard_cap, dv_max=a.dvmax)
    kinds = [dict(tag='A', launch=(0.0, 320.0, 20.0)),
             dict(tag='B', launch=(300.0, 720.0, 20.0)),
             dict(tag='C', launch=(0.0, 620.0, 20.0), perturb=0.15),
             dict(tag='D', launch=(0.0, 620.0, 20.0), first_top=20)]
    kd = kinds[k % len(kinds)]; kd = dict(kd); kd['tag'] = f'{kd["tag"]}{tag_round}'
    return PricingSpec(**base, **kd)

def append_columns(ids):
    with open(colfile, 'a') as f:
        for j in ids:
            src = store.src[j]
            if src[0] == 'inline': f.write(json.dumps(src[1]) + '\n')

# ----------------------------------------------------------------------------------------------- root column generation
Ns = [int(x) for x in a.homotopy.split(',')]; Rs = [int(x) for x in a.rounds.split(',')]
k_spec = 0; history = []
for N, R in zip(Ns, Rs):
    if a.skip_root:
        break
    say(f'=== root CG, N = {N}, up to {R} rounds')
    prize = None; dry = 0
    for r in range(R + 1):
        res = solve_twin(all_rows, N, verbose=say if r == 0 else None)
        nfrac = int(np.sum((res.y > 1e-6) & (res.y < 1 - 1e-6)))
        say(f'  N={N} round {r}: LP {res.value:.4f} (craft cost {a.craft_cost}), sum y {res.y.sum():.2f}, mu {res.mu:.3f}, #pi>0.99 {(res.pi > 0.99).sum()}, '
            f'pi median {np.median(res.pi[all_rows]):.3f}, fractional routes {nfrac}, working set {len(rmp.active)}, store {len(store)}')
        json.dump(dict(N=N, round=r, lp=res.value, mu=res.mu, pi={int(i) + 1: float(res.pi[i]) for i in all_rows}),
                  open(out / f'duals_N{N}_r{r}.json', 'w'))
        history.append(dict(N=N, round=r, lp=res.value, mu=res.mu))
        if r == R or dry >= 2:
            break
        prize = res.pi.copy() if (prize is None or dry > 0) else a.alpha * prize + (1 - a.alpha) * res.pi
        added_round = 0; rc_mins = []
        for b in range(a.beams):
            spec = spec_for(k_spec, f'{N}.{r}'); k_spec += 1
            added, st = price(eph, store, prize * a.prize_scale, res.pi, res.mu + a.craft_cost, (), spec, log=say, rnd=f'N{N}r{r}')
            append_columns(added); added_round += len(added); rc_mins.append(st['rc_min'])
            rmp.active = np.unique(np.concatenate([rmp.active, np.array(added, int)])) if added else rmp.active
        lb = res.value + N * min(0.0, min(rc_mins))
        say(f'  N={N} round {r}: added {added_round} columns, best rc {min(rc_mins):.3f}, pseudo-bound {lb:.3f}')
        dry = dry + 1 if added_round == 0 else 0
    Jm, sel, missed = rmp.mip(all_rows, N, time_limit=a.mip_time)
    if twin_recost(sel): Jm, sel, missed = rmp.mip(all_rows, N, time_limit=a.mip_time)
    summ = write_fleet(store, sel, out / f'mip_N{N}', extra=dict(lp=res.value, mip_J=Jm))
    say(f'  MIP N={N}: J {Jm:.3f}, {summ["n"]} craft, covered {summ["covered"]}, sum J_i {summ["sumJ"]:.3f}, missed {summ["missed"]}')

# ----------------------------------------------------------------------------------------------- dive
best = None
Nf = Ns[-1]
if a.dive:
    say(f'=== dive, N = {Nf}')
    rows = all_rows.copy(); fixed = []; fixed_cost = 0.0; covered = set()
    for level in range(Nf):
        N_left = Nf - level
        res = solve_twin(rows, N_left, fixed_cost=fixed_cost)
        say(f'  level {level}: rows {len(rows)}, N_left {N_left}, LP {res.value:.4f}, mu {res.mu:.3f}, #pi>0.99 {(res.pi > 0.99).sum()}')
        # re-price the residual problem before fixing (level > 0: the residual changed)
        if level > 0:
            for b in range(a.dive_beams):
                # plain full-window beam first, then perturbed-prize / hard-first variants (A and B would coincide once
                # the launch window is overridden)
                spec = spec_for([0, 2, 3][b % 3] + 4 * (k_spec // 4), f'dv{level}'); k_spec += 1
                spec.launch = (0.0, a.dive_launch_max, 20.0); spec.seed = k_spec
                # covered targets stay visitable (prize 0 = their dual outside the rows): they are the waypoints a residual
                # route needs between sparse uncovered targets (legs must end at an asteroid)
                added, st = price(eph, store, res.pi * a.prize_scale, res.pi, res.mu + a.craft_cost, sorted(covered) if a.dive_exclude else (), spec, log=say, rnd=f'dive{level}')
                append_columns(added)
                if added: rmp.active = np.unique(np.concatenate([rmp.active, np.array(added, int)]))
            res = solve_twin(rows, N_left, fixed_cost=fixed_cost)
            say(f'  level {level}: after re-pricing LP {res.value:.4f}')
        # candidates: top-3 by y, one-step look-ahead
        order = np.argsort(-res.y)[:3]
        cands = [(int(res.cols[i]), float(res.y[i])) for i in order if res.y[i] > 1e-6]
        if not cands:
            say('  no fractional route left; stop'); break
        scored = []
        for j, yj in cands:
            tj = np.array(store.targets(j)) - 1
            rows_j = np.setdiff1d(rows, tj)
            v = rmp.lp_value_fixed(rows_j, N_left - 1, fixed_cost + store.cost[j])
            scored.append((v, j, yj, len(np.intersect1d(rows, tj))))
        scored.sort()
        v, j, yj, gain = scored[0]
        fixed.append(j); fixed_cost += store.cost[j]
        tj = np.array(store.targets(j)) - 1; rows = np.setdiff1d(rows, tj); covered |= set(store.targets(j))
        say(f'  level {level}: FIX col {j} [{store.tag[j]}] y={yj:.2f}, {store.n[j]} targets ({gain} new), cost {store.cost[j]:.3f}, '
            f'look-ahead LP {v:.4f} (alternatives {[(round(s[0], 3), s[1]) for s in scored[1:]]}); covered {len(covered)}')
        if len(rows) == 0: break
    summ = write_fleet(store, fixed, out / 'dive', extra=dict(history=history))
    say(f'DIVE: {summ["n"]} craft, covered {summ["covered"]}, sum J_i {summ["sumJ"]:.3f}, planned J {summ["J"]:.3f}, missed {summ["missed"]}')
    best = (summ['J'], fixed, 'dive')

# ----------------------------------------------------------------------------------------------- final MIP over the working set
if a.final_mip > 0:
    res = solve_twin(all_rows, Nf)
    Jm, sel, missed = rmp.mip(all_rows, Nf, time_limit=a.final_mip)
    twin_recost(sel); Jm2, sel, missed = rmp.mip(all_rows, Nf, time_limit=a.final_mip); say(f'  final MIP after twin recost: {Jm:.3f} -> {Jm2:.3f}'); Jm = Jm2
    summ = write_fleet(store, sel, out / f'mip_final_N{Nf}', extra=dict(lp=res.value))
    say(f'FINAL MIP N={Nf}: LP {res.value:.3f}, J {Jm:.3f}, covered {summ["covered"]}, sum J_i {summ["sumJ"]:.3f}, missed {summ["missed"]}')
    if best is None or Jm < best[0]:
        best = (Jm, sel, 'mip')
if best:
    summ = write_fleet(store, best[1], out / 'final', extra=dict(source=best[2]))
    say(f'BEST ({best[2]}): planned J {summ["J"]:.3f}, {summ["n"]} craft, covered {summ["covered"]}; total runtime {time.time()-tic0:.0f} s')
