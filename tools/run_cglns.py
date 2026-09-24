"""Large-neighbourhood search with column-generation repair (CG-LNS).

Destroy: remove a set D of k routes from the current planned fleet. Residual rows R = targets that only routes in D cover.
Repair: small set cover over R with at most k (or k-1) new routes, chosen from the column store (all route pools +
generated columns + every route ever seen), solved exactly by HiGHS (rows = |R|, columns = store columns touching R);
then `--price` rounds: LP duals on R as prizes (waypoints: every other target visitable at prize 0.001), one or two
pricing beams (launch 0 - launch_max), new columns into the store, set cover again. Accept if the fleet J decreases.
Fleet J = sum_i c_i + misses + 2 with the colgen cost model (c_i = cost_sc(620 / (1 - 1.15 F_i / m0_search))).

Destroy sets (per round, in order): every pair among the `--small` weakest routes (fewest unique targets), every
(weak route, big route) pair whose flyby times interleave most, and `--random` random triples.
Usage: run_cglns.py fleetdir outdir --pools ... --columns results/colgen/run*/columns.jsonl [--rounds 3] [--nproc 4]
"""
import sys, json, time, pathlib, argparse, glob, itertools
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds, linprog
from scipy.sparse import csr_matrix, hstack, vstack, identity
from ctoc14.kepler import Ephemeris
from ctoc14.constants import DAY, M0_MAX
from ctoc14.colgen import CostModel, ColumnStore, PricingSpec, price, TARGETS
from ctoc14.search import UNREACHABLE

ap = argparse.ArgumentParser()
ap.add_argument('fleetdir'); ap.add_argument('outdir'); ap.add_argument('--pools', nargs='*', default=[]); ap.add_argument('--columns', nargs='*', default=[])
ap.add_argument('--rounds', type=int, default=3); ap.add_argument('--nproc', type=int, default=4); ap.add_argument('--beam', type=int, default=250)
ap.add_argument('--price', type=int, default=2, help='pricing rounds per destroy set'); ap.add_argument('--beams', type=int, default=1)
ap.add_argument('--small', type=int, default=5); ap.add_argument('--big-pairs', type=int, default=6); ap.add_argument('--random', type=int, default=4)
ap.add_argument('--launch-max', type=float, default=2000.0); ap.add_argument('--mip-time', type=float, default=60.0)
ap.add_argument('--seed', type=int, default=0)
a = ap.parse_args()
out = pathlib.Path(a.outdir); out.mkdir(parents=True, exist_ok=True); logf = open(out / 'cglns.log', 'a')
def say(s):
    line = f'[{time.strftime("%H:%M:%S")}] {s}'; print(line, flush=True); logf.write(line + '\n'); logf.flush()
rng = np.random.default_rng(a.seed)
cm = CostModel(); store = ColumnStore(cm); tic = time.time()
for pat in a.pools:
    for f in sorted(glob.glob(pat)): store.load_jsonl(f, tag=pathlib.Path(f).parent.name + '/' + pathlib.Path(f).stem)
for pat in a.columns:
    for f in sorted(glob.glob(pat)): store.load_columns_file(f)
colfile = out / 'columns.jsonl'
if colfile.exists(): store.load_columns_file(colfile)
say(f'store {len(store)} columns ({time.time()-tic:.0f} s)')
eph = Ephemeris(); TSET = set(TARGETS.tolist())

def tset(t): return {int(l['ast']) for l in t['legs']} - UNREACHABLE
def m0s(t): return float(t.get('m0_search', t['m0']))
def rcost(t): return float(cm.cost(t['fuel_est'], m0s(t)))
def fleet_J(tours):
    cov = set().union(*[tset(t) for t in tours]) if tours else set()
    return sum(rcost(t) for t in tours) + len(TSET - cov) + 2, cov
def add_tour(t, tag):
    j, _ = store.add([l['ast'] for l in t['legs']], t['fuel_est'], m0s(t), t['t_launch'], ('inline', dict(t, m0_search=m0s(t))), tag); return j
def save(tours, tag):
    for f in out.glob('tour_sc*.json'): f.unlink()
    for k, t in enumerate(sorted(tours, key=lambda t: -len(t['legs'])), 1):
        t = dict(t); t['m0'] = float(min(M0_MAX, cm.tank(t['fuel_est'], m0s(t)))); t['cost_J'] = rcost(t)
        json.dump(t, open(out / f'tour_sc{k}.json', 'w'), indent=1)
    J, cov = fleet_J(tours)
    json.dump(dict(J=J, n=len(tours), covered=len(cov), sumJ=sum(rcost(t) for t in tours), missed=sorted(TSET - cov), tag=tag),
              open(out / 'summary.json', 'w'), indent=1)

def cover(R, kmax, relax=False):
    """Set cover over rows R (0-based target idx) with <= kmax columns from the store. Returns (value, cols, duals, mu)."""
    C = store.csr(); Rr = np.array(sorted(R), int)
    sub = C[:, Rr]; touch = np.asarray(sub.sum(axis=1)).ravel()
    cand = np.where(touch >= 1)[0]
    if len(cand) > 20000:                       # keep the columns with the best (targets in R) - cost
        score = touch[cand] - store.costs()[cand]; cand = cand[np.argsort(-score)[:20000]]
    S = sub[cand].T.tocsr(); n = len(cand); nR = len(Rr); cost = store.costs()[cand]
    A = vstack([hstack([-S, -identity(nR, format='csr')]), csr_matrix((np.ones(n), (np.zeros(n, int), np.arange(n))), shape=(1, n + nR))]).tocsr()
    b = np.concatenate([-np.ones(nR), [kmax]]); c = np.concatenate([cost, np.ones(nR)])
    if relax:
        r = linprog(c, A_ub=A, b_ub=b, bounds=(0, 1), method='highs')
        pi = np.zeros(300); pi[Rr] = -r.ineqlin.marginals[:nR]; return r.fun, cand[r.x[:n] > 1e-6], pi, float(-r.ineqlin.marginals[nR])
    r = milp(c, constraints=LinearConstraint(A, np.full(nR + 1, -np.inf), b), integrality=np.concatenate([np.ones(n), np.zeros(nR)]),
             bounds=Bounds(0, 1), options=dict(time_limit=a.mip_time, mip_rel_gap=1e-4))
    if r.x is None: return np.inf, [], None, None
    return r.fun, cand[np.round(r.x[:n]) > 0], None, None

def repair(tours, D):
    """Destroy routes D (indices), repair by set cover + pricing. Returns (J_new, new tours list) or (None, None)."""
    keep = [t for i, t in enumerate(tours) if i not in D]
    kcov = set().union(*[tset(t) for t in keep]) if keep else set(); kcost = sum(rcost(t) for t in keep)
    R = {x - 1 for x in (TSET - kcov)}
    if not R: return kcost + 2, keep
    best = None
    for p in range(a.price + 1):
        val, cols, _, _ = cover(R, len(D))
        if cols is not None and len(cols):
            covR = set()
            for j in cols: covR |= set(store.tidx[j].tolist())
            J = kcost + sum(store.cost[j] for j in cols) + len(R - covR) + 2
            if best is None or J < best[0] - 1e-9: best = (J, list(cols))
        if p == a.price: break
        lv, _, pi, mu = cover(R, len(D), relax=True)
        # the residual LP is degenerate (the current routes are LP-optimal, dual weight sits on 1-2 targets), so its duals
        # mislead the pricing: price every residual target at 1 (one miss) instead, waypoints at 0.001
        pi = np.zeros(300); pi[sorted(R)] = 1.0; mu = 0.0
        prize = pi.copy()
        for x in kcov: prize[x - 1] = max(prize[x - 1], 0.001)          # waypoints
        for b in range(a.beams):
            spec = PricingSpec(tag=f'cgl{p}', launch=(0.0, a.launch_max, 20.0), beam=a.beam, nproc=a.nproc, seed=int(rng.integers(1e6)),
                               perturb=0.1 if b else 0.0)
            added, st = price(eph, store, prize, pi, mu, (), spec, log=lambda s: None, rnd='cglns')
            with open(colfile, 'a') as fh:
                for j in added: fh.write(json.dumps(store.src[j][1]) + '\n')
            say(f'      price {p}.{b}: LP {lv + kcost + 2:.3f}, added {len(added)}, rc_min {st["rc_min"]:.3f} ({st["runtime"]:.0f} s)')
    if best is None: return None, None
    new = keep + [dict(store.tour(j), m0_search=float(store.m0s[j]), tag=f'cglns:{store.tag[j]}') for j in best[1]]
    return best[0], new

src = pathlib.Path(a.fleetdir)
tours = [json.load(open(f)) for f in sorted(src.glob('tour_sc*.json'), key=lambda f: int(f.stem[7:]))]
for t in tours: t.setdefault('m0_search', 1000.0); add_tour(t, 'fleet')
J0, cov = fleet_J(tours); save(tours, 'start')
say(f'start {src}: {len(tours)} craft, covered {len(cov)}, planned J {J0:.3f}')
tried = set()
for rnd in range(a.rounds):
    improved = False
    cnt = {}
    for t in tours:
        for x in tset(t): cnt[x] = cnt.get(x, 0) + 1
    uniq = [sum(1 for x in tset(t) if cnt[x] == 1) for t in tours]
    order = sorted(range(len(tours)), key=lambda i: uniq[i])
    small = order[:a.small]; big = order[a.small:]
    sets = [tuple(sorted(p)) for p in itertools.combinations(small, 2)]
    def interleave(i, j):
        ti = [l['t_flyby'] for l in tours[i]['legs']]; tj = [l['t_flyby'] for l in tours[j]['legs']]
        return sum(1 for x in ti if min(abs(x - y) for y in tj) < 60 * DAY)
    pb = sorted(((interleave(i, j), i, j) for i in small[:3] for j in big), reverse=True)[:a.big_pairs]
    sets += [tuple(sorted((i, j))) for _, i, j in pb]
    sets += [tuple(sorted(rng.choice(len(tours), 3, replace=False).tolist())) for _ in range(a.random)]
    for D in sets:
        key = tuple(sorted(json.dumps([l['ast'] for l in tours[i]['legs']]) for i in D))
        if key in tried: continue
        tried.add(key)
        tic = time.time()
        Jn, new = repair(tours, set(D))
        if new is None:
            say(f'round {rnd+1} destroy {D}: no repair'); continue
        tag = 'ACCEPT' if Jn < J0 - 1e-6 else 'reject'
        say(f'round {rnd+1} destroy {D} ({"+".join(str(len(tours[i]["legs"])) for i in D)} flybys, unique {[uniq[i] for i in D]}): '
            f'J {J0:.3f} -> {Jn:.3f}, craft {len(tours)} -> {len(new)} [{tag}] ({time.time()-tic:.0f} s)')
        if Jn < J0 - 1e-6:
            tours = new; J0 = Jn; improved = True; save(tours, f'r{rnd+1}')
            for t in tours: add_tour(t, 'fleet')
            break                                   # indices changed: rebuild the destroy sets
    J0, cov = fleet_J(tours)
    say(f'end of round {rnd+1}: {len(tours)} craft, covered {len(cov)}, planned J {J0:.3f}; missed {sorted(TSET - cov)}')
    if not improved and rnd > 0: break
save(tours, 'final'); say(f'FINAL: {len(tours)} craft, planned J {J0:.3f}')
