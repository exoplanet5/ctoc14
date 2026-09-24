"""Fleet search on the impulsive whole-trajectory model (ctoc14/impulsive.py; ~10 s per insertion trial instead of ~7 min).

Fleet = directory of impulsive route states route_<name>.npz; J = sum cost(tank) + misses, tank = 601.5 exp(dv / ve).
Commands:
  import <exact fleetdir> <ifleetdir>      impulsive twins of run_gfleet craft states, re-optimised
  search <ifleetdir> <outdir>              rounds of: dissolve attempts (every route, weakest first) and relocate passes
                                           (move a target from route A to route B when saving(A) > insertion cost(B))
  export <ifleetdir> <exact fleetdir>      continuous-thrust conversion (to_exact + mass-free SCP + tank sizing) of
                                           every route; the result is a run_gfleet fleet (craft_<name>.npz), then
                                           `run_gfleet.py write` validates and combines
Insertion candidates: local minima of the distance between the route trajectory (Kepler arcs sampled every 2 d) and the
target asteroid below --dmax AU, nearest first; each trial = aim-point homotopy + SCP (globalopt.insert_homotopy)."""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse, multiprocessing as mp
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from ctoc14.kepler import Ephemeris, propagate_batch
from ctoc14.globalopt import (optimise, restore, insert_homotopy, insert_block, remove_flybys, problem_from_state,
                              load_state, state_of, to_thrust, tighten)
from ctoc14.impulsive import ImpulsiveProblem, from_thrust_problem, to_exact, enforce_cap
from ctoc14.constants import DAY, AU, T_MISSION, M_DRY

QUIET = lambda s: None
_E = None


def eph():
    global _E
    if _E is None:
        _E = Ephemeris()
    return _E


def cost(tank):
    x = (tank - M_DRY) / 1400.0
    return 1.0 + x + x * x


def ist(ip):
    return dict(tL=float(ip.tL), vinf=np.array(ip.vinf, float), ts=np.array(ip.ts, float), Ts=np.array(ip.Ts, float),
                tf=np.array(ip.tf, float), asts=np.array(ip.asts, int))


def ipr(st):
    return ImpulsiveProblem(eph(), st['tL'], st['vinf'], st['ts'], st['Ts'], st['tf'], [int(a) for a in st['asts']])


def settle(ip, iters=100):
    restore(ip, tol_km=100, iters=30, rho=1e3, log=QUIET)
    optimise(ip, iters=iters, rho0=100.0, log=QUIET)
    d = restore(ip, tol_km=100, iters=30, log=QUIET)
    if len(ip.Ts) and np.linalg.norm(ip.Ts, axis=1).max() > 1.02 * ip.tcap:
        # continuous-thrust capability: impulses above 0.43 N x bin / tank are not flyable (heavy routes)
        enforce_cap(ip, stages=6)
        optimise(ip, iters=40, rho0=100.0, log=QUIET)
        d = restore(ip, tol_km=100, iters=30, log=QUIET)
    return float(d.max()) if len(d) else 0.0


class IFleet:
    def __init__(self, d=None):
        self.routes = {}
        if d is not None:
            for f in sorted(pathlib.Path(d).glob('route_*.npz')):
                z = np.load(f); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL'])
                self.routes[f.stem[6:]] = dict(st=st, tank=float(ipr(st).tank()))

    def coverage(self):
        cov = {}
        for n, r in self.routes.items():
            for a in r['st']['asts']:
                cov.setdefault(int(a), set()).add(n)
        return cov

    def J(self):
        return sum(cost(r['tank']) for r in self.routes.values()) + 300 - len(self.coverage())

    def save(self, d, note=''):
        d = pathlib.Path(d); d.mkdir(parents=True, exist_ok=True)
        for f in d.glob('route_*.npz'):
            f.unlink()
        for n, r in self.routes.items():
            np.savez(d / f'route_{n}.npz', **r['st'])
        cov = self.coverage()
        json.dump(dict(J=self.J(), n=len(self.routes), covered=len(cov), note=note,
                       routes={n: dict(tank=r['tank'], J_i=cost(r['tank']), flybys=len(r['st']['asts'])) for n, r in self.routes.items()}),
                  open(d / 'fleet.json', 'w'), indent=1)

    def summary(self):
        return (f'{len(self.routes)} routes, covered {len(self.coverage())}, J {self.J():.4f}; ' +
                ' '.join(f'{n}:{len(r["st"]["asts"])}@{r["tank"]:.0f}' for n, r in sorted(self.routes.items(), key=lambda kv: -len(kv[1]['st']['asts']))))


def logger(path):
    f = open(path, 'a')
    def say(s):
        line = f'[{time.strftime("%H:%M:%S")}] {s}'; print(line, flush=True); f.write(line + '\n'); f.flush()
    return say


# ------------------------------------------------------------------ workers
def w_import(job):
    name, path = job
    gp = problem_from_state(eph(), load_state(path))
    ip = from_thrust_problem(gp)
    miss = settle(ip, 150)
    return name, ist(ip), float(ip.tank()), miss, float(gp.m0)


def w_remove(job):
    name, st, asts = job
    ip = ipr(st); remove_flybys(ip, asts)
    try:
        miss = settle(ip, 100)
    except Exception as e:
        return dict(name=name, asts=asts, ok=False, error=repr(e))
    return dict(name=name, asts=asts, ok=miss <= 150, st=ist(ip), tank=float(ip.tank()), miss=miss)


def w_insert(job):
    name, st, ast, t = job
    ip = ipr(st)
    try:
        d = insert_homotopy(ip, ast, t, stages=10, iters_stage=6, tol_km=100, log=QUIET)
        if d.max() > 1e4:
            return dict(name=name, ast=ast, t=t, ok=False, miss=float(d.max()))
        miss = settle(ip, 100)
    except Exception as e:
        return dict(name=name, ast=ast, t=t, ok=False, miss=np.inf, error=repr(e))
    return dict(name=name, ast=ast, t=t, ok=miss <= 150, st=ist(ip), tank=float(ip.tank()), miss=miss)


def trajectory(ip, step=2 * DAY):
    """Positions along the route every `step` (Kepler arcs from the post-impulse node states)."""
    Yf, Ys = ip.integrate()
    tn, kind = ip._nodes(ip.tf)
    X = np.zeros((len(tn), 6)); X[0, :3] = ip.rE; X[0, 3:] = ip.vE + ip.vinf
    for n in range(1, len(tn)):
        k = kind[n]
        X[n] = np.concatenate([Ys[k, :3], Ys[k, 3:6] + ip.Ts[k]]) if k >= 0 else Yf[-2 - k, :6]
    t_end = min(T_MISSION - DAY, tn[-1] + 400 * DAY)
    tt = np.arange(ip.tL + step, t_end, step)
    seg = np.clip(np.searchsorted(tn, tt, side='right') - 1, 0, len(tn) - 1)
    r, _ = propagate_batch(X[seg, :3], X[seg, 3:], tt - tn[seg])
    return tt, r


def w_cands(job):
    name, st, targets, dmax = job
    E = eph(); ip = ipr(st); tt, r = trajectory(ip)
    out = []
    for X in targets:
        if X in set(int(a) for a in st['asts']):
            continue
        ra, _ = E.ast_states_at(np.full(len(tt), X - 1), tt)
        d = np.linalg.norm(r - ra, axis=1) / AU
        i = np.where((d[1:-1] <= d[:-2]) & (d[1:-1] < d[2:]) & (d[1:-1] < dmax))[0] + 1
        for k in i:
            a, b, c = d[k - 1], d[k], d[k + 1]; den = a - 2 * b + c
            off = 0.5 * (a - c) / den if den > 0 else 0.0
            out.append(dict(host=name, ast=int(X), t=float(tt[k] + off * (tt[1] - tt[0])), dist=float(b)))
    return out


def w_export(job):
    """Impulsive route -> continuous thrust (mass-free SCP, tank sizing), then polish + validation with finer RK4 steps
    (0.1, 0.05, 0.025 d) until the validator accepts every flyby. Writes craft_<name>.npz and frags/frag_<name>.txt."""
    name, st, out = job
    from ctoc14.globalopt import write_craft, save_state
    out = pathlib.Path(out); (out / 'frags').mkdir(parents=True, exist_ok=True)
    ip = ipr(st); nfb = len(ip.asts); gm = to_exact(ip, eph())
    restore(gm, tol_km=300, iters=40, rho=1e3, log=QUIET)
    optimise(gm, iters=40, rho0=100.0, log=QUIET); restore(gm, tol_km=150, iters=20, log=QUIET)
    gt = to_thrust(gm, 1.5); restore(gt, tol_km=150, iters=20, log=QUIET); tighten(gt, 1.5, log=QUIET)
    ok = False; hs_used = None
    for hs in (0.1, 0.05, 0.025):
        gt = problem_from_state(eph(), state_of(gt), hstep=hs * DAY)
        restore(gt, tol_km=100, iters=30, log=QUIET); d = tighten(gt, 1.5, log=QUIET)
        frag = out / 'frags' / f'frag_{name}.txt'
        rep, fuel, pk = write_craft(gt, str(frag), eph=eph())
        if rep.ok and len(rep.flybys) >= nfb:
            ok = True; hs_used = hs; break
    save_state(state_of(gt), out / f'craft_{name}.npz')
    if not ok:
        frag.rename(out / 'frags' / f'INVALID_frag_{name}.txt')
    return name, float(gt.m0), float(d.max()), float(ip.tank()), ok, hs_used



TANK_MAX = 1900.0          # kg: keep margin below the 2000 kg launch-mass limit for the exact model


def w_block(job):
    """Insert a block of (ast, t) flybys into a host at once, after ejecting the host's own flybys inside the block
    window (+- margin). Returns the new route state, the ejected targets and the exact cost delta.
    Score of the move for the fleet objective (sum cost + misses): dJ + (#ejected - #inserted)."""
    name, st, items, margin, tag = job
    ip = ipr(st); tank0 = ip.tank(); tic = time.time()
    have = {int(a) for a in st['asts']}
    items = [(int(a), float(t)) for a, t in items if int(a) not in have]
    if not items:
        return dict(tag=tag, name=name, ok=False, why='duplicate')
    t0 = min(t for _, t in items) - margin; t1 = max(t for _, t in items) + margin
    ej = [int(a) for a, t in zip(ip.asts, ip.tf) if t0 <= t <= t1]
    if ej:
        remove_flybys(ip, ej)
    try:
        d = insert_block(ip, items, tol_km=100, ds0=0.2, iters_stage=6, log=QUIET)
        if d.max() > 1e4:
            return dict(tag=tag, name=name, ok=False, why=f'homotopy {d.max():.0f} km', sec=time.time() - tic)
        miss = settle(ip, 100)
    except Exception as e:
        return dict(tag=tag, name=name, ok=False, why=repr(e)[:60], sec=time.time() - tic)
    tank = float(ip.tank())
    if miss > 150 or tank > TANK_MAX:
        return dict(tag=tag, name=name, ok=False, why=f'miss {miss:.0f} km, tank {tank:.0f}', sec=time.time() - tic)
    dJ = cost(tank) - cost(tank0)
    return dict(tag=tag, name=name, ok=True, st=ist(ip), tank=tank, dJ=dJ, ej=ej, ins=[a for a, _ in items],
                net=len(items) - len(ej), sec=time.time() - tic)


def make_blocks(items, bmax, gap=400 * DAY):
    """Split (ast, t) pairs into time-contiguous blocks of at most bmax targets (a gap > `gap` also splits)."""
    items = sorted(items, key=lambda it: it[1]); out = []; cur = []
    for it in items:
        if cur and (len(cur) >= bmax or it[1] - cur[-1][1] > gap):
            out.append(cur); cur = []
        cur.append(it)
    if cur:
        out.append(cur)
    return out


def host_rank(fleet, items, exclude=(), E=None):
    """Hosts ranked by the mean distance between their trajectory and the block's asteroids at the block times."""
    E = E or eph(); out = []
    ta = np.array([t for _, t in items]); ia = np.array([a - 1 for a, _ in items])
    ra, _ = E.ast_states_at(ia, ta)
    for n, r in fleet.routes.items():
        if n in exclude:
            continue
        ip = ipr(r['st']); tt, rr = trajectory(ip, step=4 * DAY)
        k = np.clip(np.searchsorted(tt, ta), 0, len(tt) - 1)
        out.append((float(np.mean(np.linalg.norm(rr[k] - ra, axis=1)) / AU), n))
    return sorted(out)


def dissolve_lns(fleet, V, pool, a, say):
    """Ruin and recreate: remove route V, then place the unassigned targets in time blocks into the other routes,
    ejecting the host's flybys inside each block window (the ejected targets return to the pool). Greedy on the exact
    fleet objective: every move changes (sum cost + misses) by dJ + #ejected - #inserted. Accept if J improves."""
    J0 = fleet.J(); cov = fleet.coverage()
    trial = IFleet(); trial.routes = {n: dict(r) for n, r in fleet.routes.items() if n != V}
    vst = fleet.routes[V]['st']
    unassigned = {int(x): float(t) for x, t in zip(vst['asts'], vst['tf']) if cov[int(x)] == {V}}
    say(f'LNS dissolve {V} (cost {cost(fleet.routes[V]["tank"]):.4f}, {len(unassigned)} targets); target J < {J0:.4f}')
    tabu = {}; missed = []; cV = cost(fleet.routes[V]['tank'])
    lam_base = a.lam if a.lam > 0 else cV / max(len(unassigned), 1); lam = lam_base
    for rnd in range(a.lns_rounds):
        if not unassigned:
            break
        blocks = make_blocks(list(unassigned.items()), a.bmax)
        jobs = []
        for bi, items in enumerate(blocks):
            for _, n in host_rank(trial, items)[:a.hosts]:
                if any(tabu.get((x, n), 0) > rnd for x, _ in items):
                    continue
                for margin in a.margins:
                    jobs.append((n, trial.routes[n]['st'], items, margin * DAY, f'r{rnd}b{bi}m{margin:.0f}'))
            # also try the block split in half (cheaper moves are often smaller)
            if len(items) >= 4:
                for half in (items[:len(items) // 2], items[len(items) // 2:]):
                    for _, n in host_rank(trial, half)[:max(1, a.hosts - 1)]:
                        jobs.append((n, trial.routes[n]['st'], half, a.margins[0] * DAY, f'r{rnd}b{bi}h'))
        if not jobs:
            break
        res = [r for r in pool.map(w_block, jobs) if r.get('ok')]
        # price of one unplaced target during the dissolve: a move is worth it when its cost per NET placed target
        # is below lam. lam starts at the victim's cost per target and is raised only when the search stalls.
        lam = max(lam_base, lam * 0.7)          # relax the ratchet: each round starts near the base price
        good = [r for r in res if r['net'] > 0 and r['dJ'] - lam * r['net'] < 0]
        while not good and lam < a.lam_max:
            lam = min(lam * 1.3, a.lam_max)
            good = [r for r in res if r['net'] > 0 and r['dJ'] - lam * r['net'] < 0]
            say(f'  round {rnd}: no move under the price; lam -> {lam:.3f}')
        if not good:
            cheap = sorted((r['dJ'] / r['net'], r['name'], r['ins']) for r in res if r['net'] > 0)[:3]
            say(f'  round {rnd}: no move below lam {lam:.3f} among {len(jobs)} trials ({len(res)} feasible); '
                f'cheapest ' + ', '.join(f'{c:.3f}/target {n}{i}' for c, n, i in cheap) +
                f'; {len(unassigned)} targets stay unplaced')
            break
        best = min(good, key=lambda r: r['dJ'] / r['net'])
        pst = trial.routes[best['name']]['st']
        prev_tf = {int(x): float(t) for x, t in zip(pst['asts'], pst['tf'])}
        trial.routes[best['name']] = dict(st=best['st'], tank=best['tank'])
        for x in best['ins']:
            unassigned.pop(x, None)
        tcov = trial.coverage()
        for x in best['ej']:
            if x in tcov:                       # still covered by another route: nothing to re-place
                continue
            unassigned[x] = float(prev_tf[x])
            tabu[(x, best['name'])] = rnd + a.tabu
        spent = sum(cost(r['tank']) for r in trial.routes.values()) - (J0 - cost(fleet.routes[V]['tank']) - (300 - len(cov)))
        say(f'  round {rnd}: {best["name"]} += {best["ins"]} -= {best["ej"]} dJ {best["dJ"]:+.4f} '
            f'({best["dJ"] / best["net"]:.3f}/target, lam {lam:.3f}); unplaced {len(unassigned)}, '
            f'spent {spent:.3f}/{cV + a.anneal + a.slack:.3f}, J {trial.J():.4f}')
        if spent > cV + a.anneal + a.slack:
            say(f'  abort: placements {spent:.3f} > victim cost {cV:.3f} + anneal {a.anneal} + slack {a.slack}')
            return None
    missed = sorted(unassigned)
    say(f'LNS dissolve {V}: J {J0:.4f} -> {trial.J():.4f} (misses included) with {len(missed)} unplaced {missed}')
    if missed:
        if trial.J() >= J0 - 1e-4:
            return None
        say('  accepting a fleet with misses')
    elif trial.J() >= J0 - 1e-4:
        if not (a.anneal > 0 and trial.J() < J0 + a.anneal):
            return None
        say(f'  annealed: relocating the trial fleet (J {trial.J():.4f}, must end below {J0:.4f})')
        for k in range(a.anneal_passes):
            if not relocate_pass(trial, pool, a, say):
                break
            say(f'  annealed {V} pass {k + 1}: J {trial.J():.4f}')
            if trial.J() < J0 - 1e-4:
                break
        if trial.J() >= J0 - 1e-4:
            return None
    return trial


# ------------------------------------------------------------------ moves
def candidates(fleet, pool, targets, dmax, per_target, exclude=()):
    jobs = [(n, r['st'], targets, dmax) for n, r in fleet.routes.items() if n not in exclude]
    allc = [c for res in pool.map(w_cands, jobs) for c in res]
    by = {}
    for c in sorted(allc, key=lambda c: c['dist']):
        lst = by.setdefault(c['ast'], [])
        if len(lst) < per_target and sum(1 for x in lst if x['host'] == c['host']) < 2:
            lst.append(c)
    return by


def dissolve(fleet, V, pool, a, say):
    J0 = fleet.J(); cV = cost(fleet.routes[V]['tank']); cov = fleet.coverage()
    targets = sorted(int(x) for x in fleet.routes[V]['st']['asts'] if cov[int(x)] == {V})
    trial = IFleet(); trial.routes = {n: dict(r) for n, r in fleet.routes.items() if n != V}
    say(f'dissolve {V} (cost {cV:.4f}, {len(targets)} targets)')
    spent = 0.0; placed = []; missed = []; fails = {}
    remaining = list(targets)
    while remaining:
        dmax = {X: min(0.6, a.dmax * 1.5 ** fails.get(X, 0)) for X in remaining}
        by = {}
        for dm in sorted(set(dmax.values())):
            grp = [X for X in remaining if dmax[X] == dm]
            by.update(candidates(trial, pool, grp, dm, a.per_target))
        jobs = [(c['host'], trial.routes[c['host']]['st'], c['ast'], c['t']) for X in remaining for c in by.get(X, [])]
        res = pool.map(w_insert, jobs) if jobs else []
        best = {}
        for r in res:
            if not r['ok']:
                continue
            r['dJ'] = cost(r['tank']) - cost(trial.routes[r['name']]['tank'])
            if r['ast'] not in best or r['dJ'] < best[r['ast']]['dJ']:
                best[r['ast']] = r
        for X in [x for x in remaining if x not in best]:
            fails[X] = fails.get(X, 0) + 1
            say(f'  {X}: no insertion ({len(by.get(X, []))} candidates within {dmax[X]:.2f} AU), failure {fails[X]}')
            if fails[X] >= 3:
                missed.append(X); remaining.remove(X); say(f'  {X}: missed (+1)')
        used = set()
        for r in sorted(best.values(), key=lambda r: r['dJ']):
            if r['name'] in used:
                continue
            if r['dJ'] >= 1.0:
                say(f'  {r["ast"]}: cheapest insertion dJ {r["dJ"]:+.3f} >= 1 -> missed'); missed.append(r['ast']); remaining.remove(r['ast']); continue
            used.add(r['name']); trial.routes[r['name']] = dict(st=r['st'], tank=r['tank'])
            spent += r['dJ']; placed.append((r['ast'], r['name'], round(r['dJ'], 4))); remaining.remove(r['ast'])
            say(f'  placed {r["ast"]} on {r["name"]} at {r["t"] / DAY:.0f} d: dJ {r["dJ"]:+.4f} (placed cost {spent:.3f}, misses {len(missed)}, of {cV:.3f})')
        if spent + len(missed) > cV + a.slack + (a.anneal if not missed else 0.0):
            say(f'  placed cost {spent:.3f} + misses {len(missed)} > {cV:.3f} + slack: abort'); return None
    say(f'dissolve {V}: J {J0:.4f} -> {trial.J():.4f} {"ACCEPT" if trial.J() < J0 - 1e-4 else "reject"}; misses {missed}')
    if trial.J() < J0 - 1e-4:
        return trial
    if a.anneal > 0 and not missed and trial.J() < J0 + a.anneal:
        say(f'  annealed dissolve {V}: relocating the trial fleet (J {trial.J():.4f}, must end below {J0:.4f})')
        for k in range(a.anneal_passes):
            if not relocate_pass(trial, pool, a, say):
                break
            say(f'  annealed {V} pass {k + 1}: J {trial.J():.4f}')
            if trial.J() < J0 - 1e-4:
                break
        say(f'annealed dissolve {V}: J {J0:.4f} -> {trial.J():.4f} {"ACCEPT" if trial.J() < J0 - 1e-4 else "reject"}')
        return trial if trial.J() < J0 - 1e-4 else None
    return None


def relocate_pass(fleet, pool, a, say):
    """Removal savings of every target, then insertion trials of the most expensive ones into other routes."""
    cov = fleet.coverage()
    jobs = [(n, r['st'], [int(x)]) for n, r in fleet.routes.items() for x in r['st']['asts'] if len(cov[int(x)]) == 1]
    rem = {}
    for r in pool.map(w_remove, jobs):
        if r['ok']:
            rem[(r['name'], r['asts'][0])] = r
    sav = sorted(((cost(fleet.routes[n]['tank']) - cost(r['tank']), n, x) for (n, x), r in rem.items()), reverse=True)
    say(f'relocate: top removal savings ' + ', '.join(f'{x}@{n} {s:.3f}' for s, n, x in sav[:12]))
    top = [(s, n, x) for s, n, x in sav if s > a.min_saving][:a.relocate_top]
    if not top:
        return False
    by = candidates(fleet, pool, sorted({x for _, _, x in top}), a.dmax, a.per_target)
    jobs = [(c['host'], fleet.routes[c['host']]['st'], x, c['t']) for s, n, x in top for c in by.get(x, []) if c['host'] != n]
    if a.retime:
        # re-insert the target into its own route at another close approach (a different crossing), >= 60 d away
        for s, n, x in top:
            st_wo = rem[(n, x)]['st']; t_cur = float(fleet.routes[n]['st']['tf'][list(fleet.routes[n]['st']['asts']).index(x)])
            own = sorted(w_cands((n, st_wo, [x], a.dmax)), key=lambda c: c['dist'])
            jobs += [(n, st_wo, x, c['t']) for c in own if abs(c['t'] - t_cur) > 60 * DAY][:3]
    swap_jobs = []
    if a.swap:
        # X on A and Y on B exchanged: insert X into (B without Y) and Y into (A without X)
        topset = {x: n for s, n, x in top}
        for X, A in topset.items():
            for c in by.get(X, []):
                B = c['host']
                if B == A:
                    continue
                for Y, BY in topset.items():
                    if BY != B or frozenset((X, Y)) in {frozenset((j[0], j[2])) for j in swap_jobs}:
                        continue
                    cy = [d for d in by.get(Y, []) if d['host'] == A]
                    if not cy:
                        continue
                    swap_jobs.append((X, A, Y, B, c['t'], cy[0]['t']))
        swap_jobs = swap_jobs[:a.swap_max]
        for (X, A, Y, B, tX, tY) in swap_jobs:
            jobs.append((B + '|swap', rem[(B, Y)]['st'], X, tX))
            jobs.append((A + '|swap', rem[(A, X)]['st'], Y, tY))
    res = pool.map(w_insert, jobs)
    if a.swap:
        sw = {}
        for r in res:
            if r['name'].endswith('|swap') and r['ok']:
                host = r['name'][:-5]; sw.setdefault((host, r['ast']), []).append(r)
        res = [r for r in res if not r['name'].endswith('|swap')]
    moves = []
    src_of = {x: (s, n) for s, n, x in top}
    for r in res:
        if not r['ok']:
            continue
        s, n = src_of[r['ast']]
        if r['name'] == n:
            gain = cost(fleet.routes[n]['tank']) - cost(r['tank'])            # retime within the route
        else:
            gain = s - (cost(r['tank']) - cost(fleet.routes[r['name']]['tank']))
        if gain > 1e-3:
            moves.append((gain, r['ast'], n, r))
    if a.swap:
        for (X, A, Y, B, tX, tY) in swap_jobs:
            rb = min(sw.get((B, X), []), key=lambda r: r['tank'], default=None)
            ra = min(sw.get((A, Y), []), key=lambda r: r['tank'], default=None)
            if rb is None or ra is None:
                continue
            gain = cost(fleet.routes[A]['tank']) + cost(fleet.routes[B]['tank']) - cost(ra['tank']) - cost(rb['tank'])
            if gain > 1e-3:
                moves.append((gain, X, A, dict(name=B, st=rb['st'], tank=rb['tank'], t=rb['t'], swap=(Y, ra))))
    used = set(); applied = 0
    for gain, x, n, r in sorted(moves, key=lambda m: -m[0]):
        if n in used or r['name'] in used:
            continue
        if 'swap' in r:
            Y, ra = r['swap']; used |= {n, r['name']}
            fleet.routes[n] = dict(st=ra['st'], tank=ra['tank']); fleet.routes[r['name']] = dict(st=r['st'], tank=r['tank'])
            applied += 1
            say(f'  swapped {x}@{n} <-> {Y}@{r["name"]}: gain {gain:.4f}; J {fleet.J():.4f}')
            continue
        used |= {n, r['name']}
        if r['name'] == n:
            fleet.routes[n] = dict(st=r['st'], tank=r['tank'])
            say(f'  retimed {x} on {n} to {r["t"] / DAY:.0f} d: gain {gain:.4f}; J {fleet.J():.4f}')
        else:
            fleet.routes[n] = dict(st=rem[(n, x)]['st'], tank=rem[(n, x)]['tank'])
            fleet.routes[r['name']] = dict(st=r['st'], tank=r['tank'])
            say(f'  moved {x}: {n} -> {r["name"]} gain {gain:.4f}; J {fleet.J():.4f}')
        applied += 1
    return applied > 0




def w_retime(job):
    """Move one target to a different pass of the same asteroid: remove it, re-insert at the candidate event time.
    Uses the adaptive block homotopy (a single item), no ejection unless a margin is given."""
    name, st, ast, t_new, margin, tag = job
    ip = ipr(st); tank0 = ip.tank(); tic = time.time()
    remove_flybys(ip, [ast])
    ej = []
    if margin > 0:
        ej = [int(a) for a, t in zip(ip.asts, ip.tf) if abs(t - t_new) <= margin]
        if ej:
            remove_flybys(ip, ej)
    try:
        d = insert_block(ip, [(ast, t_new)], tol_km=100, ds0=0.2, iters_stage=6, log=QUIET)
        if d.max() > 1e4:
            return dict(tag=tag, name=name, ast=ast, ok=False, sec=time.time() - tic)
        miss = settle(ip, 100)
    except Exception as e:
        return dict(tag=tag, name=name, ast=ast, ok=False, why=repr(e)[:60], sec=time.time() - tic)
    tank = float(ip.tank())
    if miss > 150 or tank > TANK_MAX:
        return dict(tag=tag, name=name, ast=ast, ok=False, sec=time.time() - tic)
    return dict(tag=tag, name=name, ast=ast, ok=True, st=ist(ip), tank=tank, gain=cost(tank0) - cost(tank),
                t=t_new, ej=ej, sec=time.time() - tic)


def retime_pass(fleet, pool, a, say, ev):
    """For every target, try the OTHER passes of its asteroid (global event catalogue), ranked by the reduced-model
    bump cost (ctoc14/phasemodel): 8 K_DRIFT |dphi| / window for the phase plus 2 K_PLANE |dz| for the plane."""
    from ctoc14.phasemodel import bump_cost, plane_cost, phase_of
    E = eph(); YR = 365.25 * DAY
    hot = None
    if a.retime_top:
        # rank targets fleet-wide by the Delta-v of the leg that ends at them
        legs = []
        for n, r in fleet.routes.items():
            ip = ipr(r['st']); dvm = np.linalg.norm(ip.Ts, axis=1)
            edges = np.concatenate([[ip.tL], np.sort(ip.tf)])
            leg = np.clip(np.searchsorted(edges, ip.ts, side='right') - 1, 0, len(ip.tf) - 1)
            legdv = np.bincount(leg, dvm, minlength=len(ip.tf))
            for j, k in enumerate(np.argsort(ip.tf)):
                legs.append((legdv[j], n, int(ip.asts[k])))
        hot = {(n, X) for _, n, X in sorted(legs, reverse=True)[:a.retime_top]}
        say(f'retime: restricted to the {len(hot)} costliest legs')
    jobs = []
    for n, r in sorted(fleet.routes.items()):
        ip = ipr(r['st']); tt, rr = trajectory(ip, step=4 * DAY)
        th_c = phase_of(rr, tt, E); z_c = rr[:, 2] / AU
        tf = np.sort(ip.tf)
        for X, t_cur in zip(ip.asts, ip.tf):
            if hot is not None and (n, int(X)) not in hot:
                continue
            m = np.nonzero(ev['ast_id'] == X)[0]
            cands = []
            for j in m:
                te = float(ev['t'][j])
                if abs(te - t_cur) < 60 * DAY or te < ip.tL + 30 * DAY or te > T_MISSION - 30 * DAY:
                    continue
                k = int(np.clip(np.searchsorted(tt, te), 0, len(tt) - 1))
                nxt = tf[tf > te]; prv = tf[tf < te]
                win = ((nxt[0] if len(nxt) else T_MISSION) - (prv[-1] if len(prv) else ip.tL)) / YR
                dphi = abs((ev['th'][j] - th_c[k] + 180) % 360 - 180)
                est = bump_cost(dphi, win) + plane_cost(abs(ev['z'][j] - z_c[k]))
                cands.append((est, te))
            for est, te in sorted(cands)[:a.retime_cands]:
                jobs.append((n, r['st'], int(X), te, 0.0, f'{n}:{X}'))
    say(f'retime: {len(jobs)} trials over {sum(len(r["st"]["asts"]) for r in fleet.routes.values())} targets')
    res = [r for r in pool.map(w_retime, jobs) if r.get('ok') and r['gain'] > 1e-3]
    used = set(); applied = 0
    for r in sorted(res, key=lambda r: -r['gain']):
        if r['name'] in used:
            continue
        used.add(r['name']); fleet.routes[r['name']] = dict(st=r['st'], tank=r['tank']); applied += 1
        say(f'  retimed {r["ast"]} on {r["name"]} to {r["t"] / DAY:.0f} d: gain {r["gain"]:.4f}; J {fleet.J():.4f}')
    return applied


def cmd_retime(a):
    fleet = IFleet(a.src); out = pathlib.Path(a.dst); say = logger(out.with_suffix('.log'))
    ev = dict(np.load(a.events)); ev['ast_id'] = np.asarray(ev['ids'])[ev['ast']]
    say(f'retime start: {fleet.summary()}')
    with mp.get_context('fork').Pool(a.nproc, maxtasksperchild=40) as pool:
        for k in range(a.rounds):
            if not retime_pass(fleet, pool, a, say, ev):
                say('no improving retime'); break
            fleet.save(out, note=f'retime pass {k}'); say(f'fleet: {fleet.summary()}')
    fleet.save(out, note='final'); say(f'retime done: {fleet.summary()}')


# ------------------------------------------------------------------ commands
def cmd_import(a):
    src = pathlib.Path(a.src); out = pathlib.Path(a.dst); say = logger(out.with_suffix('.log'))
    jobs = [(f.stem[6:], str(f)) for f in sorted(src.glob('craft_*.npz'))]
    fleet = IFleet()
    with mp.get_context('fork').Pool(a.nproc) as pool:
        for name, st, tank, miss, m0 in pool.imap_unordered(w_import, jobs):
            fleet.routes[name] = dict(st=st, tank=tank)
            say(f'  {name}: exact m0 {m0:.1f} -> impulsive tank {tank:.1f} (miss {miss:.0f} km)')
    fleet.save(out, note=f'import {src}'); say(fleet.summary())


def w_resettle(job):
    name, st = job
    ip = ipr(st); miss = settle(ip, 150)
    return name, ist(ip), float(ip.tank()), miss


def cmd_resettle(a):
    fleet = IFleet(a.src); out = pathlib.Path(a.dst); say = logger(out.with_suffix('.log'))
    J0 = fleet.J()
    with mp.get_context('fork').Pool(a.nproc) as pool:
        for name, st, tank, miss in pool.imap_unordered(w_resettle, [(n, r['st']) for n, r in fleet.routes.items()]):
            say(f'  {name}: tank {fleet.routes[name]["tank"]:.1f} -> {tank:.1f} (miss {miss:.0f} km)')
            if miss <= 150:
                fleet.routes[name] = dict(st=st, tank=tank)
    fleet.save(out, note=f'resettle of {a.src}'); say(f'resettle: J {J0:.4f} -> {fleet.J():.4f}; {fleet.summary()}')


def cmd_search(a):
    fleet = IFleet(a.src); out = pathlib.Path(a.dst); say = logger(out.with_suffix('.log'))
    say(f'search start: {fleet.summary()}')
    tried = set()
    with mp.get_context('fork').Pool(a.nproc, maxtasksperchild=50) as pool:
        for rnd in range(a.rounds):
            improved = False
            if a.relocate and not (a.victims and rnd == 0):
                while relocate_pass(fleet, pool, a, say):
                    improved = True; fleet.save(out, note=f'round {rnd} relocate'); say(f'fleet: {fleet.summary()}')
            order = sorted((n for n in fleet.routes if n not in tried), key=lambda n: len(fleet.routes[n]['st']['asts']))
            if a.victims:
                order = [v for v in a.victims.split(',') if v in fleet.routes and v not in tried] + order
            for V in order[:a.dissolve_k]:
                tried.add(V)
                new = dissolve_lns(fleet, V, pool, a, say) if a.lns else dissolve(fleet, V, pool, a, say)
                if new is not None:
                    fleet = new; improved = True; tried = set()
                    fleet.save(out, note=f'round {rnd} dissolved {V}'); say(f'fleet: {fleet.summary()}')
                    break
            if not improved:
                say('no improving move'); break
    fleet.save(out, note='final'); say(f'search done: {fleet.summary()}')


def cmd_export(a):
    fleet = IFleet(a.src); out = pathlib.Path(a.dst); out.mkdir(parents=True, exist_ok=True); say = logger(out.with_suffix('.log'))
    with mp.get_context('fork').Pool(a.nproc, maxtasksperchild=1) as pool:
        for name, m0, miss, tank, ok, hs in pool.imap_unordered(w_export, [(n, r['st'], str(out)) for n, r in fleet.routes.items()]):
            say(f'  {name}: impulsive tank {tank:.1f} -> exact m0 {m0:.1f} (miss {miss:.0f} km), valid {ok} (RK4 step {hs} d)')
    sys.path.insert(0, str(pathlib.Path(__file__).parent))
    from run_gfleet import Fleet
    ef = Fleet(out); ef.save(out, note=f'export of {a.src}'); say(ef.summary())


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('cmd'); ap.add_argument('src'); ap.add_argument('dst')
    ap.add_argument('--nproc', type=int, default=8); ap.add_argument('--rounds', type=int, default=4)
    ap.add_argument('--dmax', type=float, default=0.3); ap.add_argument('--per-target', type=int, default=8)
    ap.add_argument('--slack', type=float, default=0.1); ap.add_argument('--dissolve-k', type=int, default=13)
    ap.add_argument('--relocate', action='store_true'); ap.add_argument('--relocate-top', type=int, default=40)
    ap.add_argument('--min-saving', type=float, default=0.01); ap.add_argument('--victims', default=None)
    ap.add_argument('--retime', action='store_true'); ap.add_argument('--swap', action='store_true'); ap.add_argument('--swap-max', type=int, default=60); ap.add_argument('--anneal', type=float, default=0.0); ap.add_argument('--anneal-passes', type=int, default=4)
    ap.add_argument('--lns', action='store_true'); ap.add_argument('--bmax', type=int, default=5)
    ap.add_argument('--hosts', type=int, default=3); ap.add_argument('--lns-rounds', type=int, default=40)
    ap.add_argument('--lam', type=float, default=0.0); ap.add_argument('--lam-max', type=float, default=0.25)
    ap.add_argument('--tabu', type=int, default=3)
    ap.add_argument('--margins', default='0,30'); ap.add_argument('--retime-cands', type=int, default=3); ap.add_argument('--retime-top', type=int, default=0)
    ap.add_argument('--events', default='results/newgen/scratch/events2.npz')
    a = ap.parse_args()
    a.margins = [float(x) for x in a.margins.split(',')]
    dict(**{'import': cmd_import, 'search': cmd_search, 'export': cmd_export, 'resettle': cmd_resettle,
            'retime': cmd_retime})[a.cmd](a)


if __name__ == '__main__':
    main()
