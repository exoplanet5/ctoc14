"""Fleet-level search on whole-trajectory-optimised craft (ctoc14/globalopt.py).

Fleet = directory of craft states craft_<name>.npz (tL, m0, ts, Ts, tf, asts, vinf). Fleet J = sum cost(m0) + misses.
Moves (every candidate evaluated with the nonlinear SCP, parallel over craft):
  init      import a submission (converted craft or go_sc<k>.npz checkpoints of tools/run_globalopt.py)
  dedupe    repeated flybys are free waypoints for the Lambert planner but pure cost here: keep each asteroid on the craft
            whose linearised removal saving is smallest, remove it from the others, re-optimise and tighten
  dissolve  remove craft V: insert its unique targets one by one (hardest first) into host trajectories, candidates =
            close approaches of each host to the target screened with the linear model, best nonlinear trial wins;
            accept when sum dJ(hosts) < cost(V)
  polish    more SCP iterations + tightening on every craft
  write     rows with the validator integrator, validation per craft, combined submission
Usage:
  run_gfleet.py init <fleetdir> --from-sub results/CTOC14_Result_TEAM.txt [--from-ckpt results/newgen/go1]
  run_gfleet.py dedupe <fleetdir> <outdir>
  run_gfleet.py dissolve <fleetdir> <outdir> [--victims 11,6 | --auto 3]
  run_gfleet.py polish <fleetdir> <outdir> [--iters 60]
  run_gfleet.py write <fleetdir> <submission.txt>"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse, shutil, multiprocessing as mp
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from ctoc14.kepler import Ephemeris, propagate_batch
from ctoc14.globalopt import (load_craft, GlobalProblem, optimise, restore, add_flyby, remove_flybys, write_craft, irls_min_fuel,
                              state_of, problem_from_state, save_state, load_state, tighten, extend_arc, to_massless, polish_ml, insert_homotopy)
from ctoc14.constants import DAY, AU, T_MISSION, M_DRY

EPS_FAST = (0.02, 0.005, 0.001)
_EPH = None


def eph():
    global _EPH
    if _EPH is None:
        _EPH = Ephemeris()
    return _EPH


def cost(m0):
    x = (m0 - M_DRY) / 1400.0
    return 1.0 + x + x * x


class Fleet:
    def __init__(self, d=None):
        self.states = {}
        if d is not None:
            for f in sorted(pathlib.Path(d).glob('craft_*.npz')):
                self.states[f.stem[6:]] = load_state(f)

    def coverage(self):
        cov = {}
        for n, st in self.states.items():
            for a in st['asts']:
                cov.setdefault(int(a), set()).add(n)
        return cov

    def J(self):
        return sum(cost(float(st['m0'])) for st in self.states.values()) + 300 - len(self.coverage())

    def save(self, d, note=''):
        d = pathlib.Path(d); d.mkdir(parents=True, exist_ok=True)
        for f in d.glob('craft_*.npz'):
            f.unlink()
        for n, st in self.states.items():
            save_state(st, d / f'craft_{n}.npz')
        cov = self.coverage()
        json.dump(dict(J=self.J(), n=len(self.states), covered=len(cov), note=note,
                       craft={n: dict(m0=float(st['m0']), J_i=cost(float(st['m0'])), flybys=len(st['asts']),
                                      unique=sum(1 for a in st['asts'] if cov[int(a)] == {n})) for n, st in self.states.items()}),
                  open(d / 'fleet.json', 'w'), indent=1)

    def summary(self):
        cov = self.coverage()
        return (f'{len(self.states)} craft, covered {len(cov)}, J {self.J():.4f}; ' +
                ' '.join(f'{n}:{len(st["asts"])}/{sum(1 for a in st["asts"] if cov[int(a)] == {n})}@{float(st["m0"]):.0f}'
                         for n, st in sorted(self.states.items(), key=lambda kv: -len(kv[1]['asts']))))


def logger(path):
    f = open(path, 'a')
    def say(s):
        line = f'[{time.strftime("%H:%M:%S")}] {s}'; print(line, flush=True); f.write(line + '\n'); f.flush()
    return say


# ------------------------------------------------------------------ workers
def w_polish(job):
    name, st, iters, margin = job
    gp = problem_from_state(eph(), st)
    gt, d = polish_ml(gp, iters, margin)
    return name, state_of(gt), float(d.max()), gt.fuel()


def w_remove(job):
    name, st, asts, iters, margin = job
    gp = problem_from_state(eph(), st); remove_flybys(gp, asts)
    gt, d = polish_ml(gp, iters, margin)
    return name, state_of(gt), float(d.max()), gt.fuel()


def close_approaches(E, gp, Ys, ast, dmax_km, t_min):
    ra, _ = E.ast_states_at(np.full(len(gp.ts), ast - 1), gp.ts)
    d = np.linalg.norm(Ys[:, 0:3] - ra, axis=1); out = []
    for i in range(1, len(d) - 1):
        if d[i] <= d[i - 1] and d[i] < d[i + 1] and d[i] < dmax_km and gp.ts[i] >= t_min:
            a, b, c = d[i - 1], d[i], d[i + 1]; den = a - 2 * b + c
            off = 0.5 * (a - c) / den if den > 0 else 0.0
            out.append((float(gp.ts[i] + off * gp.h), float(b), i))
    return out


def linear_setup(gp):
    Yf, Ys = gp.integrate(dense_samples=True); dm, vr = gp.misses(Yf)
    ch = gp.chain(Yf, Ys); kind = ch['kind']
    node_of_f = {int(-2 - kind[n]): n for n in range(len(kind)) if kind[n] <= -2}
    B, Lv = gp.rows_at(ch, gp.tf, nodes=[node_of_f[j] for j in range(len(gp.tf))])
    return Yf, Ys, dm, vr, ch, B, Lv


def w_screen(job):
    name, st, targets, dmax = job
    E = eph(); gp = to_massless(problem_from_state(E, st)); extend_arc(gp, T_MISSION - 2 * DAY)
    Yf, Ys, dm, vr, ch, B, Lv = linear_setup(gp)
    T0, _, _ = irls_min_fuel(B, Lv, vr, dm, gp.Ts, gp.vinf, gp.tcap, eps_sched=EPS_FAST, inner=3); f0 = gp.fuel(T0)
    res = []
    for X in targets:
        if X in set(int(a) for a in gp.asts):
            continue
        for (tc, dist, i) in close_approaches(E, gp, Ys, X, dmax * AU, gp.tL + 20 * DAY):
            if tc > T_MISSION - 3 * DAY:
                continue
            rs, vs = propagate_batch(Ys[i, 0:3][None], Ys[i, 3:6][None], np.array([tc - gp.ts[i]]))
            ra, va = E.ast_states_at(np.array([X - 1]), np.array([tc]))
            Bc, Lvc = gp.rows_at(ch, [tc])
            T1, _, _ = irls_min_fuel(np.concatenate([B, Bc]), np.concatenate([Lv, Lvc]), np.concatenate([vr, vs - va]),
                                     np.concatenate([dm, rs - ra]), gp.Ts, gp.vinf, gp.tcap, eps_sched=EPS_FAST, inner=3)
            res.append(dict(host=name, ast=int(X), t=tc, dist_au=dist / AU, dF=gp.fuel(T1) - f0))
    return name, res


def removal_savings(job):
    """Linear-model propellant saving of removing each listed flyby from the craft."""
    name, st, asts = job
    gp = problem_from_state(eph(), st)
    Yf, Ys, dm, vr, ch, B, Lv = linear_setup(gp)
    T0, _, _ = irls_min_fuel(B, Lv, vr, dm, gp.Ts, gp.vinf, gp.tcap, eps_sched=EPS_FAST, inner=3); f0 = gp.fuel(T0)
    out = {}
    for X in asts:
        keep = np.array([int(a) != X for a in gp.asts])
        T1, _, _ = irls_min_fuel(B[keep], Lv[keep], vr[keep], dm[keep], gp.Ts, gp.vinf, gp.tcap, eps_sched=EPS_FAST, inner=3)
        out[X] = f0 - gp.fuel(T1)
    return name, out


def w_insert(job):
    name, st, X, t, iters, margin = job
    gp = problem_from_state(eph(), st); m0_old = gp.m0
    try:
        gm = to_massless(gp)
        d = insert_homotopy(gm, X, t, stages=10, verbose=False)
        if d.max() > 300:
            return dict(host=name, ast=X, t=t, ok=False, miss=float(d.max()), fuel=gm.fuel())
        gt, d = polish_ml_from_massless(gm, iters, margin)
    except Exception as e:
        return dict(host=name, ast=X, t=t, ok=False, miss=np.inf, error=repr(e))
    ok = bool(d.max() <= 300)
    return dict(host=name, ast=X, t=t, ok=ok, miss=float(d.max()), fuel=gt.fuel(), m0=gt.m0, dJ=cost(gt.m0) - cost(m0_old),
                state=state_of(gt) if ok else None)


def polish_ml_from_massless(gm, iters, margin):
    from ctoc14.globalopt import to_thrust
    optimise(gm, iters=iters, rho0=100.0, verbose=False)
    restore(gm, tol_km=150, iters=20, verbose=False)
    gt = to_thrust(gm, margin)
    restore(gt, tol_km=150, iters=20, verbose=False)
    d = tighten(gt, margin)
    return gt, d


# ------------------------------------------------------------------ commands
def cmd_init(a):
    E = eph(); fleet = Fleet()
    from ctoc14.validator import parse
    scs = sorted({r.sc for r in parse(a.from_sub)})
    for sc in scs:
        ck = pathlib.Path(a.from_ckpt) / f'sc{sc}_s3.npz' if a.from_ckpt else None
        gp = GlobalProblem(E, load_craft(a.from_sub, sc), hstep=0.25 * DAY)
        if ck is not None and ck.exists():
            z = np.load(ck); gp.Ts, gp.tf, gp.vinf, gp.m0 = z['Ts'], z['tf'], z['vinf'], float(z['m0'])
            gp.asts = [int(x) for x in z['asts']]
        fleet.states[f'{sc:02d}'] = state_of(gp)
    fleet.save(a.fleet, note=f'init from {a.from_sub} {a.from_ckpt}')
    print(fleet.summary())


def cmd_polish(a):
    fleet = Fleet(a.fleet); say = logger(pathlib.Path(a.out).with_suffix('.log')); say(f'polish start: {fleet.summary()}')
    with mp.get_context('fork').Pool(a.nproc, maxtasksperchild=1) as pool:
        for name, st, miss, fuel in pool.imap_unordered(w_polish, [(n, st, a.iters, a.margin) for n, st in fleet.states.items()]):
            if miss <= 300:
                old = float(fleet.states[name]['m0']); fleet.states[name] = st
                say(f'  {name}: m0 {old:.1f} -> {float(st["m0"]):.1f} (miss {miss:.0f} km)')
            else:
                say(f'  {name}: polish failed (miss {miss:.0f} km), kept')
    fleet.save(a.out, note='polish'); say(f'polish done: {fleet.summary()}')


def cmd_dedupe(a):
    fleet = Fleet(a.fleet); say = logger(pathlib.Path(a.out).with_suffix('.log')); say(f'dedupe start: {fleet.summary()}')
    cov = fleet.coverage(); dups = {x: S for x, S in cov.items() if len(S) > 1}
    say(f'{len(dups)} asteroids covered more than once')
    with mp.get_context('fork').Pool(a.nproc, maxtasksperchild=1) as pool:
        per = {}
        for x, S in dups.items():
            for n in S: per.setdefault(n, []).append(x)
        sav = dict(pool.map(removal_savings, [(n, fleet.states[n], xs) for n, xs in per.items()]))
        remove = {}
        for x, S in dups.items():
            uniq = {n: sum(1 for y in fleet.states[n]['asts'] if cov[int(y)] == {n}) for n in S}
            keep = max(S, key=lambda n: (uniq[n], len(fleet.states[n]['asts']), -sav[n][x])) if a.keep_big else min(S, key=lambda n: sav[n][x])
            for n in S - {keep}:
                remove.setdefault(n, []).append(x)
            say(f'  ast {x}: keep on {keep} (savings ' + ', '.join(f'{n} {sav[n][x]:.1f}' for n in sorted(S)) + ')')
        for name, st, miss, fuel in pool.imap_unordered(w_remove, [(n, fleet.states[n], xs, a.iters, a.margin) for n, xs in remove.items()]):
            if miss <= 300:
                old = float(fleet.states[name]['m0']); fleet.states[name] = st
                say(f'  {name}: removed {remove[name]}, m0 {old:.1f} -> {float(st["m0"]):.1f}')
            else:
                say(f'  {name}: removal failed (miss {miss:.0f})')
    fleet.save(a.out, note='dedupe'); say(f'dedupe done: {fleet.summary()}')


def dissolve_one(fleet, V, pool, a, say):
    cov = fleet.coverage(); J0 = fleet.J(); cV = cost(float(fleet.states[V]['m0']))
    targets = sorted(int(x) for x in fleet.states[V]['asts'] if cov[int(x)] == {V})
    hosts = {n: st for n, st in fleet.states.items() if n != V}
    say(f'dissolve {V} (cost {cV:.4f}): {len(targets)} unique targets {targets}')
    screen = {}
    todo = set(hosts)
    spent = 0.0; placed = []
    while targets:
        for name, res in pool.imap_unordered(w_screen, [(n, hosts[n], targets, a.dmax) for n in todo]):
            screen[name] = res
        todo = set()
        by_t = {X: sorted([c for r in screen.values() for c in r if c['ast'] == X], key=lambda c: c['dF']) for X in targets}
        X = max(targets, key=lambda X: by_t[X][0]['dF'] if by_t[X] else np.inf)
        cands = by_t[X][:a.trials]
        if not cands:
            say(f'  target {X}: no close approach on any host -> abort'); return None
        say(f'  target {X}: trying ' + ', '.join(f"{c['host']}@{c['t'] / DAY:.0f}d({c['dist_au']:.3f}AU,{c['dF']:.0f}kg)" for c in cands))
        results = pool.map(w_insert, [(c['host'], hosts[c['host']], X, c['t'], a.iters, a.margin) for c in cands])
        good = [r for r in results if r['ok']]
        for r in results:
            say(f"    {r['host']}@{r['t'] / DAY:.0f}d: ok={r['ok']} miss {r['miss']:.0f} km" + (f", fuel {r['fuel']:.1f}, dJ {r['dJ']:+.4f}" if r['ok'] else ''))
        if not good:
            say(f'  target {X}: no successful insertion -> abort'); return None
        best = min(good, key=lambda r: r['dJ'])
        spent += best['dJ']; hosts[best['host']] = best['state']; placed.append((X, best['host'], best['dJ']))
        targets.remove(X); todo = {best['host']}
        say(f'  placed {X} on {best["host"]} (dJ {best["dJ"]:+.4f}); spent {spent:.4f} of {cV:.4f}')
        if spent > cV + a.slack:
            say(f'  spent more than the craft costs -> abort'); return None
    new = Fleet(); new.states = dict(hosts)
    say(f'dissolve {V}: J {J0:.4f} -> {new.J():.4f} ({"ACCEPT" if new.J() < J0 - 1e-6 else "reject"}); placements {placed}')
    return new if new.J() < J0 - 1e-6 else None


def cmd_dissolve(a):
    fleet = Fleet(a.fleet); out = pathlib.Path(a.out); say = logger(out.with_suffix('.log'))
    say(f'dissolve start: {fleet.summary()}')
    tried = set()
    with mp.get_context('fork').Pool(a.nproc, maxtasksperchild=4) as pool:
        victims = [v for v in a.victims.split(',')] if a.victims else None
        rounds = len(victims) if victims else a.auto
        for k in range(rounds):
            if victims:
                V = victims[k]
            else:
                cov = fleet.coverage()
                order = sorted((n for n in fleet.states if n not in tried),
                               key=lambda n: sum(1 for x in fleet.states[n]['asts'] if cov[int(x)] == {n}))
                if not order: break
                V = order[0]
            tried.add(V)
            new = dissolve_one(fleet, V, pool, a, say)
            if new is not None:
                fleet = new; fleet.save(out, note=f'dissolved {V}'); say(f'fleet now: {fleet.summary()}')
    fleet.save(out, note='dissolve final'); say(f'dissolve done: {fleet.summary()}')


def w_write(job):
    name, st, path = job
    gp = problem_from_state(eph(), st)
    rep, fuel, pk = write_craft(gp, path, eph=eph())
    return name, bool(rep.ok), rep.errors[:3], len(rep.flybys), fuel, pk


def cmd_write(a):
    fleet = Fleet(a.fleet); d = pathlib.Path(a.fleet) / 'frags'; d.mkdir(exist_ok=True)
    jobs = [(n, st, str(d / f'frag_{n}.txt')) for n, st in sorted(fleet.states.items())]
    with mp.get_context('fork').Pool(a.nproc, maxtasksperchild=1) as pool:
        res = pool.map(w_write, jobs)
    for r in res:
        print(f'{r[0]}: ok={r[1]} errors={r[2]} flybys={r[3]} fuel={r[4]:.2f} peak={r[5]:.3f}')
    good = [j[2] for j, r in zip(jobs, res) if r[1]]
    if len(good) < len(jobs):
        print('some craft failed validation; not combining'); return
    import subprocess
    subprocess.run([sys.executable, str(pathlib.Path(__file__).parent / 'combine_submission.py'), a.out] + good, check=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd'); ap.add_argument('fleet'); ap.add_argument('out', nargs='?')
    ap.add_argument('--from-sub'); ap.add_argument('--from-ckpt')
    ap.add_argument('--iters', type=int, default=60); ap.add_argument('--margin', type=float, default=1.5)
    ap.add_argument('--nproc', type=int, default=8); ap.add_argument('--dmax', type=float, default=0.3)
    ap.add_argument('--trials', type=int, default=4); ap.add_argument('--slack', type=float, default=0.0)
    ap.add_argument('--victims'); ap.add_argument('--auto', type=int, default=0)
    ap.add_argument('--keep-big', action='store_true', help='dedupe: keep duplicates on the craft with most unique targets')
    a = ap.parse_args()
    dict(init=cmd_init, polish=cmd_polish, dedupe=cmd_dedupe, dissolve=cmd_dissolve, write=cmd_write)[a.cmd](a)


if __name__ == '__main__':
    main()
