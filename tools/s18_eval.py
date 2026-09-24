"""Stage 18 / WP3: the FLEET REPORT and GATE VERDICTS of an RHFA run (or of any RI.IFleet directory).

Input   a run dir (results/s18/<run>: state.json + fleet/route_*.npz + fleet/fleet.json) or a bare fleet dir of
        route_*.npz files (results/s16/best/fleet works and is the reference: covered 298, sum J_i 11.7227, honest).
        A run dir without fleet/ (unfinished run) is read from its craft_*.npz; <run>/fleet also works.
Report  (REPORT dict, also written as <run>/eval.json and printed as a table)
  covered, misses [ids], misses_by_class {pin, med, fil}, n_craft, tanks, sumJi, fuel (sum tank - 600), J_raw
  (sumJi + misses + 2), per route: {name, n, tank, J_i, t_launch_d, t_end_d, cadence_d = (t_end - t_launch) / n,
  kg_per_fb, counts, honest: C.honest_check (irregular, max_window_cap, miss_km, ok)}, honest_all (bool),
  flybys_by_1461d (distinct targets flown by day 1461, for G0b), mean_fuel_by_1461d (per craft launched by day 1461:
  the propellant burnt by the impulses before day 1461, tank (1 - exp(-dv_before / ve)) -- MEASURED from the twins;
  when state.json has a history, mean_fuel_by_1461d_state = the committed honest tank - 600 after the last slice whose
  window ends <= 1461 d is reported next to it for information), slices (from state.json log when present: covered /
  fuel per slice).
Gates   (docs/stage18_rhfa.md section 3; the verdict strings are fixed here, the builders do not move the numbers)
  G0b  flybys_by_1461d >= 72 and mean_fuel_by_1461d <= 160 -> PASS; <= 60 flybys -> KILL; else WEAK
       (only for an RHFA run that has reached day 1461; a bare fleet prints N/A with the would-be verdict in the notes)
  G1   N = 8 full run: covered >= 268 and fuel <= 3900 -> PASS (portfolio); 250-267 (or >= 268 over the fuel bar) ->
       ROUND (one parameter round); < 250 -> KILL (only for a finished run; N/A otherwise, with the would-be verdict)
  G2   covered >= 292 and sumJi + lin-priced closure of the misses <= 11.50 -> EXPORT; the closure price is
       estimated with tools/s15b_close._price_host over the fleet (lin_price screen, one miss per host, greedy; an
       unplaced miss counts 1.0), marked ESTIMATE; fallback: 9 craft with sumJi + misses <= 11.64 -> EXPORT9; else NO.
Every number in the report is MEASURED from the twins on disk (one integration per route, no optimisation), except
the G2 closure price, which is an estimate and is labelled so.

usage: s18_eval.py RUN_OR_FLEET_DIR [--json OUT.json] [--no-g2] [--nproc 2] [--drop r9[,r8]] [--dmax 0.2]
--drop evaluates the fleet without the named routes (e.g. s16a minus r9 = the 8-craft view with 21 misses).
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, glob, time, pathlib, argparse, re
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import numpy as np
import run_ialns as RI
import s18_common as C
from ctoc14.constants import DAY, VE

GATES = dict(g0b_flybys=72, g0b_kill=60, g0b_fuel=160.0, g1_pass=268, g1_round=250, g1_fuel=3900.0, g2_cover=292,
             g2_sumJi=11.50, g2_fallback9=11.64, day_check=1461.0)
LIN_CAL_HI = 0.73       # true / lin above 0.2 J (docs/stage17/R1_redteam.md, 75 real insertions); ~1.0 below 0.06 J
LIN_CAL_BAND = 0.2


def _natkey(name):
    m = re.match(r'^([a-zA-Z_]*)(\d+)$', str(name))
    return (m.group(1), int(m.group(2))) if m else (str(name), 0)


def _resolve(path):
    p = pathlib.Path(path)
    return p if p.is_absolute() else (ROOT / p)


def load_fleet(path):
    """-> (routes: [{name, file, st (ist dict), ip}], state: STATE dict | None).  A run dir uses fleet/route_*.npz (or
    craft_*.npz when fleet/ is absent, i.e. an unfinished run) and state.json; a bare dir uses its route_*.npz.
    The returned state dict carries '_run_dir' (the run directory) and '_source' ('fleet' | 'craft')."""
    p = _resolve(path)
    if not p.is_dir():
        raise FileNotFoundError(f'{p} is not a directory')
    rd = None
    if (p / 'state.json').exists():
        rd = p
    elif p.name == 'fleet' and (p.parent / 'state.json').exists():
        rd = p.parent
    state = None
    if rd is not None:
        state = C.jload(rd / 'state.json')
        try:
            C.check_state(state)
        except AssertionError as e:
            state['_invalid'] = str(e)
        state['_run_dir'] = str(rd)
    if rd is not None and (rd / 'fleet').is_dir() and list((rd / 'fleet').glob('route_*.npz')):
        files = sorted((rd / 'fleet').glob('route_*.npz'), key=lambda f: _natkey(f.stem[6:])); src = 'fleet'
    elif rd is not None and not list(p.glob('route_*.npz')):
        files = sorted(rd.glob('craft_*.npz'), key=lambda f: _natkey(f.stem[6:])); src = 'craft'
    else:
        files = sorted(p.glob('route_*.npz'), key=lambda f: _natkey(f.stem[6:])); src = 'fleet'
    if state is not None:
        state['_source'] = src
    routes = []
    for f in files:
        name = f.stem[6:] if f.stem.startswith('route_') else 'c' + f.stem[6:]
        st = C.load_twin(f)
        routes.append(dict(name=name, file=str(f), st=st, ip=C.twin_of(st)))
    return routes, state


def route_report(name, st, W, ip=None):
    """Per-route facts (module docstring) from one integration: C.honest_check (which integrates once for the miss) +
    the cheap twin facts + classes + cadence."""
    ip = ip or C.twin_of(st); o = C.order_of(st)
    hc = C.honest_check(ip)
    targets = [int(st['asts'][i]) for i in o]; epochs = [float(st['tf'][i]) for i in o]
    n = len(targets); tank = float(hc['tank']); tL = float(st['tL']); t_end = max(epochs) if epochs else tL
    return dict(name=str(name), n=n, targets=targets, epochs_d=[round(e / DAY, 2) for e in epochs], tank=tank,
                J_i=C.cost(tank), dv=float(hc['dv']), t_launch_d=round(tL / DAY, 2), t_end_d=round(t_end / DAY, 2),
                cadence_d=round((t_end - tL) / DAY / n, 1) if n else None,
                kg_per_fb=round((tank - 600.0) / n, 2) if n else None, counts=C.count_classes(W, targets),
                honest=dict(irregular=int(hc['irregular']), max_imp_cap=hc['max_imp_cap'], max_window_cap=hc['max_window_cap'],
                            miss_km=hc['miss_km'], ok=bool(hc['ok'])))


def flybys_by(routes, t_d):
    """Distinct targets flown (flyby epoch <= t_d DAYS) across the routes, and the per-route fuel spent by then
    (sum of |impulses| before t_d converted to kg at the route's tank: tank (1 - exp(-dv_before / ve))).
    -> dict(t_d, n, targets, launched [names launched by t_d], fuel {name: kg}, mean_fuel (over the launched), n_by_route)."""
    ts_ = float(t_d) * DAY
    seen = set(); fuel = {}; launched = []; nby = {}
    for r in routes:
        st = r['st']; tf = np.asarray(st['tf'], float); asts = np.asarray(st['asts'], int)
        m = tf <= ts_ + 1e-6
        seen |= set(int(a) for a in asts[m]); nby[r['name']] = int(m.sum())
        if float(st['tL']) <= ts_:
            launched.append(r['name'])
            ts = np.asarray(st['ts'], float); Ts = np.asarray(st['Ts'], float).reshape(-1, 3)
            dvb = float(np.linalg.norm(Ts[ts <= ts_], axis=1).sum()) if len(ts) else 0.0
            tank = float(r['ip'].tank()) if r.get('ip') is not None else float(C.twin_of(st).tank())
            fuel[r['name']] = round(tank * (1.0 - np.exp(-dvb / VE)), 2)
    return dict(t_d=float(t_d), n=len(seen), targets=sorted(seen), launched=launched, fuel=fuel,
                mean_fuel=round(float(np.mean(list(fuel.values()))), 2) if fuel else None, n_by_route=nby)


def closure_estimate(routes, misses, nproc=2, dmax=0.2):
    """ESTIMATE of the twin cost of closing `misses` into these routes: s15b_close._price_host per host (approaches
    within dmax AU + s14_twinbeam.lin_price), then greedy one miss per host, cheapest first; returns
    dict(dJ_lin, dJ_total (= dJ_lin + 1.0 per unplaced miss), dJ_cal (lin x 0.73 above 0.2 J), placed, unplaced,
    per_miss [...], n_candidates, secs, note='ESTIMATE (lin price, calibration 0.73)')."""
    misses = sorted(int(m) for m in misses)
    note = ('ESTIMATE (lin price via s15b_close._price_host, greedy one miss per host, cheapest first; true/lin ~1.0 '
            'below 0.06 J and 0.73 above 0.2 J per docs/stage17/R1_redteam.md; an unplaced miss counts 1.0)')
    tic = time.time()
    if not misses or not routes:
        return dict(dJ_lin=0.0, dJ_total=0.0, dJ_cal=0.0, placed=[], unplaced=list(misses), per_miss=[], n_candidates=0,
                    secs=0.0, dmax=dmax, note=note)
    import s15b_close as SC       # chdir(ROOT) at import: every path here is already absolute
    jobs = [(r['name'], r['st'], list(misses), float(dmax)) for r in routes]
    if int(nproc) <= 1 or len(jobs) <= 1:
        R = [SC._price_host(j) for j in jobs]
    else:
        import multiprocessing as mp
        with mp.get_context('fork').Pool(min(int(nproc), len(jobs)), maxtasksperchild=1) as pool:
            R = pool.map(SC._price_host, jobs, chunksize=1)
    cands = sorted((c for lst in R for c in lst), key=lambda c: c['lin_dJ'])
    used = set(); placed = {}
    for c in cands:
        if c['host'] in used or c['ast'] in placed:
            continue
        placed[c['ast']] = c; used.add(c['host'])
    best = {}
    for c in cands:
        if c['ast'] not in best:
            best[c['ast']] = c
    per_miss = []
    for m in misses:
        b = best.get(m); p = placed.get(m)
        per_miss.append(dict(ast=m, cheapest_host=b['host'] if b else None, cheapest_lin_dJ=round(b['lin_dJ'], 4) if b else None,
                             cheapest_dist_au=round(b['dist'], 3) if b else None,
                             placed_host=p['host'] if p else None, placed_lin_dJ=round(p['lin_dJ'], 4) if p else None,
                             placed_t_d=round(p['t'] / DAY, 1) if p else None, n_hosts=len(set(c['host'] for c in cands if c['ast'] == m))))
    unplaced = [m for m in misses if m not in placed]
    dJ_lin = float(sum(c['lin_dJ'] for c in placed.values()))
    dJ_cal = float(sum(c['lin_dJ'] * (LIN_CAL_HI if c['lin_dJ'] > LIN_CAL_BAND else 1.0) for c in placed.values()))
    # the red team's "screen": every miss at its cheapest host, no per-host limit (ignores tank coupling; 1.0 when no host)
    dJ_screen = float(sum(best[m]['lin_dJ'] if m in best else 1.0 for m in misses))
    return dict(dJ_lin=round(dJ_lin, 4), dJ_total=round(dJ_lin + len(unplaced), 4), dJ_cal=round(dJ_cal + len(unplaced), 4),
                dJ_screen=round(dJ_screen, 4), n_without_host=sum(1 for m in misses if m not in best),
                placed=[dict(ast=m, host=c['host'], lin_dJ=round(c['lin_dJ'], 4), lin_kg=round(c['lin_kg'], 1),
                             dist_au=round(c['dist'], 3), t_d=round(c['t'] / DAY, 1)) for m, c in sorted(placed.items())],
                unplaced=unplaced, per_miss=per_miss, n_candidates=len(cands), secs=round(time.time() - tic, 1), dmax=dmax,
                note=note)


def _state_fuel_by(state, t_d):
    """Committed honest tank - 600 per launched craft after the last slice whose window (T, T + H] ends <= t_d (from
    state.json history + log); None when the state has no usable history."""
    try:
        H = float(state['params']['H']); by_slice = {int(l['slice']): float(l['T_d']) for l in state.get('log', [])}
        out = {}
        for c in state['craft']:
            hist = [h for h in c.get('history', []) if int(h['slice']) in by_slice and by_slice[int(h['slice'])] + H <= t_d + 1e-6]
            if hist:
                out[str(c['i'])] = round(float(hist[-1]['tank']) - 600.0, 2)
        return out or None
    except Exception:
        return None


def report(path, W=None, g2=True, nproc=2, drop=(), dmax=0.2):
    """The REPORT dict (module docstring) with 'gates' = gates(report)."""
    tic = time.time()
    W = W or C.load_windows()
    routes, state = load_fleet(path)
    drop = set(drop or ())
    if drop:
        routes = [r for r in routes if r['name'] not in drop]
    rr = [route_report(r['name'], r['st'], W, r['ip']) for r in routes]
    covered = set(); dup = {}
    for r in rr:
        for t in r['targets']:
            if t in covered:
                dup.setdefault(t, []).append(r['name'])
            covered.add(t)
    misses = sorted(set(C.TARGETS) - covered)
    tanks = [r['tank'] for r in rr]
    tank_note = 'honest twin tanks (C.honest_check per route)'
    if state is not None and state.get('params', {}).get('commit', 'regrid') == 'none':
        # fixture_fake run: the twins on disk are placeholders; the state's (NOT honest) tanks are the only tank numbers
        by_i = {int(c['i']): float(c['tank']) for c in state['craft'] if c['launched']}
        for r in rr:
            m = re.search(r'(\d+)$', r['name'])
            if m and int(m.group(1)) in by_i:
                r['tank_twin'] = r['tank']; r['tank'] = by_i[int(m.group(1))]; r['J_i'] = C.cost(r['tank'])
                r['kg_per_fb'] = round((r['tank'] - 600.0) / r['n'], 2) if r['n'] else None
        tanks = [r['tank'] for r in rr]; tank_note = 'tanks from state.json (commit none: NOT honest, fixture only)'
    fb = flybys_by(routes, GATES['day_check'])
    rep = dict(path=str(_resolve(path)), evaluated=C.now(), n_craft=len(rr), covered=len(covered), misses=misses,
               misses_by_class=C.count_classes(W, misses), duplicates={str(k): v for k, v in dup.items()},
               tanks=[round(t, 3) for t in tanks], sumJi=round(C.sum_Ji(tanks), 4), fuel=round(C.fuel_of(tanks), 1),
               J_raw=round(C.fleet_J(tanks, len(covered)), 4), routes=rr, honest_all=bool(rr) and all(r['honest']['ok'] for r in rr),
               flybys_by_1461d=fb['n'], mean_fuel_by_1461d=fb['mean_fuel'], fuel_by_1461d=fb['fuel'],
               launched_by_1461d=fb['launched'], dropped_routes=sorted(drop), run=None, slices=None, closure=None,
               tank_note=tank_note,
               measured='every number from one integration per twin on disk; the closure price is an ESTIMATE')
    if state is not None:
        prm = state.get('params', {})
        rep['run'] = dict(run_dir=state.get('_run_dir'), source=state.get('_source'), N=int(state['N']), T_d=round(C.s2d(state['T']), 1),
                          slice=int(state['slice']), phase=state.get('phase'), done=bool(state.get('done')), stop=state.get('stop'),
                          commit=prm.get('commit', 'regrid'), honest_commit=prm.get('commit', 'regrid') == 'regrid',
                          launched=sum(1 for c in state['craft'] if c['launched']), remaining=len(state['remaining']),
                          params={k: prm.get(k) for k in ('H', 'L', 'B', 'tries', 'beta', 'gamma', 'lam', 'qf', 'qm', 'm0max', 'proposer', 'commit')},
                          invalid=state.get('_invalid'))
        rep['slices'] = [dict(slice=l.get('slice'), T_d=l.get('T_d'), launched=l.get('launched'), covered=l.get('covered'),
                              fuel=l.get('fuel'), sumJi=l.get('sumJi'), n_columns=l.get('n_columns'), wall_s=l.get('wall_s'))
                         for l in state.get('log', [])]
        sf = _state_fuel_by(state, GATES['day_check'])
        rep['fuel_by_1461d_state'] = sf
        rep['mean_fuel_by_1461d_state'] = round(float(np.mean(list(sf.values()))), 2) if sf else None
    if g2:
        rep['closure'] = closure_estimate(routes, misses, nproc=nproc, dmax=dmax)
    rep['gates'] = gates(rep)
    rep['secs'] = round(time.time() - tic, 1)
    return rep


def gates(rep):
    """{'G0b': verdict, 'G1': verdict, 'G2': verdict, 'notes': [...]} from the REPORT numbers and GATES; a gate whose
    inputs are missing (e.g. no state.json for G0b's mean fuel) is 'N/A' with the reason in notes."""
    G = GATES; notes = []; run = rep.get('run')
    cov = int(rep['covered']); fuel = float(rep['fuel']); sumJi = float(rep['sumJi']); n = int(rep['n_craft'])
    # G0b
    fb = rep.get('flybys_by_1461d'); mf = rep.get('mean_fuel_by_1461d')
    if fb is None or mf is None:
        would = 'N/A'
    elif fb <= G['g0b_kill']:
        would = 'KILL'
    elif fb >= G['g0b_flybys'] and mf <= G['g0b_fuel']:
        would = 'PASS'
    else:
        would = 'WEAK'
    g0b_txt = f'flybys_by_1461d {fb} (>= {G["g0b_flybys"]} PASS, <= {G["g0b_kill"]} KILL), mean fuel {mf} kg (<= {G["g0b_fuel"]:.0f})'
    if run is None:
        g0b = 'N/A'; notes.append(f'G0b N/A: no state.json (not an RHFA run); {g0b_txt} would be {would}')
    elif run['T_d'] < G['day_check'] and not run['done']:
        g0b = 'N/A'; notes.append(f'G0b N/A: run has not reached day {G["day_check"]:.0f} (T = {run["T_d"]} d); {g0b_txt} so far')
    else:
        g0b = would; notes.append(f'G0b {g0b}: {g0b_txt}' + ('' if run['N'] == 8 else f' [gate defined for N = 8, this run has N = {run["N"]}]'))
        if not run.get('honest_commit', True):
            notes.append('G0b: tanks are NOT honest (commit none)')
    # G1
    if cov >= G['g1_pass'] and fuel <= G['g1_fuel']:
        w1 = 'PASS'
    elif cov >= G['g1_round']:
        w1 = 'ROUND'
    else:
        w1 = 'KILL'
    g1_txt = f'covered {cov} (>= {G["g1_pass"]} PASS, {G["g1_round"]}-{G["g1_pass"] - 1} ROUND, < {G["g1_round"]} KILL), fuel {fuel:.0f} kg (<= {G["g1_fuel"]:.0f})'
    if run is None:
        g1 = 'N/A'; notes.append(f'G1 N/A: no state.json (not an RHFA run); {g1_txt} would be {w1}')
    elif not run['done']:
        g1 = 'N/A'; notes.append(f'G1 N/A: run not finished (slice {run["slice"]}, T = {run["T_d"]} d); {g1_txt} so far ({w1})')
    else:
        g1 = w1; notes.append(f'G1 {g1}: {g1_txt}' + ('' if run['N'] == 8 else f' [gate defined for N = 8, this run has N = {run["N"]}]'))
        if cov >= G['g1_pass'] and fuel > G['g1_fuel']:
            notes.append(f'G1: coverage passes but fuel {fuel:.0f} > {G["g1_fuel"]:.0f} kg -> ROUND (fuel)')
    # G2
    cl = rep.get('closure'); n_miss = C.N_TARGETS - cov
    fallback = n == 9 and sumJi + n_miss <= G['g2_fallback9'] + 1e-9
    if cov >= G['g2_cover'] and cl is not None:
        est = sumJi + float(cl['dJ_total'])
        if est <= G['g2_sumJi'] + 1e-9:
            g2 = 'EXPORT'
        elif fallback:
            g2 = 'EXPORT9'
        else:
            g2 = 'NO'
        notes.append(f'G2 {g2}: covered {cov} (>= {G["g2_cover"]}), sumJi {sumJi:.4f} + closure ESTIMATE {cl["dJ_total"]:.4f} '
                     f'(lin {cl["dJ_lin"]:.4f}, calibrated {cl["dJ_cal"]:.4f}, {len(cl["placed"])} placed / {len(cl["unplaced"])} unplaced of '
                     f'{n_miss} misses) = {est:.4f} vs {G["g2_sumJi"]:.2f}')
    elif cov >= G['g2_cover']:
        g2 = 'EXPORT9' if fallback else 'N/A'
        notes.append(f'G2 {g2}: covered {cov}, sumJi {sumJi:.4f}; closure estimate not computed (--no-g2)')
    else:
        g2 = 'EXPORT9' if fallback else 'NO'
        notes.append(f'G2 {g2}: covered {cov} < {G["g2_cover"]} (sumJi {sumJi:.4f}, {n_miss} misses)')
    if n == 9:
        notes.append(f'G2 fallback (9 craft): sumJi + misses = {sumJi + n_miss:.4f} vs {G["g2_fallback9"]:.2f} -> ' + ('EXPORT9' if fallback else 'no'))
    if not rep.get('honest_all', False):
        notes.append('WARNING: not every route passes C.honest_check (see the honest column)')
    if rep.get('duplicates'):
        notes.append(f'WARNING: targets flown twice: {rep["duplicates"]}')
    return dict(G0b=g0b, G1=g1, G2=g2, notes=notes)


def _short(lst, n=40):
    """A list for the table: the first n ids, then '... (+k more)' (the JSON keeps the whole list)."""
    lst = list(lst)
    return str(lst) if len(lst) <= n else f'{lst[:n]} ... (+{len(lst) - n} more)'


def table(rep):
    """Printable table: one line per route (name n tank J_i launch_d end_d cadence kg/fb pin/med/fil honest) + the
    totals + gate lines."""
    L = [f'fleet {rep["path"]}' + (f'  (dropped {rep["dropped_routes"]})' if rep.get('dropped_routes') else ''),
         f'{"route":8s} {"n":>3s} {"tank":>8s} {"J_i":>7s} {"launch_d":>9s} {"end_d":>8s} {"cad_d":>6s} {"kg/fb":>6s} '
         f'{"pin/med/fil":>11s} {"irreg":>5s} {"winCap":>6s} {"miss":>6s} honest']
    for r in rep['routes']:
        h = r['honest']; c = r['counts']
        L.append(f'{r["name"]:8s} {r["n"]:3d} {r["tank"]:8.1f} {r["J_i"]:7.4f} {r["t_launch_d"]:9.1f} {r["t_end_d"]:8.1f} '
                 f'{(r["cadence_d"] if r["cadence_d"] is not None else 0):6.1f} {(r["kg_per_fb"] or 0):6.2f} '
                 f'{c["pin"]:3d}/{c["med"]:3d}/{c["fil"]:3d} {h["irregular"]:5d} {h["max_window_cap"]:6.3f} {h["miss_km"]:6.1f} '
                 f'{"ok" if h["ok"] else "FAIL"}')
    L.append(f'covered {rep["covered"]}  misses {len(rep["misses"])} {rep["misses_by_class"]}  n_craft {rep["n_craft"]}  '
             f'sumJi {rep["sumJi"]:.4f}  fuel {rep["fuel"]:.1f} kg  J_raw {rep["J_raw"]:.4f}  honest_all {rep["honest_all"]}'
             + ('' if rep.get('tank_note', '').startswith('honest') else f'  [{rep.get("tank_note")}]'))
    if rep['misses']:
        L.append(f'misses: {_short(rep["misses"])}')
    L.append(f'flybys_by_1461d {rep["flybys_by_1461d"]}  mean_fuel_by_1461d {rep["mean_fuel_by_1461d"]} kg over '
             f'{len(rep["launched_by_1461d"])} launched (MEASURED from the impulses before day 1461)'
             + (f'; from state history {rep["mean_fuel_by_1461d_state"]} kg' if rep.get('mean_fuel_by_1461d_state') is not None else ''))
    if rep.get('run'):
        r = rep['run']
        L.append(f'run {r["run_dir"]}: N {r["N"]} slice {r["slice"]} phase {r["phase"]} T {r["T_d"]} d launched {r["launched"]} '
                 f'remaining {r["remaining"]} done {r["done"]} commit {r["commit"]} source {r["source"]}'
                 + (f'  INVALID STATE: {r["invalid"]}' if r.get('invalid') else ''))
        if rep.get('slices'):
            L.append('slices: ' + '  '.join(f's{s["slice"]}@{s["T_d"]}d cov {s["covered"]} fuel {s["fuel"]}' for s in rep['slices']))
    cl = rep.get('closure')
    if cl is not None:
        L.append(f'G2 closure ESTIMATE: dJ_lin {cl["dJ_lin"]:.4f} + {len(cl["unplaced"])} unplaced = dJ_total {cl["dJ_total"]:.4f} '
                 f'(calibrated {cl["dJ_cal"]:.4f}; screen = every miss at its cheapest host {cl.get("dJ_screen", 0.0):.4f}, '
                 f'{cl.get("n_without_host", 0)} without a host within {cl["dmax"]} AU); sumJi + closure = {rep["sumJi"] + cl["dJ_total"]:.4f}; '
                 f'{len(cl["placed"])} placed, {cl["n_candidates"]} priced approaches, {cl["secs"]} s')
        L.append(f'  {cl["note"]}')
        for p in cl['placed']:
            L.append(f'  ESTIMATE place {p["ast"]:3d} -> {p["host"]:6s} lin dJ {p["lin_dJ"]:.4f} ({p["lin_kg"]:.1f} kg) at {p["t_d"]:.0f} d, {p["dist_au"]:.3f} AU')
        if cl['unplaced']:
            L.append(f'  ESTIMATE unplaced (1.0 each): {_short(cl["unplaced"])}')
    g = rep['gates']
    L.append(f'GATES  G0b {g["G0b"]}  G1 {g["G1"]}  G2 {g["G2"]}')
    for n in g['notes']:
        L.append('  ' + n)
    return '\n'.join(L)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('path'); ap.add_argument('--json', default=None)
    ap.add_argument('--no-g2', action='store_true'); ap.add_argument('--nproc', type=int, default=2)
    ap.add_argument('--drop', default=''); ap.add_argument('--dmax', type=float, default=0.2)
    a = ap.parse_args()
    tic = time.time()
    drop = [x for x in a.drop.split(',') if x]
    rep = report(a.path, g2=not a.no_g2, nproc=a.nproc, drop=drop, dmax=a.dmax)
    print(table(rep))
    out = a.json
    if out is None and rep.get('run') and not drop:
        out = str(pathlib.Path(rep['run']['run_dir']) / 'eval.json')
    if out:
        C.jdump(rep, out); print('written', out)
    print(f'[{C.now()}] eval wall {time.time() - tic:.1f} s (MEASURED)')


if __name__ == '__main__':
    main()
