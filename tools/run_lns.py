"""Large-neighbourhood improvement of a PLANNED fleet with the planner itself (phase 2 after column generation).

Moves (all priced with the colgen cost model, fleet J = sum_i c_i + misses + 2, accepted only if J decreases):
  replan(i, k)  keep route i up to flyby k (k = -1: nothing, full re-plan with launch epoch near the old one), then run
                the prize-collecting beam search from the spacecraft state after flyby k: prize 1 (one miss) on every
                target that no OTHER route covers, targets of other routes excluded (no duplicates). Every beam state of
                every depth is a candidate route; the one minimising the fleet J wins.
  spawn         a new craft over the uncovered targets only (launch anywhere in the first 2000 d), kept iff K > c.
  dissolve(i)   remove route i, re-plan the tails of the other routes (split 1/2) to absorb its targets; the chain of
                re-plans is accepted as a whole iff the fleet J decreases.
The state at flyby k is reconstructed with the planner's hybrid leg model (tools/exp_retime.py: identify_models + Chain),
so prefixes are exactly the planned chains at the search mass.

Usage: run_lns.py fleetdir outdir [--rounds 2] [--splits 0.67,0.5,0.33,-1] [--beam 300] [--nproc 5] [--spawn] [--dissolve]
Writes outdir/tour_sc*.json (current fleet, colgen format) after every accepted move, outdir/lns.log, outdir/summary.json.
"""
import sys, json, time, pathlib, argparse, shutil
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
from ctoc14.kepler import Ephemeris
from ctoc14.constants import DAY, AU, M0_MAX
from ctoc14.search import Params, beam_search, beam_search_from_state, tour_json, UNREACHABLE
from ctoc14.colgen import CostModel, TARGETS
from exp_retime import identify_models, Chain
from ctoc14.jointsearch import joint_beam_search

ap = argparse.ArgumentParser()
ap.add_argument('fleetdir'); ap.add_argument('outdir')
ap.add_argument('--rounds', type=int, default=2); ap.add_argument('--splits', default='0.67,0.5,0.33,-1')
ap.add_argument('--beam', type=int, default=300); ap.add_argument('--nproc', type=int, default=5)
ap.add_argument('--wt', type=float, default=0.9); ap.add_argument('--no-refine', action='store_true')
ap.add_argument('--cost-s', type=float, default=1.15); ap.add_argument('--reserve', type=float, default=20.0)
ap.add_argument('--dvmax', type=float, default=1.2); ap.add_argument('--lin-tofmax', type=float, default=600.0)
ap.add_argument('--spawn', action='store_true'); ap.add_argument('--dissolve', action='store_true')
ap.add_argument('--prize-covered', type=float, default=0.0, help='prize for targets other routes already cover (0 = excluded)')
ap.add_argument('--pairs', type=int, default=0, help='joint re-plans of route pairs per round (weakest pairs first)')
ap.add_argument('--pair-beam', type=int, default=200)
ap.add_argument('--dissolve-k', type=int, default=1, help='dissolve attempts per round (weakest routes first)')
ap.add_argument('--dissolve-split', type=float, default=0.5, help='split of the absorbing re-plans (-1 = full re-plan)')
ap.add_argument('--m0-up', type=float, default=0.0, help='re-plan routes at max(own search mass, this) (bigger tank allowed)')
ap.add_argument('--group', type=int, default=0, help='group move: remove the k weakest routes, joint re-plan with k-1 and k craft')
ap.add_argument('--group-launch-max', type=float, default=2000.0)
a = ap.parse_args()
out = pathlib.Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
logf = open(out / 'lns.log', 'a')
def say(s):
    line = f'[{time.strftime("%H:%M:%S")}] {s}'; print(line, flush=True); logf.write(line + '\n'); logf.flush()

cm = CostModel(s=a.cost_s, reserve=a.reserve)
eph = Ephemeris()
TSET = set(TARGETS.tolist())
def tset(t): return {int(l['ast']) for l in t['legs']} - UNREACHABLE
def m0s(t): return float(t.get('m0_search', t['m0']))
def rcost(t): return float(cm.cost(t['fuel_est'], m0s(t)))
def fleet_J(tours):
    cov = set().union(*[tset(t) for t in tours]) if tours else set()
    return sum(rcost(t) for t in tours) + len(TSET - cov) + 2, cov

def save(tours, tag=''):
    for f in out.glob('tour_sc*.json'): f.unlink()
    for k, t in enumerate(sorted(tours, key=lambda t: -len(t['legs'])), 1):
        t = dict(t); t['m0'] = float(min(M0_MAX, cm.tank(t['fuel_est'], m0s(t)))); t['cost_J'] = rcost(t)
        json.dump(t, open(out / f'tour_sc{k}.json', 'w'), indent=1)
    J, cov = fleet_J(tours)
    json.dump(dict(J=J, n=len(tours), covered=len(cov), sumJ=sum(rcost(t) for t in tours), missed=sorted(TSET - cov), tag=tag),
              open(out / 'summary.json', 'w'), indent=1)

def params(prize, m0):
    return Params(beam=a.beam, w_fuel=cm.w_fuel(360.0 * m0 / 1000.0, m0), w_t=a.wt, m_margin=40.0, dv_max=a.dvmax,
                  lin_tofs=np.arange(20, a.lin_tofmax + 1e-9, 10) * DAY, lin_drmax=0.15 * AU, vinf_cap=2.0,
                  tof_refine=not a.no_refine, prize=prize)

def state_after(tour, k):
    """Planner state after flyby k (0-based) of `tour` at its search mass."""
    asts = [l['ast'] for l in tour['legs']]; times = [l['t_flyby'] for l in tour['legs']]
    typ, rep = identify_models(eph, tour, m0s(tour))
    C = Chain(eph, 'x', tour['t_launch'], m0s(tour), asts, times, typ)
    m_after = float(C.m_start[k] - C.dm[k])
    seq = tuple((int(l['ast']), float(l['t_flyby']), float(l['dv_est']), float(l['tof'])) for l in tour['legs'][:k + 1])
    return dict(t=float(times[k]), r=C.legs[k]['r'], v=C.legs[k]['v_arr'], m=m_after, vinf=np.asarray(tour['vinf'], float),
                seq=seq, max_ddv=rep['max_ddv'], fuel_prefix=float(m0s(tour) - m_after))

def best_candidate(states, others_cost, others_cov, m0):
    """Candidate route (state) minimising the fleet J with all other routes fixed."""
    best = None
    for s in states:
        ids = {x[0] for x in s.seq} - UNREACHABLE
        c = float(cm.cost(s.fuel, m0))
        J = others_cost + c + len(TSET - (others_cov | ids)) + 2
        if best is None or J < best[0]:
            best = (J, s)
    return best

def replan(tours, i, split, prize_extra=None):
    """Returns (J_new, new_tour or None, info str)."""
    t = tours[i]; m0 = max(m0s(t), a.m0_up) if split < 0 else m0s(t)      # a larger tank only for full re-plans
    others = [tt for j, tt in enumerate(tours) if j != i]
    ocov = set().union(*[tset(tt) for tt in others]) if others else set(); ocost = sum(rcost(tt) for tt in others)
    prize = np.zeros(300)
    for x in TSET - ocov: prize[x - 1] = 1.0
    if a.prize_covered > 0:
        for x in ocov: prize[x - 1] = a.prize_covered
    if prize_extra:
        for x, v in prize_extra.items(): prize[x - 1] = max(prize[x - 1], v)
    P = params(prize, m0); P.collect = []
    n = len(t['legs']); tic = time.time()
    if split < 0:
        tL = float(t['t_launch']); grid = np.arange(max(0.0, tL - 200 * DAY), tL + 400 * DAY + 1, 20 * DAY)
        excl = sorted(ocov) if a.prize_covered <= 0 else []
        bs, beam = beam_search(eph, excluded=excl, m0=m0, t_launch_grid=grid, P=P, n_proc=a.nproc, verbose=False)
        k = -1
    else:
        k = int(round(split * (n - 1))); k = max(0, min(k, n - 2))
        st = state_after(t, k)
        vis = (ocov if a.prize_covered <= 0 else set()) | {x[0] for x in st['seq']}
        bs, beam = beam_search_from_state(eph, st['t'], st['r'], st['v'], st['m'], sorted(vis), m0, float(t['t_launch']), st['vinf'],
                                          P=P, n_proc=a.nproc, verbose=False, seq=st['seq'])
    states = P.collect + list(beam or [])
    if not states:
        return None, None, f'no states ({time.time()-tic:.0f} s)'
    Jn, s = best_candidate(states, ocost, ocov, m0)
    nt = tour_json(s); nt['m0_search'] = m0; nt['tag'] = f'lns:{t.get("tag", "")}'
    return Jn, nt, f'k={k} ({time.time()-tic:.0f} s, {len(states)} states)'

def replan_group(tours, group, K, launch_max):
    """Remove the routes in `group`, joint re-plan K craft (shared visited mask) over everything no other route covers."""
    others = [tt for k, tt in enumerate(tours) if k not in group]
    ocov = set().union(*[tset(tt) for tt in others]) if others else set(); ocost = sum(rcost(tt) for tt in others)
    prize = np.zeros(300)
    for x in TSET - ocov: prize[x - 1] = 1.0
    if a.prize_covered > 0:
        for x in ocov: prize[x - 1] = a.prize_covered
    m0 = 1000.0
    P = params(prize, m0); P.beam = a.pair_beam; P.w_rare = 1.0; P.wait = 60 * DAY; P.collect_joint = []
    tic = time.time()
    best = joint_beam_search(eph, [m0] * K, np.arange(0, launch_max, 20) * DAY, P=P, excluded=(sorted(ocov) if a.prize_covered <= 0 else []),
                             n_proc=a.nproc, verbose=False)
    cands = P.collect_joint + ([best] if best is not None else [])
    bestJ = None
    for js in cands:
        rs = [s for s in js.craft if s is not None and s.n() > 0]
        ids = set().union(*[{x[0] for x in s.seq} for s in rs]) - UNREACHABLE if rs else set()
        J = ocost + sum(float(cm.cost(s.fuel, s.m0)) for s in rs) + len(TSET - (ocov | ids)) + 2
        if bestJ is None or J < bestJ[0]:
            bestJ = (J, rs)
    if bestJ is None:
        return None, None, f'no joint states ({time.time()-tic:.0f} s)'
    new = []
    for s in bestJ[1]:
        nt = tour_json(s); nt['m0_search'] = float(s.m0); nt['tag'] = f'lns-group{K}'; new.append(nt)
    return bestJ[0], new, f'({time.time()-tic:.0f} s, {len(cands)} joint states)'

def replan_pair(tours, i, j):
    """Joint re-plan of routes i and j from launch (shared visited mask), other routes' targets excluded, prize 1 on
    every target no other route covers. Every joint beam state (0, 1 or 2 non-empty routes) is a candidate."""
    others = [tt for k, tt in enumerate(tours) if k not in (i, j)]
    ocov = set().union(*[tset(tt) for tt in others]) if others else set(); ocost = sum(rcost(tt) for tt in others)
    prize = np.zeros(300)
    for x in TSET - ocov: prize[x - 1] = 1.0
    if a.prize_covered > 0:
        for x in ocov: prize[x - 1] = a.prize_covered
    m0a, m0b = m0s(tours[i]), m0s(tours[j])
    P = params(prize, max(m0a, m0b)); P.beam = a.pair_beam; P.w_rare = 1.0; P.wait = 60 * DAY; P.collect_joint = []
    tic = time.time()
    best = joint_beam_search(eph, [m0a, m0b], np.arange(0, 720, 20) * DAY, P=P, excluded=(sorted(ocov) if a.prize_covered <= 0 else []),
                             n_proc=a.nproc, verbose=False)
    cands = P.collect_joint + ([best] if best is not None else [])
    bestJ = None
    for js in cands:
        rs = [s for s in js.craft if s is not None and s.n() > 0]
        ids = set().union(*[{x[0] for x in s.seq} for s in rs]) - UNREACHABLE if rs else set()
        c = sum(float(cm.cost(s.fuel, s.m0)) for s in rs)
        J = ocost + c + len(TSET - (ocov | ids)) + 2
        if bestJ is None or J < bestJ[0]:
            bestJ = (J, rs)
    if bestJ is None:
        return None, None, f'no joint states ({time.time()-tic:.0f} s)'
    new = []
    for s in bestJ[1]:
        nt = tour_json(s); nt['m0_search'] = float(s.m0); nt['tag'] = 'lns-pair'; new.append(nt)
    return bestJ[0], new, f'({time.time()-tic:.0f} s, {len(cands)} joint states)'

# ------------------------------------------------------------------------------------------------ main loop
src = pathlib.Path(a.fleetdir)
tours = [json.load(open(f)) for f in sorted(src.glob('tour_sc*.json'), key=lambda f: int(f.stem[7:]))]
for t in tours: t.setdefault('m0_search', 1000.0)
J0, cov = fleet_J(tours)
say(f'start {src}: {len(tours)} craft, covered {len(cov)}, sum J_i {sum(rcost(t) for t in tours):.3f}, planned J {J0:.3f}; missed {sorted(TSET - cov)}')
save(tours, 'start')
splits = [float(x) for x in a.splits.split(',')]
for rnd in range(a.rounds):
    improved = False
    # routes with the fewest unique targets first (most likely to profit)
    def uniq(i):
        o = set().union(*[tset(tt) for j, tt in enumerate(tours) if j != i]) if len(tours) > 1 else set()
        return len(tset(tours[i]) - o)
    order = sorted(range(len(tours)), key=uniq)
    for i in order:
        for sp in splits:
            Jn, nt, info = replan(tours, i, sp)
            if nt is None:
                say(f'round {rnd+1} route {i} split {sp}: {info}'); continue
            old_n = len(tours[i]['legs']); tag = 'ACCEPT' if Jn < J0 - 1e-6 else 'reject'
            say(f'round {rnd+1} route {i} split {sp}: {old_n} -> {len(nt["legs"])} flybys, fuel {tours[i]["fuel_est"]:.0f} -> {nt["fuel_est"]:.0f} kg, '
                f'J {J0:.3f} -> {Jn:.3f} [{tag}] {info}')
            if Jn < J0 - 1e-6:
                tours[i] = nt; J0 = Jn; improved = True; save(tours, f'r{rnd+1}')
    if a.pairs > 0 and len(tours) >= 2:
        U = [uniq(i) for i in range(len(tours))]
        pairs = sorted(((U[i] + U[j], i, j) for i in range(len(tours)) for j in range(i + 1, len(tours))))[:a.pairs]
        for _, i, j in pairs:
            if max(i, j) >= len(tours): continue
            Jn, new, info = replan_pair(tours, i, j)
            if new is None:
                say(f'round {rnd+1} pair ({i},{j}): {info}'); continue
            tag = 'ACCEPT' if Jn < J0 - 1e-6 else 'reject'
            say(f'round {rnd+1} pair ({i},{j}): {len(tours[i]["legs"])}+{len(tours[j]["legs"])} -> {"+".join(str(len(t["legs"])) for t in new)} flybys, '
                f'J {J0:.3f} -> {Jn:.3f} [{tag}] {info}')
            if Jn < J0 - 1e-6:
                rest = [t for k, t in enumerate(tours) if k not in (i, j)]
                tours = rest + new; J0 = Jn; improved = True; save(tours, f'r{rnd+1}-pair')
                break                                   # indices changed; next pairs next round
    if a.group >= 2 and len(tours) > a.group:
        U = [uniq(i) for i in range(len(tours))]
        grp = sorted(range(len(tours)), key=lambda i: U[i])[:a.group]
        for K in (a.group - 1, a.group):
            Jn, new, info = replan_group(tours, grp, K, a.group_launch_max)
            if new is None:
                say(f'round {rnd+1} group {grp} K={K}: {info}'); continue
            tag = 'ACCEPT' if Jn < J0 - 1e-6 else 'reject'
            say(f'round {rnd+1} group {grp} K={K}: {"+".join(str(len(tours[i]["legs"])) for i in grp)} -> '
                f'{"+".join(str(len(t["legs"])) for t in new)} flybys, J {J0:.3f} -> {Jn:.3f} [{tag}] {info}')
            if Jn < J0 - 1e-6:
                tours = [t for k, t in enumerate(tours) if k not in grp] + new; J0 = Jn; improved = True; save(tours, f'r{rnd+1}-group')
                break
    n_before_spawn = len(tours)
    if a.spawn:
        J0, cov = fleet_J(tours); unc = sorted(TSET - cov)
        if len(unc) >= 2:
            prize = np.zeros(300)
            if a.prize_covered > 0:
                for x in cov: prize[x - 1] = a.prize_covered
            prize[np.array(unc) - 1] = 1.0
            P = params(prize, 1000.0); P.collect = []; tic = time.time()
            bs, beam = beam_search(eph, excluded=(sorted(cov) if a.prize_covered <= 0 else []), m0=1000.0,
                                   t_launch_grid=np.arange(0, 2000, 20) * DAY, P=P, n_proc=a.nproc, verbose=False)
            states = P.collect + list(beam or [])
            if states:
                Jn, s = best_candidate(states, sum(rcost(t) for t in tours), cov, 1000.0)
                nt = tour_json(s); nt['m0_search'] = 1000.0; nt['tag'] = 'spawn'
                tag = 'ACCEPT' if Jn < J0 - 1e-6 else 'reject'
                say(f'round {rnd+1} spawn over {len(unc)} uncovered: {len(nt["legs"])} flybys, fuel {nt["fuel_est"]:.0f} kg, J {J0:.3f} -> {Jn:.3f} [{tag}] ({time.time()-tic:.0f} s)')
                if Jn < J0 - 1e-6:
                    tours.append(nt); J0 = Jn; improved = True; save(tours, f'r{rnd+1}-spawn')
    if a.dissolve and len(tours) > 2:
        # try dissolving the routes with the fewest unique targets (never a craft spawned in this round)
        cand_i = sorted(range(min(n_before_spawn, len(tours))), key=uniq)[:a.dissolve_k]
        for i in cand_i:
            if i >= len(tours): continue
            trial = [t for j, t in enumerate(tours) if j != i]; Jt, _ = fleet_J(trial)
            say(f'round {rnd+1} dissolve route {i} ({len(tours[i]["legs"])} flybys, {uniq(i)} unique): J {J0:.3f} -> {Jt:.3f} before re-plans')
            cov_i = tset(tours[i])
            # absorbing routes: those whose time span overlaps most with the dissolved route's flybys first
            ti = [l['t_flyby'] for l in tours[i]['legs']]
            def near(j):
                tj = [l['t_flyby'] for l in trial[j]['legs']]
                return -sum(1 for x in ti if min(abs(x - y) for y in tj) < 120 * DAY)
            for j in sorted(range(len(trial)), key=near):
                if Jt < J0 - 1e-6: break                       # already better than before: stop re-planning
                Jn, nt, info = replan(trial, j, a.dissolve_split)
                if nt is not None and Jn < Jt - 1e-6:
                    trial[j] = nt; Jt = Jn
                say(f'   dissolve re-plan {j}: J -> {Jt:.3f} {info}')
            tag = 'ACCEPT' if Jt < J0 - 1e-6 else 'reject'
            say(f'round {rnd+1} dissolve route {i}: fleet J {J0:.3f} -> {Jt:.3f} [{tag}]')
            if Jt < J0 - 1e-6:
                tours = trial; J0 = Jt; improved = True; save(tours, f'r{rnd+1}-dissolve')
                break
    J0, cov = fleet_J(tours)
    say(f'end of round {rnd+1}: {len(tours)} craft, covered {len(cov)}, planned J {J0:.3f}; missed {sorted(TSET - cov)}')
    if not improved:
        break
save(tours, 'final')
say(f'FINAL: {len(tours)} craft, planned J {J0:.3f}')
