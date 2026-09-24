"""Skeleton-first fill (docs/stage5_skeleton_fill.md): every route of the skeleton fleet (its HARD flybys with epochs)
is filled with EASY targets by the windowed beam -- the skeleton's flybys are mandatory time windows (epoch +- half),
the beam ranks children with the waypoint lookahead, and a state that misses a window dies.  Routes are filled most
constrained first (most waypoints), each over the easy targets no earlier route took; for each route a small
portfolio of beam settings is run, the (flybys, fuel) Pareto tours are settled (impulsive twin) and the twin with the
best fleet value J_i - lam * flybys is kept.  Leftover targets are listed for grow_skeletons.py / run_ialns.py.

Usage: skel_fill.py skel_dir out_dir [--prizes results/n8/prizes_H.json] [--beam 1000] [--nproc 8]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, signal, pathlib, argparse, warnings
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from greedy_cover import BP, ALL, CM, DeepCollect, _timeout
from ctoc14.search import beam_search, tour_json
import ctoc14.impulsive as IM
from ctoc14.impulsive import settle_tour
from ctoc14.constants import DAY, AU
warnings.filterwarnings('ignore')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('src'); ap.add_argument('out')
    ap.add_argument('--prizes', default=str(ROOT / 'results/n8/prizes_H.json'), help='defines the hard set (prize>=1.5)')
    ap.add_argument('--m0', type=float, default=1600.0); ap.add_argument('--vinf', type=float, default=4.0)
    ap.add_argument('--dvmax', type=float, default=1.2); ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--portfolio', default='1000:0.15:45,1000:0.15:90,1000:0.20:60',
                    help='comma list of beam:lin_drmax(AU):window half-width(d)')
    ap.add_argument('--wp-prize', type=float, default=2.0); ap.add_argument('--max-depth', type=int, default=60)
    ap.add_argument('--try', dest='ntry', type=int, default=4, help='Pareto depths to settle per route (deepest first)')
    ap.add_argument('--lam', type=float, default=0.05, help='fleet value of one flyby (J) when choosing among settled twins')
    ap.add_argument('--settle-timeout', type=float, default=240.0); ap.add_argument('--iters', type=int, default=100)
    ap.add_argument('--max-tank', type=float, default=1100.0); ap.add_argument('--lin-fallback', type=float, default=2.5)
    ap.add_argument('--order', default='', help='explicit route order (comma list); default most waypoints first')
    ap.add_argument('--wp-dvmax', type=float, default=2.5, help='per-leg dv cap for legs ENDING at a waypoint (they are '
                    'mandatory; the skeleton itself may need a 1.9 km/s leg that the ordinary --dvmax forbids)')
    ap.add_argument('--retry', default='1000:0.15:120,600:0.20:150', help='portfolio used when no state passes all windows')
    ap.add_argument('--mode', default='windowed', choices=['windowed', 'segmented'],
                    help='windowed: one beam with all waypoint windows; segmented: waypoint-to-waypoint beams (fill_segmented)')
    ap.add_argument('--seg-k', type=int, default=30, help='segmented: states kept at each waypoint')
    ap.add_argument('--seg-depth', type=int, default=14, help='segmented: max flybys per segment')
    ap.add_argument('--refill', default='', help='REFILL mode: an already filled fleet dir; every route is re-filled over its '
                    'own easy targets + the fleet leftovers (prize --left-prize) and replaced when the net value improves')
    ap.add_argument('--left-prize', type=float, default=2.5); ap.add_argument('--left-lam', type=float, default=0.15,
                    help='fleet value (J) of covering one leftover (an uninsertable one costs 1 J)')
    ap.add_argument('--refill-rounds', type=int, default=2)
    ap.add_argument('--easy-prizes', default='', help='JSON list of 300 prizes for the easy targets (dual loop); default 1.0')
    ap.add_argument('--assign', action='store_true', help='pre-assign every easy target to the skeleton whose trajectory '
                    'passes closest (balanced, capacity --assign-cap x mean); each route is filled over its group + a pad')
    ap.add_argument('--assign-cap', type=float, default=1.25); ap.add_argument('--pad', type=int, default=24)
    ap.add_argument('--columns', default='', help='COLUMN mode: write diverse Pareto tours per skeleton to this jsonl '
                    '(no settling); --ncol columns per skeleton, each built with the earlier columns\' easy targets at --used-prize')
    ap.add_argument('--ncol', type=int, default=5); ap.add_argument('--used-prize', type=float, default=0.3)
    ap.add_argument('--col-min', type=int, default=24, help='keep Pareto tours with at least this many flybys')
    ap.add_argument('--pad-prize', type=float, default=0.5); ap.add_argument('--assign-dmax', type=float, default=2.0)
    ap.add_argument('--covered-prize', type=float, default=0.0, help='STEPPING STONES (stage 6 L1): keep the easy targets '
                    'already covered by earlier routes in the pool at this prize instead of excluding them -- every leg must '
                    'end at an asteroid, so a depleted pool has no relays.  Route value then counts NEW targets only.')
    ap.add_argument('--start', default='', help='an existing filled fleet: routes not in --only are kept as they are and '
                    'their targets are stepping stones for the refilled ones')
    ap.add_argument('--only', default='', help='comma list of route names to (re)fill; default all')
    ap.add_argument('--pool-frac', type=float, default=1.0, help='keep this fraction of the UNCOVERED easy pool at random '
                    '(waypoints and covered stepping stones are always kept): decorrelates the columns of stage 6 B')
    ap.add_argument('--pool-seed', type=int, default=0)
    ap.add_argument('--snapshots', action='store_true', help='save the fleet after each route to <out>/snap_<k>, so a '
                    'partition search can restart the sequential fill from the first route it changed')
    ap.add_argument('--drop-wp', default='', help='DIAGNOSTIC: waypoint ids (or "all") removed from the skeleton for this '
                    'run; a dropped hard target stays in the pool as an ordinary prize-1 target, so the fill may still '
                    'pick it up.  The depth difference is that waypoint\'s capacity cost on this route.')
    ap.add_argument('--pool-fracs', default='1.0,0.7,0.5,0.35', help='column mode: one fill per (fraction, seed)')
    ap.add_argument('--pool-seeds', type=int, default=6, help='column mode: random pools per fraction (< 1.0)')
    ap.add_argument('--pools', default='', help='JSON {route: [easy ids]} giving each route its own pool (stage 6 C): the '
                    'assignment comes from the MEASURED reach matrix of the column library, not from geometry.  Targets '
                    'outside a route\'s pool are still available at --pad-prize so the route is never pool-starved.')
    a = ap.parse_args()
    IM.LIN_FALLBACK = a.lin_fallback if a.lin_fallback > 0 else None
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt')
    RI.eph()
    H = np.asarray(json.load(open(a.prizes)), float)
    hard = set(t for t in ALL if H[t - 1] >= 1.5)
    skel = RI.IFleet(a.src)
    say(f'skeletons: {skel.summary()}')
    if a.drop_wp:                     # diagnostic: the dropped hard targets fall back into the ordinary easy pool
        drop = set(range(1, 301)) if a.drop_wp == 'all' else set(int(x) for x in a.drop_wp.split(',') if x)
        for n, r in skel.routes.items():
            st = r['st']; asts = [int(x) for x in st['asts']]
            keep = [i for i, x in enumerate(asts) if x not in drop]
            if len(keep) == len(asts): continue
            say(f'  {n}: dropping waypoints {[x for x in asts if x in drop]} -> {len(keep)} left')
            st['asts'] = np.asarray(asts, float)[keep]; st['tf'] = np.asarray(st['tf'], float)[keep]
    EP = np.ones(300)
    if a.easy_prizes:
        EP = np.asarray(json.load(open(a.easy_prizes)), float)
        say(f'easy prizes from {a.easy_prizes}: {int((EP > 1.0).sum())} targets above 1.0, max {EP.max():.2f}')
    wf = CM.w_fuel(0.5 * (a.m0 - 640.0), a.m0); grid = np.arange(0, 800 + 1e-9, 20) * DAY
    port = [tuple(float(x) for x in p.split(':')) for p in a.portfolio.split(',')]
    fleet = RI.IFleet()
    used = set(); hard_all = set()
    for n, r in skel.routes.items():
        hard_all |= set(int(t) for t in r['st']['asts'])
    used |= hard_all
    order = a.order.split(',') if a.order else sorted(skel.routes, key=lambda n: -len(skel.routes[n]['st']['asts']))
    if a.start:
        fleet = RI.IFleet(a.start); used |= set(fleet.coverage())
        say(f'start from {a.start}: {fleet.summary()}')
    if a.only:
        order = [n for n in a.only.split(',')]
    rng = np.random.default_rng(a.pool_seed)
    POOLS = json.load(open(a.pools)) if a.pools else None
    grp = pad = None
    if a.assign:
        import multiprocessing as mp
        easy0 = [t for t in ALL if t not in used]; keys = sorted(skel.routes)
        with mp.get_context('fork').Pool(a.nproc) as pool:
            res = pool.map(RI.w_cands, [(n, skel.routes[n]['st'], easy0, a.assign_dmax) for n in keys])
        D = {n: {} for n in keys}
        for n, lst in zip(keys, res):
            for c in lst: D[n][c['ast']] = min(D[n].get(c['ast'], 9e9), c['dist'])
        cap = int(np.ceil(a.assign_cap * len(easy0) / len(keys)))
        grp = {n: set() for n in keys}; done = set()
        for d, t, n in sorted((D[n].get(t, 9.0), t, n) for t in easy0 for n in keys):
            if t in done or len(grp[n]) >= cap: continue
            grp[n].add(t); done.add(t)
        pad = {n: [t for t in sorted(D[n], key=lambda t: D[n][t]) if t not in grp[n]][:a.pad] for n in keys}
        say('assignment (cap %d): ' % cap + ' '.join(f'{n}:{len(grp[n])}' for n in keys) +
            '; mean closest approach ' + ' '.join(f'{n}:{np.mean([D[n].get(t, 9.0) for t in grp[n]]):.2f}' for n in keys) + ' AU')
    if a.refill:
        fleet = RI.IFleet(a.refill); say(f'refill from {a.refill}: {fleet.summary()}')
        order = order * a.refill_rounds
    if a.columns:
        return gen_columns(a, skel, order, used, EP, wf, grid, say)
    for name in order:
        st = skel.routes[name]['st']; wps = sorted((float(t), int(x)) for x, t in zip(st['asts'], st['tf']))
        wpset = set(x for _, x in wps)
        others = set()
        for m, r in fleet.routes.items():
            if m != name: others |= set(int(x) for x in r['st']['asts'])
        own_now = set(int(x) for x in fleet.routes[name]['st']['asts']) if name in fleet.routes else set()
        own_old = len(own_now - others); covered = others - hard_all
        if a.refill:
            cov = fleet.coverage(); own = set(int(x) for x in fleet.routes[name]['st']['asts']) - wpset
            left = [t for t in ALL if t not in cov]
            easy = sorted(own | set(left)); c_old = RI.cost(fleet.routes[name]['tank'])
            say(f'--- refill {name}: {len(own)} own easy + {len(left)} leftovers {left}')
        elif a.assign:
            easy = sorted((grp[name] | set(pad[name])) - used)
        elif POOLS is not None:
            mine = set(POOLS.get(name, [])) - used
            easy = sorted(mine | set(t for t in ALL if t not in used and t not in hard_all))
        elif a.covered_prize > 0.0:
            # stepping stones: the whole easy pool stays available; what the OTHER routes cover is worth ~nothing but
            # is still a legal leg end.  A route being REFILLED releases its own targets (they are new again for it).
            easy = [t for t in ALL if t not in hard_all]
            if a.pool_frac < 1.0:
                fresh = [t for t in easy if t not in covered]
                keep = set(rng.choice(fresh, size=max(1, int(round(a.pool_frac * len(fresh)))), replace=False).tolist())
                easy = [t for t in easy if t in covered or t in keep]
                say(f'--- {name}: pool_frac {a.pool_frac} seed {a.pool_seed}: {len(fresh)} fresh -> {len(keep)}')
        else:
            easy = [t for t in ALL if t not in used]
        avail = easy + sorted(wpset)
        prize = np.zeros(300)
        for t in avail: prize[t - 1] = EP[t - 1]
        if POOLS is not None:
            for t in easy:
                if t not in mine: prize[t - 1] = a.pad_prize * EP[t - 1]
            say(f'--- {name}: pool {len(mine)} assigned + {len(easy)-len(mine)} shared @ {a.pad_prize}')
        if a.covered_prize > 0.0 and not a.refill and not a.assign:
            for t in easy:
                if t in covered: prize[t - 1] = a.covered_prize
            ns = len(covered & set(easy))
            say(f'--- {name}: pool {len(easy)} easy = {len(easy)-ns} fresh + {ns} stepping stones @ {a.covered_prize}')
        if a.assign:
            for t in pad[name]:
                if t not in grp[name]: prize[t - 1] = a.pad_prize * EP[t - 1]
        if a.refill:
            for t in left: prize[t - 1] = a.left_prize
        for x in wpset: prize[x - 1] = a.wp_prize
        stone = a.covered_prize > 0.0 and not a.refill and not a.assign
        nkey = (lambda s: len(set(x[0] for x in s.seq) - covered)) if stone else (lambda s: len(s.seq))
        byn = {}; tic = time.time()
        dvt = np.full(300, a.dvmax)
        for x in wpset: dvt[x - 1] = a.wp_dvmax
        retry = [tuple(float(x) for x in p.split(':')) for p in a.retry.split(',')] if a.retry else []
        for beam_w, drmax, half in port + retry:
            if byn and (beam_w, drmax, half) in retry:
                break                                             # the retry portfolio is only for routes with no full state
            if a.mode == 'segmented':
                t0 = time.time()
                full = fill_segmented(name, wps, easy, a, wf, grid, say, half, beam_w, drmax, K=a.seg_k, seg_depth=a.seg_depth,
                                      prize=prize)
                for s in full:
                    k = nkey(s)
                    if k not in byn or s.fuel < byn[k].fuel: byn[k] = s
                say(f'  {name}: segmented beam {beam_w:.0f} drmax {drmax:.2f} half {half:.0f}: {len(full)} final states, deepest '
                    f'{max((len(s.seq) for s in full), default=0)}  ({time.time()-t0:.0f} s)')
                continue
            lo = np.full(300, -np.inf); hi = np.full(300, np.inf)
            for t, x in wps: lo[x - 1] = t - half * DAY; hi[x - 1] = t + half * DAY
            P = BP(beam=int(beam_w), w_fuel=wf, m_margin=40.0, dv_max=a.dvmax, lin_tofs=np.arange(20, 600 + 1e-9, 10) * DAY,
                   tofs=np.arange(15, 400 + 1e-9, 5) * DAY, lin_drmax=drmax * AU, vinf_cap=a.vinf, tof_refine=True,
                   max_depth=a.max_depth)
            P.win_lo, P.win_hi = lo, hi; P.prize = prize; P.collect = DeepCollect(8); P.dv_max_t = dvt; t0 = time.time()
            best, beam = beam_search(RI.eph(), excluded=sorted(set(ALL) - set(avail)), m0=a.m0, t_launch_grid=grid,
                                     P=P, n_proc=a.nproc, verbose=False)
            full = [s for s in list(P.collect) + list(beam) if wpset <= set(x[0] for x in s.seq)]
            for s in full:
                k = nkey(s)
                if k not in byn or s.fuel < byn[k].fuel: byn[k] = s
            say(f'  {name}: beam {beam_w:.0f} drmax {drmax:.2f} half {half:.0f}: {len(full)} full states, deepest '
                f'{max((len(s.seq) for s in full), default=0)}  ({time.time()-t0:.0f} s)')
        if not byn:
            say(f'{name}: NO full state -- keeping the bare skeleton'); fleet.routes[name] = skel.routes[name]; continue
        say(f'  {name}: Pareto ' + ' '.join(f'{k}:{CM.tank(byn[k].fuel, a.m0):.0f}' for k in sorted(byn)))
        took = None
        for k in sorted(byn, reverse=True)[:a.ntry]:
            s = byn[k]; tour = tour_json(s); t1 = time.time()
            try:
                signal.signal(signal.SIGALRM, _timeout); signal.setitimer(signal.ITIMER_REAL, a.settle_timeout)
                ip, miss, lag = settle_tour(RI.eph(), tour, lambda ip: RI.settle(ip, a.iters))
            except Exception as e:
                say(f'    {k} fb: settle raised {type(e).__name__} ({time.time()-t1:.0f} s)'); continue
            finally:
                signal.setitimer(signal.ITIMER_REAL, 0)
            if ip is None or miss > 150.0:
                say(f'    {k} fb: settle failed (miss {miss:.0f}, {time.time()-t1:.0f} s)'); continue
            tank = float(ip.tank()); val = RI.cost(tank) - a.lam * k
            if a.refill:
                tg = set(int(x) for x in ip.asts); got_left = len(tg & set(left)); lost_own = len(own - tg)
                val = RI.cost(tank) - a.lam * len(tg & own) - a.left_lam * got_left + a.left_lam * lost_own
            ntot = len(ip.asts); nnew = len(set(int(x) for x in ip.asts) - covered)
            val = RI.cost(tank) - a.lam * nnew if not a.refill else val
            say(f'    {k} fb: planner {CM.tank(s.fuel, a.m0):.0f} -> twin {tank:.0f} kg ({ntot} fb, {nnew} new), '
                f'{ip.dv()/ntot:.3f} km/s/fb, J_i {RI.cost(tank):.3f}, value {val:.3f} ({time.time()-t1:.0f} s)')
            if tank <= a.max_tank and (took is None or val < took[0]):
                took = (val, ip, tank, nnew, ntot)
        if took is None:
            if a.refill:
                say(f'{name}: nothing settled -- unchanged'); continue
            say(f'{name}: nothing settled -- keeping the bare skeleton'); fleet.routes[name] = skel.routes[name]; continue
        val, ip, tank, k, ntot = took
        if name in fleet.routes and a.start:                 # refilling an existing route: keep the better one
            tank_old = fleet.routes[name].get('tank', 9e9)
            if k < own_old or (k == own_old and tank >= tank_old):
                say(f'[{name}] refill REJECTED: {k} new vs {own_old} at tank {tank:.0f} vs {tank_old:.0f}'); continue
            say(f'[{name}] refill ACCEPTED: {own_old} -> {k} new targets, tank {tank_old:.0f} -> {tank:.0f}')
        if a.refill:
            val_old = c_old - a.lam * len(own)
            if val >= val_old - 1e-6:
                say(f'[{name}] refill rejected: value {val:.3f} vs current {val_old:.3f}'); continue
            tg = set(int(x) for x in ip.asts)
            say(f'[{name}] refill ACCEPTED: {len(own) + len(wpset)} -> {k} flybys, +{len(tg & set(left))} leftovers, '
                f'-{len(own - tg)} own released, tank {fleet.routes[name]["tank"]:.0f} -> {tank:.0f}, value {val_old:.3f} -> {val:.3f}')
        fleet.routes[name] = dict(st=RI.ist(ip), tank=tank)
        used |= set(int(x) for x in ip.asts)
        fleet.save(out, note=f'after {name}')
        if a.snapshots:
            snap = out / f'snap_{order.index(name):02d}'; snap.mkdir(parents=True, exist_ok=True)
            fleet.save(snap, note=f'prefix through {name}')
        cov = fleet.coverage()
        say(f'[{name}] {len(wps)} waypoints -> {ntot} flybys ({k} new), tank {tank:.0f}, J_i {RI.cost(tank):.3f}; '
            f'fleet covered {len(cov)}, uncovered {len(ALL) - len(cov)}, fleet J {fleet.J():.3f}  ({time.time()-tic:.0f} s)')
    cov = fleet.coverage(); left = sorted(t for t in ALL if t not in cov)
    sumJ = sum(RI.cost(r['tank']) for r in fleet.routes.values())
    json.dump(dict(covered=len(cov), leftovers=left, sum_Ji=sumJ, J_if_all_inserted_free=sumJ + 2,
                   J_now=sumJ + 2 + len(left), tanks={n: r['tank'] for n, r in fleet.routes.items()}),
              open(out / 'result.json', 'w'), indent=1)
    fleet.save(out, note='filled')
    say(f'FINAL: {fleet.summary()}; leftovers {len(left)}: {left}')
    say(f'sum J_i {sumJ:.3f} -> J {sumJ+2:.3f} if the leftovers insert for free, {sumJ+2+len(left):.3f} as misses')



# ----------------------------------------------------------------------------------------------- segmented fill
def fill_segmented(name, wps, avail_easy, a, wf, grid, say, half, beam_w, drmax, K=30, seg_depth=14, prize=None):
    """Waypoint-to-waypoint beam: segment k packs easy targets and must END at waypoint k (its window), the K best
    states at that waypoint start segment k+1; after the last waypoint a free run to the end of the mission.  Every
    segment uses the whole beam width on a small problem instead of losing most of the beam at each deadline.
    Returns the list of final states (all waypoints visited)."""
    from ctoc14.search import _beam_loop, roots, State, T_MISSION
    P = BP(beam=int(beam_w), w_fuel=wf, m_margin=40.0, dv_max=a.dvmax, lin_tofs=np.arange(20, 600 + 1e-9, 10) * DAY,
           tofs=np.arange(15, 400 + 1e-9, 5) * DAY, lin_drmax=drmax * AU, vinf_cap=a.vinf, tof_refine=True, max_depth=seg_depth)
    wpset = [x for _, x in wps]
    dvt = np.full(300, a.dvmax)
    for x in wpset: dvt[x - 1] = a.wp_dvmax
    P.dv_max_t = dvt
    if prize is None:
        prize = np.zeros(300)
        for t in avail_easy: prize[t - 1] = 1.0
        for x in wpset: prize[x - 1] = a.wp_prize
    P.prize = prize
    excluded_base = sorted(set(ALL) - set(avail_easy) - set(wpset))
    starts = None; fronts = []
    for k, (tw, x) in enumerate(wps + [(None, None)]):
        lo = np.full(300, -np.inf); hi = np.full(300, np.inf)
        if x is not None:
            t_end = tw + half * DAY
            for t in avail_easy: hi[t - 1] = t_end          # easy targets only before the segment end
            for y in wpset: hi[y - 1] = -np.inf             # other waypoints are not visitable in this segment
            lo[x - 1] = tw - half * DAY; hi[x - 1] = t_end
            P.mandatory = [x]
        else:
            P.mandatory = []                                # free tail to the end of the mission
            for y in wpset: hi[y - 1] = -np.inf
        P.win_lo, P.win_hi = lo, hi
        ends = []
        for relaxed in (False, True):
            if relaxed:                                       # a skeleton leg the ordinary leg models cannot do
                Q = BP(**{kk: vv for kk, vv in P.__dict__.items() if kk not in ('collect',)})
                Q.tofs = np.arange(15, 700 + 1e-9, 5) * DAY; Q.lin_tofs = np.arange(20, 1000 + 1e-9, 10) * DAY
                Q.lin_drmax = 0.30 * AU; Q.dv_max_t = dvt.copy()
                for y in wpset: Q.dv_max_t[y - 1] = max(a.wp_dvmax, 4.0)
            else:
                Q = P
            Q.collect = DeepCollect(1)
            beam0 = roots(RI.eph(), grid, a.m0, excluded_base, Q) if starts is None else starts
            _beam_loop(RI.eph(), beam0, Q, n_proc=a.nproc, verbose=False)
            pool = list(Q.collect) + (list(beam0) if starts is not None else [])
            ends += [s for s in pool if s.seq and s.seq[-1][0] == x] if x is not None else pool
            if len(ends) >= 5 or x is None:
                break
            say(f'    {name} seg {k}: {len(ends)} states reach waypoint {x} -- {"relaxed retry" if not relaxed else "giving up"}')
        if not ends:
            if x is None:
                return []
            say(f'    {name} seg {k}: DROPPING waypoint {x} ([{(tw-half*DAY)/DAY:.0f},{(tw+half*DAY)/DAY:.0f}] d unreachable); '
                f'it goes to the leftovers')
            fronts.append((k, x, 0, 0, 0)); continue
        uniq = {}
        for s in sorted(ends, key=lambda s: s.score(P)):
            key = s.visited
            if key not in uniq: uniq[key] = s
        ends = sorted(uniq.values(), key=lambda s: s.score(P))
        starts = ends[:K]
        fronts.append((k, x, len(ends), max(len(s.seq) for s in ends), min(len(s.seq) for s in starts)))
        if x is None:
            say('    ' + name + ' segments: ' + ' '.join(f'{k}:{x}->{n}st/{dmax}fb' for k, x, n, dmax, _ in fronts))
            return ends
    return starts


def gen_columns(a, skel, order, hardset, EP, wf, grid, say):
    """Column mode (stage 6 B): for every skeleton run one segmented fill per (pool fraction, seed) pair.  Sub-sampling
    the easy pool is what DECORRELATES the columns -- at full pool every skeleton's best column takes the same popular
    targets, which is why the stage-5 MIP stalled at 223.  A fraction f simulates "half of the pool was taken by
    somebody else", which is the situation the later routes of a fleet actually face.
    EVERY final state of every segmented beam is a column (deduplicated by target set, cheapest fuel kept), not just
    the Pareto front: one 100 s fill then yields tens of columns instead of one.  No settling (planner tank = twin
    within ~3-7 %; tools/skel_select.py applies --tank-cal)."""
    port = [tuple(float(x) for x in p.split(':')) for p in a.portfolio.split(',')]
    cf = open(a.columns, 'a'); ncols = 0
    easy_all = [t for t in ALL if t not in hardset]
    fracs = [float(x) for x in a.pool_fracs.split(',')]
    for name in order:
        st = skel.routes[name]['st']; wps = sorted((float(t), int(x)) for x, t in zip(st['asts'], st['tf']))
        wpset = set(x for _, x in wps); seen = {}; tic = time.time(); j = -1
        for frac in fracs:
            for seed in range(a.pool_seeds if frac < 1.0 else 1):
                j += 1
                rng = np.random.default_rng(hash((name, frac, seed)) % (2 ** 31))
                pool = easy_all if frac >= 1.0 else sorted(
                    rng.choice(easy_all, size=max(8, int(round(frac * len(easy_all)))), replace=False).tolist())
                prize = np.zeros(300)
                for t in pool: prize[t - 1] = EP[t - 1]
                for x in wpset: prize[x - 1] = a.wp_prize
                got = 0; deep = 0
                beam_w, drmax, half = port[j % len(port)]    # cycle the beam settings over the pool configs
                for s in fill_segmented(name, wps, pool, a, wf, grid, say, half, beam_w, drmax, K=a.seg_k,
                                        seg_depth=a.seg_depth, prize=prize):
                    if len(s.seq) < a.col_min: continue
                    key = frozenset(x[0] for x in s.seq)
                    if key in seen and seen[key].fuel <= s.fuel: continue
                    seen[key] = s; got += 1; deep = max(deep, len(s.seq))
                    tg = sorted(x[0] for x in s.seq)         # write as we go: a killed job keeps its columns
                    cf.write(json.dumps(dict(skel=name, col=len(tg), n=len(tg), targets=tg, fuel=float(s.fuel),
                                             tank=float(CM.tank(s.fuel, a.m0)), tour=tour_json(s))) + '\n'); ncols += 1
                cf.flush()
                say(f'  {name} f{frac:.2f}/s{seed} half{half:.0f}: pool {len(pool)}, +{got} columns (deepest {deep}), '
                    f'{len(seen)} distinct')
        reach = set()
        for s in seen.values(): reach |= set(x[0] for x in s.seq)
        say(f'[{name}] {len(seen)} columns, depth {min((len(s.seq) for s in seen.values()), default=0)}-'
            f'{max((len(s.seq) for s in seen.values()), default=0)}, {len(reach - wpset)} distinct easy targets reached '
            f'({time.time()-tic:.0f} s)')
    cf.close(); say(f'{ncols} columns written to {a.columns}')


if __name__ == '__main__':
    main()
