"""SR probe driver = tools/skel_fill.py main() with ONE fix in fill_segmented (copied here, shared tool untouched):
segment 0 also collects the ROOT states themselves, so a host whose first flyby is a waypoint (every host, when all its
own flybys are pinned) can start.  Used to price the segment re-plan closer: host flybys pinned +-half, leftovers prize 1."""
import os, sys, pathlib
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import skel_fill as SF
import run_ialns as RI
from greedy_cover import BP, ALL, DeepCollect
from ctoc14.constants import DAY, AU


def fill_segmented(name, wps, avail_easy, a, wf, grid, say, half, beam_w, drmax, K=30, seg_depth=14, prize=None):
    from ctoc14.search import _beam_loop, roots
    P = BP(beam=int(beam_w), w_fuel=wf, m_margin=40.0, dv_max=a.dvmax, lin_tofs=np.arange(20, 600 + 1e-9, 10) * DAY,
           tofs=np.arange(15, 400 + 1e-9, 5) * DAY, lin_drmax=drmax * AU, vinf_cap=a.vinf, tof_refine=True, max_depth=seg_depth)
    P.n_per_target = int(os.environ.get('SR_NPT', '2')); P.max_children = int(os.environ.get('SR_MC', '60'))
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
            for t in avail_easy: hi[t - 1] = t_end
            for y in wpset: hi[y - 1] = -np.inf
            lo[x - 1] = tw - half * DAY; hi[x - 1] = t_end
            P.mandatory = [x]
        else:
            P.mandatory = []
            for y in wpset: hi[y - 1] = -np.inf
        P.win_lo, P.win_hi = lo, hi
        ends = []
        for relaxed in (False, True):
            if relaxed:
                Q = BP(**{kk: vv for kk, vv in P.__dict__.items() if kk not in ('collect',)})
                Q.tofs = np.arange(15, 700 + 1e-9, 5) * DAY; Q.lin_tofs = np.arange(20, 1000 + 1e-9, 10) * DAY
                Q.lin_drmax = 0.30 * AU; Q.dv_max_t = dvt.copy()
                for y in wpset: Q.dv_max_t[y - 1] = max(a.wp_dvmax, 4.0)
            else:
                Q = P
            Q.collect = DeepCollect(1)
            beam0 = roots(RI.eph(), grid, a.m0, excluded_base, Q) if starts is None else starts
            _beam_loop(RI.eph(), beam0, Q, n_proc=a.nproc, verbose=False)
            pool = list(Q.collect) + list(beam0)                      # FIX: roots count as segment-0 ends too
            ends += [s for s in pool if s.seq and s.seq[-1][0] == x] if x is not None else pool
            if len(ends) >= 5 or x is None:
                break
            say(f'    {name} seg {k}: {len(ends)} states reach waypoint {x} -- {"relaxed retry" if not relaxed else "giving up"}')
        if not ends:
            if x is None:
                return []
            say(f'    {name} seg {k}: DROPPING waypoint {x}')
            fronts.append((k, x, 0, 0, 0)); continue
        uniq = {}
        TB = float(os.environ.get('SR_TBIN', '0'))            # >0: keep timing variants (visited, t bin) at a waypoint
        for s in sorted(ends, key=lambda s: s.score(P)):
            key = (s.visited, int(s.t // (TB * DAY))) if TB > 0 else s.visited
            if key not in uniq: uniq[key] = s
        ends = sorted(uniq.values(), key=lambda s: s.score(P))
        starts = ends[:K]
        fronts.append((k, x, len(ends), max(len(s.seq) for s in ends), min(len(s.seq) for s in starts)))
        if x is None:
            say('    ' + name + ' segments: ' + ' '.join(f'{k}:{x}->{n}st/{dmax}fb' for k, x, n, dmax, _ in fronts))
            return ends
    return starts


SF.fill_segmented = fill_segmented
if __name__ == '__main__':
    SF.main()
