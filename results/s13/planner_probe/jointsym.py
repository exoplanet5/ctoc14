"""Copy of ctoc14.jointsearch.joint_beam_search (the shared module is NOT modified) with two probe options:
  canon=True : the beam dedup key ignores craft ORDER -- (visited, sorted per-craft 10-d time bins) instead of
               (visited, per-craft time bins in craft order).  With K identical craft the original key keeps every
               permutation of the same fleet (up to K! copies with identical score), which eats the beam.
  stats      : per level, how many of the kept states are permutation duplicates of another kept state.
Everything else (expansion of the earliest craft, wait/retire, scoring) is ctoc14.jointsearch unchanged."""
import numpy as np
import ctoc14.jointsearch as JS
from ctoc14.jointsearch import Joint, _batch, _init, _CTX
from ctoc14.search import Params, mask_of, UNREACHABLE
from ctoc14.constants import DAY


def _canon(c):
    return (c.visited, tuple(sorted(int((s.t if s is not None else 0) // (10 * DAY)) for s in c.craft)))


def _perm_free(c):
    """Order-free identity of a joint state: the multiset of per-craft target sequences."""
    return tuple(sorted(tuple(x[0] for x in s.seq) if s is not None else () for s in c.craft))


def joint_beam_search(eph, m0s, launch_grid, P=None, excluded=(), n_proc=1, verbose=True, canon=False, stats=None):
    P = P or Params(); K = len(m0s)
    exm = mask_of(set(excluded) | UNREACHABLE)
    j0 = Joint(craft=(None,) * K, active=(True,) * K, visited=0, launch_t=(0.0,) * K)
    beam = [j0]; best = j0; pool = None
    _init(eph, P, m0s, launch_grid, exm)
    if n_proc > 1:
        import multiprocessing as mp
        pool = mp.get_context('fork').Pool(n_proc, initializer=_init, initargs=(eph, P, m0s, launch_grid, exm, _CTX['roots']))
    try:
        for depth in range(1, 400):
            if pool is not None and len(beam) >= 2 * n_proc:
                chunks = [beam[i::n_proc] for i in range(n_proc)]
                children = [c for r in pool.map(_batch, chunks) for c in r]
            else:
                children = _batch(beam)
            if not children:
                break
            children.sort(key=lambda j: j.score(P))
            seen = set(); nb = []
            for c in children:
                key = _canon(c) if canon else (c.visited, tuple(int((s.t if s is not None else 0) // (10 * DAY)) for s in c.craft))
                if key in seen: continue
                seen.add(key); nb.append(c)
                if len(nb) >= P.beam: break
            beam = nb; top = beam[0]
            coll = getattr(P, 'collect_joint', None)
            if coll is not None: coll.extend(beam)
            if stats is not None:
                stats.append(dict(depth=depth, children=len(children), kept=len(nb), perm_free=len(set(_perm_free(c) for c in nb))))
            if top.n() > best.n() or (top.n() == best.n() and top.score(P) < best.score(P)): best = top
            if verbose and (depth % 10 == 0 or depth < 3):
                pf = stats[-1]['perm_free'] if stats else -1
                print(f'depth {depth:3d}: {len(children):6d} children -> beam {len(beam):4d} ({pf} order-free distinct); n={top.n()} '
                      f'fuel={top.fuel():.0f} per craft {[s.n() if s else 0 for s in top.craft]}', flush=True)
    finally:
        if pool is not None: pool.close(); pool.join()
    return best
