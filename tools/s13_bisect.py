"""Stage 13 (PCC-8G) settle machinery and bisection of failing planner tours.  docs/stage13_solver_plan.md 4.2.

Settled (plan section 4): miss <= 150 km AND twin tank <= 1.15 x planner tank AND twin tank <= the route's cap.
Every settle runs `impulsive.settle_tour` (or `twin_wait.settle_tour_wait` for joint-planner tours) with
`IM.LIN_FALLBACK = 2.5`, `RI.settle(ip, 100)`, under a 240 s SIGALRM (`greedy_cover._timeout`).

  settle_one(tour, kind)           -> dict(st, tank, miss, sec, err)      (st = RI.ist(ip) or None)
  failing_leg(tour, ...)           -> (k_ok, k_fail, leg, info)            binary (or k-ary, with a pool) prefix search
  settled_prefix(tour, ...)        -> st | None                            the longest settling prefix found
  prefix_planner_tanks(tour)       -> planner tank estimate of every prefix (used for the 1.15x guard on prefixes)

A small content-addressed disk cache (Cache) makes every beam / settle result replayable, so a restarted pass continues
from disk instead of recomputing (and two passes over the same pool share work).
Built from results/s13/seed_passes/scripts/bisect_settle.py (the shared modules are imported, never edited).

CLI (measurement): s13_bisect.py TOUR.json [--nproc 4]  -> prints (k_ok, k_fail, failing leg)."""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, signal, pathlib, hashlib, pickle
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools'), str(ROOT / 'results/s13/planner_probe')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import numpy as np
import run_ialns as RI
import ctoc14.impulsive as IM
from greedy_cover import CM, _timeout
from ctoc14.constants import VE, TMAX

IM.LIN_FALLBACK = 2.5
SETTLE_TIMEOUT = 240.0
MISS_OK = 150.0
GUARD = 1.15
CODE_TAG = 's13v1'


# ----------------------------------------------------------------------------------------------------------- cache
class Cache:
    """Content-addressed pickle store: key = sha1 of a canonical JSON of the job inputs.  Writes are atomic (tmp +
    os.replace), so concurrent passes may share one directory."""

    def __init__(self, d):
        self.d = pathlib.Path(d) if d else None
        if self.d is not None:
            self.d.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def key(obj):
        return hashlib.sha1(json.dumps(obj, sort_keys=True, default=_jdefault).encode()).hexdigest()

    def _p(self, kind, k):
        return self.d / kind / k[:2] / f'{k}.pkl'

    def get(self, kind, k):
        if self.d is None:
            return None
        p = self._p(kind, k)
        if not p.exists():
            return None
        try:
            with open(p, 'rb') as f:
                return pickle.load(f)
        except Exception:
            return None

    def put(self, kind, k, val):
        if self.d is None:
            return
        p = self._p(kind, k); p.parent.mkdir(parents=True, exist_ok=True)
        tmp = p.with_suffix(f'.tmp{os.getpid()}')
        with open(tmp, 'wb') as f:
            pickle.dump(val, f, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp, p)


def _jdefault(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (set, frozenset)):
        return sorted(o)
    raise TypeError(type(o))


# ---------------------------------------------------------------------------------------------------------- settle
def settle_one(tour, kind='stock', timeout=SETTLE_TIMEOUT, iters=100):
    """Settle one planner tour.  kind 'stock' = impulsive.settle_tour, 'wait' = twin_wait.settle_tour_wait (joint
    tours with coasts).  Returns dict(st, tank, miss, sec, err); st is None when settle_tour gave no ip."""
    tic = time.time(); ip = None; miss = float('inf'); err = None
    try:
        signal.signal(signal.SIGALRM, _timeout); signal.setitimer(signal.ITIMER_REAL, timeout)
        if kind == 'wait':
            import twin_wait as TW
            ip, miss, lag = TW.settle_tour_wait(RI.eph(), tour, lambda ip: RI.settle(ip, iters))
        else:
            ip, miss, lag = IM.settle_tour(RI.eph(), tour, lambda ip: RI.settle(ip, iters))
    except Exception as e:                       # singular linearisation / timeout
        ip = None; miss = float('inf'); err = type(e).__name__
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    if ip is None:
        return dict(st=None, tank=None, miss=float(miss) if np.isfinite(miss) else float('inf'), sec=time.time() - tic, err=err)
    return dict(st=RI.ist(ip), tank=float(ip.tank()), miss=float(miss), sec=time.time() - tic, err=err)


def settle_key(tour, kind, timeout=SETTLE_TIMEOUT, iters=100):
    t = dict(t_launch=tour['t_launch'], vinf=list(tour['vinf']),
             legs=[[int(l['ast']), float(l['t_flyby']), float(l['tof'])] for l in tour['legs']])
    return Cache.key(dict(v=CODE_TAG, kind=kind, tour=t, timeout=timeout, iters=iters, lf=IM.LIN_FALLBACK))


def is_ok(res, ptank, cap=None, guard=GUARD):
    """The plan's 'settled': miss <= 150 km, twin <= guard x planner, twin <= cap."""
    if res is None or res.get('st') is None or res.get('tank') is None:
        return False
    if not (res['miss'] <= MISS_OK):
        return False
    if res['tank'] > guard * ptank + 1e-9:
        return False
    if cap is not None and res['tank'] > cap + 1e-9:
        return False
    return True


def failed_settle(res, ptank, guard=GUARD):
    """A real settle FAILURE (not merely over the route cap): no ip, miss > 150 km, or the 1.15x guard violated."""
    return not is_ok(res, ptank, cap=None, guard=guard)


def w_settle(job):
    """Pool worker: job = dict(tour, kind, timeout).  Returns the settle_one dict plus the job's 'id'."""
    r = settle_one(job['tour'], job.get('kind', 'stock'), job.get('timeout', SETTLE_TIMEOUT))
    r['id'] = job.get('id')
    return r


def settle_many(jobs, pmap, cache=None):
    """Settle a list of jobs (dict(tour, kind, id)) through `pmap`, with the cache.  Returns results in job order."""
    out = [None] * len(jobs); todo = []
    for i, j in enumerate(jobs):
        k = settle_key(j['tour'], j.get('kind', 'stock'), j.get('timeout', SETTLE_TIMEOUT))
        j['_key'] = k
        v = cache.get('settle', k) if cache is not None else None
        if v is not None:
            v = dict(v); v['id'] = j.get('id'); v['cached'] = True; out[i] = v
        else:
            todo.append(i)
    if todo:
        res = pmap(w_settle, [jobs[i] for i in todo])
        for i, r in zip(todo, res):
            if cache is not None:
                cache.put('settle', jobs[i]['_key'], {k: v for k, v in r.items() if k != 'id'})
            r['cached'] = False; out[i] = r
    return out


# ----------------------------------------------------------------------------------------------- prefix planner tanks
def prefix_planner_tanks(tour, m0=None):
    """Planner tank of every prefix legs[:k], k = 1..n: the planner's per-leg propellant rebuilt from the leg dv
    estimates (Lambert-junction formula, dm = m (1 - exp(-kappa dv / ve)), kappa = 1 + dv / (2 a tof)), scaled so the
    total equals the tour's fuel_est, then CostModel.tank.  Used only for the 1.15x guard on prefixes."""
    m0 = float(tour.get('m0', 1600.0) if m0 is None else m0)
    m = m0; dms = []
    for l in tour['legs']:
        dv = float(l['dv_est']); tof = max(float(l['tof']), 1.0)
        a = TMAX / m * 1e-3
        kappa = 1.0 + dv / (2 * a * tof)
        dm = m * (1.0 - np.exp(-kappa * dv / VE)); dms.append(dm); m -= dm
    cum = np.cumsum(dms)
    tot = float(tour.get('fuel_est', cum[-1] if len(cum) else 0.0))
    scale = tot / cum[-1] if len(cum) and cum[-1] > 1e-9 else 1.0
    return np.asarray(CM.tank(cum * scale, m0), float)


def prefix_tour(tour, k, ptanks=None):
    t = dict(tour); t['legs'] = list(tour['legs'][:k]); t['n_flybys'] = k
    if ptanks is not None:
        t['_ptank'] = float(ptanks[k - 1])
    return t


# ------------------------------------------------------------------------------------------------------ bisection
def failing_leg(tour, max_settles=6, pmap=None, width=1, cache=None, guard=GUARD, say=None, full_failed=True,
                timeout=SETTLE_TIMEOUT):
    """Find the first prefix length at which a planner tour stops settling (settle-ability assumed monotone in k).
    Serial binary search (width 1), or k-ary with `width` parallel probes per round through `pmap`.
    A prefix 'settles' iff is_ok(res, planner prefix tank) (miss <= 150 km, 1.15x guard; no cap).
    Returns (k_ok, k_fail, leg, info); info = dict(probes=[(k, ok, tank, ptank, miss, sec)], st_ok=state of the longest
    settled prefix or None, settles=#settles).  full_failed: the caller already knows the full tour fails."""
    legs = tour['legs']; n = len(legs); pt = prefix_planner_tanks(tour)
    lo, hi = 0, n
    probes = []; st_ok = None; used = 0
    pm = pmap if pmap is not None else (lambda f, js: [f(j) for j in js])
    if not full_failed:
        r = settle_many([dict(tour=tour, kind='stock', id=n, timeout=timeout)], pm, cache)[0]; used += 1
        ok = is_ok(r, pt[-1], guard=guard)
        probes.append((n, ok, r['tank'], float(pt[-1]), r['miss'], round(r['sec'], 1)))
        if ok:
            return n, None, None, dict(probes=probes, st_ok=r['st'], settles=used)
    while hi - lo > 1 and used < max_settles:
        q = max(1, min(width, hi - lo - 1, max_settles - used))
        ks = sorted(set(int(round(lo + (hi - lo) * (i + 1) / (q + 1))) for i in range(q)))
        ks = [k for k in ks if lo < k < hi] or [(lo + hi) // 2]
        res = settle_many([dict(tour=prefix_tour(tour, k), kind='stock', id=k, timeout=timeout) for k in ks], pm, cache)
        used += len(ks)
        oks = []
        for k, r in zip(ks, res):
            ok = is_ok(r, pt[k - 1], guard=guard); oks.append(ok)
            probes.append((k, ok, r['tank'], float(pt[k - 1]), r['miss'], round(r['sec'], 1)))
            if say:
                say(f'      bisect prefix {k}/{n}: {"OK" if ok else "fail"} (twin {r["tank"] if r["tank"] is None else round(r["tank"])}'
                    f' vs planner {pt[k - 1]:.0f}, miss {r["miss"]:.0f} km, {r["sec"]:.0f} s)')
        first_fail = next((k for k, ok in zip(ks, oks) if not ok), None)
        if first_fail is None:
            lo = ks[-1]; st_ok = res[-1]['st']
        else:
            below = [(k, r) for k, ok, r in zip(ks, oks, res) if ok and k < first_fail]
            if below:
                lo = below[-1][0]; st_ok = below[-1][1]['st']
            hi = first_fail
    leg = legs[hi - 1] if hi is not None and hi >= 1 else None
    return lo, hi, leg, dict(probes=probes, st_ok=st_ok, settles=used)


def settled_prefix(tour, **kw):
    """State (RI.ist) of the longest settling prefix found by failing_leg, or None (kept as a library column)."""
    lo, hi, leg, info = failing_leg(tour, **kw)
    return info['st_ok'] if lo and lo > 0 else None


def main():
    import argparse, multiprocessing as mp
    ap = argparse.ArgumentParser(); ap.add_argument('tour'); ap.add_argument('--nproc', type=int, default=1)
    ap.add_argument('--max-settles', type=int, default=6); ap.add_argument('--cache', default='')
    a = ap.parse_args()
    tour = json.load(open(a.tour))
    cache = Cache(a.cache) if a.cache else None
    RI.eph()
    if a.nproc > 1:
        with mp.get_context('fork').Pool(a.nproc, maxtasksperchild=2) as pool:
            r = failing_leg(tour, max_settles=a.max_settles, pmap=lambda f, js: pool.map(f, js, chunksize=1),
                            width=a.nproc, cache=cache, say=print, full_failed=False)
    else:
        r = failing_leg(tour, max_settles=a.max_settles, cache=cache, say=print, full_failed=False)
    print('last settling prefix', r[0], 'first failing', r[1], 'failing leg', r[2], 'settles', r[3]['settles'])


if __name__ == '__main__':
    main()
