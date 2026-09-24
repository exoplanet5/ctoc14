"""s15b: KEEP-OWN production-planner re-plan of one route (tools/s15_isogen.run_column with patched prizes).

The production beam (ctoc14.search, m0 1600, Lambert + linear legs; the planner that built the heads) re-plans a route
over pool = targets(r) + misses M, with the route's own targets prized own_prize (default 2.0) and M prized 1 + p, so a
miss enters only as an ADDITION; legs into M get the dvS cap (s15_isogen.params).  The front is settled in the twin by
s15_isogen's walk (deepest first) and every settled state is saved as a column; then s15_isogen's absorb step offers M
to the best settled hosts (twin insertion).  One job ~10-20 min (vs ~1 h for a twin-beam re-plan).

Job (JSON): tag, out, pool, S, p, own, own_prize, member ('H'|'T1'|'T2'|'T3'|'T5'|'T6'), seed, beam, ... (s15_isogen.DEFAULTS)
usage: s15b_isokeep.py run JOB.json"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, pathlib
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import numpy as np
import s15_isogen as IG

_OWN = set(); _OWN_PRIZE = 1.0
_orig = IG.prize_vector


def prize_vector(pool, S, p, seed, jit):
    pz = _orig(pool, S, p, seed, jit)
    for t in _OWN:
        if pz[t - 1] > 0:
            pz[t - 1] = max(pz[t - 1], _OWN_PRIZE)
    return pz


IG.prize_vector = prize_vector


def run(J):
    global _OWN, _OWN_PRIZE
    _OWN = set(int(x) for x in J.get('own', [])); _OWN_PRIZE = float(J.get('own_prize', 2.0))
    J = dict(J); J.pop('own', None); J.pop('own_prize', None); J.pop('runner', None)
    for k in ('fleet_pre', 'replaced', 'src_fleet', 'src_route'):
        J.pop(k, None)
    J['kind'] = 'column'
    return IG.run_column(J)


if __name__ == '__main__':
    if len(sys.argv) != 3 or sys.argv[1] != 'run':
        raise SystemExit(__doc__)
    J = json.load(open(sys.argv[2]))
    m = run(J)
    print(json.dumps(dict(cols=len(m.get('cols', [])), wall_s=m.get('wall_s'))))
