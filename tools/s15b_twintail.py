"""s15b: KEEP-OWN variant of tools/s15_twintail.py (imports it; the shared tool is not edited).

s15_twintail prizes every pool target 1 and the premium set S 1 + p, so a re-plan of route r over targets(r) + M TRADES
own targets for misses (EF296a_h04: 33 and 280 taken, 5-6 of h04's own dropped).  Here the job may also carry
    own        the route's own targets, prized own_prize (e.g. 2.0) -- dropping one costs more than any miss gains,
    own_prize  (default 2.0),
so the beam keeps the whole own set whenever it can and takes a miss only as an ADDITION (designed-in absorption).
Everything else (beam, front, columns, checkpoints, chain) is s15_twintail.run.

usage: s15b_twintail.py run JOB.json"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, pathlib
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import numpy as np
import s15_twintail as TT

_OWN = set(); _OWN_PRIZE = 1.0
_orig_prize_vec = TT.prize_vec


def prize_vec(pool, S, p):
    pz = _orig_prize_vec(pool, S, p)
    for t in _OWN:
        if pz[t - 1] > 0:
            pz[t - 1] = max(pz[t - 1], _OWN_PRIZE)
    return pz


TT.prize_vec = prize_vec


def run(J):
    global _OWN, _OWN_PRIZE
    _OWN = set(int(x) for x in J.get('own', [])); _OWN_PRIZE = float(J.get('own_prize', 2.0))
    return TT.run(J)


if __name__ == '__main__':
    if len(sys.argv) != 3 or sys.argv[1] != 'run':
        raise SystemExit(__doc__)
    m = run(json.load(open(sys.argv[2])))
    print(json.dumps(dict(covered=m.get('covered'), cols=len(m.get('cols', [])), wall_s=m.get('wall_s'))))
