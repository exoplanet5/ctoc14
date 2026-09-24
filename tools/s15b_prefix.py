"""s15b: PREFIX-SEEDED keep-own twin beam -- re-plan only the TAIL of a route, from its k-th flyby on.

The launch-rooted twin beam cannot always re-fly a production route's early part (K2Bt1x64 skipped t1's 2nd target at
level 8).  Here the first k flybys (time order) of the route are kept as a settled twin PREFIX (impulse nodes after the
k-th flyby dropped, re-settled with run_ialns.settle), written as the level-k checkpoint of s15_twintail's beam, and the
beam (keep-own prizes, tools/s15b_twintail.py) continues over the route's remaining targets + the extra targets S.

Job (JSON) = s15_twintail job + route (file), k (flybys kept), own (the route's targets; prize own_prize), S (extra).
usage: s15b_prefix.py run JOB.json"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, signal, pickle, pathlib
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)
import numpy as np
import run_ialns as RI
import s14_twinbeam as TB
import s15b_twintail as KW
from greedy_cover import _timeout
from ctoc14.impulsive import ImpulsiveProblem
from ctoc14.constants import DAY


def make_prefix(route_file, k, timeout=600.0):
    z = np.load(ROOT / route_file); st = {kk: z[kk] for kk in z.files}; st['tL'] = float(st['tL'])
    ip = RI.ipr(st)
    o = np.argsort(ip.tf); keep = o[:k]; tk = float(ip.tf[o[k - 1]])
    m = ip.ts <= tk
    ip2 = ImpulsiveProblem(RI.eph(), ip.tL, ip.vinf, ip.ts[m], ip.Ts[m], ip.tf[keep], [ip.asts[i] for i in keep])
    try:
        signal.signal(signal.SIGALRM, _timeout); signal.setitimer(signal.ITIMER_REAL, timeout)
        miss = float(RI.settle(ip2, 100))
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    nreg = len(TB.centres(float(ip.tL), tk))
    return TB.pack(ip2, nreg, miss, 'prefix'), miss


def run(J):
    out = ROOT / J['out']; out.mkdir(parents=True, exist_ok=True)
    if (out / 'meta.json').exists():
        return json.load(open(out / 'meta.json'))
    sd = out / 's1'; sd.mkdir(parents=True, exist_ok=True)
    if not (sd / 'ckpt.pkl').exists():
        say = RI.logger(out / 'log.txt'); tic = time.time()
        pre, miss = make_prefix(J['route'], int(J['k']))
        say(f'prefix {J["route"]} k={J["k"]}: {len(pre["asts"])} flybys to {pre["t_end"] / DAY:.0f} d, tank {pre["tank"]:.1f} '
            f'kg, miss {miss:.0f} km ({time.time() - tic:.0f} s)')
        if miss > 150.0:
            json.dump(dict(job=J, error=f'prefix did not settle (miss {miss:.0f} km)', cols=[], chain=dict(steps=[], cols=[])),
                      open(out / 'meta.json', 'w'), indent=1)
            return json.load(open(out / 'meta.json'))
        C = dict(level=len(pre['asts']), beam=[pre], front={}, hist=[], end=None, settles=0, settle_ok=0, wall=0.0)
        TB.save_ckpt(sd, C)
    return KW.run(J)


if __name__ == '__main__':
    if len(sys.argv) != 3 or sys.argv[1] != 'run':
        raise SystemExit(__doc__)
    m = run(json.load(open(sys.argv[2])))
    print(json.dumps(dict(covered=m.get('covered'), cols=len(m.get('cols', [])), wall_s=m.get('wall_s'))))
