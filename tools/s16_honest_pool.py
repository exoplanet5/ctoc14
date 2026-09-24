"""Stage 16: re-price EVERY candidate route with an honest thrust cap, so the fleet selector compares real costs.

Why.  ImpulsiveProblem.tcap (what 0.43 N delivers over one dtau bin) is enforced PER IMPULSE.  Planner-born twins
(impulsive.from_tour / settle_tour) carry extra irregular nodes a fraction of a day after other nodes, so one 20 d window
can hold two impulses at the cap -- up to 1.6x the real 0.5 N capability.  Their tanks are optimistic by 0-30 kg, and by
DIFFERENT amounts per route, so every fleet selection made over this pool (stages 13-15) ranked routes on biased prices.
On the submitted 9-craft fleet the honest re-pricing cost +0.056 J before relocation (docs: stage-15b memory).

What.  For every distinct route file: if it already has only regular nodes and no 20 d window exceeds the cap, it is
honest as is (copied, tank unchanged).  Otherwise tools/s15b_regrid.regrid (sum impulses into regular dtau bins) +
run_ialns.settle, kept only if the settle converges (miss <= 150 km) and the tank stays <= --max-tank.

Output: OUT/route_<key>.npz (run_ialns ist format) + OUT/index.jsonl (src, key, n, tank_in, tank_out, miss, status).
Resumable: keys already in the index are skipped.

Usage: s16_honest_pool.py OUT PAT [PAT ...] [--nproc 10] [--timeout 240]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, glob, time, signal, hashlib, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
os.chdir(ROOT)
import numpy as np
import run_ialns as RI
from s15b_regrid import regrid
from ctoc14.constants import DAY

A = None


class _TO(Exception):
    pass


def _alarm(signum, frame):
    raise _TO()


def load(f):
    z = np.load(f); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL'])
    return st


def key_of(st):
    asts = tuple(int(x) for x in st['asts'])
    tf = tuple(int(round(float(t) / DAY)) for t in st['tf'])
    h = hashlib.sha1(repr((round(st['tL'] / DAY, 1), asts, tf)).encode()).hexdigest()[:16]
    return h


def honest_already(ip):
    ts = np.asarray(ip.ts)
    if len(ts) > 1 and np.any(np.abs(np.diff(ts) - ip.h) > 0.01 * DAY):
        return False
    return bool((np.linalg.norm(ip.Ts, axis=1) / ip.tcap).max() <= 1.001) if len(ip.Ts) else True


def job(args):
    f, key = args
    out = pathlib.Path(A.out)
    t0 = time.time()
    rec = dict(src=f, key=key)
    signal.signal(signal.SIGALRM, _alarm); signal.alarm(int(A.timeout))
    try:
        st = load(f); ip = RI.ipr(st); tank_in = float(ip.tank())
        rec.update(n=len(ip.asts), tank_in=tank_in)
        if honest_already(ip):
            Yf, _ = ip.integrate(); dm, _ = ip.misses(Yf); miss = float(np.linalg.norm(dm, axis=1).max()) if len(dm) else 0.0
            if miss > 150.0:
                miss = RI.settle(ip, A.iters)
            g = ip; status = 'regular'
        else:
            g = regrid(ip, A.dtau * DAY); miss = RI.settle(g, A.iters); status = 'regridded'
        tank = float(g.tank())
        rec.update(tank_out=tank, miss=float(miss))
        if miss <= 150.0 and tank <= A.max_tank:
            np.savez(out / f'route_{key}.npz', **RI.ist(g)); rec['status'] = status
        else:
            rec['status'] = 'rejected'
    except _TO:
        rec['status'] = 'timeout'
    except Exception as e:
        rec['status'] = 'error'; rec['err'] = repr(e)[:120]
    finally:
        signal.alarm(0)
    rec['dt'] = round(time.time() - t0, 1)
    return rec


def main():
    global A
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('pats', nargs='+')
    ap.add_argument('--nproc', type=int, default=10); ap.add_argument('--timeout', type=float, default=240.0)
    ap.add_argument('--dtau', type=float, default=20.0); ap.add_argument('--iters', type=int, default=150)
    ap.add_argument('--max-tank', type=float, default=1400.0)
    A = ap.parse_args()
    out = pathlib.Path(A.out); out.mkdir(parents=True, exist_ok=True)
    idx = out / 'index.jsonl'
    done = set()
    if idx.exists():
        for l in open(idx):
            try: done.add(json.loads(l)['key'])
            except Exception: pass
    files = sorted(set(f for p in A.pats for f in glob.glob(p, recursive=True)))
    seen = {}
    for f in files:
        try:
            k = key_of(load(f))
        except Exception:
            continue
        seen.setdefault(k, f)
    todo = [(f, k) for k, f in seen.items() if k not in done]
    print(f'{len(files)} files -> {len(seen)} distinct routes; {len(done)} already done; {len(todo)} to do on {A.nproc} procs', flush=True)
    t0 = time.time(); n = 0; stat = {}
    # longest routes first so the tail of the run is short jobs
    todo.sort(key=lambda fk: -os.path.getsize(fk[0]))
    with mp.get_context('fork').Pool(A.nproc, maxtasksperchild=50) as pl, open(idx, 'a') as fo:
        for rec in pl.imap_unordered(job, todo, chunksize=1):
            fo.write(json.dumps(rec) + '\n'); fo.flush(); n += 1
            stat[rec['status']] = stat.get(rec['status'], 0) + 1
            if n % 50 == 0 or n == len(todo):
                el = time.time() - t0
                print(f'[{time.strftime("%H:%M:%S")}] {n}/{len(todo)} ({el/60:.0f} min, eta {el/n*(len(todo)-n)/60:.0f} min) {stat}', flush=True)
    print('done', stat, flush=True)


if __name__ == '__main__':
    main()
