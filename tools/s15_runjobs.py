"""Stage 15: run a JSON list of s15_isogen jobs (column / chain / absorb / close / grow) on a fork Pool(--nproc <= 8),
longest first, skipping jobs whose meta.json exists (resumable, like s15_loop.run_jobs).

usage: s15_runjobs.py JOBS.json [--nproc 2]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)
import s15_isogen as IG
import run_ialns as RI


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('jobs'); ap.add_argument('--nproc', type=int, default=2)
    a = ap.parse_args()
    if a.nproc > 8:
        raise SystemExit('CPU budget: --nproc <= 8')
    jobs = json.load(open(a.jobs))
    todo = [j for j in jobs if not (ROOT / j['out'] / 'meta.json').exists()]
    todo.sort(key=lambda j: -j.get('_est', 1.0))
    print(f'[{time.strftime("%H:%M:%S")}] {len(jobs)} jobs, {len(todo)} to run on {a.nproc} procs', flush=True)
    RI.eph()
    with mp.get_context('fork').Pool(a.nproc, maxtasksperchild=1) as pool:
        for res in pool.imap_unordered(IG.run_job, [{k: v for k, v in j.items() if k != '_est'} for j in todo], chunksize=1):
            if res.get('error'):
                print(f'[{time.strftime("%H:%M:%S")}] job {res.get("tag")} ERROR {res["error"][-300:]}', flush=True); continue
            jb = res.get('job', {}); cols = res.get('cols', [])
            ch = res.get('chain', {})
            print(f'[{time.strftime("%H:%M:%S")}] {jb.get("tag")}: {len(cols)} columns, deepest {max([c["n"] for c in cols] + [0])}'
                  + (f', chain covered {len(ch.get("covered", []))}/{res.get("pool_n")} taken {[t["n"] for t in ch.get("taken", [])]}'
                     if ch else '') + f' ({res.get("wall_s")} s)', flush=True)
    print(f'[{time.strftime("%H:%M:%S")}] done', flush=True)


if __name__ == '__main__':
    main()
