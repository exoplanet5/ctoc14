"""Stage 16: EXACT conversion (continuous thrust + validation) of individual route files, for exact-aware variant choice.

Why.  The twin sum J_i is only a proxy: the export (run_ialns.w_export) turns each impulsive route into a continuous-thrust
craft whose m0 differs from the twin tank by -1.2..+4.3 kg (s15a: +8.7 kg over 9 craft; fleet E: +22.5 kg).  Many
twin variants of the same target set (polish columns) sit within 1 kg of each other, so the variant that CONVERTS best
is worth picking, and the mass-free optimise of the export (40 iterations) may not be converged.

What.  Exactly the steps of run_ialns.w_export (to_exact -> restore -> mass-free optimise -> restore -> to_thrust ->
restore -> tighten -> RK4 hstep 0.1/0.05/0.025 d restore + tighten + write_craft validation), with the mass-free optimise
iteration count as a parameter (--opt-iters, 40 = the export).  Every route file ->
OUT/crafts/craft_<key>.npz + OUT/frags/frag_<key>.txt (INVALID_frag_<key>.txt when the validator rejects it) and one line
in OUT/index.jsonl {key, src, n, tank, m0, miss, ok, hs, opt_iters, n_opt, sec}; key = s16_iter.chash + '_o<iters>'.
Resumable: keys already in the index are skipped.  Never writes results/CTOC14_Result_*.txt.

usage: s16_xconv.py OUT FILE_OR_SPEC [...] [--opt-iters 40] [--min-rel-gain 1e-4] [--nproc 10] [--timeout 3600]
FILE_OR_SPEC = route .npz file | fleet dir | FRONTIER_JSON:KEY (tools/s15b_fleet.fleet_files)."""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, signal, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)
import numpy as np
import run_ialns as RI
import s15b_fleet as FA
from s16_iter import chash
from ctoc14.globalopt import optimise, restore, problem_from_state, state_of, to_thrust, tighten, write_craft, save_state
from ctoc14.impulsive import to_exact
from ctoc14.constants import DAY

Q = lambda s: None


class _TO(Exception):
    pass


def _alarm(signum, frame):
    raise _TO()


def convert(st, out, key, opt_iters=40, mrg=1e-4):
    """run_ialns.w_export with a free optimise iteration count; returns the index record (without src/key)."""
    out = pathlib.Path(out)
    ip = RI.ipr(st); nfb = len(ip.asts); gm = to_exact(ip, RI.eph())
    restore(gm, tol_km=300, iters=40, rho=1e3, log=Q)
    _, hist = optimise(gm, iters=opt_iters, rho0=100.0, log=Q, min_rel_gain=mrg); restore(gm, tol_km=150, iters=20, log=Q)
    gt = to_thrust(gm, 1.5); restore(gt, tol_km=150, iters=20, log=Q); tighten(gt, 1.5, log=Q)
    ok = False; hs_used = None; d = np.zeros(1)
    frag = out / 'frags' / f'frag_{key}.txt'
    for hs in (0.1, 0.05, 0.025):
        gt = problem_from_state(RI.eph(), state_of(gt), hstep=hs * DAY)
        restore(gt, tol_km=100, iters=30, log=Q); d = tighten(gt, 1.5, log=Q)
        rep, fuel, pk = write_craft(gt, str(frag), eph=RI.eph())
        if rep.ok and len(rep.flybys) >= nfb:
            ok = True; hs_used = hs; break
    save_state(state_of(gt), out / 'crafts' / f'craft_{key}.npz')
    if not ok:
        frag.rename(out / 'frags' / f'INVALID_frag_{key}.txt')
    return dict(n=nfb, tank=float(ip.tank()), m0=float(gt.m0), miss=float(np.max(d)) if len(d) else 0.0, ok=ok, hs=hs_used,
                n_opt=len(hist))


def job(args):
    f, key, opt_iters, timeout, out, mrg = args
    tic = time.time(); rec = dict(key=key, src=f, opt_iters=opt_iters, mrg=mrg)
    signal.signal(signal.SIGALRM, _alarm); signal.alarm(int(timeout))
    try:
        rec.update(convert(FA.load_st(f), out, key, opt_iters, mrg))
    except _TO:
        rec.update(ok=False, err='timeout')
    except Exception as e:
        rec.update(ok=False, err=repr(e)[:160])
    finally:
        signal.alarm(0)
    rec['sec'] = round(time.time() - tic)
    return rec


def expand(specs):
    files = []
    for s in specs:
        if s.endswith('.npz') and os.path.isfile(s):
            files.append(s)
        else:
            files += list(FA.fleet_files(s).values())
    return files


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('specs', nargs='+')
    ap.add_argument('--opt-iters', type=int, default=40); ap.add_argument('--nproc', type=int, default=10)
    ap.add_argument('--timeout', type=float, default=3600.0)
    ap.add_argument('--min-rel-gain', type=float, default=1e-4, help='optimise stall test (1e-4 = the export)')
    a = ap.parse_args()
    out = ROOT / a.out; (out / 'crafts').mkdir(parents=True, exist_ok=True); (out / 'frags').mkdir(parents=True, exist_ok=True)
    idx = out / 'index.jsonl'; done = set()
    if idx.exists():
        for l in open(idx):
            try:
                done.add(json.loads(l)['key'])
            except Exception:
                pass
    RI.eph()
    todo = []; seen = set()
    for f in expand(a.specs):
        key = f'{chash(FA.load_st(f))}_o{a.opt_iters}' + ('' if a.min_rel_gain == 1e-4 else f'g{a.min_rel_gain:.0e}')
        if key in done or key in seen:
            continue
        seen.add(key); todo.append((f, key, a.opt_iters, a.timeout, str(out), a.min_rel_gain))
    todo.sort(key=lambda j: -len(FA.load_st(j[0])['asts']))            # deepest (slowest) first
    print(f'[{time.strftime("%H:%M:%S")}] {len(todo)} conversions (opt_iters {a.opt_iters}) on {a.nproc} procs; {len(done)} in index', flush=True)
    with mp.get_context('fork').Pool(a.nproc, maxtasksperchild=1) as pool, open(idx, 'a') as fo:
        for r in pool.imap_unordered(job, todo, chunksize=1):
            fo.write(json.dumps(r) + '\n'); fo.flush()
            if 'm0' in r:
                print(f'[{time.strftime("%H:%M:%S")}] {r["key"]} {r["n"]} fb: tank {r["tank"]:.2f} -> m0 {r["m0"]:.2f} '
                      f'({r["m0"] - r["tank"]:+.2f} kg) miss {r["miss"]:.0f} km valid {r["ok"]} hs {r["hs"]} n_opt {r["n_opt"]} '
                      f'({r["sec"]} s) {r["src"]}', flush=True)
            else:
                print(f'[{time.strftime("%H:%M:%S")}] {r["key"]} FAILED {r.get("err")} ({r["sec"]} s) {r["src"]}', flush=True)
    print(f'[{time.strftime("%H:%M:%S")}] done', flush=True)


if __name__ == '__main__':
    main()
