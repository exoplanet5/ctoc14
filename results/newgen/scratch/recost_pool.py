"""Gate F: true (impulsive-twin) cost of the planner-born pool columns.
For every column with >= --min-targets targets: rebuild the tour as a Lambert chain (junction impulses at the departure
of each leg, launch v_inf from the first leg), put the impulses on the 20-d node grid of ctoc14/impulsive.py, then the
standard settle (restore, L1 optimise, restore, enforce_cap). twin_from_tour below is the first (1 s lag) version, kept for
the record; the run uses ctoc14.impulsive.from_tour (1 d lag). Writes one JSON line per column.
Usage: recost_pool.py out.jsonl [--min-targets 35] [--max-cols 0] [--nproc 8] [--iters 100]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, glob, json, time, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
from ctoc14.colgen import CostModel, ColumnStore
from ctoc14.lambert import lambert_best
from ctoc14.impulsive import ImpulsiveProblem, from_tour, settle_tour
from ctoc14.constants import DAY, VINF_MAX, M_DRY, VE
import run_ialns as RI

DTAU = 20 * DAY
EXACT_NODES = True


def twin_from_tour(eph, tour):
    tL = float(tour['t_launch'])
    asts = [int(l['ast']) for l in tour['legs']]; tf = np.array([float(l['t_flyby']) for l in tour['legs']])
    rE, vE = eph.earth_state(tL)
    R, V = eph.ast_states_at(np.array(asts) - 1, tf)
    # Lambert chain: leg k from node k-1 (Earth or previous asteroid) to asteroid k
    r_prev = rE; v_prev = vE + np.asarray(tour['vinf'], float); t_prev = tL
    imp_t = []; imp_dv = []
    for k in range(len(asts)):
        tof = tf[k] - t_prev
        v1, v2, nrev, br = lambert_best(r_prev[None], R[k][None], np.array([tof]), v_prev[None], nrev_max=2)
        v1 = v1[0]; v2 = v2[0]
        if not np.all(np.isfinite(v1)):
            return None
        if k == 0:
            vinf = v1 - vE; nv = np.linalg.norm(vinf)
            if nv > VINF_MAX:
                imp_t.append(tL + 1.0); imp_dv.append(vinf * (1 - VINF_MAX / nv)); vinf *= VINF_MAX / nv
        else:
            imp_t.append(t_prev + 1.0); imp_dv.append(v1 - v_prev)
        r_prev = R[k]; v_prev = v2; t_prev = tf[k]
    t_end = tf.max()
    edges = np.arange(tL, t_end + DTAU, DTAU); ts = 0.5 * (edges[:-1] + edges[1:])
    Ts = np.zeros((len(ts), 3))
    if EXACT_NODES:      # junction impulses at their exact departure times (extra nodes): initial miss ~ 0
        ts = np.concatenate([ts, np.array(imp_t)]); Ts = np.concatenate([Ts, np.array(imp_dv).reshape(-1, 3)])
        o = np.argsort(ts); ts = ts[o]; Ts = Ts[o]
    else:
        for t, dv in zip(imp_t, imp_dv):
            b = int(np.clip(np.searchsorted(edges, t, side='right') - 1, 0, len(ts) - 1)); Ts[b] += dv
    return ImpulsiveProblem(eph, tL, vinf, ts, Ts, tf, asts)


def work(job):
    j, tour, iters = job
    tic = time.time()
    try:
        ip, miss, lag = settle_tour(RI.eph(), tour, lambda ip: RI.settle(ip, iters))
        if ip is None:
            return dict(col=j, ok=False, why=f'miss {miss:.0f} km after the lag ladder', secs=time.time() - tic)
        dv0 = float('nan'); tank0 = float('nan'); d0 = lag
        return dict(col=j, ok=bool(miss <= 150.0), miss=float(miss), dv=float(ip.dv()), tank=float(ip.tank()),
                    dv_lambert=dv0, tank_lambert=tank0, miss0_km=d0, n=len(ip.asts), secs=time.time() - tic,
                    st={k: (v.tolist() if hasattr(v, 'tolist') else v) for k, v in RI.ist(ip).items()})
    except Exception as e:  # noqa
        return dict(col=j, ok=False, why=repr(e)[:200], secs=time.time() - tic)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--min-targets', type=int, default=35)
    ap.add_argument('--max-cols', type=int, default=0); ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--iters', type=int, default=100); ap.add_argument('--no-state', action='store_true')
    a = ap.parse_args()
    cm = CostModel(s=0.70, reserve=2.0); store = ColumnStore(cm)
    for pat in ['results/phasing/camp*/pool.jsonl', 'results/phasing/pool/*.jsonl']:
        for f in sorted(glob.glob(str(ROOT / pat))): store.load_jsonl(f, tag=pathlib.Path(f).parent.name)
    for f in sorted(glob.glob(str(ROOT / 'results/colgen/run*/columns.jsonl'))) + sorted(glob.glob(str(ROOT / 'results/colgen/*/columns.jsonl'))):
        store.load_columns_file(f)
    n = np.array(store.n); c = np.array(store.cost, float)
    sel = np.where((n >= a.min_targets) & (c < 9))[0]
    sel = sel[np.argsort(c[sel] / n[sel])]          # best planner ratio first
    if a.max_cols:
        sel = sel[:a.max_cols]
    print(f'{len(store)} columns; recosting {len(sel)} with >= {a.min_targets} targets on {a.nproc} procs', flush=True)
    RI.eph()
    jobs = [(int(j), store.tour(int(j)), a.iters) for j in sel]
    out = open(a.out, 'a'); tic = time.time(); k = 0
    with mp.get_context('fork').Pool(a.nproc) as pool:
        for r in pool.imap_unordered(work, jobs):
            j = r['col']; r['n_targets'] = int(n[j]); r['cost_planner'] = float(c[j]); r['fuel_planner'] = float(store.fuel[j])
            r['targets'] = store.targets(j); r['tag'] = store.tag[j]
            if a.no_state: r.pop('st', None)
            out.write(json.dumps(r) + '\n'); out.flush(); k += 1
            if r['ok']:
                print(f'[{k:5d}/{len(sel)}] col {j:6d} n {r["n"]:2d}: planner {r["cost_planner"]:.3f} (fuel {r["fuel_planner"]:.0f}) '
                      f'-> twin tank {r["tank"]:.1f} J_i {RI.cost(r["tank"]):.3f} dv {r["dv"]:.2f} ({r["dv"]/r["n"]:.3f}/flyby) '
                      f'lag {r["miss0_km"]/86400:.2f} d  {r["secs"]:.0f} s  [{(time.time()-tic)/60:.1f} min]', flush=True)
            else:
                print(f'[{k:5d}/{len(sel)}] col {j:6d} n {r["n_targets"]:2d}: FAIL {r.get("why", "miss %.0f km" % r.get("miss", -1))}  {r["secs"]:.0f} s', flush=True)


if __name__ == '__main__':
    main()
