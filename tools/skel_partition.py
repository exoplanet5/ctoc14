"""Search the HARD-TARGET PARTITION (stage 6, the last mechanism standing).

Everything else was measured and killed on 2026-09-20 (docs/stage6_results.md): per-route capacity is ample
(35-41 flybys, sum 304 >= 298) but the eight routes want the same targets, and forcing a route onto a target it
did not choose costs ~2 free flybys.  Every fill therefore saturates near 240 of 298.

What decides which targets a route wants is its HARD chain: the 8-9 mandatory waypoints fix where and when the
craft has to be, and the easy sequence is whatever fits around them.  skel8b's chains came from t10d, where the
hard targets were INSERTED into finished routes -- nobody ever chose the partition.  drop-one showed how much it
matters: hard target 199 alone costs route 10 four flybys.

So: local search over "which craft flies which hard target", scored by the thing we actually need, the coverage of
a full sequential fill.  A move keeps the target's t10d epoch (a real crossing) and only changes the craft, subject
to --min-gap days from that craft's other waypoints.  Because the fill is sequential, a move that first touches
position p only invalidates routes p.. ; the driver restarts the fill from the cached snapshot of positions < p,
which roughly halves the cost of an evaluation.

Usage: skel_partition.py skel_dir out_dir [--iters 200] [--beam 300] [--nproc 8] [--anneal 1]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, shutil, pathlib, argparse, subprocess
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from ctoc14.search import UNREACHABLE
from ctoc14.constants import DAY
PY = sys.executable
ALL = [t for t in range(1, 301) if t not in UNREACHABLE]


def write_skel(assign, out):
    """assign: {route: [(t_epoch, target), ...]} -> a skeleton IFleet dir (waypoint ids + epochs only)."""
    out = pathlib.Path(out)
    if out.exists(): shutil.rmtree(out)
    out.mkdir(parents=True)
    fl = RI.IFleet()
    for n, wps in assign.items():
        wps = sorted(wps)
        st = dict(BASE_ST[n])
        st['asts'] = np.array([x for _, x in wps], float)
        st['tf'] = np.array([t for t, _ in wps], float)
        fl.routes[n] = dict(st=st, tank=700.0)
    fl.save(out, note='partition trial')
    return out


def run_fill(skel, out, a, only=None, start=None):
    cmd = [PY, str(ROOT / 'tools/skel_fill.py'), str(skel), str(out), '--mode', 'segmented',
           '--portfolio', f'{a.beam}:0.15:45', '--m0', '1600', '--vinf', '4.0', '--dvmax', '1.2',
           '--wp-dvmax', '2.5', '--seg-k', str(a.segk), '--seg-depth', '20', '--max-tank', str(a.max_tank),
           '--try', str(a.ntry), '--lam', str(a.lam), '--nproc', str(a.nproc), '--order', ','.join(a.order),
           '--snapshots']
    if only: cmd += ['--only', ','.join(only)]
    if start: cmd += ['--start', str(start)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=a.timeout)
    if r.returncode != 0:
        return None, r.stderr[-400:]
    try:
        return RI.IFleet(out), ''
    except Exception as e:
        return None, repr(e)


def score(fleet):
    cov = set(fleet.coverage())
    sumJ = sum(RI.cost(r['tank']) for r in fleet.routes.values())
    return len(cov), sumJ


def main():
    global BASE_ST
    ap = argparse.ArgumentParser(); ap.add_argument('src'); ap.add_argument('out')
    ap.add_argument('--iters', type=int, default=200); ap.add_argument('--beam', type=int, default=300)
    ap.add_argument('--segk', type=int, default=30); ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--lam', type=float, default=0.4); ap.add_argument('--max-tank', type=float, default=1150.0)
    ap.add_argument('--ntry', type=int, default=3); ap.add_argument('--timeout', type=float, default=2400.0)
    ap.add_argument('--min-gap', type=float, default=40.0, help='days between two waypoints of one craft')
    ap.add_argument('--anneal', type=int, default=1, help='accept a solution this many targets worse (coverage)')
    ap.add_argument('--seed', type=int, default=0); ap.add_argument('--max-hard', type=int, default=11)
    ap.add_argument('--swap-prob', type=float, default=0.5)
    a = ap.parse_args()
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True); say = RI.logger(out / 'log.txt')
    rng = np.random.default_rng(a.seed)
    skel = RI.IFleet(a.src)
    BASE_ST = {n: r['st'] for n, r in skel.routes.items()}
    assign = {n: [(float(t), int(x)) for x, t in zip(r['st']['asts'], r['st']['tf'])] for n, r in skel.routes.items()}
    a.order = sorted(assign, key=lambda n: -len(assign[n]))
    say(f'partition search from {a.src}: ' + ' '.join(f'{n}:{len(assign[n])}' for n in a.order))
    say(f'fill order {a.order}, beam {a.beam}, seg-k {a.segk}, lam {a.lam}, anneal {a.anneal}')

    t0 = time.time()
    sk = write_skel(assign, out / 'skel_cur')
    fleet, err = run_fill(sk, out / 'fill_cur', a)
    if fleet is None: say(f'baseline fill FAILED: {err}'); return
    cov, sumJ = score(fleet)
    best = (cov, -sumJ); best_assign = {n: list(v) for n, v in assign.items()}
    say(f'baseline: covered {cov}, sum J_i {sumJ:.3f}  ({time.time()-t0:.0f} s)')
    shutil.rmtree(out / 'best', ignore_errors=True); shutil.copytree(out / 'fill_cur', out / 'best')
    json.dump({n: sorted(x for _, x in v) for n, v in assign.items()}, open(out / 'best' / 'partition.json', 'w'), indent=1)

    hist = []
    for it in range(1, a.iters + 1):
        # per-route coverage of the CURRENT fill, to bias moves away from the starved routes
        percov = {n: len(set(int(x) for x in fleet.routes[n]['st']['asts'])) if n in fleet.routes else 0
                  for n in a.order}
        weak = sorted(a.order, key=lambda n: percov[n])
        pA = np.array([1.0 / (1 + percov[n]) for n in a.order]); pA /= pA.sum()
        trial = {n: list(v) for n, v in assign.items()}
        A = str(rng.choice(a.order, p=pA)); B = str(rng.choice([n for n in a.order if n != A]))
        swap = rng.random() < a.swap_prob and len(trial[B]) > 1
        if not trial[A]: continue
        i = int(rng.integers(len(trial[A]))); mv = trial[A][i]
        moves = [(A, B, mv)]
        if swap:
            j = int(rng.integers(len(trial[B]))); mv2 = trial[B][j]
            moves.append((B, A, mv2))
        elif len(trial[B]) >= a.max_hard:
            continue
        ok = True
        for src, dst, (te, x) in moves:
            others = [u for u, _ in trial[dst] if not any(abs(u - m[2][0]) < 1e-9 for m in moves if m[0] == dst)]
            if any(abs(te - u) < a.min_gap * DAY for u in others): ok = False; break
        if not ok: continue
        for src, dst, m in moves:
            trial[src] = [w for w in trial[src] if w != m]; trial[dst] = trial[dst] + [m]
        tag = ' + '.join(f'{m[1]}<-{m[2][1]}({m[0]})' for m in moves)
        # restart the sequential fill from the first position the move touched
        p = min(a.order.index(A), a.order.index(B))
        start = (out / 'fill_cur' / f'snap_{p-1:02d}') if p > 0 and (out / 'fill_cur' / f'snap_{p-1:02d}').exists() else None
        only = a.order[p:] if start else None
        tic = time.time()
        sk = write_skel(trial, out / 'skel_try')
        f2, err = run_fill(sk, out / 'fill_try', a, only=only, start=start)
        if f2 is None:
            say(f'[{it:3d}] {tag}: fill failed ({err[:80]})'); continue
        c2, j2 = score(f2)
        key = (c2, -j2)
        acc = key > best or (c2 >= best[0] - a.anneal and rng.random() < 0.25)
        say(f'[{it:3d}] {tag}: covered {c2} (best {best[0]}), sum J_i {j2:.3f}, from pos {p}'
            f'  {"ACCEPT" if acc else "reject"}  ({time.time()-tic:.0f} s)')
        hist.append(dict(it=it, move=tag, covered=c2, sumJ=j2, acc=bool(acc)))
        if acc:
            assign = trial; fleet = f2
            shutil.rmtree(out / 'fill_cur', ignore_errors=True); shutil.copytree(out / 'fill_try', out / 'fill_cur')
            if key > best:
                best = key
                shutil.rmtree(out / 'best', ignore_errors=True); shutil.copytree(out / 'fill_try', out / 'best')
                json.dump({n: sorted(x for _, x in v) for n, v in trial.items()},
                          open(out / 'best' / 'partition.json', 'w'), indent=1)
                say(f'      NEW BEST: covered {c2}, sum J_i {j2:.3f} -> J {j2 + 2 + (298 - c2):.3f}')
        json.dump(hist, open(out / 'history.json', 'w'), indent=1)
    say(f'done: best covered {best[0]}, sum J_i {-best[1]:.3f}  ({time.time()-t0:.0f} s)')


if __name__ == '__main__':
    main()
