"""Subgradient (column-generation) loop around tools/greedy_cover.py.

A beam with a uniform prize of 1 J per target is INDIFFERENT between targets, so it takes the cheapest ones and
leaves the orphans: after seven greedy routes the 67 remaining targets chain into a 6-flyby route, i.e. they can
only be covered woven INTO a deep route, never on their own. So raise the price of the targets the fleet misses and
run the greedy again: pi_t <- pi_t * up for every uncovered t (and decay towards 1 for the covered ones), which is the
Lagrangian dual of the set-cover. Every pass is scored with the TRUE objective (sum cost + misses) and the best fleet
is kept.
Usage: dual_greedy.py out_dir [--iters 6] [--n 8] [--m0 1600] [--up 1.7] [--nproc 8]
"""
import sys, json, time, shutil, pathlib, argparse, subprocess
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]
PY = str(pathlib.Path.home() / '.venvs/astro313/bin/python')

ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--iters', type=int, default=6)
ap.add_argument('--n', type=int, default=8); ap.add_argument('--m0', type=float, default=1600.0)
ap.add_argument('--up', type=float, default=1.7); ap.add_argument('--decay', type=float, default=0.85)
ap.add_argument('--nproc', type=int, default=8); ap.add_argument('--beam', type=int, default=300)
ap.add_argument('--max-depth', type=int, default=120); ap.add_argument('--columns', default='')
ap.add_argument('--scan', type=int, default=40); ap.add_argument('--settle-timeout', type=float, default=180.0)
ap.add_argument('--settle-budget', type=float, default=600.0)
ap.add_argument('--cap', type=float, default=6.0, help='maximum prize')
ap.add_argument('--vinf', type=float, default=2.0); ap.add_argument('--dvmax', type=float, default=1.2)
ap.add_argument('--tofmax', type=float, default=400.0); ap.add_argument('--lintofmax', type=float, default=600.0)
a = ap.parse_args()
out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
log = open(out / 'dual.log', 'a')
def say(s):
    line = f'[{time.strftime("%H:%M:%S")}] {s}'; print(line, flush=True); log.write(line + '\n'); log.flush()

pz = np.ones(300)
best = (1e9, None)
for it in range(a.iters):
    pf = out / f'prizes_{it:02d}.json'; json.dump(list(pz), open(pf, 'w'))
    d = out / f'it{it:02d}'
    if d.exists(): shutil.rmtree(d)
    cmd = [PY, str(ROOT / 'tools/greedy_cover.py'), str(d), '--n', str(a.n), '--m0', str(a.m0),
           '--beam', str(a.beam), '--nproc', str(a.nproc), '--prizes', str(pf), '--max-depth', str(a.max_depth),
           '--vinf', str(a.vinf), '--dvmax', str(a.dvmax), '--tofmax', str(a.tofmax), '--lintofmax', str(a.lintofmax),
           '--scan', str(a.scan), '--settle-timeout', str(a.settle_timeout), '--settle-budget', str(a.settle_budget)]
    if a.columns: cmd += ['--columns', a.columns]
    tic = time.time()
    subprocess.run(cmd, stdout=open(d.with_suffix('.out'), 'w'), stderr=subprocess.STDOUT)
    rf = d / 'result.json'
    if not rf.exists():
        say(f'iter {it}: greedy pass FAILED (no result.json, see {d}.out); keeping prizes and going on')
        continue
    r = json.load(open(rf))
    J = float(r['J']); nun = len(r['uncovered'])
    sizes = {n: len(v) for n, v in r['routes'].items()}
    say(f'iter {it}: J {J:.4f}, covered {len(r["covered"])}, uncovered {nun}, routes '
        + ' '.join(f'{n}:{sizes[n]}@{r["tanks"][n]:.0f}' for n in sorted(sizes)) + f'  ({time.time()-tic:.0f} s)')
    if J < best[0]:
        best = (J, it); say(f'   new best fleet: it{it:02d}')
    # subgradient step
    unc = set(r['uncovered'])
    for t in range(1, 301):
        if t in unc: pz[t - 1] = min(a.cap, pz[t - 1] * a.up)
        else: pz[t - 1] = max(1.0, 1.0 + (pz[t - 1] - 1.0) * a.decay)
    say(f'   prizes: {int((pz > 1.01).sum())} above 1, max {pz.max():.2f}')
say(f'best: iter {best[1]} with J {best[0]:.4f} -> {out}/it{best[1]:02d}')
