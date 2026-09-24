"""Dual-price loop around tools/skel_fill.py (segmented mode): the leftovers of one full fill pass get their prize
multiplied by --up (capped at --cap) and the fill is re-run from the same skeletons, so the early routes -- which have
depth slack -- route through the targets the greedy tail keeps missing.  Keeps the pass with the fewest leftovers
(ties: lowest sum J_i) in out_dir/best.
Usage: dual_fill.py skel_dir out_dir [--passes 3] [--up 1.8] [--cap 4.0] [--nproc 8] [--lam 0.1]"""
import sys, json, shutil, pathlib, argparse, subprocess
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]
PY = str(pathlib.Path.home() / '.venvs/astro313/bin/python')
ap = argparse.ArgumentParser(); ap.add_argument('skel'); ap.add_argument('out')
ap.add_argument('--passes', type=int, default=3); ap.add_argument('--up', type=float, default=1.8)
ap.add_argument('--cap', type=float, default=4.0); ap.add_argument('--nproc', type=int, default=8)
ap.add_argument('--lam', type=float, default=0.1); ap.add_argument('--portfolio', default='300:0.15:45,300:0.15:90')
ap.add_argument('--start-prizes', default='', help='JSON prizes to start from (e.g. the leftovers of an earlier pass)')
a = ap.parse_args()
out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
pz = np.asarray(json.load(open(a.start_prizes)), float) if a.start_prizes else np.ones(300)
best = None
for k in range(a.passes):
    pf = out / f'prizes_{k}.json'; json.dump([float(x) for x in pz], open(pf, 'w'))
    d = out / f'pass{k}'
    cmd = [PY, str(ROOT / 'tools/skel_fill.py'), a.skel, str(d), '--mode', 'segmented', '--portfolio', a.portfolio,
           '--retry', '', '--nproc', str(a.nproc), '--lam', str(a.lam), '--easy-prizes', str(pf)]
    print(f'=== pass {k}: {int((pz > 1).sum())} priced targets (max {pz.max():.2f})', flush=True)
    subprocess.run(cmd, check=True)
    r = json.load(open(d / 'result.json')); left = r['leftovers']
    print(f'=== pass {k}: covered {r["covered"]}, leftovers {len(left)}, sum J_i {r["sum_Ji"]:.3f}, '
          f'J(free insert) {r["J_if_all_inserted_free"]:.3f}', flush=True)
    key = (len(left), r['sum_Ji'])
    if best is None or key < best[0]:
        best = (key, k)
        if (out / 'best').exists(): shutil.rmtree(out / 'best')
        shutil.copytree(d, out / 'best')
    for t in left: pz[t - 1] = min(a.cap, pz[t - 1] * a.up)
print(f'best: pass {best[1]} leftovers {best[0][0]} sum J_i {best[0][1]:.3f} -> {out}/best', flush=True)
