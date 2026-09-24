"""Column generation over the skeleton fill: LP relaxation of the one-column-per-skeleton cover, coverage duals as
the prizes of the next skel_fill --columns round (the beam's value cost - sum prize IS the reduced cost), repeat;
finally the MIP (skel_select.py) over all columns.
Usage: skel_cg.py skel_dir columns.jsonl out_dir [--rounds 3] [--nproc 8] [--portfolio 300:0.15:45,300:0.15:90]"""
import sys, json, subprocess, pathlib, argparse
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from ctoc14.search import UNREACHABLE
PY = str(pathlib.Path.home() / '.venvs/astro313/bin/python')
M_DRY = 600.0
cost = lambda tank: 1.0 + (tank - M_DRY) / 1400.0 + ((tank - M_DRY) / 1400.0) ** 2
ALL = [t for t in range(1, 301) if t not in UNREACHABLE]; tid = {t: i for i, t in enumerate(ALL)}


def lp_duals(cols, miss=1.0):
    skels = sorted(set(c['skel'] for c in cols)); nc, nt = len(cols), len(ALL)
    c_obj = np.concatenate([[cost(c['tank']) for c in cols], np.full(nt, -miss)])
    A_ub = lil_matrix((nt, nc + nt)); A_eq = lil_matrix((len(skels), nc + nt))
    for t in ALL:
        A_ub[tid[t], nc + tid[t]] = 1.0
    for j, c in enumerate(cols):
        A_eq[skels.index(c['skel']), j] = 1.0
        for t in c['targets']: A_ub[tid[t], j] = -1.0
    res = linprog(c_obj, A_ub=A_ub.tocsr(), b_ub=np.zeros(nt), A_eq=A_eq.tocsr(), b_eq=np.ones(len(skels)),
                  bounds=[(0, 1)] * (nc + nt), method='highs')
    pi = -np.asarray(res.ineqlin.marginals)              # value of covering target t (0..miss)
    x = res.x[:nc]; y = res.x[nc:]
    return res.fun, pi, x, y, skels


ap = argparse.ArgumentParser(); ap.add_argument('skel'); ap.add_argument('cols'); ap.add_argument('out')
ap.add_argument('--rounds', type=int, default=3); ap.add_argument('--nproc', type=int, default=8)
ap.add_argument('--portfolio', default='300:0.15:45,300:0.15:90'); ap.add_argument('--floor', type=float, default=0.05)
ap.add_argument('--ncol', type=int, default=1)
a = ap.parse_args()
out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
for r in range(a.rounds):
    cols = [json.loads(l) for l in open(a.cols)]
    obj, pi, x, y, skels = lp_duals(cols)
    cov_lp = float(y.sum())
    print(f'=== CG round {r}: {len(cols)} columns, LP objective {obj:.3f} (J ~ {obj + 2 + len(ALL):.3f}), LP coverage {cov_lp:.1f}, '
          f'duals: {int((pi > 0.5).sum())} targets > 0.5, {int((pi < 0.05).sum())} < 0.05', flush=True)
    EP = np.zeros(300)
    for t in ALL: EP[t - 1] = max(a.floor, min(1.0, pi[tid[t]]))
    pf = out / f'prizes_r{r}.json'; json.dump([float(v) for v in EP], open(pf, 'w'))
    cmd = [PY, str(ROOT / 'tools/skel_fill.py'), a.skel, str(out / f'gen_r{r}'), '--mode', 'segmented', '--columns', a.cols,
           '--ncol', str(a.ncol), '--portfolio', a.portfolio, '--retry', '', '--nproc', str(a.nproc), '--easy-prizes', str(pf),
           '--col-min', '20']
    subprocess.run(cmd, check=True)
cols = [json.loads(l) for l in open(a.cols)]
obj, pi, x, y, skels = lp_duals(cols)
print(f'=== after {a.rounds} rounds: {len(cols)} columns, LP objective {obj:.3f}, LP coverage {float(y.sum()):.1f}', flush=True)
subprocess.run([PY, str(ROOT / 'tools/skel_select.py'), a.cols, str(out / 'chosen.json'), '--time', '600'], check=True)
