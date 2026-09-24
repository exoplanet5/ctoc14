"""Set-cover MILP over CONVERTED fragments: min sum_j J_j y_j + sum_a (1 - z_a), z_a <= sum_{j covers a} y_j.
Duplicated flybys across selected spacecraft are allowed (only the first detection counts, no penalty).
Usage: select_frags_milp.py out.txt dir_or_file ...   (fragments frag_*_sc*.txt / *_frag.txt)"""
import sys, glob, os, pathlib, subprocess, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix
from ctoc14.validator import parse
from ctoc14.constants import cost_sc
dry = '--dry' in sys.argv
if dry: sys.argv.remove('--dry')                 # selection only: no combined file, no validation
out = sys.argv[1]; args = sys.argv[2:]
files = []
for a in args:
    files += [a] if os.path.isfile(a) else sorted(glob.glob(f'{a}/frag_*_sc*.txt'))
frags = []
for f in files:
    rows = parse(f); S = {r.ast for r in rows if r.event == 3}
    if S: frags.append((f, S, rows[0].m))
n = len(frags); print(f'{n} fragments')
# variables: y_j (n), z_a (300)
c = np.concatenate([np.array([cost_sc(m0) for _, _, m0 in frags]), -np.ones(300)])   # minimise sum J_j y_j - sum z_a  (+300 const)
A = lil_matrix((300, n + 300))
for j, (_, S, _) in enumerate(frags):
    for a in S: A[a - 1, j] = -1.0
for a in range(300): A[a, n + a] = 1.0
cons = LinearConstraint(A.tocsr(), -np.inf, 0.0)           # z_a - sum y_j <= 0
res = milp(c, constraints=cons, integrality=np.concatenate([np.ones(n), np.zeros(300)]), bounds=Bounds(0, 1))
y = np.round(res.x[:n]).astype(int); sel = [j for j in range(n) if y[j]]
cov = set().union(*[frags[j][1] for j in sel]) if sel else set()
J = sum(cost_sc(frags[j][2]) for j in sel) + 300 - len(cov)
for j in sel: print(f'  {frags[j][0]}  m0={frags[j][2]:.0f} J_i={cost_sc(frags[j][2]):.3f} targets={len(frags[j][1])}')
print(f'selected {len(sel)} spacecraft, covered {len(cov)}, J = {J:.3f}; missing {sorted(set(range(1,301))-cov)}', flush=True)
if dry:
    open(out + '.selection', 'w').write('\n'.join(frags[j][0] for j in sel) + '\n')
else:
    subprocess.run([sys.executable, 'tools/combine_submission.py', out] + [frags[j][0] for j in sel])
