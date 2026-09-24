"""Pick, per spacecraft, the best available converted fragment (max distinct new targets, then lowest J_i), then combine.
Usage: assemble_best.py out.txt dir1 dir2 ...   (fragments frag_*_scN.txt in each dir; SC numbering per dir)"""
import sys, glob, re, pathlib, subprocess, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ctoc14.validator import parse
from ctoc14.constants import cost_sc
out = sys.argv[1]; dirs = sys.argv[2:]
cands = {}   # (dir, sc) -> list of (frag, targets, m0); a file argument forms its own group
import os
for d in dirs:
    files = [d] if os.path.isfile(d) else glob.glob(f'{d}/frag_*_sc*.txt')
    for f in files:
        m = re.search(r'frag_(\w+)_sc(\d+)\.txt$', f)
        key = (d, int(m.group(2))) if (m and not os.path.isfile(d)) else (f, 0)
        rows = parse(f); m0 = rows[0].m; targets = {r.ast for r in rows if r.event == 3}
        cands.setdefault(key, []).append((f, targets, m0))
# greedy: order candidates groups by best size, pick fragment maximizing (new targets - cost)
chosen = []; covered = set()
groups = sorted(cands.items(), key=lambda kv: -max(len(t) for _, t, _ in kv[1]))
for key, lst in groups:
    best = max(lst, key=lambda c: (len(c[1] - covered) - cost_sc(c[2]), -c[2]))
    gain = len(best[1] - covered) - cost_sc(best[2])
    if gain <= 0: print(f'skip {key}: best gain {gain:.2f}'); continue
    chosen.append(best); covered |= best[1]
    print(f'{key}: {best[0]}  m0={best[2]:.1f} J_i={cost_sc(best[2]):.3f} new={len(best[1]-covered|best[1])} targets={len(best[1])}')
J = sum(cost_sc(c[2]) for c in chosen) + 300 - len(covered)
print(f'assembled {len(chosen)} spacecraft, covered {len(covered)}, J = {J:.3f}; missing {sorted(set(range(1,301))-covered)}')
subprocess.run([sys.executable, 'tools/combine_submission.py', out] + [c[0] for c in chosen])
