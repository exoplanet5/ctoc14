"""Gate G0 (docs/stage13_solver_plan.md section 6): does the PCC-8G build reproduce the tail-probe fleet T-N6a?
  (1) s13_pass on N6a's head (routes 01-05 kept) -> 273 +- 2 covered at sum J_i 10.48 +- 0.05
  (2) s13_master given N6a's own routes returns N6a.
usage: g0_check.py PASS_DIR [PASS_DIR ...]   (prints a JSON verdict per pass dir, re-integrating every route)"""
import os, sys, json, pathlib
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
os.chdir(ROOT)
import numpy as np
import run_ialns as RI


def fleet_numbers(d):
    fl = RI.IFleet(d); rows = []
    for n, r in sorted(fl.routes.items()):
        ip = RI.ipr(r['st']); Yf, _ = ip.integrate(); dm, _ = ip.misses(Yf)
        rows.append(dict(route=n, n=len(r['st']['asts']), tank=round(r['tank'], 1), miss_km=round(float(np.linalg.norm(dm, axis=1).max()), 1)))
    cov = fl.coverage(); sj = sum(RI.cost(r['tank']) for r in fl.routes.values())
    return dict(covered=len(cov), sumJi=round(sj, 4), J=round(sj + 298 - len(cov) + 2, 4), dup=sum(1 for h in cov.values() if len(h) > 1),
                routes=rows, max_miss_km=max(x['miss_km'] for x in rows))


ref = fleet_numbers(ROOT / 'results/s13/plan13/tail_N6a_k5/fleet')
print('T-N6a reference', json.dumps({k: v for k, v in ref.items() if k != 'routes'}), [(x['n'], x['tank']) for x in ref['routes']])
for d in sys.argv[1:]:
    x = fleet_numbers(ROOT / d)
    ok = abs(x['covered'] - 273) <= 2 and abs(x['sumJi'] - 10.48) <= 0.05 and x['max_miss_km'] <= 150 and x['dup'] == 0
    pj = json.load(open(ROOT / d / 'pass.json'))
    print(d, 'G0 pass-part', 'PASS' if ok else 'FAIL', json.dumps({k: v for k, v in x.items() if k != 'routes'}),
          [(r['n'], r['tank']) for r in x['routes']], 'wall', pj.get('summary', {}).get('wall_s'), 's')
m = json.load(open(ROOT / 'results/s13/g0/master_N6a/master.json'))
picks = sorted(p['f'] for p in m['N']['8']['picks'])
want = sorted(f'results/s13/seed_passes/N6a/route_{k:02d}.npz' for k in range(1, 9))
print('master on N6a own routes:', 'PASS' if picks == want else 'FAIL', m['N']['8']['covered'], m['N']['8']['sumJi'], picks == want)
