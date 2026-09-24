"""Part-A link scoring applied to the part-B (open-pool) columns of L2r6k17O. Same rules as chain_screen.py; targets taken
from r5/r7/r8 are counted at 0 kg (their removal refund is 1.7-3.2 kg median, REDTEAM 1.1)."""
import sys, json, glob, pathlib
sys.path.insert(0, '/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0, '/Users/mickey/solarsystem/ctoc14/tools')
import numpy as np
import run_ialns as RI
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
an = json.load(open(ROOT / 'results/s17/lns/anatomy.json'))
cs = {v: json.load(open(ROOT / f'results/s17/lns/census_{v}.json')) for v in ('r6', 'r9')}
HARD = set(json.load(open(ROOT / 'results/s17/lns/census_summary.json'))['r9']['hard'])
def best(v, X, excl):
    c = [x for x in cs[v]['cands'] if x['ast'] == X and x['host'] not in excl]
    return (min(min(x['lin_kg'] for x in c), 400.0) if c else 400.0)
own = set(json.load(open(ROOT / 'results/s17/lns/job_L1r6k17S.json'))['own']); r9 = set(an['r9']['asts']); tank0 = an['r6']['tank']
out = []
for f in sorted(glob.glob(str(ROOT / 'results/s17/lns/r2/L2r6k17O/route_*.npz'))):
    st = dict(np.load(f, allow_pickle=True)); a = [int(x) for x in st['asts']]
    if max(st['tf']) / 86400 < 5200: continue          # only (near-)complete suffixes
    tank = float(RI.ipr(st).tank()); drop = sorted(own - set(a)); add = sorted(set(a) & r9)
    rh = sum(best('r6', d, ('r6', 'r9')) for d in drop); alt = sum(best('r9', x, ('r6',)) for x in add)
    out.append(dict(col=pathlib.Path(f).name, n=len(a), tank=round(tank, 1), dtank=round(tank - tank0, 1), add=add,
                    hard=sorted(set(add) & HARD), drop=drop, rehome=round(rh), alt=round(alt), saving=round(alt - (tank - tank0 + rh))))
for r in sorted(out, key=lambda r: -r['saving']):
    print(r)
json.dump(out, open(ROOT / 'results/s17/lns/r2/chain_screen_B.json', 'w'), indent=1)
