"""R2 LNS quick check, part A (offline, seconds): score every saved suffix column of the two F1 jobs as an ejection-chain
link.  For a column of host h: dtank vs original route; dropped own targets priced for re-homing by census_h (best lin
host, excluding h and r9 which is dissolved; cap 400 kg if no approach); added r9 targets priced by census_r9 (their
forced alternative, excluding h).  net = dtank + rehome(dropped);  saving = alt(added) - net.
A link HELPS closure only if saving > 0 AND it adds a hard-core target at net <= 46 kg per hard-core target."""
import sys, json, glob, pathlib
sys.path.insert(0, '/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0, '/Users/mickey/solarsystem/ctoc14/tools')
import numpy as np
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
import run_ialns as RI
an = json.load(open(ROOT / 'results/s17/lns/anatomy.json'))
cs = {v: json.load(open(ROOT / f'results/s17/lns/census_{v}.json')) for v in ('r6', 'r8', 'r9')}
HARD = set(json.load(open(ROOT / 'results/s17/lns/census_summary.json'))['r9']['hard'])
def best(v, X, excl):
    c = [x for x in cs[v]['cands'] if x['ast'] == X and x['host'] not in excl]
    if not c: return (400.0, None, None)
    b = min(c, key=lambda x: x['lin_kg']); return (min(b['lin_kg'], 400.0), b['host'], b['dist'])
out = []
for tag, host, k in (('L1r6k17S', 'r6', 17), ('L1r8k13S', 'r8', 13)):
    own = set(json.load(open(ROOT / f'results/s17/lns/job_{tag}.json'))['own'])
    r9 = set(an['r9']['asts']); tank0 = an[host]['tank']
    for f in sorted(glob.glob(str(ROOT / f'results/s17/lns/pre/{tag}/route_*.npz'))):
        st = dict(np.load(f, allow_pickle=True)); a = [int(x) for x in st['asts']]
        tank = float(RI.ipr(st).tank())
        drop = sorted(own - set(a)); add = sorted(set(a) & r9)
        rh = {d: best(host, d, (host, 'r9')) for d in drop}
        al = {x: best('r9', x, (host,)) for x in add}
        net = tank - tank0 + sum(v[0] for v in rh.values()); alt = sum(v[0] for v in al.values())
        nh = len(set(add) & HARD)
        out.append(dict(col=pathlib.Path(f).name, host=host, n=len(a), tank=round(tank, 1), dtank=round(tank - tank0, 1),
                        add=add, hard_added=sorted(set(add) & HARD), drop=drop,
                        rehome={d: [round(v[0], 1), v[1], None if v[2] is None else round(v[2], 3)] for d, v in rh.items()},
                        alt={x: [round(v[0], 1), v[1]] for x, v in al.items()},
                        net=round(net, 1), saving=round(alt - net, 1),
                        net_per_hard=None if nh == 0 else round((net - sum(al[x][0] for x in add if x not in HARD)) / nh, 1)))
out.sort(key=lambda r: -r['saving'])
for r in out:
    if r['add'] or r['n'] >= (30 if r['host'] == 'r6' else 23):
        print(f"{r['col']:32s} n{r['n']} dtank {r['dtank']:+7.1f} add {r['add']} drop {r['drop']} rehome "
              f"{ {d: v[:2] for d, v in r['rehome'].items()} } net {r['net']:+.1f} alt {sum(v[0] for v in r['alt'].values()):.0f} "
              f"saving {r['saving']:+.1f} net/hard {r['net_per_hard']}")
json.dump(out, open(ROOT / 'results/s17/lns/r2/chain_screen.json', 'w'), indent=1)
