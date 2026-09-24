"""R2 part B analysis: per depth n of the open-pool job L2r6k17O -> min tank, min t_end, composition of the deepest
states (own / other-tail / r9); the pre-registered metric n* (max n with tank <= 1065 kg) and the depth-26 projection."""
import sys, json, pickle, pathlib
sys.path.insert(0, '/Users/mickey/solarsystem/ctoc14'); sys.path.insert(0, '/Users/mickey/solarsystem/ctoc14/tools')
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
an = json.load(open(ROOT / 'results/s17/lns/anatomy.json'))
J = json.load(open(ROOT / 'results/s17/lns/r2/job_L2r6k17O.json'))
own = set(J['own']); r9 = set(an['r9']['asts']); other = set(J['pool']) - own - r9
HARD = set(json.load(open(ROOT / 'results/s17/lns/census_summary.json'))['r9']['hard'])
C = pickle.load(open(ROOT / J['out'] / 's1' / 'ckpt.pkl', 'rb'))
states = [ts for f in C['front'].values() for (v, ts) in f.values()] + list(C['beam'])
rows = {}
for ts in states:
    a = set(int(x) for x in ts['asts']); n = len(a)
    r = dict(n=n, tank=round(float(ts['tank']), 1), t_end=round(float(ts['t_end']) / 86400), own=len(a & own),
             other=sorted(a & other), r9=sorted(a & r9), hard=sorted(a & r9 & HARD), own_missing=sorted(own - a))
    rows.setdefault(n, []).append(r)
out = dict(level=C['level'], end=C['end'], settles=C['settles'], settle_ok=C['settle_ok'], by_depth={})
for n in sorted(rows):
    L = rows[n]
    out['by_depth'][n] = dict(min_tank=min(r['tank'] for r in L), min_t_end=min(r['t_end'] for r in L),
                              max_r9=max(len(r['r9']) for r in L), best=sorted(L, key=lambda r: r['tank'])[:3])
    b = out['by_depth'][n]
    print(f"n {n}: states {len(L)} min tank {b['min_tank']} min t_end {b['min_t_end']} max r9 {b['max_r9']}; cheapest: "
          f"own {b['best'][0]['own']}/13 other {b['best'][0]['other']} r9 {b['best'][0]['r9']}")
ok = [n for n in rows if any(r['tank'] <= 1065 for r in rows[n])]
out['n_star'] = max(ok) if ok else None
out['t26'] = out['by_depth'].get(26, {}).get('min_t_end')
print('n* =', out['n_star'], ' min t_end at depth 26 =', out['t26'], ' end', C['end'], ' level', C['level'])
json.dump(out, open(ROOT / 'results/s17/lns/r2/open_r6.json', 'w'), indent=1)
