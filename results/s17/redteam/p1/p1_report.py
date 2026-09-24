"""Print P1 frontier tables + the arc list of chosen solutions (markdown) from p1/*.json frontier files."""
import json, sys, pathlib
P1 = pathlib.Path('/Users/mickey/solarsystem/ctoc14/results/s17/redteam/p1')
R = []
for f in sys.argv[1:]:
    R += json.load(open(P1 / f))
R.sort(key=lambda r: (r['thr'], r['N'], r['K']))
print('| arcs thr | N | K | status | opt ΣJ_i (additive) | dual bound | convex ΣJ_i | cols | arcs used (>0.06) | arcs dJ | overlap | sec |')
print('|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|')
for r in R:
    if 'obj' in r:
        nu = len(r['arcs']); nun = sum(a['untrusted'] for a in r['arcs'])
        print(f"| {r['thr']} | {r['N']} | {r['K']} | {'opt' if r['status'] == 0 else r['message'][:18]} | {r['obj']:.4f} | {r['dual_bound']} | "
              f"{r['convex_sumJi']:.4f} | {len(r['columns'])} | {nu} ({nun}) | {r['arcs_dJ']:.3f} | {r['overlap']} | {r['sec']} |")
    else:
        print(f"| {r['thr']} | {r['N']} | {r['K']} | {r['message'][:30]} | — | {r['dual_bound']} | — | — | — | — | — | {r['sec']} |")
for r in R:
    if 'obj' in r and r['K'] == 298 and r['arcs']:
        print(f"\n**thr {r['thr']} N {r['N']} K 298**: columns")
        for c in r['columns']:
            print(f"- {c['key']} ({c['file']}) n {c['n']} tank {c['tank']} J_i {c['J_i']} +{c['add_kg']} kg arcs")
        print('arcs (host, ast, t_day, dist AU, lin kg, lin dJ, untrusted):')
        for a in r['arcs']:
            print(f"  {a['host']} {a['ast']} t{a['t_day']} d{a['dist']:.3f} {a['lin_kg']:.1f} kg dJ {a['lin_dJ']:.4f}{' UNTRUSTED' if a['untrusted'] else ''}")
