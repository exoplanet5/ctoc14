"""Stage 14 track A stage 2: tables for the look-back planner (s14_lookback.py).

  guided  results/s14/A/s2/guided/gRR_L{L}: admissibility of t10d's own legs from the planner's own label, per L
  refly   results/s14/A/s2/<dir>/rRR: free re-fly of t10d's routes over their own targets (pool = own)
  pools   results/s14/A/s2/<dir>/pRR: the linchpin 10 -> 8 pools (fleet covered, sum J_i, J)
Usage: s14_s2summary.py guided | refly DIR | pools DIR   (prints a table and writes DIR/summary.json)
"""
import sys, json, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np

S2 = ROOT / 'results/s14/A/s2'
T10 = {'00': 918.2, '01': 904.4, '02': 903.7, '03': 795.3, '04': 886.7, '05': 976.3, '06': 791.9, '07': 815.1,
       '08': 870.6, '09': 876.5}
POOL10 = {'00': 37, '01': 36, '02': 33, '03': 30, '04': 30, '05': 28, '06': 27, '07': 27, '08': 26, '09': 24}
S1_LT20 = {'00': 35, '01': 33, '02': 28, '03': 30, '04': 30, '05': 28, '06': 25, '07': 24, '08': 23, '09': 21}
S1_BEST = {'00': 35, '01': 34, '02': 28, '03': 30, '04': 30, '05': 28, '06': 26, '07': 24, '08': 24, '09': 24}


def Ji(tank):
    x = (tank - 600.0) / 1400.0
    return 1 + x + x * x


def guided():
    out = {}
    for L in (1, 2, 3, 4):
        rows = []
        for r in sorted(T10):
            f = S2 / f'guided/g{r}_L{L}/result.json'
            if not f.exists():
                continue
            d = json.load(open(f))
            rows.append(dict(route=r, legs=d['legs'], **{k: d[k] for k in d if k.startswith('adm_')},
                             fails=len(d['fails']), planned_tank=d['planned_tank'],
                             settled=d.get('settle', {}).get('tank'), src=d['src_tank'], wall=d['wall_s'],
                             settle_s=d.get('settle', {}).get('dt_s')))
        if not rows:
            continue
        tot = {k: sum(r[k] for r in rows) for k in rows[0] if k.startswith('adm_') or k in ('legs', 'fails')}
        tot['wall'] = round(sum(r['wall'] for r in rows), 1)
        tot['settle_s'] = round(sum(r['settle_s'] or 0 for r in rows), 1)
        tot['settled_sum'] = round(sum(r['settled'] or 0 for r in rows), 1)
        tot['src_sum'] = round(sum(r['src'] for r in rows), 1)
        tot['planned_sum'] = round(sum(r['planned_tank'] for r in rows), 1)
        tot['n_routes'] = len(rows)
        out[f'L{L}'] = dict(rows=rows, total=tot)
        print(f'--- L = {L}: {len(rows)} routes')
        print(' route legs | src-epoch <=1.2 2.0 3.0 | best+-15d <=1.2 2.0 3.0 | fails | planned kg | settled kg | t10d kg | wall s')
        for x in rows:
            print(f'  {x["route"]}   {x["legs"]:3d} |  {x["adm_src_1.2"]:3d} {x["adm_src_2.0"]:3d} {x["adm_src_3.0"]:3d}  |  '
                  f'{x["adm_best_1.2"]:3d} {x["adm_best_2.0"]:3d} {x["adm_best_3.0"]:3d}  | {x["fails"]:3d}  | {x["planned_tank"]:7.1f} | '
                  f'{(x["settled"] or float("nan")):7.1f} | {x["src"]:6.1f} | {x["wall"]:5.1f}')
        print(f'  tot  {tot["legs"]:3d} |  {tot["adm_src_1.2"]:3d} {tot["adm_src_2.0"]:3d} {tot["adm_src_3.0"]:3d}  |  '
              f'{tot["adm_best_1.2"]:3d} {tot["adm_best_2.0"]:3d} {tot["adm_best_3.0"]:3d}  | {tot["fails"]:3d}  | {tot["planned_sum"]:7.1f} | '
              f'{tot["settled_sum"]:7.1f} | {tot["src_sum"]:6.1f} | {tot["wall"]:5.1f} (+ settles {tot["settle_s"]} s)')
    json.dump(out, open(S2 / 'guided/summary.json', 'w'), indent=1)


def refly(d):
    d = S2 / d; rows = []
    for r in sorted(T10):
        f = d / f'r{r}/result.json'
        if not f.exists():
            continue
        x = json.load(open(f))
        if 'wall_s' not in x:                       # planning done, final settle still running
            continue
        s = x.get('settled') or {}
        rows.append(dict(route=r, pool=x['pool_size'], planned=x['planned_n'], planned_tank=x['planned_tank'],
                         settled=s.get('n'), tank=s.get('tank'), missing=s.get('missing'), t10d=T10[r], s1_lt20=S1_LT20[r],
                         s1_best=S1_BEST[r], plan_wall=x['plan_wall_s'], wall=x['wall_s'], end=x['end']))
    print(' route pool | planned (kg)  | settled (twin kg) | t10d kg | stage1 lt20 / best | plan s  total s | missing')
    for x in rows:
        print(f'  {x["route"]}   {x["pool"]:3d} | {x["planned"]:3d} ({x["planned_tank"]:6.1f}) | {x["settled"] or 0:3d} '
              f'({(x["tank"] or float("nan")):6.1f})    | {x["t10d"]:6.1f} | {x["s1_lt20"]:3d} / {x["s1_best"]:3d} | '
              f'{x["plan_wall"]:6.0f} {x["wall"]:7.0f} | {x["missing"]}')
    tot = dict(pool=sum(x['pool'] for x in rows), planned=sum(x['planned'] for x in rows),
               settled=sum(x['settled'] or 0 for x in rows), plan_cpu_s=round(sum(x['plan_wall'] for x in rows), 1),
               cpu_s=round(sum(x['wall'] for x in rows), 1),
               sum_Ji=round(sum(Ji(x['tank']) for x in rows if x['tank']), 4), n=len(rows))
    print(f'  TOTAL pool {tot["pool"]}: planned {tot["planned"]}, settled {tot["settled"]} ({tot["settled"] / max(tot["pool"], 1):.1%}); '
          f'sum J_i {tot["sum_Ji"]}; planning CPU {tot["plan_cpu_s"] / 3600:.2f} h, total CPU {tot["cpu_s"] / 3600:.2f} h')
    json.dump(dict(rows=rows, total=tot), open(d / 'summary.json', 'w'), indent=1)


def pools(d, pools_file=ROOT / 'results/s13/linchpin/pools_A_c38.json'):
    d = S2 / d; P = json.load(open(pools_file))['pools']; rows = []; cov = set()
    for r in sorted(P):
        f = d / f'p{r}/result.json'
        if not f.exists():
            continue
        x = json.load(open(f)); s = x.get('settled') or {}
        if 'wall_s' not in x:
            continue
        if s:
            cov |= set(s['targets'])
        rows.append(dict(route=r, pool=len(P[r]), planned=x['planned_n'], planned_tank=x['planned_tank'], flown=s.get('n'),
                         tank=s.get('tank'), kg_fb=s.get('kg_per_fb'), J_i=s.get('J_i'), wall=x['wall_s']))
    print(' route pool | planned (kg)  | flown  tank  kg/fb  J_i | wall s')
    for x in rows:
        print(f'  {x["route"]}   {x["pool"]:3d} | {x["planned"]:3d} ({x["planned_tank"]:6.1f}) | {x["flown"] or 0:3d} '
              f'{(x["tank"] or float("nan")):7.1f} {(x["kg_fb"] or float("nan")):5.2f} {(x["J_i"] or float("nan")):.4f} | {x["wall"]:6.0f}')
    sJ = sum(x['J_i'] for x in rows if x['J_i']); n = len(cov)
    fl = [x['flown'] for x in rows if x['flown']]
    tot = dict(routes=len(rows), covered=n, sum_Ji=round(sJ, 4), J=round(sJ + 300 - n, 4),
               mean_flown=round(float(np.mean(fl)), 2) if fl else None,
               kg_per_fb=round(sum((x['tank'] - 600) for x in rows if x['tank']) / max(sum(fl), 1), 3),
               cpu_s=round(sum(x['wall'] for x in rows), 1))
    print(f'  FLEET: {tot["routes"]} routes, covered {n}/298, sum J_i {tot["sum_Ji"]}, J {tot["J"]}; mean flown {tot["mean_flown"]}, '
          f'{tot["kg_per_fb"]} kg/fb; CPU {tot["cpu_s"] / 3600:.2f} h')
    json.dump(dict(rows=rows, total=tot), open(d / 'summary.json', 'w'), indent=1)


if __name__ == '__main__':
    if sys.argv[1] == 'guided':
        guided()
    elif sys.argv[1] == 'refly':
        refly(sys.argv[2])
    else:
        pools(sys.argv[2])
