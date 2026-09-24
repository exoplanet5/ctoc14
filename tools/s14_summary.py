"""Stage 14 track A: aggregate the ceiling-test outputs under results/s14/A into one table + summary.json.

Reads   results/s14/A/free/<variant>_<route>/result.json   (s14_twinbeam free; variants prod, rel, lt, ...)
        results/s14/A/guided*/g<route>/result.json          (s14_twinbeam guided)
        results/s14/A/chaincheck.json                        (s14_chaincheck, optional)
Usage   s14_summary.py [--root results/s14/A]
"""
import sys, json, pathlib, argparse
from collections import defaultdict

ROUTES = ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09']
BASE = {'00': 28, '01': 32, '02': 26, '03': 28, '04': 25, '05': 24, '06': 23, '07': 25, '08': 22, '09': 19}  # best old builder


def src_legs(d):
    """t10d legs reproduced: consecutive pairs (a -> b) of the beam route that are consecutive in the source route
    with both epochs within 5 d of the source epochs."""
    so = d.get('src_order') or []
    if not so:
        return None
    st = d['src_epochs_d']; pos = {a: i for i, a in enumerate(so)}
    n = 0
    for (a, ta), (b, tb) in zip(zip(d['order'], d['epochs_d']), list(zip(d['order'], d['epochs_d']))[1:]):
        i = pos.get(a); j = pos.get(b)
        if i is not None and j == i + 1 and abs(ta - st[i]) <= 5 and abs(tb - st[j]) <= 5:
            n += 1
    return n

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--root', default='results/s14/A')
    a = ap.parse_args(); R = pathlib.Path(a.root)
    out = dict(free={}, guided={}, chain=None)
    var = defaultdict(dict)
    for f in sorted(R.glob('free/*/result.json')):
        d = json.load(open(f)); v = f.parent.name.rsplit('_', 1)[0]
        g = R / 'guided' / f'g{d["route"]}' / 'result.json'
        if g.exists():
            d['src_order'] = [st['ast'] for st in json.load(open(g))['steps']]
        d['src_legs'] = src_legs(d)
        var[v][d['route']] = d
    print('FREE twin-state beams (re-flown / pool, settled twin tank vs t10d, end, wall):')
    for v, rr in sorted(var.items()):
        tot = sum(d['reflown'] for d in rr.values()); npool = sum(d['pool_size'] for d in rr.values())
        out['free'][v] = dict(routes={r: dict(pool=d['pool_size'], reflown=d['reflown'], tank=d['tank'], src_tank=d['src_tank'],
                                              end=d['end'], wall_s=d['wall_s'], missing=d['missing'],
                                              pace=d.get('pace_d_per_fb'), src_pace=d.get('src_pace_d_per_fb'),
                                              src_legs_reproduced=d.get('src_legs'))
                                      for r, d in rr.items()}, reflown=tot, pool=npool)
        print(f'  {v}: {tot}/{npool} over {len(rr)} routes')
        for r in ROUTES:
            if r in rr:
                d = rr[r]
                print(f'    {r}: {d["reflown"]:2d}/{d["pool_size"]:2d} (old best {BASE[r]:2d}) tank {d["tank"]:7.1f} vs {d["src_tank"]:7.1f}  '
                      f'end {d["end"]["kind"]}@{d["end"].get("depth")} left {d["end"].get("time_left_d", "")} d  '
                      f'pace {d.get("pace_d_per_fb")} vs {d.get("src_pace_d_per_fb")} d/fb  t10d legs reproduced '
                      f'{d.get("src_legs")}  wall {d["wall_s"]:.0f} s')
    for gdir in sorted(R.glob('guided*')):
        rows = {}
        for r in ROUTES:
            f = gdir / f'g{r}' / 'result.json'
            if f.exists():
                rows[r] = json.load(open(f))
        if not rows:
            continue
        n = sum(d['legs'] for d in rows.values())
        agg = {k: sum(d[k] for d in rows.values()) for k in ('adm_prod', 'adm_rel', 'adm_full_prod', 'adm_full_rel', 'settled')}
        lt = [s for d in rows.values() for s in d['steps'][1:] if 'lt_price_src' in s]
        extra = {}
        if lt:
            near = [s for s in lt if s.get('lt_epoch_dt_d') is not None and abs(s['lt_epoch_dt_d']) <= 15]
            err = [s['lt_price_src'] - s['dv_marg'] for s in lt if s.get('dv_marg') is not None]
            extra = dict(lt_legs=len(lt), lt_epoch_within_15d=len(near),
                         lt_epoch_within_20d=sum(1 for s in lt if s.get('lt_epoch_dt_d') is not None and abs(s['lt_epoch_dt_d']) <= 20),
                         lt_price_src_le_3=sum(1 for s in lt if s['lt_price_src'] <= 3.0),
                         lt_price_src_le_2=sum(1 for s in lt if s['lt_price_src'] <= 2.0),
                         lt_price_src_le_1p2=sum(1 for s in lt if s['lt_price_src'] <= 1.2),
                         lt_price_err_median=round(sorted(err)[len(err) // 2], 4) if err else None,
                         lt_price_mae=round(sum(abs(e) for e in err) / len(err), 4) if err else None)
        out['guided'][gdir.name] = dict(legs=n, **agg, **extra)
        print(f'GUIDED {gdir.name}: {n} legs; myopic-twin admissible prod {agg["adm_prod"]} rel {agg["adm_rel"]}; '
              f'full-route prod {agg["adm_full_prod"]} rel {agg["adm_full_rel"]}; twin settled {agg["settled"]} {extra}')
    cf = R / 'chaincheck.json'
    if cf.exists():
        c = json.load(open(cf)); out['chain'] = c['total']
        print(f'PLANNER-STATE chain along t10d: {c["total"]}')
    json.dump(out, open(R / 'summary.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
