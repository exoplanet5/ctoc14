"""Per-generation trajectory of the s13 generation loop (reads results/s13/gen/{gen.json, gen2.json, gNN/...}).

usage: trajectory.py [GEN_DIR]   -> prints the table and writes GEN_DIR/trajectory.json"""
import sys, json, pathlib, collections
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
D = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / 'results/s13/gen'
if not D.is_absolute():
    D = ROOT / D


def src_of(f):
    """route file -> short provenance label (slot of this loop, or archive family)."""
    p = pathlib.Path(f)
    s = str(p.relative_to(ROOT)) if str(p).startswith(str(ROOT)) else str(p)
    parts = s.split('/')
    if s.startswith('results/s13/gen/g') and len(parts) > 4:
        return f'{parts[3]}/{parts[4]}' + ('(lib)' if 'lib' in parts else '')
    if s.startswith('results/s13/g0/'):
        return 'g0/' + parts[3]
    if s.startswith('results/s13/plan13/'):
        return 'plan13/' + parts[3]
    if s.startswith('results/s13/seed_passes/'):
        return 'seed/' + parts[3]
    return s


def main():
    G = json.load(open(D / 'gen.json')); S = json.load(open(D / 'gen2.json')) if (D / 'gen2.json').exists() else {}
    out = dict(best0=G.get('best0'), gens={}, stopped=G.get('stopped'), stopped2=S.get('stopped'))
    print(f'best0 (seed master, frontier 10.75): {G.get("best0")}')
    for key in sorted(G['gens']):
        rec = G['gens'][key]
        row = dict(slots={})
        for name in rec.get('slots', []):
            d = D / f'g{key}' / name.split('/')[-1]
            pj = d / 'pass.json'
            if pj.exists():
                P = json.load(open(pj)); sm = P.get('summary')
                row['slots'][name] = dict(done=P.get('done'), choice=P['args'].get('choice'), summary=sm,
                                          routes={k: (r.get('member'), r.get('n'), r.get('twin'), r.get('pool'))
                                                  for k, r in P.get('routes', {}).items()})
        if 'master' in rec:
            m = rec['master']; row['master'] = m
            md = ROOT / m['dir']
            prov = {}
            for sub in sorted(p.name for p in md.iterdir() if p.is_dir()):
                pk = md / sub / 'picked.json'
                if pk.exists():
                    prov[sub] = [(src_of(x['f']), x['n'], round(x['tank'])) for x in json.load(open(pk))]
            row['picked'] = prov
        if 'probe' in rec:
            row['probe'] = rec['probe']
            pf = D / f'g{key}' / 'probe.json'
            if pf.exists():
                P = json.load(open(pf)); row['probe_best'] = P.get('best'); row['probe_sec'] = P.get('sec')
        row['wall_s'] = S.get('gens', {}).get(key, {}).get('wall_s'); row['g05'] = S.get('gens', {}).get(key, {}).get('g05')
        out['gens'][key] = row
        # ---- print
        print(f'\n=== generation {key}  (wall {row["wall_s"]} s)')
        for name, s in row['slots'].items():
            sm = s['summary'] or {}
            print(f'  {name:14s} {s["choice"]:5s} ' + (f'{sm.get("covered")} @ {sm.get("sumJi")}  depths {sm.get("depths")} '
                                                       f'tanks {sm.get("tanks")}' if sm else f'running/partial {s["routes"]}'))
        if 'master' in row:
            m = row['master']; n8 = m['N8']; n9 = m.get('N9')
            print(f'  master N<=8 J-opt: {n8["covered"]} @ {n8["sumJi"]} (J {n8["J"]}, LP {n8["lp"]})' +
                  (f';  N<=9 J-opt: {n9["covered"]} @ {n9["sumJi"]} (J {n9["J"]}, LP {n9["lp"]})' if n9 else ''))
            print('  frontier N<=8: ' + ', '.join(f'<= {B}: {v["covered"]} @ {v["sumJi"]}' for B, v in m['frontier'].items()))
            for sub, picks in row['picked'].items():
                print(f'    {sub:10s}: ' + ', '.join(f'{s}:{n}@{t}' for s, n, t in picks))
        if 'probe' in row:
            p = row['probe']
            print(f'  probe ({row.get("probe_sec")} s) on {p["covered"]} @ {p["sumJi"]}: cheap {len(p["cheap"])} {p["cheap"]}, '
                  f'dear {len(p["dear"])} {p["dear"]}, unplaceable {len(p["unplaceable"])} {p["unplaceable"]}; '
                  f'projected closed J {p["projected_J"]}')
            if row.get('probe_best'):
                print('    prices: ' + ', '.join(f'{u}:{v["dJ"]}@{v["host"]}' for u, v in
                                                sorted(row['probe_best'].items(), key=lambda kv: kv[1]['dJ'])))
        if row.get('g05'):
            print(f'  G0.5: {row["g05"]}')
    print(f'\nstopped: {out["stopped"]} | driver: {out["stopped2"]}')
    json.dump(out, open(D / 'trajectory.json', 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
