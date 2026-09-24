"""Print, for every route of the given fleet directories, its best existing converted fragment: among
frag_<sig>_p1/_p2(+ _t, _t_t, ...) in the fragment directory (tools/convert_pool.py) and the conv_a-style fragments
registered there, the valid one with the smallest initial mass that keeps the most flybys.
Usage: fleet_frags.py fragdir fleetdir [fleetdir ...]   (one path per line on stdout; missing routes on stderr)"""
import sys, json, glob, pathlib, re
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1])); sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from convert_pool import signature
from ctoc14.validator import parse

fd = pathlib.Path(sys.argv[1]); reg = {}
for line in open(fd / 'registry.jsonl'):
    try:
        r = json.loads(line); reg[r['sig']] = r
    except Exception:
        pass

def versions(r):
    """All fragment files derived from a registry record (p1, p2 and their tightened chains)."""
    out = []
    for key in ('p1', 'p2'):
        p = (r.get(key) or {}).get('out')
        if not p: continue
        base = p[:-4]
        pat = re.compile(re.escape(base) + r'(_t)*\.txt$')
        out += [f for f in glob.glob(base + '*.txt') if pat.match(f)]
    return sorted(set(out))

seen = set()
for d in sys.argv[2:]:
    for f in sorted(pathlib.Path(d).glob('tour_sc*.json')):
        s = signature(json.load(open(f)))
        if s in seen: continue
        seen.add(s)
        r = reg.get(s)
        if r is None:
            print(f'# no fragment for {f}', file=sys.stderr); continue
        best = None
        for v in versions(r):
            try:
                rows = parse(v)
            except Exception:
                continue
            info = pathlib.Path(v[:-4] + '_info.json')
            if info.exists() and not json.load(open(info)).get('ok', True): continue
            n = sum(1 for x in rows if x.event == 3); m0 = rows[0].m
            key = (-n, m0)
            if best is None or key < best[0]: best = (key, v)
        if best: print(best[1])
        else: print(f'# no valid fragment for {f}', file=sys.stderr)
