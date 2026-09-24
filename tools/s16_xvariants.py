"""Stage 16: list the twin VARIANTS (same target set, different settled state) of every route of a fleet, from the column
cache of tools/s16_select.py (results/s16/colcache.pkl) -- the candidates for exact-aware variant choice (s16_xconv.py).

Per route: every cached column whose target set equals the route's, deduped by content hash (s16_iter.chash), lightest
twin first; keep at most --per variants with twin tank <= lightest + --dkg, skipping variants whose twin tank is within
--sep kg of one already kept (near-identical states convert alike).
Writes OUT_JSON {route: [[file, tank], ...]} and prints the file list (one per line) with --files.

usage: s16_xvariants.py FLEET OUT_JSON [--per 6] [--dkg 2.5] [--sep 0.05] [--files]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, pickle, pathlib, argparse, collections
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)
import s15b_fleet as FA
from s16_iter import chash


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('fleet'); ap.add_argument('out')
    ap.add_argument('--per', type=int, default=6); ap.add_argument('--dkg', type=float, default=2.5)
    ap.add_argument('--sep', type=float, default=0.05); ap.add_argument('--files', action='store_true')
    a = ap.parse_args()
    cache = pickle.load(open(ROOT / 'results/s16/colcache.pkl', 'rb'))
    byset = collections.defaultdict(dict)
    for k, (c, why) in cache.items():
        if c is None or not os.path.exists(ROOT / k[0]):
            continue
        byset[frozenset(int(x) for x in c['asts'])][k[0]] = float(c['tank'])
    res = {}
    for name, f in FA.fleet_files(a.fleet).items():
        st = FA.load_st(f); s = frozenset(int(x) for x in st['asts'])
        cands = sorted((t, g) for g, t in byset[s].items())
        t0 = min([t for t, _ in cands] + [99999.0])
        keep, hs = [], set()
        for t, g in [(None, f)] + cands:                 # the fleet's own route first
            if len(keep) >= a.per:
                break
            h = chash(FA.load_st(g))
            if h in hs:
                continue
            tt = t if t is not None else float(__import__('run_ialns').ipr(st).tank())
            if tt > t0 + a.dkg or any(abs(tt - k2) < a.sep for _, k2 in keep):
                hs.add(h); continue
            hs.add(h); keep.append((g, round(tt, 3)))
        res[name] = keep
    json.dump(res, open(ROOT / a.out, 'w'), indent=1)
    if a.files:
        for name, keep in res.items():
            for g, t in keep:
                print(g)
    else:
        for name, keep in res.items():
            print(name, ' '.join(f'{t:.2f}' for _, t in keep))


if __name__ == '__main__':
    main()
