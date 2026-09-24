"""Incrementally convert extra tours (e.g. mop-up chains) given the set of already-covered targets, then merge all valid fragments.
Usage: convert_extra.py out.txt covered.json frag_existing1.txt ... -- tour1.json[:m0] tour2.json[:m0] ..."""
import sys, json, pathlib, multiprocessing as mp
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ctoc14.kepler import Ephemeris
from ctoc14.lowthrust import convert_tour
from ctoc14.submission import write_submission
from ctoc14.validator import validate
from ctoc14.constants import cost_sc, M_DRY, M0_MAX
_EPH = None
def _init():
    global _EPH; _EPH = Ephemeris()
def _convert(args):
    tour, out, excluded = args
    try:
        traj, info = convert_tour(_EPH, tour, verbose=False, global_excluded=excluded, allowed_extra=())
    except Exception as e:
        return dict(out=out, ok=False, error=repr(e))
    write_submission(out, [traj], header=f'm0={tour["m0"]}')
    rep = validate(out, eph=_EPH, verbose=False)
    return dict(out=out, ok=rep.ok, n=rep.n_covered, fuel=info['fuel'], m0=tour['m0'], dropped=info['dropped'], errors=rep.errors[:3], ids=sorted(rep.flybys))
def main():
    out = sys.argv[1]; covered = set(json.load(open(sys.argv[2])))
    k = sys.argv.index('--'); frags = sys.argv[3:k]; specs = sys.argv[k + 1:]
    jobs = []
    for sp in specs:
        path, _, m0 = sp.partition(':'); t = json.load(open(path))
        if m0: t['m0'] = float(m0)
        jobs.append((t, path.replace('.json', '_frag.txt'), covered))
    with mp.get_context('fork').Pool(min(8, len(jobs)), initializer=_init) as pool:
        res = pool.map(_convert, jobs)
    for r in res: print('extra:', r, flush=True)
    keep = [r for r in res if r.get('ok') and r['n'] >= 2]
    # right-size tanks for kept extras (second pass)
    jobs2 = []
    for r, (t, path, _) in zip(res, jobs):
        if r.get('ok') and r['n'] >= 2:
            m0n = min(M0_MAX, M_DRY + 1.03 * r['fuel'] + 10.0)
            if m0n < t['m0'] - 5:
                t2 = dict(t); t2['m0'] = m0n; jobs2.append((t2, path.replace('.json', '_frag2.txt'), covered))
    if jobs2:
        with mp.get_context('fork').Pool(min(8, len(jobs2)), initializer=_init) as pool:
            res2 = pool.map(_convert, jobs2)
        for r2 in res2:
            print('extra pass2:', r2, flush=True)
            if r2.get('ok') and r2['n'] >= 2:
                base = r2['out'].replace('_frag2.txt', '_frag.txt')
                for i, r in enumerate(keep):
                    if r['out'] == base and r2['n'] >= r['n']:
                        keep[i] = r2
    allfrags = list(frags) + [r['out'] for r in keep]
    line = 1
    with open(out, 'w', encoding='utf-8') as f:
        f.write('# CTOC14 Problem A submission (ctoc14 solver)\n')
        for sc, fr in enumerate(allfrags, 1):
            for row in open(fr, encoding='utf-8'):
                if row.startswith('#') or not row.strip(): continue
                p = row.split(); p[0] = str(line); p[1] = str(sc); f.write(' '.join(p) + '\n'); line += 1
    rep = validate(out, verbose=True)
    print(f'FINAL: {len(allfrags)} spacecraft, covered {rep.n_covered}, J = {rep.J:.3f}, missed {sorted(set(range(1,301)) - set(rep.flybys))}')
if __name__ == '__main__':
    main()
