"""Convert all tours in a directory into one validated submission.
Usage: convert_fleet.py tourdir out.txt [--nproc 8] [--no-resize] [--tag name]
Pass 1: convert each tour with its planned m0 (re-planning may pick up leftover targets not in any tour).
Pass 2: right-size m0 = 600 + 1.03*fuel + 10 (if smaller than planned) and re-convert; keep whichever pass is valid with lower cost.
"""
import sys, json, time, pathlib, argparse, multiprocessing as mp
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from ctoc14.kepler import Ephemeris
from ctoc14.lowthrust import convert_tour
from ctoc14.submission import write_submission
from ctoc14.validator import validate
from ctoc14.constants import cost_sc, M_DRY, FUEL_MAX, M0_MAX

_EPH = None
def _init():
    global _EPH; _EPH = Ephemeris()

def _convert(args):
    tour, out, global_excluded, allowed_extra, m0 = args
    tour = dict(tour); tour['m0'] = m0
    tic = time.time()
    try:
        traj, info = convert_tour(_EPH, tour, verbose=False, global_excluded=global_excluded, allowed_extra=allowed_extra)
    except Exception as e:
        return dict(out=out, ok=False, error=repr(e), m0=m0)
    write_submission(out, [traj], header=f'm0={m0}')
    rep = validate(out, eph=_EPH, verbose=False)
    info.update(out=out, ok=rep.ok, errors=rep.errors[:5], m0=m0, n_covered=rep.n_covered, flyby_ids=sorted(rep.flybys), runtime=time.time() - tic)
    info['legs'] = [dict(ast=l['ast'], t_flyby=float(l['t_flyby']), dist=float(l['dist'])) for l in info['legs']]
    return info

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('tourdir'); ap.add_argument('out'); ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--no-resize', action='store_true'); ap.add_argument('--tag', default='p')
    a = ap.parse_args()
    td = pathlib.Path(a.tourdir)
    tours = [json.load(open(f)) for f in sorted([f for f in td.glob('tour_sc*.json') if f.stem.split('sc')[1].isdigit()], key=lambda f: int(f.stem.split('sc')[1]))]
    sets = [{l['ast'] for l in t['legs']} for t in tours]
    covered = set().union(*sets)
    leftovers = set(range(1, 301)) - covered - {131, 144}
    print(f'{len(tours)} tours, planned coverage {len(covered)}, leftovers {sorted(leftovers)}', flush=True)
    jobs = []
    for k, t in enumerate(tours):
        others = set().union(*[s for j, s in enumerate(sets) if j != k])
        jobs.append((t, str(td / f'frag_{a.tag}1_sc{k+1}.txt'), others, leftovers, float(t['m0'])))
    with mp.get_context('fork').Pool(a.nproc, initializer=_init) as pool:
        res1 = pool.map(_convert, jobs)
    for k, r in enumerate(res1):
        print(f"pass1 SC{k+1}: ok={r.get('ok')} flybys={r.get('n_flybys')} fuel={r.get('fuel', float('nan')):.1f} dropped={r.get('dropped')} err={r.get('error') or r.get('errors')}", flush=True)
    results = list(res1)
    if not a.no_resize:
        jobs2 = []; idx2 = []
        for k, r in enumerate(res1):
            if not r.get('ok'): continue
            m0n = min(M0_MAX, M_DRY + 1.03 * r['fuel'] + 10.0)
            if m0n < r['m0'] - 5:
                jobs2.append((tours[k], str(td / f'frag_{a.tag}2_sc{k+1}.txt'), jobs[k][2], leftovers, m0n)); idx2.append(k)
        if jobs2:
            with mp.get_context('fork').Pool(min(a.nproc, len(jobs2)), initializer=_init) as pool:
                res2 = pool.map(_convert, jobs2)
            for k, r in zip(idx2, res2):
                print(f"pass2 SC{k+1}: m0={r['m0']:.1f} ok={r.get('ok')} flybys={r.get('n_flybys')} fuel={r.get('fuel', float('nan')):.1f} dropped={r.get('dropped')} err={r.get('error') or r.get('errors')}", flush=True)
                if r.get('ok') and r['n_flybys'] >= res1[k]['n_flybys'] and r['fuel'] <= r['m0'] - M_DRY:
                    results[k] = r
    # combine valid fragments
    frags = [r for r in results if r.get('ok')]
    line = 1
    with open(a.out, 'w', encoding='utf-8') as f:
        f.write('# CTOC14 Problem A submission (ctoc14 solver)\n')
        for sc, r in enumerate(frags, 1):
            for row in open(r['out'], encoding='utf-8'):
                if row.startswith('#') or not row.strip(): continue
                p = row.split(); p[0] = str(line); p[1] = str(sc); f.write(' '.join(p) + '\n'); line += 1
    rep = validate(a.out, verbose=True)
    json.dump(dict(results=results, J=rep.J, n_covered=rep.n_covered, missed=sorted(set(range(1, 301)) - set(rep.flybys))), open(a.out.replace('.txt', '_report.json'), 'w'), indent=1, default=float)
    print(f'FINAL: {len(frags)} spacecraft, covered {rep.n_covered}, J = {rep.J:.3f}, missed {sorted(set(range(1,301)) - set(rep.flybys))}')

if __name__ == '__main__':
    main()
