"""Convert a few tours with a chosen set of convert_tour options, in parallel (<= 3 processes), validate, and log.
Usage: convert_review2.py TAG [--tours 3,8,7] [--nproc 3] [--opt key=value ...] [--verbose]
Outputs: results/review2/TAG_scN.txt (fragment), TAG_scN_info.json, TAG_scN.log (per-tour log), TAG_summary.json
"""
import sys, json, time, pathlib, argparse, contextlib, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from ctoc14.kepler import Ephemeris
from ctoc14.lowthrust import convert_tour
from ctoc14.submission import write_submission
from ctoc14.validator import validate

_EPH = None
def _init():
    global _EPH; _EPH = Ephemeris()

def _parse_val(v):
    try:
        return json.loads(v)
    except Exception:
        return v

def _convert(args):
    sc, tag, opts, verbose, outdir = args
    tour = json.load(open(ROOT / 'results/fleet_b300' / f'tour_sc{sc}.json'))
    out = outdir / f'{tag}_sc{sc}.txt'
    log = outdir / f'{tag}_sc{sc}.log'
    tic = time.time()
    with open(log, 'w') as lf, contextlib.redirect_stdout(lf):
        print(f'tour_sc{sc} m0={tour["m0"]} planned={tour["n_flybys"]} opts={opts}', flush=True)
        try:
            traj, info = convert_tour(_EPH, tour, verbose=verbose, **opts)
        except Exception as e:
            import traceback; traceback.print_exc()
            return dict(sc=sc, ok=False, error=repr(e), runtime=time.time() - tic)
        rt = time.time() - tic
        write_submission(str(out), [traj], header=f'tour_sc{sc} m0={tour["m0"]} opts={opts}')
        rep = validate(str(out), eph=_EPH, verbose=True)
        res = dict(sc=sc, ok=rep.ok, errors=rep.errors[:5], n_flybys=info['n_flybys'], n_covered=rep.n_covered,
                   fuel=float(info['fuel']), planned=tour['n_flybys'], dropped=info.get('dropped'), runtime=rt,
                   n_replans=info.get('n_replans', 0), events=info.get('events', []), warn=info.get('warn', []),
                   flyby_ids=sorted(rep.flybys), max_thrust=rep.max_thrust, legs=info['legs'])
        json.dump(res, open(outdir / f'{tag}_sc{sc}_info.json', 'w'), indent=1, default=float)
        print(f'DONE sc{sc}: flybys {info["n_flybys"]}/{tour["n_flybys"]} fuel {info["fuel"]:.1f} valid={rep.ok} runtime {rt:.0f}s', flush=True)
    return res

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('tag'); ap.add_argument('--tours', default='3,8,7'); ap.add_argument('--nproc', type=int, default=3)
    ap.add_argument('--opt', action='append', default=[]); ap.add_argument('--verbose', action='store_true')
    a = ap.parse_args()
    opts = {}
    for kv in a.opt:
        k, v = kv.split('=', 1); opts[k] = _parse_val(v)
    outdir = ROOT / 'results/review2'; outdir.mkdir(parents=True, exist_ok=True)
    scs = [int(x) for x in a.tours.split(',')]
    jobs = [(sc, a.tag, opts, a.verbose, outdir) for sc in scs]
    tic = time.time()
    with mp.get_context('fork').Pool(min(a.nproc, len(jobs)), initializer=_init) as pool:
        res = pool.map(_convert, jobs)
    summ = dict(tag=a.tag, opts=opts, results=[{k: v for k, v in r.items() if k != 'legs'} for r in res], wall=time.time() - tic)
    json.dump(summ, open(outdir / f'{a.tag}_summary.json', 'w'), indent=1, default=float)
    for r in res:
        print(f"{a.tag} sc{r['sc']}: ok={r.get('ok')} flybys={r.get('n_flybys')}/{r.get('planned')} fuel={r.get('fuel', float('nan')):.1f} "
              f"runtime={r.get('runtime', 0):.0f}s replans={r.get('n_replans')} dropped={r.get('dropped')} err={r.get('error') or r.get('errors')}", flush=True)

if __name__ == '__main__':
    main()
