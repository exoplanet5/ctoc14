"""Convert every route of one or more planned fleets that has no converted fragment yet (routes keyed by a signature of
launch epoch + flyby sequence), into a shared fragment directory; then right-size the tank (pass 2).
Pass 1: m0 = 600 + k * fuel_est + c (default 1.3, 20). Pass 2: m0 = 600 + 1.03 * actual + 10, kept if valid with >= as many
flybys. Fragments: <fragdir>/frag_<sig>_p1.txt / _p2.txt (+ _info.json); registry <fragdir>/registry.jsonl.
Usage: convert_pool.py fragdir fleetdir [fleetdir ...] [--nproc 3] [--k 1.3] [--c 20]"""
import sys, json, time, pathlib, argparse, hashlib, multiprocessing as mp
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ctoc14.kepler import Ephemeris
from ctoc14.lowthrust import convert_tour
from ctoc14.submission import write_submission
from ctoc14.validator import validate
from ctoc14.constants import M_DRY, M0_MAX

def signature(t):
    s = f"{round(float(t['t_launch']))}|" + ','.join(f"{int(l['ast'])}@{round(float(l['t_flyby']) / 3600)}" for l in t['legs'])
    return hashlib.sha1(s.encode()).hexdigest()[:16]

_EPH = None
def _init():
    global _EPH; _EPH = Ephemeris()

def _convert(args):
    tour, out, m0 = args
    tour = dict(tour); tour['m0'] = float(m0); tic = time.time()
    try:
        traj, info = convert_tour(_EPH, tour, verbose=False)
    except Exception as e:
        return dict(out=out, ok=False, error=repr(e), m0=m0)
    write_submission(out, [traj], header=f'm0={m0}')
    rep = validate(out, eph=_EPH, verbose=False)
    res = dict(out=out, ok=bool(rep.ok), m0=float(m0), fuel=float(info['fuel']), n_flybys=int(info['n_flybys']),
               dropped=list(info['dropped']), flyby_ids=sorted(rep.flybys), runtime=time.time() - tic)
    json.dump(dict(info, **res), open(out.replace('.txt', '_info.json'), 'w'), indent=1, default=float)
    return res

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('fragdir'); ap.add_argument('fleets', nargs='+')
    ap.add_argument('--nproc', type=int, default=3); ap.add_argument('--k', type=float, default=1.3); ap.add_argument('--c', type=float, default=20.0)
    a = ap.parse_args()
    fd = pathlib.Path(a.fragdir); fd.mkdir(parents=True, exist_ok=True); reg = fd / 'registry.jsonl'
    done = set()
    if reg.exists():
        for line in open(reg):
            try: done.add(json.loads(line)['sig'])
            except Exception: pass
    todo = {}
    for d in a.fleets:
        for f in sorted(pathlib.Path(d).glob('tour_sc*.json')):
            t = json.load(open(f)); s = signature(t)
            if s not in done and s not in todo: todo[s] = (t, str(f))
    print(f'{len(todo)} routes to convert ({len(done)} already converted)', flush=True)
    if not todo: return
    jobs = [(t, str(fd / f'frag_{s}_p1.txt'), min(M0_MAX, M_DRY + a.k * float(t['fuel_est']) + a.c)) for s, (t, _) in todo.items()]
    with mp.get_context('fork').Pool(a.nproc, initializer=_init) as pool:
        res1 = pool.map(_convert, jobs, chunksize=1)
    jobs2 = []; idx2 = []
    for (s, (t, src)), r in zip(todo.items(), res1):
        print(f"pass1 {s} ({len(t['legs'])} planned, {src}): ok={r.get('ok')} flybys={r.get('n_flybys')} fuel={r.get('fuel', float('nan')):.1f} "
              f"m0={r['m0']:.1f} dropped={r.get('dropped')} err={r.get('error')}", flush=True)
        if r.get('ok'):
            m0n = min(M0_MAX, M_DRY + 1.03 * r['fuel'] + 10.0)
            if m0n < r['m0'] - 5:
                jobs2.append((t, str(fd / f'frag_{s}_p2.txt'), m0n)); idx2.append(len(jobs2) - 1)
    res2 = []
    if jobs2:
        with mp.get_context('fork').Pool(a.nproc, initializer=_init) as pool:
            res2 = pool.map(_convert, jobs2, chunksize=1)
    by_out = {r['out']: r for r in res2}
    with open(reg, 'a') as fh:
        for (s, (t, src)), r1 in zip(todo.items(), res1):
            r2 = by_out.get(str(fd / f'frag_{s}_p2.txt'))
            if r2 is not None:
                print(f"pass2 {s}: m0={r2['m0']:.1f} ok={r2.get('ok')} flybys={r2.get('n_flybys')} (pass1 {r1.get('n_flybys')}) fuel={r2.get('fuel', float('nan')):.1f}", flush=True)
            fh.write(json.dumps(dict(sig=s, src=src, n_plan=len(t['legs']), fuel_plan=float(t['fuel_est']), p1=r1, p2=r2), default=float) + '\n')

if __name__ == '__main__':
    main()
