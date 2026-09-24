"""Tank tightening pass: re-convert fragments with m0 = 600 + (propellant actually used) + margin.
Right-sizing (convert_fleet / convert_pool pass 2: 600 + 1.03 x pass-1 fuel + 10) leaves 25-50 kg unused because the lighter
craft burns less than the heavier pass-1 craft. The tightened fragment is kept only if it validates with at least as many
flybys (the fragment MILP chooses among all versions anyway).
Tour lookup: conv_a-style fragments frag_<tag><pass>_sc<K>.txt -> <dir>/tour_sc<K>.json; convert_pool fragments
frag_<sig>_p<N>.txt -> registry.jsonl in the same directory. Output: <frag>_t.txt next to the input.
Usage: tighten_frags.py fraglist.txt [--margin 8] [--nproc 3] [--min-left 12]"""
import sys, json, re, time, pathlib, argparse, multiprocessing as mp
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from convert_pool import _convert, _init
from ctoc14.validator import parse
from ctoc14.constants import M_DRY

def tour_for(frag):
    p = pathlib.Path(frag)
    m = re.match(r'frag_([0-9a-f]{16})_p\d', p.stem)
    if m:
        for line in open(p.parent / 'registry.jsonl'):
            r = json.loads(line)
            if r['sig'] == m.group(1): return json.load(open(r['src']))
        return None
    m = re.match(r'frag_[a-z]+\d*_sc(\d+)', p.stem)
    if m and (p.parent / f'tour_sc{m.group(1)}.json').exists():
        return json.load(open(p.parent / f'tour_sc{m.group(1)}.json'))
    return None

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('fraglist'); ap.add_argument('--margin', type=float, default=8.0)
    ap.add_argument('--nproc', type=int, default=3); ap.add_argument('--min-left', type=float, default=12.0)
    ap.add_argument('--rounds', type=int, default=1, help='re-tighten the kept outputs while they still have > min-left kg')
    a = ap.parse_args()
    frags = [l.strip() for l in open(a.fraglist) if l.strip()]
    for rnd in range(a.rounds):
        kept = tighten_round(frags, a)
        print(f'round {rnd + 1}: {len(kept)} kept', flush=True)
        if not kept: break
        frags = kept

def tighten_round(frags, a):
    jobs = []; meta = []
    for frag in frags:
        rows = parse(frag); m0 = rows[0].m; mend = rows[-1].m; n = sum(1 for r in rows if r.event == 3)
        left = mend - M_DRY
        if left < a.min_left:
            print(f'skip {frag}: only {left:.1f} kg left'); continue
        t = tour_for(frag)
        if t is None:
            print(f'skip {frag}: tour not found'); continue
        m0n = M_DRY + (m0 - mend) + a.margin
        jobs.append((t, frag.replace('.txt', '_t.txt'), m0n)); meta.append((frag, m0, n, left))
    print(f'{len(jobs)} fragments to tighten', flush=True)
    if not jobs: return []
    with mp.get_context('fork').Pool(a.nproc, initializer=_init) as pool:
        res = pool.map(_convert, jobs, chunksize=1)
    kept = []
    for (frag, m0, n, left), r in zip(meta, res):
        ok = r.get('ok') and r.get('n_flybys', 0) >= n
        print(f"{frag}: m0 {m0:.1f} -> {r['m0']:.1f}, ok={r.get('ok')} flybys {r.get('n_flybys')} (was {n}), fuel {r.get('fuel', float('nan')):.1f} "
              f"-> {'KEEP' if ok else 'discard'}", flush=True)
        if ok: kept.append(r['out'])
        else: pathlib.Path(r['out']).unlink(missing_ok=True)
    return kept

if __name__ == '__main__':
    main()
