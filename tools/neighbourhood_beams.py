"""Neighbourhood beams: build a pool of DEEP, HARD-RICH columns for the 8/9-craft partition MIP.

A plain (prize-free) beam returns 33-36-flyby routes only when it is confined to a subset of ~50-70 mutually
compatible targets (results/newgen/pbeam8.out); on the full 298 it drifts to the easy targets and on subsets of
<= 37 it starves.  So: take every route we trust as a SEED, define its neighbourhood = its own targets + the K
nearest other targets (closest approach of the asteroid to the seed trajectory, run_ialns.w_cands), optionally
discounting the distance of HARD targets so the subset is hard-rich, and run the plain beam inside it.
Every beam state with >= --min-keep flybys is written as a tour_json record (columns.jsonl) for the twin recost
and the partition MIP.

Usage: neighbourhood_beams.py out_dir [--seeds dir,dir,...] [--sizes 55,70] [--hard-w 1.0,0.5] [--nproc 8]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, glob, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from ctoc14.search import Params, beam_search, tour_json, UNREACHABLE
from ctoc14.colgen import CostModel
from ctoc14.constants import DAY, AU

ALL = [t for t in range(1, 301) if t not in UNREACHABLE]


def hard_set(path=ROOT / 'results/newgen/recost35.jsonl', thr=5):
    cnt = np.zeros(301, int)
    for line in open(path):
        r = json.loads(line)
        if r.get('ok'):
            for t in r['targets']: cnt[t] += 1
    return set(t for t in ALL if cnt[t] <= thr)


class DeepCollect(list):
    """Beam collector that keeps only states worth costing (P.collect is extended with the whole beam at every depth)."""
    def __init__(self, minn): super().__init__(); self.minn = minn
    def extend(self, states):
        for s in states:
            if len(s.seq) >= self.minn: list.append(self, s)


class BP(Params):
    def __getstate__(self):            # never ship the collector to the beam workers
        d = dict(self.__dict__); d.pop('collect', None); return d


def load_seeds(specs, max_per_dir=0):
    """-> list of (name, state dict). Directories of route_*.npz / col_*.npz, or explicit npz paths."""
    seeds = []
    for spec in specs:
        p = pathlib.Path(spec)
        if not p.is_absolute(): p = ROOT / spec
        files = sorted(glob.glob(str(p / '*.npz'))) if p.is_dir() else sorted(glob.glob(str(p)))
        if max_per_dir: files = files[:max_per_dir]
        for f in files:
            z = np.load(f); st = {k: z[k] for k in z.files}; st['tL'] = float(st['tL'])
            seeds.append((f'{pathlib.Path(f).parent.name}/{pathlib.Path(f).stem}', st))
    out = []; seen = set()
    for name, st in seeds:
        k = frozenset(int(a) for a in st['asts'])
        if len(k) < 8 or k in seen: continue
        seen.add(k); out.append((name, st))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out')
    ap.add_argument('--seeds', default='results/newgen/skel8,results/newgen/skel8b,results/newgen/ifleet_t10d')
    ap.add_argument('--sizes', default='55,70'); ap.add_argument('--hard-w', default='1.0,0.5')
    ap.add_argument('--hard-frac', default='0.0', help='fraction of each subset reserved for the nearest HARD targets')
    ap.add_argument('--own', default='keep', choices=('keep', 'drop'), help='force the seed own targets into the subset')
    ap.add_argument('--dmax', type=float, default=1.5)
    ap.add_argument('--beam', type=int, default=400); ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--m0', type=float, default=1000.0)
    ap.add_argument('--wfuel', type=float, default=-1.0, help='<0: from CostModel(0.70, 2) at this m0')
    ap.add_argument('--prize', type=float, default=0.0, help='uniform prize (J) per target in the subset; 0 = plain beam')
    ap.add_argument('--wt', type=float, default=1.0)
    ap.add_argument('--vinf', type=float, default=2.0); ap.add_argument('--dvmax', type=float, default=1.2)
    ap.add_argument('--tofmax', type=float, default=400.0); ap.add_argument('--lintofmax', type=float, default=600.0)
    ap.add_argument('--grid', default='0,800,20'); ap.add_argument('--min-keep', type=int, default=28)
    ap.add_argument('--keep', type=int, default=300); ap.add_argument('--max-seeds', type=int, default=0)
    ap.add_argument('--max-runs', type=int, default=0); ap.add_argument('--hard-thr', type=int, default=5)
    a = ap.parse_args()
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    say = RI.logger(out / 'log.txt')
    hard = hard_set(thr=a.hard_thr)
    sizes = [int(x) for x in a.sizes.split(',')]; hws = [float(x) for x in a.hard_w.split(',')]
    hfs = [float(x) for x in a.hard_frac.split(',')]
    seeds = load_seeds(a.seeds.split(','))
    if a.max_seeds: seeds = seeds[:a.max_seeds]
    say(f'{len(seeds)} distinct seeds, {len(hard)} hard targets, sizes {sizes}, hard weights {hws}')
    RI.eph()
    tic = time.time()
    with mp.get_context('fork').Pool(a.nproc) as pool:
        res = pool.map(RI.w_cands, [(n, st, ALL, a.dmax) for n, st in seeds])
    rank = {}
    for (n, st), lst in zip(seeds, res):
        d = {}
        for c in lst: d[c['ast']] = min(d.get(c['ast'], 9e9), c['dist'])
        rank[n] = d
    say(f'neighbourhood ranking: {time.time()-tic:.0f} s; '
        f'median reachable targets per seed {np.median([len(v) for v in rank.values()]):.0f}')

    g = [float(x) for x in a.grid.split(',')]; grid = np.arange(g[0], g[1] + 1e-9, g[2]) * DAY
    CM = CostModel(s=0.70, reserve=2.0)
    wf = a.wfuel if a.wfuel >= 0 else CM.w_fuel(0.5 * (a.m0 - 640.0), a.m0)
    P = BP(beam=a.beam, w_fuel=wf, m_margin=40.0, dv_max=a.dvmax, lin_tofs=np.arange(20, a.lintofmax + 1e-9, 10) * DAY,
           tofs=np.arange(15, a.tofmax + 1e-9, 5) * DAY,
           lin_drmax=0.15 * AU, vinf_cap=a.vinf, tof_refine=True, w_t=a.wt)
    say(f'm0 {a.m0:.0f}, w_fuel {wf:.3f}, prize {a.prize}, tank ceiling {CM.tank(a.m0-640.0, a.m0):.0f} kg')
    done = set()
    dfile = out / 'done.json'
    if dfile.exists(): done = set(tuple(x) for x in json.load(open(dfile)))
    cf = open(out / 'columns.jsonl', 'a')
    jobs = [(n, size, hw, hf) for n, _ in seeds for size in sizes for hw in hws for hf in hfs]
    if a.max_runs: jobs = jobs[:a.max_runs]
    say(f'{len(jobs)} beam runs ({len(done)} already done)')
    stt = {n: set(int(x) for x in st['asts']) for n, st in seeds}
    for i, (name, size, hw, hf) in enumerate(jobs, 1):
        if (name, size, hw, hf) in done: continue
        own = stt[name] if a.own == 'keep' else set()
        d_of = dict(rank[name])
        if a.own == 'keep':
            for t in stt[name]: d_of.setdefault(t, 0.0)
        else:                                  # the seed is only a distance reference: its own targets are excluded
            for t in stt[name]: d_of.pop(t, None)
        cand = sorted(((d * (hw if t in hard else 1.0), t) for t, d in d_of.items() if t not in own))
        sub = set(own)
        nh = int(round(hf * size)) - len(sub & hard)
        if nh > 0:
            for _, t in [c for c in cand if c[1] in hard][:nh]: sub.add(t)
        for _, t in cand:
            if len(sub) >= size: break
            sub.add(t)
        sub = sorted(sub)
        t0 = time.time()
        if a.prize:
            pz = np.zeros(300)
            for t in sub: pz[t - 1] = a.prize
            P.prize = pz
        P.collect = DeepCollect(a.min_keep)
        best, beam = beam_search(RI.eph(), excluded=sorted(set(ALL) - set(sub)), m0=a.m0,
                                 t_launch_grid=grid, P=P, n_proc=a.nproc, verbose=False)
        states = list(P.collect) + [s for s in beam if len(s.seq) >= a.min_keep]
        uniq = {}
        for s in states:
            k = frozenset(x[0] for x in s.seq)
            if k not in uniq or s.fuel < uniq[k].fuel: uniq[k] = s
        keep = sorted(uniq.values(), key=lambda s: (-len(s.seq), s.fuel))[:a.keep]
        tag = f'{name}:{size}:{hw}:{hf}'
        nw = 0
        for s in keep:
            r = tour_json(s); r['tag'] = tag; r['targets'] = sorted(x[0] for x in s.seq)
            r['cost_planner'] = float(CM.cost(s.fuel, s.m0)); r['tank_planner'] = float(CM.tank(s.fuel, s.m0))
            r['n_hard'] = len(set(r['targets']) & hard)
            cf.write(json.dumps(r) + '\n'); nw += 1
        cf.flush()
        deep = [s for s in keep if len(s.seq) >= 33 and len(set(x[0] for x in s.seq) & hard) >= 8]
        bn = len(best.seq) if best else 0
        say(f'[{i:4d}/{len(jobs)}] {tag}: subset {len(sub)} ({len(set(sub)&hard)} hard) -> best {bn} flybys '
            f'({len(set(x[0] for x in best.seq)&hard) if best else 0} hard, fuel {best.fuel:.0f} kg, '
            f'launch {best.t_launch/DAY:.0f} d, end {best.t/DAY/365.25:.1f} yr), kept {nw} cols, '
            f'{len(deep)} deep+hard  ({time.time()-t0:.0f} s, {(time.time()-tic)/60:.0f} min total)')
        done.add((name, size, hw, hf)); json.dump(sorted(done), open(dfile, 'w'))
    cf.close()


if __name__ == '__main__':
    main()
