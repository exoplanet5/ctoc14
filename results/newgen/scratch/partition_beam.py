"""Partition first, sequence second: assign every target to the nearest skeleton route (closest approach), then run the
PLAIN beam (no prizes) inside each subset (all other targets excluded). Reports depth / fuel / hard count per subset.
Usage: partition_beam.py skel_dir out_dir [--nproc 8] [--beam 400]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import run_ialns as RI
from ctoc14.search import Params, beam_search, tour_json, UNREACHABLE
from ctoc14.constants import DAY, AU

ap = argparse.ArgumentParser(); ap.add_argument('src'); ap.add_argument('dst'); ap.add_argument('--nproc', type=int, default=8)
ap.add_argument('--beam', type=int, default=400); ap.add_argument('--wfuel', type=float, default=1.0); ap.add_argument('--m0', type=float, default=1000.0)
a = ap.parse_args()
out = pathlib.Path(a.dst); out.mkdir(parents=True, exist_ok=True)
fleet = RI.IFleet(a.src); eph = RI.eph()
cnt = np.zeros(301, int)
for line in open(ROOT / 'results/newgen/recost35.jsonl'):
    r = json.loads(line)
    if r['ok']:
        for t in r['targets']: cnt[t] += 1
hard = set(t for t in range(1, 301) if cnt[t] <= 5) - set(UNREACHABLE)
allt = [t for t in range(1, 301) if t not in UNREACHABLE]
own = {n: set(int(x) for x in r['st']['asts']) for n, r in fleet.routes.items()}
covered = set().union(*own.values())
free = [t for t in allt if t not in covered]
with mp.get_context('fork').Pool(a.nproc) as pool:
    by = RI.candidates(fleet, pool, free, 1.5, 1)
subset = {n: set(s) for n, s in own.items()}
for t in free:
    if t in by and by[t]:
        subset[by[t][0]['host']].add(t)
    else:
        n = min(subset, key=lambda n: len(subset[n])); subset[n].add(t)
# overlaps between skeletons: keep in the first (dedupe later)
seen = set()
for n in sorted(subset):
    subset[n] -= seen; seen |= subset[n]
print({n: (len(s), len(s & hard)) for n, s in sorted(subset.items())}, flush=True)
json.dump({n: sorted(s) for n, s in subset.items()}, open(out / 'subsets.json', 'w'))
P = Params(beam=a.beam, w_fuel=a.wfuel, m_margin=40.0, dv_max=1.2, lin_tofs=np.arange(20, 600 + 1e-9, 10) * DAY,
           lin_drmax=0.15 * AU, vinf_cap=2.0, tof_refine=True)
grid = np.arange(0.0, 600.0, 20.0) * DAY
res = {}
for n in sorted(subset):
    tic = time.time()
    excluded = sorted(set(allt) - subset[n])
    best, beam = beam_search(eph, excluded=excluded, m0=a.m0, t_launch_grid=grid, P=P, n_proc=a.nproc, verbose=False)
    ids = [x[0] for x in best.seq]
    print(f'subset {n}: {len(subset[n])} targets ({len(subset[n] & hard)} hard) -> best route {len(ids)} flybys ({len(set(ids) & hard)} hard), '
          f'fuel {best.fuel:.0f} kg, launch {best.t_launch/DAY:.0f} d, end {best.t/DAY/365.25:.1f} yr  ({time.time()-tic:.0f} s)', flush=True)
    t = tour_json(best); t['subset'] = sorted(subset[n]); json.dump(t, open(out / f'tour_{n}.json', 'w'))
    res[n] = ids
cov = set().union(*[set(v) for v in res.values()])
print(f'partition beam: {sum(len(v) for v in res.values())} flybys, {len(cov)} distinct covered of 298; hard covered {len(cov & hard)}/{len(hard)}')
