"""End-to-end campaign: planned fleet dir -> conversion -> greedy mop-up chains for leftovers -> incremental conversion -> merged,
validated submission. Usage: run_campaign.py plandir out.txt [--nproc 6] [--m0s 800,1000,1300,1600,2000] [--max-chains 4]"""
import sys, json, pathlib, subprocess, argparse, time
ROOT = pathlib.Path(__file__).resolve().parents[1]; PY = sys.executable
ap = argparse.ArgumentParser(); ap.add_argument('plandir'); ap.add_argument('out'); ap.add_argument('--nproc', type=int, default=6)
ap.add_argument('--m0s', default='800,1000,1300,1600,2000'); ap.add_argument('--max-chains', type=int, default=4); ap.add_argument('--skip-fleet', action='store_true')
a = ap.parse_args(); pd_ = pathlib.Path(a.plandir); out = pathlib.Path(a.out)
def run(cmd, log):
    print('>>', ' '.join(cmd), flush=True); t = time.time()
    with open(log, 'w') as f: subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, cwd=ROOT)
    print(f'   done in {time.time()-t:.0f} s', flush=True)
sys.path.insert(0, str(ROOT))
from ctoc14.constants import cost_sc
# 1. fleet conversion
merged0 = str(out).replace('.txt', '_fleet.txt')
if not a.skip_fleet:
    run([PY, 'tools/convert_fleet.py', str(pd_), merged0, '--nproc', str(a.nproc), '--tag', 'p'], pd_ / 'convert_fleet.log')
rep = json.load(open(merged0.replace('.txt', '_report.json')))
frags = [r['out'] for r in rep['results'] if r.get('ok') and r.get('n_flybys', 0) >= 2]
covered = set()
for r in rep['results']:
    if r.get('ok') and r.get('n_flybys', 0) >= 2: covered |= set(r['flyby_ids'])
print(f'fleet: {len(frags)} usable spacecraft, covered {len(covered)}', flush=True)
# 2. greedy mop-up chains
extras = []
for ci in range(a.max_chains):
    left = sorted(set(range(1, 301)) - covered - {131, 144})
    if len(left) < 2: break
    tj = pd_ / f'mop{ci}_targets.json'; json.dump(left, open(tj, 'w'))
    run([PY, 'tools/run_tail_search.py', str(pd_ / f'mop{ci}'), str(tj)] + a.m0s.split(','), pd_ / f'mop{ci}_search.log')
    best = None
    for m0 in a.m0s.split(','):
        f = pd_ / f'mop{ci}_m{int(float(m0))}.json'
        if not f.exists(): continue
        t = json.load(open(f)); gain = t['n_flybys'] - cost_sc(t['m0'])
        if gain > 0 and (best is None or gain > best[0]): best = (gain, f, t)
    if best is None: print('no worthwhile chain'); break
    gain, f, t = best; print(f'mop-up chain {ci}: {t["n_flybys"]} targets, m0 {t["m0"]}, planned gain {gain:.2f}', flush=True)
    ef = pd_ / f'extra_mop{ci}.json'; json.dump(t, open(ef, 'w'))
    run([PY, 'tools/convert_one.py', str(ef), str(pd_ / f'extra_mop{ci}_frag.txt')], pd_ / f'convert_mop{ci}.log')
    info = json.load(open(pd_ / f'extra_mop{ci}_frag_info.json'))
    ids = set()
    for line in open(pd_ / f'extra_mop{ci}_frag.txt'):
        p = line.split()
        if len(p) == 15 and p[2] == '3': ids.add(int(p[14]))
    newc = ids - covered
    if len(newc) - cost_sc(t['m0']) > 0 and 'valid=True' in open(pd_ / f'convert_mop{ci}.log').read():
        extras.append(str(pd_ / f'extra_mop{ci}_frag.txt')); covered |= newc
        print(f'   converted: {len(ids)} flybys, {len(newc)} new; fuel {info["fuel"]:.0f} kg', flush=True)
    else:
        print(f'   chain not worthwhile after conversion ({len(newc)} new)', flush=True)
# 3. merge + validate
run([PY, 'tools/combine_submission.py', str(out)] + frags + extras, str(out).replace('.txt', '_combine.log'))
print(open(str(out).replace('.txt', '_combine.log')).read()[-800:])
