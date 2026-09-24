"""Stage 15 [H]: HYBRID chains -- hand the small-pool end of finished production-beam chains (tools/s15_isogen.py
run_chain, results/s15/cols/rNN/[CK]*) to the lintwin twin beam (tools/s15_twintail.py).

Round 1/2 chains starve at their last steps (pools 16-38 give the production beam 4-18 flybys, while the twin beam
takes 12 of 16 and 22 of 38 there, run TB1).  For every chain: cut at the first step whose pool (chain pool minus
the targets taken before it) is <= --cut; the twin job re-plans steps k..L over that pool with the chain's premium set
and carries the chain's heads + taken routes before k as 'fleet_pre' (tools/s15_twinfleet.py assembles the fleet).

usage: s15_hybrid.py OUT_JOBS.json --round 2 [--cut 45] [--nproc 1] [--tag h02]"""
import os, sys, json, pathlib, argparse
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
sys.path.insert(0, str(ROOT / 'tools'))
os.chdir(ROOT)


def taken_files(C):
    out = []
    for t in C.get('taken', []):
        tg = sorted(t['targets'])
        m = [c for c in C['cols'] if sorted(c['targets']) == tg]
        if not m:
            return None
        out.append(min(m, key=lambda c: abs(c['tank'] - t['tank']))['f'])
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--round', type=int, default=2)
    ap.add_argument('--cut', type=int, default=45); ap.add_argument('--nproc', type=int, default=1)
    ap.add_argument('--tag', default=None); ap.add_argument('--only', default='')
    ap.add_argument('--p', type=float, default=None); ap.add_argument('--w-fuel', type=float, default=1.0)
    ap.add_argument('--partial', action='store_true', help='also chains that are still running')
    a = ap.parse_args()
    L = json.load(open(ROOT / 'results/s15/loop.json'))
    R = L['rounds'][f'{a.round:02d}']
    tag = a.tag or f'h{a.round:02d}'
    jobs = []
    for j in R['jobs']:
        if j['kind'] != 'chain' or (a.only and j['tag'] not in a.only.split(',')):
            continue
        d = ROOT / j['out']
        ck = d / 'chain.json'
        if not ck.exists():
            print(f'{j["tag"]}: no chain.json yet'); continue
        if not (d / 'meta.json').exists() and not a.partial:
            print(f'{j["tag"]}: chain still running'); continue
        C = json.load(open(ck))
        tf = taken_files(C)
        if tf is None:
            print(f'{j["tag"]}: taken route file not found'); continue
        pool0 = sorted(j['pool']); cov = set(); k_cut = None
        for i, t in enumerate(C.get('taken', [])):
            if len(set(pool0) - cov) <= a.cut:
                k_cut = i; break
            cov |= set(t['targets'])
        if k_cut is None:
            k_cut = len(C.get('taken', []))
            if len(set(pool0) - cov) > a.cut + 25:
                print(f'{j["tag"]}: residual {len(set(pool0) - cov)} after all {k_cut} taken steps > cut+25: skipped'); continue
        pool = sorted(set(pool0) - cov)
        steps = j['steps'] - k_cut
        if steps < 1 or len(pool) < 3:
            print(f'{j["tag"]}: nothing left to re-plan (steps {steps}, pool {len(pool)})'); continue
        pre = [str(pathlib.Path(h).relative_to(ROOT)) if h.startswith('/') else h for h in j.get('heads', [])] + tf[:k_cut]
        jt = f'{j["tag"]}{tag}'
        jobs.append(dict(tag=jt, out=f'results/s15/cols/{tag}/{jt}', pool=pool, S=sorted(set(j.get('S', [])) & set(pool)),
                         p=a.p if a.p is not None else j.get('p', 0.015), steps=steps, beam=20, tries=60, nproc=a.nproc,
                         min_save=3, w_fuel=a.w_fuel, fleet_pre=pre, src_chain=j['tag'], k_cut=k_cut + 1))
        print(f'{jt}: chain {j["tag"]} cut before step {k_cut + 1} (production took {[t["n"] for t in C["taken"]]}); '
              f'twin pool {len(pool)}, steps {steps}, S {len(jobs[-1]["S"])}; fleet_pre {len(pre)} routes')
    json.dump(jobs, open(a.out, 'w'), indent=1)
    print(f'{len(jobs)} hybrid jobs -> {a.out}')


if __name__ == '__main__':
    main()
