"""Stage 15 [TT]: grow a TAIL TREE from finished twin-tail chains (tools/s15_twintail.py) behind fixed heads.

For each source chain (2 twin steps behind a pre-taken tail-1 or a fleet_pre), with M = the leftovers of its fleet:
  C-jobs  step-2 alternatives: for the top step-1 front states (score-best at the deepest two depths + the lightest
          deepest), one twin step over the residual that state leaves (beam --beam2, tries --tries2);
  D-jobs  premium on the leftovers: both steps re-planned with S = ISO-in-pool + M at premium --pd.
Fleet assembly uses the chain's tail1 / fleet_pre (+ 'pre' = the step-1 state for C-jobs): tools/s15_twinfleet.py.

usage: s15_tailtree.py OUT_JOBS.json SRC_JOBDIR [SRC_JOBDIR ...] [--heads results/s14/probe_cover/N9/fleet] [--tag t04]"""
import os, sys, json, pathlib, argparse
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
sys.path.insert(0, str(ROOT / 'tools'))
os.chdir(ROOT)
import numpy as np

ISO = [64, 119, 137, 138, 157, 158, 168, 216]


def targets(f):
    return set(int(x) for x in np.load(ROOT / f)['asts'])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('src', nargs='+')
    ap.add_argument('--heads', default='results/s14/probe_cover/N9/fleet'); ap.add_argument('--head-names', default='01,02,03,04,05,06')
    ap.add_argument('--tag', default='t04'); ap.add_argument('--w-fuel', type=float, default=1.5)
    ap.add_argument('--beam2', type=int, default=30); ap.add_argument('--tries2', type=int, default=90)
    ap.add_argument('--pd', type=float, default=0.05); ap.add_argument('--nc', type=int, default=3)
    ap.add_argument('--no-d', action='store_true'); ap.add_argument('--no-c', action='store_true')
    ap.add_argument('--stones', action='store_true', help='also an S-job: last step with stepping stones')
    a = ap.parse_args()
    from greedy_cover import ALL
    jobs = []
    for sd in a.src:
        m = json.load(open(ROOT / sd / 'meta.json')); J = m['job']; nm = pathlib.Path(sd).name
        if J.get('fleet_pre'):
            fixed = list(J['fleet_pre'])
        else:
            fixed = [str(pathlib.Path(a.heads) / f'route_{n}.npz') for n in a.head_names.split(',')]
            if J.get('tail1'):
                fixed.append(J['tail1'])
            fixed += list(J.get('pre') or [])
        cov = set()
        for f in fixed:
            cov |= targets(f)
        pool = sorted(set(ALL) - cov)                     # the chain's own pool (what the fixed part leaves)
        taken = m['chain'].get('taken', [])
        cov_all = set(cov)
        for t in taken:
            cov_all |= set(t['targets'])
        M = sorted(set(ALL) - cov_all)
        S0 = sorted(set(ISO) & set(pool))
        base = dict(p=0.03, beam=20, tries=60, nproc=1, min_save=3, w_fuel=a.w_fuel)
        extra = {k: J[k] for k in ('tail1', 'fleet_pre') if J.get(k)}
        pre0 = list(J.get('pre') or [])
        # C-jobs: step-2 alternatives from the step-1 front
        s1 = sorted([c for c in m['cols'] if c.get('step') == 1], key=lambda c: (-c['n'], c['tank']))
        if s1 and not a.no_c:
            dmax = s1[0]['n']
            cand = [c for c in s1 if c['kind'] == 'score' and c['n'] >= dmax - 1] + [c for c in s1 if c['kind'] == 'light' and c['n'] == dmax]
            seen = set(); k = 0
            for c in cand:
                key = frozenset(c['targets'])
                if key in seen or k >= a.nc:
                    continue
                seen.add(key)
                p2 = sorted(set(pool) - set(c['targets']))
                jt = f'{nm}C{k}'
                jobs.append(dict(base, tag=jt, out=f'results/s15/cols/{a.tag}/{jt}', pool=p2, S=sorted(set(ISO) & set(p2)), steps=1,
                                 beam=a.beam2, tries=a.tries2, pre=pre0 + [c['f']], src=nm, **extra))
                print(f'{jt}: step-1 state {c["n"]} fb @ {c["tank"]:.0f} ({c["kind"]}); step-2 pool {len(p2)}')
                k += 1
        # S-job: the LAST step again with stepping stones (the other tails' targets, prize eps) and a wider launch grid
        # (TB2f: a 16-target last pool has 36 launch legs -> 10 roots -> the beam dies at depth 3)
        if a.stones and taken:
            last = taken[-1] if len(taken) >= 2 else None
            prev = taken[:-1] if len(taken) >= 2 else taken
            cov_prev = set(cov)
            for t in prev:
                cov_prev |= set(t['targets'])
            p3 = sorted(set(ALL) - cov_prev)
            st_ = set()
            for f in fixed[len(a.head_names.split(',')):]:
                st_ |= targets(f)
            for t in prev:
                st_ |= set(t['targets'])
            jt = f'{nm}S'
            prev_f = [t['f'] for t in prev]
            jobs.append(dict(base, tag=jt, out=f'results/s15/cols/{a.tag}/{jt}', pool=p3, S=sorted((set(ISO) | set(M)) & set(p3)), p=0.05,
                             steps=1, pre=pre0 + prev_f, stones=sorted(st_ - set(p3)), stone_prize=0.01, grid='0,1500,20',
                             src=nm, **extra))
            print(f'{jt}: last step over {len(p3)} with {len(st_ - set(p3))} stepping stones (replaces {None if last is None else last["n"]} fb)')
        # D-job: premium on the leftovers
        if not a.no_d and M:
            jt = f'{nm}D'
            jobs.append(dict(base, tag=jt, out=f'results/s15/cols/{a.tag}/{jt}', pool=pool, S=sorted(set(S0) | set(M)), p=a.pd,
                             steps=len(taken) or 2, pre=pre0, src=nm, **extra))
            print(f'{jt}: pool {len(pool)}, premium on leftovers {M} + iso {S0} at p {a.pd}')
    json.dump(jobs, open(a.out, 'w'), indent=1)
    print(f'{len(jobs)} jobs -> {a.out}')


if __name__ == '__main__':
    main()
