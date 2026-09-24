"""Stage-13 linchpin, step 2: dissolve routes V1,V2 of a settled fleet and distribute their targets among the survivors'
pools by trajectory proximity, balanced: min sum of closest-approach distance (AU) subject to pool size <= --cap
(survivors keep all their own targets).  Exact capacitated assignment (Hungarian over capacity slots).
Output JSON: {"dissolved": [...], "cap": C, "pools": {route: sorted ids}, "extra": {route: [[id, d_au, t_s], ...]},
              "own": {route: [...]}}
Usage: linchpin_pools.py dist.json out.json --dissolve 05,09 [--cap 38]
"""
import sys, json, argparse
import numpy as np
from scipy.optimize import linear_sum_assignment


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('dist'); ap.add_argument('out')
    ap.add_argument('--dissolve', required=True); ap.add_argument('--cap', type=int, default=38)
    a = ap.parse_args()
    D = json.load(open(a.dist)); R = D['routes']; M = D['minima']
    V = a.dissolve.split(','); surv = sorted(n for n in R if n not in V)
    X = sorted(set(x for v in V for x in R[v]['asts']))
    own = {n: sorted(R[n]['asts']) for n in surv}
    slots = []
    for n in surv:
        slots += [n] * max(0, a.cap - len(own[n]))
    if len(slots) < len(X):
        sys.exit(f'capacity {len(slots)} < {len(X)} targets: raise --cap')
    best = {}
    C = np.zeros((len(X), len(slots)))
    for i, x in enumerate(X):
        for n in surv:
            lst = M[n][str(x)]
            k = int(np.argmin([d for t, d in lst])) if lst else -1
            best[(x, n)] = (lst[k][1], lst[k][0]) if lst else (9.0, float('nan'))
        for j, n in enumerate(slots):
            C[i, j] = best[(x, n)][0]
    ri, cj = linear_sum_assignment(C)
    extra = {n: [] for n in surv}
    for i, j in zip(ri, cj):
        n = slots[j]; d, t = best[(X[i], n)]; extra[n].append([X[i], round(d, 4), t])
    pools = {n: sorted(own[n] + [e[0] for e in extra[n]]) for n in surv}
    tot = sum(len(p) for p in pools.values())
    nearest = {x: min(surv, key=lambda n: best[(x, n)][0]) for x in X}
    json.dump(dict(dissolved=V, cap=a.cap, pools=pools, extra=extra, own=own,
                   nearest_unconstrained={str(x): [nearest[x], round(best[(x, nearest[x])][0], 4)] for x in X}),
              open(a.out, 'w'), indent=1)
    print(f'dissolve {V}: {len(X)} targets -> {len(surv)} survivors, cap {a.cap}, total pool {tot}')
    for n in surv:
        ds = [e[1] for e in extra[n]]
        print(f'  {n}: own {len(own[n])} + {len(extra[n])} = {len(pools[n])}   extra dist ' +
              (f'mean {np.mean(ds):.3f} max {np.max(ds):.3f} AU  ' if ds else '') + str([e[0] for e in extra[n]]))
    got = np.array([C[i, j] for i, j in zip(ri, cj)])
    print(f'assigned distance: mean {got.mean():.3f} median {np.median(got):.3f} max {got.max():.3f} AU; '
          f'<0.03 {int((got < 0.03).sum())}, <0.05 {int((got < 0.05).sum())}, <0.10 {int((got < 0.10).sum())}, '
          f'>=0.15 {int((got >= 0.15).sum())}; moved off their nearest survivor by the cap: '
          f'{sum(1 for i, j in zip(ri, cj) if slots[j] != nearest[X[i]])}')


if __name__ == '__main__':
    main()
