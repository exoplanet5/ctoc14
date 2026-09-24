"""s15b ledger (results/s15b/ledger.json): best CLOSED (298 covered, 9 craft) and best PARTIAL 9-craft fleets + frontier.

usage: s15b_ledger.py partial DIR COVERED SUMJI "MISSES" "NOTE"
       s15b_ledger.py closed  DIR SUMJI "NOTE"            (298 covered; SUMJI = twin sum J_i, re-settled when stated)
       s15b_ledger.py frontier "K:SUMJI,K:SUMJI,..." "SRC"
       s15b_ledger.py event "TEXT"
Partial fleets are kept per coverage level (best sum J_i per covered count); 'best_partial' = the one with the smallest
sum J_i + misses (J up to the constant 2)."""
import sys, json, time, pathlib
L = pathlib.Path('/Users/mickey/solarsystem/ctoc14/results/s15b/ledger.json')


def load():
    return json.load(open(L)) if L.exists() else dict(history=[])


def save(D):
    D['updated'] = time.strftime('%m-%d %H:%M')
    tmp = L.with_suffix('.tmp'); json.dump(D, open(tmp, 'w'), indent=1); tmp.replace(L)


def main():
    D = load(); cmd = sys.argv[1]; now = time.strftime('%m-%d %H:%M')
    if cmd == 'partial':
        d, cov, sj, ms, note = sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), sys.argv[5], sys.argv[6]
        misses = [int(x) for x in ms.replace(',', ' ').split()]
        P = D.setdefault('partials', {})
        cur = P.get(str(cov))
        if cur is None or sj < cur['sumJi'] - 1e-9:
            P[str(cov)] = dict(dir=d, covered=cov, sumJi=sj, misses=misses, note=note, when=now)
            D['history'].append(dict(when=now, event=f'partial {cov} @ {sj:.4f} ({d}; misses {misses}) {note}'))
        best = min(P.values(), key=lambda p: p['sumJi'] + (298 - p['covered']))
        D['best_partial'] = best
    elif cmd == 'closed':
        d, sj, note = sys.argv[2], float(sys.argv[3]), sys.argv[4]
        cur = D.get('best_closed')
        if cur is None or sj < cur['sumJi'] - 1e-9:
            D['best_closed'] = dict(dir=d, covered=298, sumJi=sj, note=note, when=now)
        D['history'].append(dict(when=now, event=f'closed 298 @ {sj:.4f} ({d}) {note}'))
    elif cmd == 'frontier':
        fr = {k: float(v) for k, v in (x.split(':') for x in sys.argv[2].split(','))}
        D['frontier'] = dict(when=now, src=sys.argv[3], K=fr)
        D['history'].append(dict(when=now, event=f'frontier {sys.argv[3]}: ' + ' | '.join(f'{k} @ {v:.4f}' for k, v in fr.items())))
    elif cmd == 'event':
        D['history'].append(dict(when=now, event=sys.argv[2]))
    save(D)
    print(json.dumps({k: D.get(k) for k in ('best_closed', 'best_partial')}, indent=1))


if __name__ == '__main__':
    main()
