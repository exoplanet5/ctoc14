"""Stage 16: EXACT-aware best keeper.  results/s16/best/fleet (twin) + results/s16/efleet (its validated exact crafts).

The twin sum J_i is a proxy; the submitted J is the exact one.  Given 9 twin route files, this tool
  1. finds for every route the lightest VALID exact conversion of that exact state (tools/s16_xconv.py indexes
     results/s16/xconv/*/index.jsonl, matched by content hash s16_iter.chash; any --opt-iters);
  2. assembles the twin fleet OUT_FLEET (r1..r9 by depth, then tank; RI.IFleet) and runs the honest check
     (tools/s16_pipe.honest_check: 0 irregular nodes, 20 d windows <= cap, twin miss <= 150 km, coverage);
  3. if every route has a valid conversion, the fleet is honest with 9 craft / 298 covered, and its EXACT sum J_i is below
     the stored best's (best.json exact_sumJi, or the matching efleet), it replaces results/s16/best/fleet and
     results/s16/efleet (old copies kept under results/s16/hist/<time>/), writes efleet/SOURCE_MD5 (md5 of
     cat best/fleet/route_*.npz, as results/s16/export_s16a.sh checks it) and best.json (twin + exact sums).
Never writes results/CTOC14_Result_*.txt.

usage: s16_xbest.py OUT_FLEET ROUTE_FILE x9 [--note TEXT] [--dry]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, glob, time, shutil, pathlib, argparse, subprocess
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)
import run_ialns as RI
import s15b_fleet as FA
from s16_iter import chash
from s16_pipe import honest_check
from run_gfleet import Fleet

BEST = ROOT / 'results/s16/best'; EF = ROOT / 'results/s16/efleet'; HIST = ROOT / 'results/s16/hist'


def conversions():
    """{chash: [record, ...]} of every valid exact conversion in results/s16/xconv/*/index.jsonl (record + 'dir')."""
    out = {}
    for idx in glob.glob(str(ROOT / 'results/s16/xconv/*/index.jsonl')):
        d = pathlib.Path(idx).parent
        for l in open(idx):
            try:
                r = json.loads(l)
            except Exception:
                continue
            if not r.get('ok') or 'm0' not in r:
                continue
            if not (d / 'frags' / f'frag_{r["key"]}.txt').exists() or not (d / 'crafts' / f'craft_{r["key"]}.npz').exists():
                continue
            r['dir'] = str(d); out.setdefault(r['key'].split('_o')[0], []).append(r)
    return out


def md5_of(fleet_dir):
    return subprocess.run(['zsh', '-c', f'setopt nonomatch; cat {fleet_dir}/route_*.npz | md5 -q'], capture_output=True,
                          text=True).stdout.strip()


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('routes', nargs='+')
    ap.add_argument('--note', default=''); ap.add_argument('--dry', action='store_true')
    a = ap.parse_args()
    BEST.mkdir(parents=True, exist_ok=True); say = RI.logger(BEST / 'log.txt'); RI.eph()
    conv = conversions()
    sts = [(f, FA.load_st(f)) for f in a.routes]
    sts.sort(key=lambda fs: (-len(fs[1]['asts']), float(RI.ipr(fs[1]).tank())))
    F = RI.IFleet(); pick = {}
    for i, (f, st) in enumerate(sts):
        name = f'r{i + 1}'; F.routes[name] = dict(st=st, tank=float(RI.ipr(st).tank()))
        cs = sorted(conv.get(chash(st), []), key=lambda r: r['m0'])
        pick[name] = dict(src=f, tank=F.routes[name]['tank'], conv=cs[0] if cs else None)
    F.save(ROOT / a.out, note=a.note or 's16_xbest')
    rep = honest_check(ROOT / a.out, 9, say)
    missing = [n for n, p in pick.items() if p['conv'] is None]
    ex = None if missing else sum(RI.cost(p['conv']['m0']) for p in pick.values())
    for n, p in sorted(pick.items()):
        c = p['conv']
        say(f'   {n}: twin {p["tank"]:.2f}' + (f' -> exact {c["m0"]:.2f} ({c["m0"] - p["tank"]:+.2f} kg, o{c["opt_iters"]}, hs {c["hs"]}) {c["key"]}'
                                              if c else ' -> NO valid conversion') + f'  [{p["src"]}]')
    old = json.load(open(BEST / 'best.json')) if (BEST / 'best.json').exists() else {}
    old_ex = old.get('exact_sumJi')
    if old_ex is None and (EF / 'fleet.json').exists() and (EF / 'SOURCE_MD5').exists() \
            and open(EF / 'SOURCE_MD5').read().strip() == md5_of(BEST / 'fleet'):
        j = json.load(open(EF / 'fleet.json')); old_ex = j['J'] - (300 - j['covered'])
    ok = rep['honest'] and rep['covered'] == 298 and rep['craft'] == 9 and ex is not None
    say(f'  exact-aware {a.out}: twin {rep["sumJi"]:.5f}, exact {ex if ex is None else round(ex, 5)}, honest {rep["honest"]}, '
        f'covered {rep["covered"]}; stored best exact {old_ex}')
    if not ok or (old_ex is not None and ex >= old_ex - 1e-6) or a.dry:
        say('  not stored' + (' (dry run)' if a.dry else '')); return
    stamp = time.strftime('%m%d_%H%M%S'); (HIST / stamp).mkdir(parents=True, exist_ok=True)
    for d in (BEST / 'fleet', EF):
        if d.exists():
            shutil.copytree(d, HIST / stamp / d.name)
    if (BEST / 'best.json').exists():
        shutil.copy(BEST / 'best.json', HIST / stamp / 'best.json')
    tmp = BEST / 'fleet.tmp'
    if tmp.exists():
        shutil.rmtree(tmp)
    shutil.copytree(ROOT / a.out, tmp)
    if (BEST / 'fleet').exists():
        shutil.rmtree(BEST / 'fleet')
    tmp.rename(BEST / 'fleet')
    etmp = ROOT / 'results/s16/efleet.tmp'
    if etmp.exists():
        shutil.rmtree(etmp)
    (etmp / 'frags').mkdir(parents=True)
    for n, p in pick.items():
        c = p['conv']; d = pathlib.Path(c['dir'])
        shutil.copy(d / 'crafts' / f'craft_{c["key"]}.npz', etmp / f'craft_{n}.npz')
        shutil.copy(d / 'frags' / f'frag_{c["key"]}.txt', etmp / 'frags' / f'frag_{n}.txt')
    ef = Fleet(etmp); ef.save(etmp, note=f'exact-aware conversions of results/s16/best/fleet ({a.out})')
    open(etmp / 'SOURCE_MD5', 'w').write(md5_of(BEST / 'fleet') + '\n')
    if EF.exists():
        shutil.rmtree(EF)
    etmp.rename(EF)
    rep.update(source=a.out, note=a.note, when=time.strftime('%m-%d %H:%M'), exact_sumJi=round(ex, 6),
               exact_J=round(ef.J(), 6), conversions={n: dict(src=p['src'], key=p['conv']['key'], dir=p['conv']['dir'],
                                                              m0=p['conv']['m0'], hs=p['conv']['hs']) for n, p in pick.items()})
    json.dump(rep, open(BEST / 'best.json', 'w'), indent=1)
    say(f'NEW BEST (exact) {ex:.5f} (twin {rep["sumJi"]:.5f}) from {a.out}: {ef.summary()}; previous kept in {HIST / stamp}')


if __name__ == '__main__':
    main()
