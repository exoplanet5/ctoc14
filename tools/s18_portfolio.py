"""Stage 18 / WP3: the PORTFOLIO runner -- parameter sets of s18_rhfa run one after the other on the whole machine.

Each entry of the list is one run: {"name": "p01_base", "args": {"N": 8, "H": 540, "L": 360, "beta": 1.5, ...}}; the
runner starts `s18_rhfa.py run results/s18/<name> --<k> <v> ... --nproc <nproc>` as a child process (nice -n 5,
threads 1, stdout to results/s18/<name>/nohup.out), waits, then s18_eval.report on the finished run and appends one
line to results/s18/portfolio.jsonl: {name, args, covered, sumJi, fuel, J_raw, gates, wall_s, finished_at, status}.
Runs whose fleet/fleet.json exists are skipped (status 'skipped'; a finished run without a row is evaluated and
recorded first); a run with a state.json but no fleet is RESUMED (s18_rhfa resumes by itself).  The runner itself is
resumable: kill it and restart with the same list.  Status of a row: 'done' (fleet.json exists after the child
returned), 'failed' (no fleet.json; rc recorded), 'skipped' (already recorded), 'dry'.
Default portfolio (PORTFOLIO below): the base parameters and the one-factor variations named in docs/stage18_rhfa.md
G1 "one parameter round" (H 450/630, beta 1.0/2.5, quotas 8/14 and 12/18, gamma 0.15/0.5).

usage: s18_portfolio.py run [--list results/s18/portfolio.json] [--nproc 8] [--out results/s18/portfolio.jsonl] [--dry]
       s18_portfolio.py write-default results/s18/portfolio.json
       s18_portfolio.py table [--out results/s18/portfolio.jsonl]
--dry prints the commands and writes nothing (the WP3 acceptance test).
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse, subprocess
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import s18_common as C

BASE = dict(N=8, H=540, L=360, B=8, tries=16, beta=1.5, gamma=0.3, lam=1.0, qf=12, qm=18, m0max=1070, wall_job=2400)
# Fixer 09-24: BASE follows the new s18_common.DEF (m0max 1070 = the bar, lam 1.0 as the starting dual, quotas 12/18,
# wall_job 2400 as a safety only -- the beam now stops by itself after 2 lookahead-only levels, finding 6).
# Integrator 09-24 (G0b, results/s18/g0b*): the default beam (w_t 1.0) + auction (lam 0.15) spends 188 kg per craft by
# day 1461 for 72 flybys (s16a: 94 kg for 70); w_t 0.25 spends 98 kg but reaches 55 flybys.  Fuel is the binding
# axis, so the portfolio varies the beam time weight and the auction fuel price first (p01-p03), then the design's
# one-factor round.  wall_job 720 s caps a job at ~12 min so a 10-slice run stays near 2 h on 8 cores (MEASURED
# slices 290 / 673 / 818 s at wall_job 2400 with the deepest jobs at 800 s by slice 2).
PORTFOLIO = [dict(name='p01_wt05_lam1', args=dict(BASE, w_t=0.5, lam=1.0)),
             dict(name='p02_wt05_lam1_beta1', args=dict(BASE, w_t=0.5, lam=1.0, beta=1.0)),
             dict(name='p03_wt025_lam1', args=dict(BASE, w_t=0.25, lam=1.0)),
             dict(name='p04_base', args=dict(BASE)),
             dict(name='p05_H450', args=dict(BASE, w_t=0.5, lam=1.0, H=450, L=300)),
             dict(name='p06_H630', args=dict(BASE, w_t=0.5, lam=1.0, H=630, L=420)),
             dict(name='p07_beta25', args=dict(BASE, w_t=0.5, lam=1.0, beta=2.5)),
             dict(name='p08_q8_14', args=dict(BASE, w_t=0.5, lam=1.0, qf=8, qm=14)),
             dict(name='p09_q12_18', args=dict(BASE, w_t=0.5, lam=1.0, qf=12, qm=18)),
             dict(name='p10_gamma015', args=dict(BASE, w_t=0.5, lam=1.0, gamma=0.15)),
             dict(name='p11_gamma05', args=dict(BASE, w_t=0.5, lam=1.0, gamma=0.5))]
RHFA = ROOT / 'tools' / 's18_rhfa.py'
THREAD_VARS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS')


def _fmt(v):
    """CLI text of a value: booleans as 0/1, everything else as str (540 -> '540', 1.5 -> '1.5')."""
    return str(int(v)) if isinstance(v, bool) else str(v)


def command(entry, nproc):
    """The argv list of one run: [nice -n 5, PY, tools/s18_rhfa.py, run, results/s18/<name>, --k v ..., --nproc n]
    (underscores in keys become dashes; booleans as 0/1)."""
    argv = list(C.NICE) + [C.PY, str(RHFA), 'run', f'results/s18/{entry["name"]}']
    for k, v in entry.get('args', {}).items():
        if k == 'nproc' or v is None:
            continue
        argv += ['--' + str(k).replace('_', '-'), _fmt(v)]
    argv += ['--nproc', str(int(nproc))]
    return argv


def finished(name):
    """True when the run is complete: state.json says done with stop kind 'horizon' or 'covered' and fleet/fleet.json
    exists (fixer 09-24, reviewer minor finding: a --stop-T run also has a fleet.json but must be RESUMED)."""
    rd = C.run_dir(name)
    if not (rd / 'fleet' / 'fleet.json').exists() or not (rd / 'state.json').exists():
        return False
    try:
        S = C.jload(rd / 'state.json')
    except Exception:
        return False
    return bool(S.get('done')) and (S.get('stop') or {}).get('kind') in ('horizon', 'covered')


def _rows(out_jsonl):
    p = pathlib.Path(out_jsonl); rows = []
    if p.exists():
        for line in p.read_text().splitlines():
            line = line.strip()
            if line:
                try:
                    rows.append(json.loads(line))
                except json.JSONDecodeError:
                    pass
    return rows


def _append(out_jsonl, row):
    p = pathlib.Path(out_jsonl); p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, 'a') as f:
        f.write(json.dumps(row, default=float) + '\n')


def _evaluate(entry, say):
    """s18_eval.report of the run (G2 estimate included, 2 workers); never raises: {covered, sumJi, fuel, J_raw, gates,
    flybys_by_1461d, eval_error}."""
    rd = C.run_dir(entry['name'])
    try:
        import s18_eval as EV
        rep = EV.report(rd, g2=True, nproc=2)
        C.jdump(rep, rd / 'eval.json')
        for line in EV.table(rep).splitlines():
            say('    ' + line)
        return dict(covered=rep['covered'], sumJi=rep['sumJi'], fuel=rep['fuel'], J_raw=rep['J_raw'],
                    gates={k: rep['gates'][k] for k in ('G0b', 'G1', 'G2')}, flybys_by_1461d=rep['flybys_by_1461d'],
                    closure=(rep['closure'] or {}).get('dJ_total'), n_craft=rep['n_craft'], honest_all=rep['honest_all'],
                    eval_error=None)
    except Exception as e:
        say(f'    eval failed: {repr(e)[:200]}')
        return dict(covered=None, sumJi=None, fuel=None, J_raw=None, gates=None, flybys_by_1461d=None, closure=None,
                    n_craft=None, honest_all=None, eval_error=repr(e)[:300])


def run_one(entry, nproc, out_jsonl, say, dry=False):
    """Start (or skip / resume) one run, wait, evaluate, append the row.  Returns the row dict.  dry: print the command,
    return {name, status 'dry'}.  The child's environment gets OMP/OPENBLAS/VECLIB/MKL threads = 1."""
    name = entry['name']; argv = command(entry, nproc)
    if dry:
        print(' '.join(argv))
        return dict(name=name, status='dry', args=entry.get('args', {}))
    rd = C.run_dir(name)
    recorded = {r.get('name'): r for r in _rows(out_jsonl) if r.get('status') in ('done', 'finished')}
    if finished(name):
        if name in recorded:
            say(f'{name}: finished and recorded -> skipped')
            return dict(recorded[name], status='skipped')
        say(f'{name}: finished (fleet/fleet.json exists) but not recorded -> evaluating')
        row = dict(name=name, args=entry.get('args', {}), **_evaluate(entry, say), wall_s=None, finished_at=C.now(),
                   status='finished', rc=None)
        _append(out_jsonl, row)
        return row
    resume = (rd / 'state.json').exists()
    rd.mkdir(parents=True, exist_ok=True)
    say(f'{name}: {"RESUME" if resume else "START"}  ' + ' '.join(argv))
    env = dict(os.environ); env.update({v: '1' for v in THREAD_VARS})
    tic = time.time()
    with open(rd / 'nohup.out', 'a') as fo:
        fo.write(f'\n===== portfolio {C.now()} {"resume" if resume else "start"}: {" ".join(argv)}\n'); fo.flush()
        rc = subprocess.call(argv, cwd=str(ROOT), env=env, stdout=fo, stderr=subprocess.STDOUT)
    wall = time.time() - tic
    ok = finished(name)
    say(f'{name}: child rc {rc} after {wall / 60:.1f} min; fleet.json {"present" if ok else "MISSING"}')
    row = dict(name=name, args=entry.get('args', {}), **(_evaluate(entry, say) if ok else dict(
        covered=None, sumJi=None, fuel=None, J_raw=None, gates=None, flybys_by_1461d=None, closure=None, n_craft=None,
        honest_all=None, eval_error='no fleet.json')), wall_s=round(wall, 1), finished_at=C.now(),
        status='done' if ok else 'failed', rc=rc)
    _append(out_jsonl, row)
    return row


def run(lst, nproc, out_jsonl, dry=False):
    """Sequential over the list; skips finished; logs to results/s18/portfolio.log."""
    if dry:
        for e in lst:
            run_one(e, nproc, out_jsonl, print, dry=True)
        return []
    say = C.logger(C.S18 / 'portfolio.log')
    say(f'portfolio: {len(lst)} entries, nproc {nproc}, rows -> {out_jsonl}')
    rows = []
    for e in lst:
        rows.append(run_one(e, nproc, out_jsonl, say))
        r = rows[-1]
        say(f'  row {r["name"]}: status {r["status"]} covered {r.get("covered")} sumJi {r.get("sumJi")} fuel {r.get("fuel")} '
            f'J_raw {r.get("J_raw")} gates {r.get("gates")}')
    say('portfolio done\n' + table(out_jsonl))
    return rows


def table(out_jsonl):
    """Printable table of portfolio.jsonl sorted by (covered desc, sumJi asc)."""
    rows = _rows(out_jsonl)
    if not rows:
        return f'(no rows in {out_jsonl})'
    latest = {}
    for r in rows:
        latest[r.get('name')] = r
    rows = sorted(latest.values(), key=lambda r: (-(r.get('covered') or -1), r.get('sumJi') if r.get('sumJi') is not None else 1e9))
    L = [f'{"name":14s} {"status":8s} {"cov":>4s} {"sumJi":>8s} {"fuel":>7s} {"J_raw":>8s} {"fb1461":>6s} {"G0b":>5s} {"G1":>5s} {"G2":>7s} {"wall_min":>8s}  args']
    for r in rows:
        g = r.get('gates') or {}
        args = ' '.join(f'{k}={v}' for k, v in (r.get('args') or {}).items() if BASE.get(k) != v) or 'base'
        L.append(f'{str(r.get("name")):14s} {str(r.get("status")):8s} {str(r.get("covered", "")):>4s} '
                 f'{(f"{r["sumJi"]:.4f}" if r.get("sumJi") is not None else ""):>8s} '
                 f'{(f"{r["fuel"]:.0f}" if r.get("fuel") is not None else ""):>7s} '
                 f'{(f"{r["J_raw"]:.4f}" if r.get("J_raw") is not None else ""):>8s} {str(r.get("flybys_by_1461d", "")):>6s} '
                 f'{str(g.get("G0b", "")):>5s} {str(g.get("G1", "")):>5s} {str(g.get("G2", "")):>7s} '
                 f'{(f"{r["wall_s"] / 60:.0f}" if r.get("wall_s") else ""):>8s}  {args}')
    return '\n'.join(L)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('cmd', choices=['run', 'write-default', 'table']); ap.add_argument('path', nargs='?')
    ap.add_argument('--list', default=str(C.S18 / 'portfolio.json')); ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--out', default=str(C.S18 / 'portfolio.jsonl')); ap.add_argument('--dry', action='store_true')
    a = ap.parse_args()
    if a.cmd == 'write-default':
        C.jdump(PORTFOLIO, a.path or a.list); print('written', a.path or a.list)
    elif a.cmd == 'run':
        lst = C.jload(a.list) if pathlib.Path(a.list).exists() else PORTFOLIO
        run(lst, a.nproc, a.out, a.dry)
    else:
        print(table(a.out))


if __name__ == '__main__':
    main()
