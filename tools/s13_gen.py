"""Stage 13 (PCC-8G) [G]: the generation loop.  docs/stage13_solver_plan.md sections 3.1, 3.2 #11, 4.3 and 6.

usage: s13_gen.py results/s13/gen --gens 6 --slots 8 \
           --seed-fleets results/s13/plan13/tail_N6a_k5/fleet,results/s13/seed_passes/N6a [--budget 4] [--slot-nproc 1]
           [--jitter 0.004] [--dry-run]

Generation g (state in OUT/gen.json; every step is skipped on restart once recorded, slots resume themselves):
  slots   4 full passes (premiums jittered +-0.004, seeded), 2 tail re-runs of the current best N=8 fleet (--keep 01..05
          and --keep 01..04), 2 pilot passes (--pilot 4) from generation 2 on (generation 1: 2 more full passes).
          One tools/s13_pass.py process per slot under `nice -n 5`; slots run in waves so that the sum of their
          --nproc never exceeds --budget (4 cores: 4 slots x --nproc 1 by default).
  (a) master   tools/s13_master.py over the archive (--N 8,9 --exact);
  (b) probe    price probe of the best N=8 fleet's leftovers: one RI.candidates (radius 0.30, 6 per target) +
               RI.w_insert round, NOTHING applied; each leftover is cheap (<= 0.08 J), dear, or unplaceable;
  (c) premiums closability-aware update: +0.006 unplaceable, +0.002 dear, 0 cheap, -0.002 multiply covered, lambda
               clipped to [-0.01, 0.02] (s13_pass then clips prizes to [1.00, 1.02]);
  (d) log      best N=8 (covered, sum J_i), projected closed J (ledger = sum J_i + sum of best leftover prices, 1 per
               unplaceable, + 2), LP bound.
Gates (section 6): after generation 3, G0.5 = fleet >= 278 at <= 10.55 or projected < 13.9, else STOP;
G1 (a) >= 290 at <= 10.75 with every leftover <= 0.08 J, (b) >= 282 at <= 10.6 with no unplaceable -> stop (handoff)."""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, signal, pathlib, argparse, subprocess, multiprocessing as mp
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
_CWD0 = os.getcwd()
os.chdir(ROOT)
import numpy as np
import run_ialns as RI
from greedy_cover import ALL, _timeout

PY = str(pathlib.Path.home() / '.venvs/astro313/bin/python')
ENV = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1', MKL_NUM_THREADS='1')
CHEAP = 0.08
BUDGETS = '10.55,10.6,10.75'           # capacity frontier of the gates (G0.5, G1(b), G1(a))
ARCHIVE = ['results/s13/gen/**/route_*.npz', 'results/s13/g0/**/route_*.npz', 'results/s13/plan13/**/route_*.npz',
           'results/s13/seed_passes/**/route_*.npz']


def rel(p):
    p = pathlib.Path(p)
    return str(p.relative_to(ROOT)) if str(p).startswith(str(ROOT)) else str(p)


# ------------------------------------------------------------------------------------------------ premiums / probe
def prem_initial(delta=0.01):
    import s13_pass as SP
    return {str(t): delta for t in sorted(SP.hard_set())}


def prem_update(lam, probe, fleet_cov):
    """Closability-aware rule (plan 3.2 #11)."""
    lam = {int(k): float(v) for k, v in lam.items()}
    for t in probe['unplaceable']:
        lam[t] = lam.get(t, 0.0) + 0.006
    for t in probe['dear']:
        lam[t] = lam.get(t, 0.0) + 0.002
    for t, hosts in fleet_cov.items():
        if len(hosts) > 1:
            lam[t] = lam.get(t, 0.0) - 0.002
    return {str(t): round(float(np.clip(v, -0.01, 0.02)), 6) for t, v in sorted(lam.items()) if abs(v) > 1e-12}


def w_insert_t(job):
    """RI.w_insert under the plan's 240 s SIGALRM."""
    try:
        signal.signal(signal.SIGALRM, _timeout); signal.setitimer(signal.ITIMER_REAL, 240.0)
        r = RI.w_insert(job)
    except Exception as e:
        r = dict(name=job[0], ast=job[2], t=job[3], ok=False, miss=float('inf'), error=type(e).__name__)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    r.pop('st', None)
    return r


def price_probe(fleet_dir, nproc=4, radius=0.30, per_target=6, say=print, max_targets=0):
    """One candidates + w_insert round over the fleet's leftovers; nothing is applied.  (max_targets > 0: smoke test on
    the first few leftovers only.)"""
    fl = RI.IFleet(fleet_dir); cov = fl.coverage()
    left = sorted(set(ALL) - set(cov)); tic = time.time()
    if max_targets > 0:
        left = left[:max_targets]
    with mp.get_context('fork').Pool(nproc, maxtasksperchild=20) as pool:
        by = RI.candidates(fl, pool, left, radius, per_target)
        jobs = [(c['host'], fl.routes[c['host']]['st'], c['ast'], c['t']) for u in left for c in by.get(u, [])]
        res = pool.map(w_insert_t, jobs, chunksize=1) if jobs else []
    best = {}
    for r in res:
        if not r.get('ok'):
            continue
        dJ = RI.cost(r['tank']) - RI.cost(fl.routes[r['name']]['tank'])
        if r['ast'] not in best or dJ < best[r['ast']]['dJ']:
            best[r['ast']] = dict(dJ=round(dJ, 5), host=r['name'], t=r['t'], tank=round(r['tank'], 1))
    cheap = sorted(u for u in left if u in best and best[u]['dJ'] <= CHEAP)
    dear = sorted(u for u in left if u in best and best[u]['dJ'] > CHEAP)
    unpl = sorted(u for u in left if u not in best)
    sj = sum(RI.cost(r['tank']) for r in fl.routes.values())
    ledger = sj + sum(best[u]['dJ'] for u in cheap + dear) + len(unpl) + 2.0
    out = dict(fleet=rel(fleet_dir), leftovers=left, n_cand=len(jobs), n_ok=sum(1 for r in res if r.get('ok')),
               best={str(u): v for u, v in best.items()}, cheap=cheap, dear=dear, unplaceable=unpl,
               sumJi=round(sj, 4), covered=len(cov), projected_J=round(ledger, 4), sec=round(time.time() - tic))
    say(f'  probe: {len(left)} leftovers, {len(jobs)} insertion trials ({out["n_ok"]} settled): cheap {len(cheap)}, '
        f'dear {len(dear)}, unplaceable {len(unpl)}; projected closed J {ledger:.3f} ({out["sec"]} s)')
    return out


# ---------------------------------------------------------------------------------------------------- scheduler
def _alive(pid):
    try:
        os.kill(pid, 0); return True
    except OSError:
        return False


def run_slots(slots, budget, say, dry=False, poll=30.0):
    """Run slot commands (dicts: name, dir, cmd, nproc) with sum(nproc of running) <= budget.  A slot is skipped when
    its pass.json says done; a slot whose recorded pid is still alive (orphan of a killed generation) is waited for."""
    pending = []; running = {}
    for s in slots:
        pj = ROOT / s['dir'] / 'pass.json'
        if pj.exists() and json.load(open(pj)).get('done'):
            say(f'  slot {s["name"]}: done'); continue
        pidf = ROOT / s['dir'] / '.pid'
        if pidf.exists():
            try:
                pid = int(pidf.read_text().split()[0])
            except Exception:
                pid = -1
            if pid > 0 and _alive(pid):
                say(f'  slot {s["name"]}: still running as pid {pid} (orphan) -- waiting'); running[s['name']] = (None, s, pid); continue
        pending.append(s)
    if dry:
        for s in pending:
            say('  DRY ' + ' '.join(s['cmd']))
        return
    while pending or running:
        used = sum(s['nproc'] for _, s, _ in running.values())
        while pending and used + pending[0]['nproc'] <= budget:
            s = pending.pop(0); d = ROOT / s['dir']; d.mkdir(parents=True, exist_ok=True)
            fo = open(d.parent / f'{d.name}.out', 'a')
            p = subprocess.Popen(['nice', '-n', '5'] + s['cmd'], cwd=ROOT, env=ENV, stdout=fo, stderr=subprocess.STDOUT)
            (d / '.pid').write_text(f'{p.pid} {time.time():.0f}\n')
            running[s['name']] = (p, s, p.pid); used += s['nproc']
            say(f'  slot {s["name"]}: started pid {p.pid} ({s["nproc"]} proc): {" ".join(s["cmd"][2:])}')
        time.sleep(poll)
        for name in list(running):
            p, s, pid = running[name]
            rc = p.poll() if p is not None else (None if _alive(pid) else 0)
            if rc is not None:
                pj = ROOT / s['dir'] / 'pass.json'
                summ = json.load(open(pj)).get('summary') if pj.exists() else None
                say(f'  slot {name}: finished rc {rc}: {summ}')
                del running[name]


# ----------------------------------------------------------------------------------------------------------- loop
def best_n8(master_dir):
    """The 'best N=8 fleet' of a master run: max coverage within the J<13 budget (frontier 10.75), else J-optimal."""
    m = json.load(open(ROOT / master_dir / 'master.json'))
    return (m.get('frontier', {}).get('8', {}).get('10.75')) or m['N'].get('8')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out')
    ap.add_argument('--gens', type=int, default=6); ap.add_argument('--slots', type=int, default=8)
    ap.add_argument('--seed-fleets', default='results/s13/plan13/tail_N6a_k5/fleet,results/s13/seed_passes/N6a')
    ap.add_argument('--budget', type=int, default=4); ap.add_argument('--slot-nproc', type=int, default=1)
    ap.add_argument('--probe-nproc', type=int, default=4); ap.add_argument('--jitter', type=float, default=0.004)
    ap.add_argument('--delta', type=float, default=0.01); ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--pass-args', default='', help='extra s13_pass.py arguments for every slot')
    ap.add_argument('--probe-radius', type=float, default=0.30); ap.add_argument('--probe-per-target', type=int, default=6)
    ap.add_argument('--probe-max-targets', type=int, default=0, help='smoke tests only: probe the first K leftovers')
    a = ap.parse_args()
    if a.budget > 4 or a.probe_nproc > 4 or a.slot_nproc > a.budget:
        raise SystemExit('CPU budget: at most 4 cores')
    out = pathlib.Path(a.out) if os.path.isabs(a.out) else pathlib.Path(_CWD0) / a.out
    out.mkdir(parents=True, exist_ok=True); (out / 'prem').mkdir(exist_ok=True)
    say = RI.logger(out / 'log.txt')
    gj = out / 'gen.json'
    G = json.load(open(gj)) if gj.exists() else dict(gens={}, started=time.strftime('%Y-%m-%d %H:%M:%S'))
    def save():
        if a.dry_run:
            return
        json.dump(G, open(gj.with_suffix('.tmp'), 'w'), indent=1); os.replace(gj.with_suffix('.tmp'), gj)
    say(f's13_gen: {vars(a)}; resume {sorted(G["gens"])}')
    if G.get('stopped'):
        say(f'already stopped: {G["stopped"]}'); return
    pf = out / 'prem' / 'PREM_01.json'
    if not pf.exists():
        json.dump({'lambda': prem_initial(a.delta), 'from': f'delta {a.delta} on the gc16 hard set'}, open(pf, 'w'), indent=1)
    # ---- generation 0: the master over the seeds gives the first "current best" fleet
    if 'best0' not in G:
        m0 = out / 'g00' / 'master'
        pats = [f'{s}/route_*.npz' for s in a.seed_fleets.split(',')]
        cmd = [PY, 'tools/s13_master.py', rel(m0)] + pats + ['--N', '8', '--exact', '--budgets', BUDGETS,
                                                             '--nproc', str(min(4, a.budget))]
        if a.dry_run:
            say('  DRY ' + ' '.join(cmd)); G['best0'] = None
        else:
            subprocess.run(['nice', '-n', '5'] + cmd, cwd=ROOT, env=ENV, check=True, stdout=subprocess.DEVNULL)
            b = best_n8(rel(m0)); G['best0'] = dict(fleet=b['fleet'], covered=b['covered'], sumJi=b['sumJi'])
            save()
    best = G['best0']
    for g in range(1, a.gens + 1):
        key = f'{g:02d}'; rec = G['gens'].setdefault(key, {})
        if rec.get('done'):
            best = rec['best']; continue
        pfile = out / 'prem' / f'PREM_{key}.json'
        gd = out / f'g{key}'
        base = [PY, 'tools/s13_pass.py']
        extra = a.pass_args.split() if a.pass_args else []
        common = ['--prem', rel(pfile), '--jitter', str(a.jitter), '--nproc', str(a.slot_nproc)] + extra
        slots = []
        nfull = 4 if g >= 2 else 6
        for i in range(nfull):
            d = rel(gd / f'full{i + 1}')
            slots.append(dict(name=f'g{key}/full{i + 1}', dir=d, nproc=a.slot_nproc,
                              cmd=base + [d, '--seed', str(1000 * g + i + 1)] + common))
        if best:
            for i, keep in enumerate(('01,02,03,04,05', '01,02,03,04')):
                d = rel(gd / f'tail{len(keep.split(","))}')
                slots.append(dict(name=f'g{key}/tail{len(keep.split(","))}', dir=d, nproc=a.slot_nproc,
                                  cmd=base + [d, '--start', best['fleet'], '--keep', keep, '--seed', str(1000 * g + 11 + i)] + common))
        if g >= 2:
            for i in range(2):
                d = rel(gd / f'pilot{i + 1}')
                slots.append(dict(name=f'g{key}/pilot{i + 1}', dir=d, nproc=a.slot_nproc,
                                  cmd=base + [d, '--pilot', '4', '--seed', str(1000 * g + 21 + i)] + common))
        slots = slots[:a.slots]
        say(f'=== generation {g}: {len(slots)} slots; premiums {rel(pfile)}; best so far {best}')
        rec['slots'] = [s['name'] for s in slots]; save()
        run_slots(slots, a.budget, say, dry=a.dry_run)
        if a.dry_run:
            say('dry run: stop after listing generation 1'); return
        # (a) master
        if 'master' not in rec:
            md = gd / 'master'
            pats = [f'{rel(out)}/**/route_*.npz'] + [p for p in ARCHIVE if not p.startswith(rel(out) + '/')]
            cmd = [PY, 'tools/s13_master.py', rel(md)] + pats + ['--N', '8,9', '--exact', '--strip', '--budgets', BUDGETS,
                                                                 '--nproc', str(min(4, a.budget))]
            subprocess.run(['nice', '-n', '5'] + cmd, cwd=ROOT, env=ENV, check=True, stdout=subprocess.DEVNULL)
            m = json.load(open(md / 'master.json'))
            keys = ('covered', 'sumJi', 'J', 'lp', 'fleet', 'dup_targets', 'craft')
            b8 = m['N'].get('8'); b9 = m['N'].get('9')
            rec['master'] = dict(dir=rel(md), N8={k: b8.get(k) for k in keys}, N9=None if not b9 else {k: b9.get(k) for k in keys},
                                 frontier={B: {k: v.get(k) for k in keys} for B, v in m.get('frontier', {}).get('8', {}).items()})
            save()
            say(f'  master: J-optimal N8 {rec["master"]["N8"]["covered"]} @ {rec["master"]["N8"]["sumJi"]} (J {rec["master"]["N8"]["J"]}, '
                f'LP {rec["master"]["N8"]["lp"]}); frontier ' +
                ', '.join(f'<= {B}: {v["covered"]} @ {v["sumJi"]}' for B, v in rec['master']['frontier'].items()) +
                (f'; N9 {rec["master"]["N9"]["covered"]} @ {rec["master"]["N9"]["sumJi"]}' if rec['master']['N9'] else ''))
        fr = rec['master'].get('frontier') or {}
        b8 = fr.get('10.75') or rec['master']['N8']          # the "best N=8 fleet": max coverage within the J<13 budget
        # (b) price probe of the best N=8 fleet
        if 'probe' not in rec:
            pr = price_probe(ROOT / b8['fleet'], nproc=a.probe_nproc, radius=a.probe_radius, per_target=a.probe_per_target,
                             say=say, max_targets=a.probe_max_targets)
            json.dump(pr, open(gd / 'probe.json', 'w'), indent=1)
            rec['probe'] = {k: pr[k] for k in ('covered', 'sumJi', 'projected_J', 'cheap', 'dear', 'unplaceable')}
            save()
        pr = rec['probe']
        # (c) premium update
        nfile = out / 'prem' / f'PREM_{g + 1:02d}.json'
        if not nfile.exists():
            lam = json.load(open(pfile))['lambda']
            if b8.get('dup_targets') is not None:           # multiply covered in the master's pick (before any strip)
                cov = {int(t): {'a', 'b'} for t in b8['dup_targets']}
            else:
                cov = RI.IFleet(ROOT / b8['fleet']).coverage()
            new = prem_update(lam, pr, cov)
            json.dump({'lambda': new, 'from': rel(pfile), 'probe': rel(gd / 'probe.json')}, open(nfile, 'w'), indent=1)
        # (d) log + gates
        best = dict(fleet=b8['fleet'], covered=b8['covered'], sumJi=b8['sumJi'])
        rec.update(best=best, done=True, finished=time.strftime('%Y-%m-%d %H:%M:%S')); save()
        lp = rec['master']['N8'].get('lp')
        say(f'=== generation {g} done: best N=8 (frontier <= 10.75) {b8["covered"]} covered at sum J_i {b8["sumJi"]:.4f} '
            f'(J {b8["J"]:.4f}; LP bound of the J-optimal N8 {lp}); projected closed J {pr["projected_J"]:.3f}; leftovers cheap '
            f'{len(pr["cheap"])}, dear {len(pr["dear"])}, unplaceable {len(pr["unplaceable"])}')
        c, s = b8['covered'], b8['sumJi']
        f106 = fr.get('10.6') or b8; f1055 = fr.get('10.55') or b8
        if c >= 290 and s <= 10.75 and not pr['dear'] and not pr['unplaceable']:
            G['stopped'] = f'G1(a) passed at generation {g}: J<13 attempt'; save(); say(G['stopped']); return
        if f106['covered'] >= 282 and f106['sumJi'] <= 10.6:
            same = f106['covered'] == b8['covered'] and abs(f106['sumJi'] - b8['sumJi']) < 1e-6
            pr6 = pr if same else None
            if pr6 is None:
                if 'probe106' not in rec:
                    p6 = price_probe(ROOT / f106['fleet'], nproc=a.probe_nproc, radius=a.probe_radius,
                                     per_target=a.probe_per_target, say=say, max_targets=a.probe_max_targets)
                    json.dump(p6, open(gd / 'probe106.json', 'w'), indent=1)
                    rec['probe106'] = {k: p6[k] for k in ('covered', 'sumJi', 'projected_J', 'cheap', 'dear', 'unplaceable')}
                    save()
                pr6 = rec['probe106']
            if not pr6['unplaceable']:
                G['stopped'] = f'G1(b) passed at generation {g}: improvement attempt'; save(); say(G['stopped']); return
        if g == 3 and not ((f1055['covered'] >= 278 and f1055['sumJi'] <= 10.55) or pr['projected_J'] < 13.9):
            G['stopped'] = (f'G0.5 FAILED after generation 3: {f1055["covered"]} at {f1055["sumJi"]:.4f} (<= 10.55), '
                            f'projected {pr["projected_J"]:.3f}')
            save(); say(G['stopped']); return
    G['stopped'] = f'{a.gens} generations done without G1'; save(); say(G['stopped'])


if __name__ == '__main__':
    main()
