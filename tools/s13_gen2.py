"""Stage 13 track B, stage 2: generation driver around tools/s13_gen.py (imported and run unchanged).

usage: s13_gen2.py OUT --gens 3 [--alt full2,full4,full6,tail4,pilot2 --alt-args="--choice value"]
                   [--deadline "2026-09-23 06:00"] [s13_gen.py options, e.g. --pass-args="--members ..."]

What it adds to s13_gen.py, without editing it:
  * --alt / --alt-args: extra s13_pass.py arguments for the named slots of every generation (slot base names such as
    full2, tail4, pilot2).  Used to run half of each generation with the kappa-value choice rule and half with the
    plan-literal deepest-first (lex) rule, so the master sees both kinds of tail.
  * generation-by-generation driving with a G0.5 decision stop (plan section 6): after each generation the G0.5 test is
    evaluated on s13_gen's own record -- PASS when the master's capacity frontier at sum J_i <= 10.55 reaches >= 278
    covered, or the price-probe ledger projects a closed J < 13.9.  PASS -> stop (handoff to G1);
    s13_gen itself stops with 'G0.5 FAILED' after generation 3, and on G1(a)/G1(b).
  * --deadline: no new generation is started when the mean wall time of the finished generations would overrun it.
  * --pass-tool (default tools/s13_pass2.py): the pass tool of every slot; s13_pass2 = s13_pass with premium-only
    jitter and a 24-level candidate front (see its docstring for the measurements behind both fixes).
State: OUT/gen.json (s13_gen's, untouched except the benign 'N generations done without G1' stop marker, which is
removed before the next generation is driven) and OUT/gen2.json (per-generation wall times and gate evaluations).
Restart-safe: re-running the same command continues where it stopped (s13_gen resumes slots and steps)."""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pathlib, argparse, datetime
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
sys.path.insert(0, str(ROOT / 'tools'))
_CWD0 = os.getcwd()
import s13_gen as SG                      # chdirs to ROOT

BENIGN = 'generations done without G1'


def g05(rec):
    """G0.5 on one finished generation record of s13_gen -> (passed, text)."""
    m = rec['master']; fr = m.get('frontier') or {}
    f = fr.get('10.55') or m['N8']; pr = rec['probe']
    cap = f['covered'] >= 278 and f['sumJi'] <= 10.55
    led = pr['projected_J'] < 13.9
    txt = (f'frontier <= 10.55: {f["covered"]} @ {f["sumJi"]:.4f} (need >= 278); projected closed J '
           f'{pr["projected_J"]:.3f} (need < 13.9)')
    return (cap or led), txt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out'); ap.add_argument('--gens', type=int, default=3)
    ap.add_argument('--alt', default=''); ap.add_argument('--alt-args', default='')
    ap.add_argument('--deadline', default='')
    ap.add_argument('--pass-tool', default='tools/s13_pass2.py', help='pass tool run in every slot (s13_gen: tools/s13_pass.py)')
    a, rest = ap.parse_known_args()
    out = pathlib.Path(a.out) if os.path.isabs(a.out) else pathlib.Path(_CWD0) / a.out
    out.mkdir(parents=True, exist_ok=True)
    say = SG.RI.logger(out / 'log.txt')
    alt = {x.strip() for x in a.alt.split(',') if x.strip()}; alt_args = a.alt_args.split()
    orig = SG.run_slots

    def run_slots(slots, budget, say_, dry=False, poll=30.0):
        for s in slots:
            s['cmd'] = [a.pass_tool if x == 'tools/s13_pass.py' else x for x in s['cmd']]
            if s['name'].split('/')[-1] in alt and alt_args and s['cmd'][-len(alt_args):] != alt_args:
                s['cmd'] = s['cmd'] + alt_args
        return orig(slots, budget, say_, dry=dry, poll=poll)
    SG.run_slots = run_slots

    gj = out / 'gen.json'; g2 = out / 'gen2.json'
    S = json.load(open(g2)) if g2.exists() else dict(args=vars(a), rest=rest, gens={})

    def save2():
        json.dump(S, open(g2.with_suffix('.tmp'), 'w'), indent=1); os.replace(g2.with_suffix('.tmp'), g2)
    say(f's13_gen2: {vars(a)}; s13_gen args {rest}')
    dl = datetime.datetime.strptime(a.deadline, '%Y-%m-%d %H:%M').timestamp() if a.deadline else None
    for g in range(1, a.gens + 1):
        key = f'{g:02d}'
        G = json.load(open(gj)) if gj.exists() else {}
        st = G.get('stopped', '')
        if st and not st.endswith(BENIGN):
            say(f's13_gen2: stopped: {st}'); break
        if S.get('stopped'):
            say(f's13_gen2: stopped: {S["stopped"]}'); break
        done = G.get('gens', {}).get(key, {}).get('done')
        if not done:
            walls = [v['wall_s'] for v in S['gens'].values() if v.get('wall_s')]
            if dl is not None and walls and time.time() + sum(walls) / len(walls) > dl:
                S['stopped'] = f'deadline {a.deadline}: generation {g} would overrun (mean {sum(walls) / len(walls) / 3600:.2f} h)'
                save2(); say(S['stopped']); break
            if st:
                G.pop('stopped'); json.dump(G, open(gj.with_suffix('.tmp'), 'w'), indent=1); os.replace(gj.with_suffix('.tmp'), gj)
            S['gens'].setdefault(key, {})['start'] = S['gens'].get(key, {}).get('start') or time.time(); save2()
            sys.argv = ['s13_gen.py', str(out), '--gens', str(g)] + rest
            SG.main()
            G = json.load(open(gj)) if gj.exists() else {}
            if not G.get('gens', {}).get(key, {}).get('done'):
                say(f's13_gen2: generation {g} did not finish ({G.get("stopped")})'); break
            S['gens'][key]['wall_s'] = round(time.time() - S['gens'][key]['start'])
        rec = G['gens'][key]
        ok, txt = g05(rec)
        S['gens'].setdefault(key, {})['g05'] = dict(passed=ok, text=txt); save2()
        say(f's13_gen2: generation {g}: G0.5 {"PASS" if ok else "not yet"}: {txt}')
        if ok:
            S['stopped'] = f'G0.5 PASSED at generation {g}: {txt}'; save2(); say(S['stopped']); break
        st = G.get('stopped', '')
        if st and not st.endswith(BENIGN):
            say(f's13_gen2: s13_gen stopped: {st}'); break
    say('s13_gen2: exit')


if __name__ == '__main__':
    main()
