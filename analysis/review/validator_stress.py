"""Adversarial tests of ctoc14.validator: parsing (tabs, comments, blank lines, BOM), NaN/inf values, Event=3 rows inside
arcs, arcs split by Event=2, thrust check between samples, tolerance edges, time window on Event=2/4 rows, duplicate flybys,
Asteroid_ID range, structural checks (Event=0 in the middle, SC order, Line column)."""
import sys, pathlib, io, contextlib, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from ctoc14.validator import validate, parse, Row
from ctoc14.submission import SCTrajectory, write_submission
from ctoc14.kepler import Ephemeris
from ctoc14.constants import DAY, T_MISSION

SCR = pathlib.Path('/private/tmp/claude-502/-Users-mickey-solarsystem-ctoc14/ee06b008-a7f5-4987-a150-04e70812e22d/scratchpad/valstress')
SCR.mkdir(parents=True, exist_ok=True)
eph = Ephemeris()
fails = []
def check(name, cond, detail=''):
    (print(f'ok   {name}: {detail}') if cond else (fails.append(name), print(f'FAIL {name}: {detail}')))

def run(path, expect_ok, name, expect_msg=None):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        try:
            rep = validate(path, eph=eph)
        except Exception as ex:
            rep = None; err = f'EXCEPTION {type(ex).__name__}: {ex}'
    if rep is None:
        check(name, expect_ok is None, err); return None
    msg = '; '.join(rep.errors)[:300]
    ok = (rep.ok == expect_ok) and (expect_msg is None or any(expect_msg in e for e in rep.errors))
    check(name, ok, f'ok={rep.ok} covered={rep.n_covered} J={rep.J:.3f} errors=[{msg}]')
    return rep

def rows_to_lines(rows):
    return [f'{i+1} {sc} {ev} {t:.17g} {r[0]:.17g} {r[1]:.17g} {r[2]:.17g} {v[0]:.17g} {v[1]:.17g} {v[2]:.17g} {m:.17g} {T[0]:.17g} {T[1]:.17g} {T[2]:.17g} {a}'
            for i, (sc, ev, t, r, v, m, T, a) in enumerate(rows)]

def write(path, lines):
    with open(path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')

# ---------------------------------------------------------------- base trajectory: PDF example (ballistic flyby of 174 + arc)
r1 = np.array([-1.9500328647e+07, 1.4581848347e+08, -7.6372048832e+03]); v1 = np.array([-2.9389477847e+01, -4.3166912970e+00, -3.9417153392e+00])
def base_traj(m0=2000.0, arc=True, flyby_in_arc=False, n_samples=41, mag=0.5):
    sc = SCTrajectory(0.0, r1, v1, m0)
    sc.coast_to(1.0008479152e+07); sc.flyby(174)
    if arc:
        t3 = 1.0008489152e+07
        ts = t3 + np.arange(n_samples) * DAY
        mg = mag * np.sin(np.pi * np.arange(n_samples) / (n_samples - 1.0))
        az = np.radians(40.107 + 9.0 * np.arange(n_samples)); el = np.radians(17.189)
        Ts = mg[:, None] * np.stack([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el) * np.ones(n_samples)], axis=1)
        sc.thrust_arc(ts, Ts, flybys=[(t3 + 5.5 * DAY, 174)] if flyby_in_arc else ())
    sc.coast(30 * DAY, node=True); sc.coast(100 * DAY); sc.end()
    return sc

base = base_traj()
p0 = SCR / 'base.txt'; write_submission(p0, [base], header='base'); rep0 = run(p0, True, 'base PDF-like file passes')
lines0 = [l for l in open(p0).read().splitlines() if l and not l.startswith('#')]

# ---------------------------------------------------------------- 1. parsing: tabs, comments, blank lines, mixed separators, CRLF, BOM
p = SCR / 'tabs.txt'
write(p, ['# comment', '', '   ', '\t'] + [l.replace(' ', '\t') for l in lines0[:3]] + ['# mid comment'] + [l.replace(' ', ' \t ') for l in lines0[3:]] + ['', '#end'])
run(p, True, 'tabs / comments / blank lines parse')
p = SCR / 'crlf.txt'
open(p, 'w', newline='').write('\r\n'.join(lines0) + '\r\n'); run(p, True, 'CRLF line endings parse')
p = SCR / 'bom.txt'
open(p, 'w', encoding='utf-8-sig').write('\n'.join(lines0) + '\n'); run(p, None, 'UTF-8 BOM (expect exception or pass; informational)')
p = SCR / 'inline_comment.txt'
write(p, [lines0[0] + '  # trailing comment'] + lines0[1:]); run(p, None, 'inline trailing comment -> column-count exception (informational)')
p = SCR / 'fortran_exp.txt'
write(p, [lines0[0].replace('2000', '2.0D3')] + lines0[1:]); run(p, None, 'Fortran D exponent -> exception (informational)')

# ---------------------------------------------------------------- 2. NaN / inf values
for field, name in [(4, 'x'), (7, 'vx'), (10, 'm'), (3, 't')]:
    p = SCR / f'nan_{name}.txt'
    ls = lines0.copy(); parts = ls[5].split(); parts[field] = 'nan'; ls[5] = ' '.join(parts); write(p, ls)
    run(p, False, f'NaN in {name} column of an Event=1 row must FAIL')
p = SCR / 'inf_x.txt'; ls = lines0.copy(); parts = ls[5].split(); parts[4] = 'inf'; ls[5] = ' '.join(parts); write(p, ls)
run(p, False, 'inf in x column must FAIL')
p = SCR / 'nan_T.txt'; ls = lines0.copy(); parts = ls[5].split(); parts[11] = 'nan'; ls[5] = ' '.join(parts); write(p, ls)
run(p, False, 'NaN thrust component must FAIL')

# ---------------------------------------------------------------- 3. Event=3 rows inside arcs and arcs separated by Event=2
sc = base_traj(flyby_in_arc=True)
p = SCR / 'flyby_in_arc.txt'; write_submission(p, [sc]); rep = run(p, True, 'Event=3 row inside arc (writer-generated) passes')
check('flyby-in-arc row order', [r[0] for r in sc.rows[:12]] == [0, 3, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1], str([r[0] for r in sc.rows[:12]]))
# the Event=3 row inside the arc with NON-zero thrust columns must fail
ls = [l for l in open(p).read().splitlines() if l and not l.startswith('#')]
k = [i for i, l in enumerate(ls) if l.split()[2] == '3'][1]
parts = ls[k].split(); parts[11] = '0.1'; ls2 = ls.copy(); ls2[k] = ' '.join(parts); p2 = SCR / 'flyby_in_arc_T.txt'; write(p2, ls2)
run(p2, False, 'Event=3 inside arc with nonzero thrust must FAIL', 'non-thrust row must have zero thrust')
# same file but the Event=3 row's mass made stale (as if written with the pre-arc mass)
parts = ls[k].split(); parts[10] = '2000'; ls2 = ls.copy(); ls2[k] = ' '.join(parts); p2 = SCR / 'flyby_in_arc_m.txt'; write(p2, ls2)
run(p2, False, 'Event=3 inside arc with stale mass must FAIL', 'mass error')
# an Event=2 inserted between two samples must split the arc -> dynamics mismatch (thrust zero on the split) -> FAIL
k1 = [i for i, l in enumerate(ls) if l.split()[2] == '1'][10]
parts = ls[k1].split(); parts[2] = '2'; parts[11] = parts[12] = parts[13] = '0'; ls2 = ls.copy(); ls2[k1] = ' '.join(parts); p2 = SCR / 'event2_in_arc.txt'; write(p2, ls2)
run(p2, False, 'Event=2 between samples splits the arc -> must FAIL dynamics')

# ---------------------------------------------------------------- 4. thrust checks: |T| between samples, spacing, sample magnitude
# plateau [0, .5, .5, 0] overshoots to 0.5625 between samples: samples are legal, interpolant is not
sc = SCTrajectory(0.0, r1, v1, 2000.0); sc.coast(10 * DAY)
ts = sc.t + np.arange(4) * DAY; Ts = np.array([[0, 0, 0], [.5, 0, 0], [.5, 0, 0], [0, 0, 0]], float)
sc.thrust_arc(ts, Ts); sc.coast(10 * DAY); sc.end()
p = SCR / 'overshoot.txt'; write_submission(p, [sc]); run(p, False, 'interpolation overshoot 0.5625 N must FAIL', 'interpolated |T|')
# constant 0.5 N arc exactly at the limit: legal
sc = SCTrajectory(0.0, r1, v1, 2000.0); sc.coast(10 * DAY)
ts = sc.t + np.arange(6) * DAY; Ts = np.tile([0.3, 0.4, 0.0], (6, 1))
sc.thrust_arc(ts, Ts); sc.coast(10 * DAY); sc.end()
p = SCR / 'limit.txt'; write_submission(p, [sc]); run(p, True, 'constant |T| = 0.5 N exactly passes')
# rounded components: (0.5/sqrt3)*3 written with 10 digits may exceed 0.5 by 1e-11
sc = SCTrajectory(0.0, r1, v1, 2000.0); sc.coast(10 * DAY)
ts = sc.t + np.arange(6) * DAY; Ts = np.tile([0.5 / np.sqrt(3)] * 3, (6, 1))
sc.thrust_arc(ts, Ts); sc.coast(10 * DAY); sc.end()
p = SCR / 'limit_sqrt3.txt'; write_submission(p, [sc]); rep = run(p, None, 'constant 0.5 N with irrational components (informational)')
ls = [l for l in open(p).read().splitlines() if l and not l.startswith('#')]
print('   max |T| reported:', rep.max_thrust if rep else None, ' 0.5 - max|T| =', 0.5 - rep.max_thrust if rep else None)
# sample spacing 8639 s must fail; 8640 s must pass
for dt, exp in [(8639.0, False), (8640.0, True)]:
    sc = SCTrajectory(0.0, r1, v1, 2000.0); sc.coast(10 * DAY)
    ts = sc.t + np.arange(5) * dt; Ts = np.tile([0.1, 0.1, 0.1], (5, 1))
    sc.thrust_arc(ts, Ts); sc.coast(10 * DAY); sc.end()
    p = SCR / f'spacing_{int(dt)}.txt'; write_submission(p, [sc]); run(p, exp, f'sample spacing {dt} s -> {"pass" if exp else "fail"}')
# spacing check must ignore an Event=3 row between samples (PDF 6.1-4)
sc = SCTrajectory(0.0, r1, v1, 2000.0); sc.coast(10 * DAY)
ts = sc.t + np.arange(5) * 8640.0; Ts = np.tile([0.1, 0.1, 0.1], (5, 1))
sc.thrust_arc(ts, Ts, flybys=[(ts[1] + 4000.0, 174)]); sc.coast(10 * DAY); sc.end()
p = SCR / 'spacing_flyby.txt'; write_submission(p, [sc]); rep = run(p, False, 'flyby row 4000 s after a sample: spacing must NOT be flagged (only the 1000 km miss)')
check('spacing not flagged with Event=3 in between', rep is not None and not any('spaced' in e for e in rep.errors), '; '.join(rep.errors)[:200] if rep else '')
# sample |T| = 0.5 + 1e-11 (rounding) -> currently fails? tolerance 1e-12
ls = [l for l in open(SCR / 'limit.txt').read().splitlines() if l and not l.startswith('#')]
k = [i for i, l in enumerate(ls) if l.split()[2] == '1'][2]; parts = ls[k].split(); parts[11] = '0.30000000001'; ls2 = ls.copy(); ls2[k] = ' '.join(parts)
p = SCR / 'limit_plus.txt'; write(p, ls2); run(p, False, 'sample |T| = 0.5 + 6e-12 fails (hard bound)')

# ---------------------------------------------------------------- 5. tolerance edges on dynamics consistency
def perturb(lines, idx, col, delta, name, expect):
    ls = lines.copy(); parts = ls[idx].split(); parts[col] = f'{float(parts[col]) + delta:.17g}'; ls[idx] = ' '.join(parts)
    p = SCR / f'tol_{name}.txt'; write(p, ls); return run(p, expect, name)
perturb(lines0, 10, 4, 0.99, 'pos +0.99 km passes', True)
perturb(lines0, 10, 4, 1.01, 'pos +1.01 km fails', False)
perturb(lines0, 10, 7, 0.00099, 'vel +0.99 m/s passes (start-of-segment error only)', True)
perturb(lines0, 10, 7, 0.00101, 'vel +1.01 m/s fails', False)
perturb(lines0, 10, 10, 0.0099, 'mass +0.0099 kg passes', True)
perturb(lines0, 10, 10, 0.0101, 'mass +0.0101 kg fails', False)
# mass below 600 on a row
perturb(lines0, 0, 10, -1400.5, 'm0 = 599.5 fails', False)
perturb(lines0, 0, 10, 0.5, 'm0 = 2000.5 fails', False)

# ---------------------------------------------------------------- 6. time window on Event=2 / Event=4 rows, launch time
sc = SCTrajectory(0.0, r1, v1, 2000.0); sc.coast(T_MISSION - 1.0, node=True); sc.coast(2.0); sc.end()
p = SCR / 'end_late.txt'; write_submission(p, [sc]); run(p, False, 'Event=4 at T+1 s fails', 'outside')
sc = SCTrajectory(0.0, r1, v1, 2000.0); sc.coast(T_MISSION, node=True); sc.end()   # Event=2 and Event=4 same time -> also strictness
p = SCR / 'end_exact.txt'; write_submission(p, [sc]); run(p, False, 'Event=2 and Event=4 at the same time fails (strictly increasing)')
sc = SCTrajectory(0.0, r1, v1, 2000.0); sc.coast(T_MISSION - 10.0, node=True); sc.coast(10.0); sc.end()
p = SCR / 'end_at_T.txt'; write_submission(p, [sc]); run(p, True, 'Event=4 exactly at T passes')
sc = SCTrajectory(-1.0, r1, v1, 2000.0); sc.coast(10 * DAY); sc.end()
p = SCR / 'launch_neg.txt'; write_submission(p, [sc]); run(p, False, 'launch at t=-1 fails')
# flyby after T (with end after T): both flagged
sc = SCTrajectory(0.0, r1, v1, 2000.0); sc.coast(T_MISSION + 5.0); sc.flyby(1); sc.coast(1.0); sc.end()
p = SCR / 'flyby_late.txt'; write_submission(p, [sc]); run(p, False, 'flyby after T fails')

# ---------------------------------------------------------------- 7. duplicate flybys, Asteroid_ID range, bogus flyby
sc = SCTrajectory(0.0, r1, v1, 2000.0); sc.coast_to(1.0008479152e+07); sc.flyby(174); sc.coast(1.0); sc.flyby(174); sc.coast(10 * DAY); sc.end()
p = SCR / 'dup_flyby.txt'; write_submission(p, [sc]); rep = run(p, True, 'duplicate flyby (both within 1000 km) passes, counted once')
check('duplicate counted once', rep is not None and rep.n_covered == 1, f'covered={rep.n_covered if rep else None}')
sc = SCTrajectory(0.0, r1, v1, 2000.0); sc.coast_to(1.0008479152e+07); sc.flyby(174); sc.coast(20 * DAY); sc.flyby(174); sc.coast(10 * DAY); sc.end()
p = SCR / 'dup_flyby_far.txt'; write_submission(p, [sc]); run(p, False, 'second (far) declaration of an already-covered asteroid fails', 'flyby distance')
for aid, exp in [(0, False), (301, False), (-5, False), (300, False)]:
    sc = SCTrajectory(0.0, r1, v1, 2000.0); sc.coast_to(1.0008479152e+07); sc.flyby(aid); sc.coast(10 * DAY); sc.end()
    p = SCR / f'aid_{aid}.txt'; write_submission(p, [sc]); run(p, exp, f'Asteroid_ID {aid} on Event=3 fails')
# Asteroid_ID nonzero on a non-flyby row (PDF: must be 0)
ls = lines0.copy(); parts = ls[3].split(); parts[14] = '5'; ls[3] = ' '.join(parts); p = SCR / 'aid_on_event1.txt'; write(p, ls)
run(p, None, 'Asteroid_ID=5 on an Event=1 row (PDF says 0) - informational')
# two spacecraft declare the same asteroid: counted once
sc2 = base_traj(arc=False)
p = SCR / 'two_sc_same.txt'; write_submission(p, [base_traj(arc=False), sc2]); rep = run(p, True, 'two SC same asteroid passes')
check('two SC same asteroid counted once, J = 2 + 299', rep is not None and rep.n_covered == 1 and abs(rep.J - 301) < 1e-9, f'covered={rep.n_covered if rep else None} J={rep.J if rep else None}')

# ---------------------------------------------------------------- 8. structural: Event=0 in the middle, Event=4 in the middle, SC order, SC gap, Line column
ls = lines0.copy(); parts = ls[44 if len(ls) > 44 else -2].split()
k = [i for i, l in enumerate(ls) if l.split()[2] == '2'][0]; parts = ls[k].split(); parts[2] = '0'; ls[k] = ' '.join(parts); p = SCR / 'event0_mid.txt'; write(p, ls)
run(p, False, 'second Event=0 in the middle of a spacecraft must FAIL (informational: expected by PDF structure)')
ls = lines0.copy(); parts = ls[k].split(); parts[2] = '4'; ls[k] = ' '.join(parts); p = SCR / 'event4_mid.txt'; write(p, ls)
run(p, False, 'Event=4 in the middle must FAIL (informational)')
ls = lines0.copy(); ls = [l.replace(' 1 ', ' 2 ', 1) for l in ls]; p = SCR / 'sc_starts_at_2.txt'; write(p, ls); run(p, False, 'SC_ID starting at 2 fails')
# SC blocks interleaved / out of order (PDF: sorted by SC_ID)
two = SCR / 'two_sc.txt'; write_submission(two, [base_traj(arc=False), base_traj(arc=False, m0=1000.0)])
lt = [l for l in open(two).read().splitlines() if l and not l.startswith('#')]
n1 = sum(1 for l in lt if l.split()[1] == '1'); p = SCR / 'sc_swapped.txt'; write(p, lt[n1:] + lt[:n1]); run(p, None, 'SC2 block before SC1 block (PDF requires SC_ID ascending) - informational')
p = SCR / 'lines_zero.txt'; write(p, ['0' + l[1:] if i == 0 else l for i, l in enumerate(lines0)]); run(p, None, 'Line column = 0 on first row (PDF: positive) - informational')

# ---------------------------------------------------------------- 9. launch checks
sc = SCTrajectory(0.0, r1 + np.array([1.5, 0, 0]), v1, 2000.0); sc.coast(10 * DAY); sc.end()
p = SCR / 'launch_pos.txt'; write_submission(p, [sc]); run(p, False, 'launch position 1.5 km off fails', 'launch position')
rE, vE = eph.earth_state(0.0)
sc = SCTrajectory(0.0, rE, vE + np.array([4.0105, 0, 0]), 2000.0); sc.coast(10 * DAY); sc.end()
p = SCR / 'vinf_4011.txt'; write_submission(p, [sc]); run(p, False, 'v_inf 4.0105 fails')
sc = SCTrajectory(0.0, rE, vE + np.array([4.0095, 0, 0]), 2000.0); sc.coast(10 * DAY); sc.end()
p = SCR / 'vinf_4009.txt'; write_submission(p, [sc]); run(p, True, 'v_inf 4.0095 passes (tolerance)')
# zero-fuel spacecraft with a thrust arc -> mass < 600 must fail
sc = SCTrajectory(0.0, rE, vE, 600.0); sc.coast(10 * DAY)
ts = sc.t + np.arange(4) * DAY; sc.thrust_arc(ts, np.tile([0.1, 0, 0], (4, 1))); sc.coast(DAY); sc.end()
p = SCR / 'dry_thrust.txt'; write_submission(p, [sc]); run(p, False, 'm0 = 600 with a thrust arc fails (mass < 600)')
# mass 599.995 (within 0.01 kg of dry): passes here, may fail a strict official checker
sc = SCTrajectory(0.0, rE, vE, 600.0 + 0.004); sc.coast(10 * DAY)
ts = sc.t + np.arange(4) * DAY; sc.thrust_arc(ts, np.tile([0.01, 0, 0], (4, 1))); sc.coast(DAY); sc.end()
p = SCR / 'dry_margin.txt'; write_submission(p, [sc]); rep = run(p, None, f'final mass {sc.m:.4f} kg (< 600 by less than 0.01) - informational')

print('\nFAILURES:', fails if fails else 'none')
