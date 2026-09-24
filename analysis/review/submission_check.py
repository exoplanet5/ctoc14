"""Submission writer edge cases: flyby coinciding with a sample time, flyby at the arc ends, coast with negative dt,
%.17g round trip, row ordering, thrust_arc starting before the current time, multi-spacecraft numbering."""
import sys, pathlib, io, contextlib, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from ctoc14.submission import SCTrajectory, write_submission
from ctoc14.validator import validate, parse
from ctoc14.kepler import Ephemeris
from ctoc14.constants import DAY

SCR = pathlib.Path('/private/tmp/claude-502/-Users-mickey-solarsystem-ctoc14/ee06b008-a7f5-4987-a150-04e70812e22d/scratchpad/subcheck'); SCR.mkdir(parents=True, exist_ok=True)
eph = Ephemeris(); rE, vE = eph.earth_state(0.0)
fails = []
def check(name, cond, detail=''):
    (print(f'ok   {name}: {detail}') if cond else (fails.append(name), print(f'FAIL {name}: {detail}')))
def val(path):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        return validate(path, eph=eph)

# 1. flyby exactly at a sample time -> duplicate Time rows
sc = SCTrajectory(0.0, rE, vE + np.array([3.0, 0, 0]), 2000.0); sc.coast(10 * DAY)
ts = sc.t + np.arange(6) * DAY; Ts = np.tile([0.1, 0.1, 0.0], (6, 1))
sc.thrust_arc(ts, Ts, flybys=[(ts[2], 7)]); sc.coast(DAY); sc.end()
times = [r[1] for r in sc.rows]; evs = [r[0] for r in sc.rows]
p = SCR / 'flyby_at_sample.txt'; write_submission(p, [sc]); rep = val(p)
check('flyby at a sample time yields strictly increasing times', np.all(np.diff(times) > 0), f'events={evs} dup={[t for t in times if times.count(t) > 1]}; validator errors={rep.errors[:2]}')
# 2. flyby at the arc start / end sample
for which in ['start', 'end']:
    sc = SCTrajectory(0.0, rE, vE + np.array([3.0, 0, 0]), 2000.0); sc.coast(10 * DAY)
    ts = sc.t + np.arange(6) * DAY
    sc.thrust_arc(ts, Ts, flybys=[(ts[0] if which == 'start' else ts[-1], 7)]); sc.coast(DAY); sc.end()
    times = [r[1] for r in sc.rows]
    check(f'flyby at arc {which} sample yields strictly increasing times', np.all(np.diff(times) > 0), f'events={[r[0] for r in sc.rows]}')
# 3. flyby outside the arc passed through thrust_arc (after last sample)
sc = SCTrajectory(0.0, rE, vE + np.array([3.0, 0, 0]), 2000.0); sc.coast(10 * DAY)
ts = sc.t + np.arange(6) * DAY
sc.thrust_arc(ts, Ts, flybys=[(ts[-1] + 3600.0, 7)]); sc.coast(DAY); sc.end()
p = SCR / 'flyby_after_arc.txt'; write_submission(p, [sc]); rep = val(p)
check('flyby 1 h after the last sample (via thrust_arc): rows consistent', all('error' not in e for e in rep.errors), f'errors={rep.errors[:3]}')
# 4. coast with negative dt / coast_to into the past silently ignored
sc = SCTrajectory(0.0, rE, vE, 2000.0); sc.coast(10 * DAY, node=True); sc.coast_to(5 * DAY, node=True); sc.end()
times = [r[1] for r in sc.rows]
check('coast_to into the past is rejected or produces increasing times', np.all(np.diff(times) > 0), f'times={times}')
# 5. thrust_arc whose first sample is before the current time
sc = SCTrajectory(0.0, rE, vE, 2000.0); sc.coast(10 * DAY)
try:
    sc.thrust_arc(sc.t - DAY + np.arange(4) * DAY, np.tile([0.1, 0, 0], (4, 1))); res = 'no error'
except AssertionError as ex:
    res = f'AssertionError: {ex}'
check('thrust_arc starting in the past raises', res.startswith('AssertionError'), res)
# 6. %.17g round trip of every field
sc = SCTrajectory(1234.56789012345678, rE, vE + np.array([1 / 3, 0, 0]), 2000.0 / 3); sc.coast(DAY / 7, node=True); sc.end()
p = SCR / 'roundtrip.txt'; write_submission(p, [sc]); rows = parse(p)
exact = all(rows[i].t == sc.rows[i][1] and np.array_equal(rows[i].r, sc.rows[i][2]) and np.array_equal(rows[i].v, sc.rows[i][3]) and rows[i].m == sc.rows[i][4] for i in range(len(rows)))
check('%.17g round-trips every double exactly', exact, open(p).read().splitlines()[0][:120])
# 7. line numbering across spacecraft and SC ids
sc1 = SCTrajectory(0.0, rE, vE, 700.0); sc1.coast(DAY); sc1.end()
sc2 = SCTrajectory(DAY, *eph.earth_state(DAY), 800.0); sc2.coast(DAY); sc2.end()
p = SCR / 'two.txt'; write_submission(p, [sc1, sc2], header='h1\nh2'); rows = parse(p)
check('global consecutive Line numbers and SC ids', [r.line for r in rows] == [1, 2, 3, 4] and [r.sc for r in rows] == [1, 1, 2, 2], str([(r.line, r.sc) for r in rows]))
txt = open(p).read(); check('multi-line header is commented', txt.startswith('# h1\n# h2\n'), txt[:20].replace('\n', '|'))
# 8. thrust rows carry the sample vector; non-thrust rows exact zeros; Asteroid_ID 0 elsewhere
sc = SCTrajectory(0.0, rE, vE + np.array([3.0, 0, 0]), 2000.0); sc.coast(10 * DAY)
ts = sc.t + np.arange(4) * DAY; Ts = np.array([[0.1, -0.2, 0.3], [0.0, 0.0, 0.0], [0.4, 0.0, 0.1], [0.0, 0.2, 0.0]])
sc.thrust_arc(ts, Ts, flybys=[(ts[1] + 100.0, 9)]); sc.coast(DAY, node=True); sc.end()
p = SCR / 'fields.txt'; write_submission(p, [sc]); rows = parse(p)
ok = all((np.array_equal(r.T, Ts[[i for i, x in enumerate(ts) if x == r.t][0]]) if r.event == 1 else np.all(r.T == 0)) for r in rows) and all((r.ast == 9) == (r.event == 3) for r in rows)
check('row fields (thrust / asteroid id) correct', ok, str([(r.event, r.ast, list(r.T)) for r in rows]))
# 9. writer output re-validated after parse/rewrite (combine_submission path)
rep = val(p); check('fields file validates', rep.ok, '; '.join(rep.errors)[:200])
print('\nFAILURES:', fails if fails else 'none')
