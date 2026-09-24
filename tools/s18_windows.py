"""Stage 18 / WP1: the WINDOW and CLASS table  results/s18/windows.json.

A target's windows are its node events (minima of the distance between the asteroid and the 1 AU ecliptic circle, the
"torus distance") closer than D_WIN = 0.05 AU; its class is pin (<= 4 windows in 15 yr), fil (>= 10) or med.  The
reference builder is s18_common.build_windows (it reads results/s17/physics/events.json and writes the file on first
use); this module is the CLI around it plus an independent RECOMPUTE from MEA.txt (the events file is itself derived
data: results/s17/physics/events.py, 0.5 d sampling, minima < 0.15 AU) and a consistency report.

File format (s18_common.build_windows):
  {"meta": {source, d_win, d_fallback, pin_max, fill_min, built, n_fallback},
   "targets": {"<id>": {"t_d": [DAYS, ascending], "n": <windows at d_win>, "cls": "pin"|"med"|"fil", "fallback": bool}}}
Runtime table (s18_common.load_windows): {id (int): {"t": np.array SECONDS, "n", "cls", "fallback"}, "_meta": {...}}.

usage: s18_windows.py build  [--events results/s17/physics/events.json] [--out results/s18/windows.json] [--force]
       s18_windows.py recompute [--mea MEA.txt] [--step-d 0.5] [--dmax 0.15] [--out results/s18/events_recomputed.json]
       s18_windows.py report [--windows results/s18/windows.json]
       s18_windows.py check  [--windows ...] [--events ...]     (exit 1 on any mismatch with history_classes)
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, pathlib, argparse
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import numpy as np
import s18_common as C
from ctoc14.search import UNREACHABLE
from ctoc14.constants import DAY, AU, T_MISSION

# expected class sizes (results/s17/physics/history_classes.txt: catalogue {'pin': 95, 'med': 123, 'fill': 80})
EXPECTED = dict(pin=95, med=123, fil=80)
# the 19 targets without any node event < 0.05 AU (docs/stage18/BUILD.md 1.1); their windows are the events < 0.08 AU
FALLBACK_IDS = [1, 31, 39, 50, 66, 70, 74, 75, 96, 107, 118, 140, 149, 179, 186, 269, 270, 280, 295]
A_EARTH = 1.0009175           # AU, semi-major axis of the Earth ellipse (constants.EARTH_ELEMENTS[0])


def build(events=C.EVENTS_SRC, out=C.WINDOWS_FILE, force=False):
    """Build (or reuse, unless force) the window table.  Returns the runtime table (s18_common.load_windows format)."""
    out = pathlib.Path(out)
    if out.exists() and not force:
        return C.load_windows(out)
    return C.build_windows(events=events, out=out)


def recompute_events(mea=ROOT / 'MEA.txt', step_d=0.5, d_max=0.15):
    """Recompute the node events from the ephemeris (independent of results/s17): sample every asteroid on a step_d grid
    over [0, T_MISSION), torus distance d = hypot(rho - a_E, z) [AU] with rho the ecliptic cylindrical radius and
    a_E = 1.0009175, keep the strict local minima with d < d_max.  Returns [{"ast": id, "t": DAYS, "d": AU}]
    (results/s17/physics/events.py, minus the lag / speed diagnostics)."""
    import run_ialns as RI
    mea = pathlib.Path(mea)
    if mea.resolve() == (ROOT / 'MEA.txt').resolve():
        E = RI.eph()
    else:
        from ctoc14.kepler import Ephemeris
        E = Ephemeris(str(mea))
    tt = np.arange(0.0, T_MISSION, float(step_d) * DAY)
    ev = []
    for X in range(1, E.n_ast + 1):
        r, _ = E.ast_states_at(np.full(len(tt), X - 1), tt)
        rho = np.hypot(r[:, 0], r[:, 1]) / AU; z = r[:, 2] / AU
        d = np.hypot(rho - A_EARTH, z)
        k = np.where((d[1:-1] <= d[:-2]) & (d[1:-1] < d[2:]) & (d[1:-1] < d_max))[0] + 1
        for j in k:
            ev.append(dict(ast=int(X), t=float(tt[j] / DAY), d=float(d[j])))
    return ev


def _quartiles(v):
    if not len(v):
        return [0, 0, 0, 0, 0]
    return [float(x) for x in np.percentile(np.asarray(v, float), [0, 25, 50, 75, 100])]


def report(W):
    """Text report: class histogram, number of fallback targets, per-class quartiles of the window count, the targets
    whose LAST window is before 2 yr / after 13 yr (the "now or never" ends), and a table of the pins with <= 2 windows."""
    ids = sorted(k for k in W if isinstance(k, int))
    meta = W.get('_meta', {})
    L = [f'windows table: {len(ids)} targets; meta {json.dumps(meta, default=str)}']
    hist = {c: [t for t in ids if W[t]['cls'] == c] for c in C.CLASSES}
    L.append('classes: ' + ', '.join(f'{c} {len(hist[c])} (expected {EXPECTED[c]})' for c in C.CLASSES))
    fb = [t for t in ids if W[t]['fallback']]
    L.append(f'fallback targets (no event < {meta.get("d_win", C.D_WIN)} AU; windows at < {meta.get("d_fallback", C.D_WIN_FALLBACK)} AU): '
             f'{len(fb)}: {fb}')
    for c in C.CLASSES:
        q = _quartiles([W[t]['n'] for t in hist[c]])
        qw = _quartiles([len(W[t]['t']) for t in hist[c]])
        L.append(f'  {c}: windows(n at d_win) min/q1/med/q3/max = {q[0]:.0f}/{q[1]:.1f}/{q[2]:.1f}/{q[3]:.1f}/{q[4]:.0f}; '
                 f'usable windows {qw[0]:.0f}/{qw[1]:.1f}/{qw[2]:.1f}/{qw[3]:.1f}/{qw[4]:.0f}')
    last = {t: (float(W[t]['t'][-1]) / DAY if len(W[t]['t']) else None) for t in ids}
    early = sorted(t for t in ids if last[t] is not None and last[t] < 2 * 365.25)
    late = sorted(t for t in ids if last[t] is not None and last[t] > 13 * 365.25)
    none = sorted(t for t in ids if last[t] is None)
    L.append(f'last window before 2 yr ({len(early)}): ' + ', '.join(f'{t}({W[t]["cls"]},{last[t]:.0f}d)' for t in early))
    L.append(f'last window after 13 yr ({len(late)}): ' + ', '.join(f'{t}({W[t]["cls"]})' for t in late))
    if none:
        L.append(f'NO window at all ({len(none)}): {none}')
    pins2 = sorted((t for t in ids if W[t]['cls'] == 'pin' and len(W[t]['t']) <= 2), key=lambda t: (len(W[t]['t']), t))
    L.append(f'pins with <= 2 windows ({len(pins2)}):')
    L.append('   id  n  fb  windows [d]')
    for t in pins2:
        L.append(f'  {t:3d}  {len(W[t]["t"])}  {"y" if W[t]["fallback"] else "-"}   ' + ', '.join(f'{x / DAY:.0f}' for x in W[t]['t']))
    return '\n'.join(L)


def check(W, events=C.EVENTS_SRC):
    """Consistency: class sizes == EXPECTED; every target has >= 1 window (fallback rule); every window epoch is inside
    (0, T_MISSION); the 19 fallback ids are exactly those with no event < 0.05 AU in `events`; 131/144 absent.
    Returns (ok: bool, problems: [str])."""
    probs = []
    ids = sorted(k for k in W if isinstance(k, int))
    if set(ids) != set(C.TARGETS):
        probs.append(f'target set != TARGETS: missing {sorted(set(C.TARGETS) - set(ids))}, extra {sorted(set(ids) - set(C.TARGETS))}')
    for u in UNREACHABLE:
        if u in W:
            probs.append(f'unreachable {u} present')
    from collections import Counter
    cnt = Counter(W[t]['cls'] for t in ids)
    for c in C.CLASSES:
        if cnt.get(c, 0) != EXPECTED[c]:
            probs.append(f'class {c}: {cnt.get(c, 0)} != expected {EXPECTED[c]}')
    for t in ids:
        w = W[t]['t']
        if len(w) < 1:
            probs.append(f'target {t}: no window')
        elif not (np.all(w > 0) and np.all(w < T_MISSION)):
            probs.append(f'target {t}: window epoch outside (0, T_MISSION)')
        if len(w) > 1 and not np.all(np.diff(w) > 0):
            probs.append(f'target {t}: windows not ascending')
        n_at = int(W[t]['n']); cls = W[t]['cls']
        exp = 'pin' if n_at <= C.PIN_MAX else ('fil' if n_at >= C.FILL_MIN else 'med')
        if cls != exp:
            probs.append(f'target {t}: class {cls} but n = {n_at} -> {exp}')
        if W[t]['fallback'] != (n_at == 0):
            probs.append(f'target {t}: fallback flag {W[t]["fallback"]} but n = {n_at}')
    fb = sorted(t for t in ids if W[t]['fallback'])
    if fb != FALLBACK_IDS:
        probs.append(f'fallback ids {fb} != BUILD.md list {FALLBACK_IDS}')
    events = pathlib.Path(events)
    if events.exists():
        ev = C.jload(events); n05 = {}
        for e in ev:
            if float(e['d']) < C.D_WIN:
                n05[int(e['ast'])] = n05.get(int(e['ast']), 0) + 1
        fb_ev = sorted(t for t in C.TARGETS if n05.get(t, 0) == 0)
        if fb_ev != fb:
            probs.append(f'fallback ids from events {fb_ev} != table {fb}')
        for t in ids:
            if n05.get(t, 0) != int(W[t]['n']):
                probs.append(f'target {t}: n {W[t]["n"]} != events count {n05.get(t, 0)}')
    else:
        probs.append(f'events file {events} missing (fallback cross-check skipped)')
    return (len(probs) == 0), probs


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('cmd', choices=['build', 'recompute', 'report', 'check'])
    ap.add_argument('--events', default=str(C.EVENTS_SRC)); ap.add_argument('--out', default=None)
    ap.add_argument('--windows', default=str(C.WINDOWS_FILE)); ap.add_argument('--force', action='store_true')
    ap.add_argument('--mea', default=str(ROOT / 'MEA.txt')); ap.add_argument('--step-d', type=float, default=0.5)
    ap.add_argument('--dmax', type=float, default=0.15)
    a = ap.parse_args()
    if a.cmd == 'build':
        W = build(a.events, a.out or str(C.WINDOWS_FILE), a.force); print(report(W))
    elif a.cmd == 'recompute':
        ev = recompute_events(a.mea, a.step_d, a.dmax)
        C.jdump(ev, a.out or str(C.S18 / 'events_recomputed.json')); print(len(ev), 'events')
    elif a.cmd == 'report':
        print(report(C.load_windows(a.windows)))
    else:
        ok, probs = check(C.load_windows(a.windows), a.events)
        print('OK' if ok else 'FAIL'); [print(' ', p) for p in probs]
        sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
