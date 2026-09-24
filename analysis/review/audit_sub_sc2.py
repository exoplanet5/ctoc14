#!/usr/bin/env python
"""Independent audit of a CTOC14 submission file (default: results/fleet_b300/sub_sc2.txt).

Everything here is implemented from the PDF rules (Sections 3-7) WITHOUT importing the ctoc14 package:
own Kepler ephemeris (Earth Table 3 / MEA.txt, MJD-epoch mean-motion rule), own sliding-window cubic Lagrange
interpolant (Reading A, 0-based centred window; Reading B forward window as a sensitivity test), exact per-interval
maximisation of |T(t)| (degree-6 polynomial, roots of its derivative), own vectorised fixed-step RK4 dynamics check
at several step sizes (to model an organisers' checker with a different integrator), plus text-level format checks
(15 columns, integer columns, %.17g round-trip, NaN/inf, CR/LF, tabs, BOM, non-ASCII).
Sensitivity sections quantify what checker variants that deviate from the PDF-confirmed model would see: forward
Lagrange window (Reading B, dynamics + exact peak thrust), Event=3 rows breaking the arc / used as zero samples,
continuous non-restarting integration, and the 8640 s rule misapplied to Event=3 rows. A fleet-context section looks
(read-only) at sibling sub_sc*.txt files for duplicate flybys. Exit status 1 iff a FAIL finding was produced.

Usage:  ~/.venvs/astro313/bin/python analysis/review/audit_sub_sc2.py [file]
"""
import sys, re, math
import numpy as np
from numpy.polynomial import polynomial as P
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SUB = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / 'results/fleet_b300/sub_sc2.txt'
MEA = ROOT / 'MEA.txt'

# ---------------------------------------------------------------- constants (PDF Table 2 / 3)
MU = 1.32712440018e11; AU = 149597870.7; G0 = 9.80665; ISP = 4000.0; TMAX = 0.5
MDRY = 600.0; M0MAX = 2000.0; VINFMAX = 4.0; DFLY = 1000.0; DTMIN = 8640.0; TMISSION = 473364000.0
T0_MJD = 62502.0; EPH_AST_MJD = 61200.0; EPH_EARTH_MJD = 60676.0; DAY = 86400.0
EARTH = (1.0009175020, 0.017566762041, 0.002976847126, 189.953211282428, 273.196254000254, 357.4135031077)
VEX = ISP * G0  # m/s

FINDINGS = []   # (severity, text)
def finding(sev, txt):
    FINDINGS.append((sev, txt)); print(f'  [{sev}] {txt}')
def ok(txt):
    print(f'  [ok] {txt}')

# ---------------------------------------------------------------- independent Kepler ephemeris
def kepler_E(M, e):
    M = math.fmod(M, 2 * math.pi)
    if M < 0: M += 2 * math.pi
    E = M if e < 0.8 else math.pi
    for _ in range(200):
        f = E - e * math.sin(E) - M
        d = -f / (1.0 - e * math.cos(E))
        E += d
        if abs(d) < 1e-15: break
    return E

def elem_state(elem, epoch_mjd, t):
    """Heliocentric ecliptic state (km, km/s) of a Keplerian body at mission time t [s since MJD 62502.0]."""
    a_au, e, i, Om, w, M0 = elem
    a = a_au * AU
    i, Om, w, M0 = (math.radians(x) for x in (i, Om, w, M0))
    n = math.sqrt(MU / a ** 3)                                  # eq (5), rad/s
    M = M0 + n * ((T0_MJD - epoch_mjd) * DAY + t)               # eq (4), seconds
    E = kepler_E(M, e)
    r = a * (1 - e * math.cos(E))
    xp = a * (math.cos(E) - e); yp = a * math.sqrt(1 - e * e) * math.sin(E)
    f = math.sqrt(MU * a) / r
    vxp = -f * math.sin(E); vyp = f * math.sqrt(1 - e * e) * math.cos(E)
    cO, sO, ci, si, cw, sw = math.cos(Om), math.sin(Om), math.cos(i), math.sin(i), math.cos(w), math.sin(w)
    R = np.array([[cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si],
                  [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si],
                  [sw * si, cw * si, ci]])
    return R @ np.array([xp, yp, 0.0]), R @ np.array([vxp, vyp, 0.0])

def load_mea(path):
    el = {}
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith('#'): continue
        p = s.split()
        el[int(p[0])] = tuple(float(x) for x in p[1:7])
    return el

# ---------------------------------------------------------------- 1. text-level parse
print(f'=== AUDIT of {SUB.relative_to(ROOT)} ===')
raw = SUB.read_bytes()
print('--- file / text format')
if raw.startswith(b'\xef\xbb\xbf'): finding('WARN', 'file starts with a UTF-8 BOM')
else: ok('no BOM')
if b'\r' in raw: finding('WARN', f'{raw.count(b"\r")} carriage returns (CRLF line endings)')
else: ok('LF line endings only')
if b'\t' in raw: finding('INFO', f'{raw.count(b"\t")} tab characters (allowed by rule 6.1-8)')
else: ok('no tabs (single-space separated)')
try:
    text = raw.decode('utf-8')
except UnicodeDecodeError as ex:
    finding('FAIL', f'not valid UTF-8: {ex}'); sys.exit(1)
if any(ord(c) > 127 for c in text): finding('WARN', 'non-ASCII characters present')
else: ok('pure ASCII')
if not raw.endswith(b'\n'): finding('WARN', 'file does not end with a newline')
else: ok('file ends with newline')

lines = text.split('\n')
if lines and lines[-1] == '': lines = lines[:-1]
n_comment = sum(1 for l in lines if l.strip().startswith('#')); n_blank = sum(1 for l in lines if not l.strip())
rows = []          # dicts
int_re = re.compile(r'^[+-]?\d+$')
fmt_bad = 0; roundtrip_bad = 0; nonfinite = 0; ncol_bad = 0; trailing_ws = 0; double_space = 0
for ln, line in enumerate(lines, 1):
    s = line.strip()
    if not s or s.startswith('#'): continue
    if line != s: trailing_ws += 1
    if '  ' in s: double_space += 1
    p = s.replace('\t', ' ').split()
    if len(p) != 15:
        ncol_bad += 1; finding('FAIL', f'text line {ln}: {len(p)} columns (need 15)'); continue
    for k in (0, 1, 2, 14):
        if not int_re.match(p[k]): fmt_bad += 1; finding('FAIL', f'text line {ln}: column {k+1} not an integer literal: {p[k]!r}')
    vals = []
    for k in range(3, 14):
        try: v = float(p[k])
        except ValueError: v = float('nan'); fmt_bad += 1; finding('FAIL', f'text line {ln}: column {k+1} not a float: {p[k]!r}')
        if not math.isfinite(v): nonfinite += 1; finding('FAIL', f'text line {ln}: column {k+1} non-finite {p[k]!r}')
        if ('%.17g' % v) != p[k]: roundtrip_bad += 1
        vals.append(v)
    rows.append(dict(ln=ln, line=int(p[0]), sc=int(p[1]), ev=int(p[2]), t=vals[0], r=np.array(vals[1:4]), v=np.array(vals[4:7]),
                     m=vals[7], T=np.array(vals[8:11]), Ttok=p[11:14], ast=int(p[14]), tok=p))
N = len(rows)
print(f'  data rows {N}, comment lines {n_comment}, blank lines {n_blank}')
if ncol_bad == 0: ok('every data row has exactly 15 whitespace-separated columns')
if fmt_bad == 0: ok('Line/SC_ID/Event/Asteroid_ID are integer literals; all 11 float columns parse')
if nonfinite == 0: ok('no NaN/inf tokens')
if roundtrip_bad == 0: ok("every float token is the canonical '%.17g' rendering of its double (round-trips exactly)")
else: finding('INFO', f'{roundtrip_bad} float tokens are not canonical %.17g renderings (they still parse; only cosmetic)')
negzero = sum(1 for R in rows for tok in R['tok'][3:14] if tok.startswith('-0') and float(tok) == 0.0)
intlit = sum(1 for R in rows for tok in R['tok'][3:14] if int_re.match(tok))
intlit_state = sum(1 for R in rows for tok in R['tok'][4:11] if int_re.match(tok))
if negzero: finding('INFO', f'{negzero} float tokens are negative zero ("-0"); numerically 0 but a literal string compare to "0" would fail')
else: ok('no negative-zero ("-0") tokens')
finding('INFO', f'{intlit} float-column tokens are bare integer literals (e.g. Time "3456000", mass "2000", thrust "0"; state columns: {intlit_state}); the PDF example prints every float as %.10e. strtod/float()/sscanf %lf/textscan all accept bare integers, so this is a robustness note only (a %.16e format would remove even that doubt).')
if trailing_ws: finding('INFO', f'{trailing_ws} rows with leading/trailing whitespace')
if double_space: finding('INFO', f'{double_space} rows contain double spaces')
# digit-count statistics (PDF recommends >= 10 significant digits)
sig = []
for R in rows:
    for tok in R['tok'][3:11]:
        m = re.match(r'^[+-]?0*(\d[\d.]*?)0*(?:[eE].*)?$', tok)
        digits = re.sub(r'[^\d]', '', tok.split('e')[0].split('E')[0]).lstrip('0')
        sig.append(len(digits))
print(f'  significant digits in Time/state columns: min {min(sig)}, max {max(sig)} (PDF recommends >= 10; integers like 3456000 / 2000 are exact)')
few = [(R['ln'], R['tok'][3]) for R in rows if len(re.sub(r"[^\d]", "", R['tok'][3]).lstrip('0')) < 10]
print(f'  Time tokens with < 10 significant digits: {len(few)} (all exact integers => no loss)')

# ---------------------------------------------------------------- 2. structural checks
print('--- structure')
lines_no = np.array([R['line'] for R in rows]); scs = np.array([R['sc'] for R in rows]); evs = np.array([R['ev'] for R in rows])
ts = np.array([R['t'] for R in rows]); ms = np.array([R['m'] for R in rows]); asts = np.array([R['ast'] for R in rows])
if np.all(lines_no > 0) and np.all(np.diff(lines_no) > 0): ok(f'Line numbers positive and strictly increasing ({lines_no[0]}..{lines_no[-1]})')
else: finding('FAIL', 'Line numbers not positive/strictly increasing')
if np.array_equal(lines_no, np.arange(1, N + 1)): ok('Line numbers consecutive 1..N (NOTE: must be renumbered globally when merged into the fleet file)')
else: finding('WARN', 'Line numbers not consecutive from 1')
if np.all(scs == 1): ok('SC_ID == 1 on every row (renumber when merging with other spacecraft)')
else: finding('FAIL', f'SC_ID values {sorted(set(scs))} (expected all 1 for a single-spacecraft file)')
if set(evs) <= {0, 1, 2, 3, 4}: ok(f'Event codes in {{0..4}}; counts: ' + ', '.join(f'{e}:{int((evs == e).sum())}' for e in range(5)))
else: finding('FAIL', f'illegal Event codes {sorted(set(evs) - {0,1,2,3,4})}')
if evs[0] == 0: ok('first row Event=0')
else: finding('FAIL', f'first row Event={evs[0]}')
if evs[-1] == 4: ok('last row Event=4')
else: finding('FAIL', f'last row Event={evs[-1]}')
if int((evs == 0).sum()) == 1 and int((evs == 4).sum()) == 1: ok('exactly one Event=0 and one Event=4 row')
else: finding('FAIL', f'{int((evs==0).sum())} Event=0 rows, {int((evs==4).sum())} Event=4 rows')
dts = np.diff(ts)
if np.all(dts > 0): ok(f'Time strictly increasing on all {N} rows; min row gap {dts.min():.6g} s, max row gap {dts.max():.6g} s')
else: finding('FAIL', f'Time not strictly increasing at rows {np.where(dts <= 0)[0] + 1}')
if ts.min() >= 0 and ts.max() <= TMISSION:
    ok(f'ALL rows 0 <= Time <= 473364000: first {ts[0]:.0f} s ({ts[0]/DAY:.2f} d), last {ts[-1]:.6f} s ({ts[-1]/DAY:.3f} d); margin to window end {TMISSION - ts[-1]:.1f} s = {(TMISSION - ts[-1])/DAY:.3f} d')
else: finding('FAIL', f'rows outside [0, T_max]: {int(((ts < 0) | (ts > TMISSION)).sum())}')
# thrust columns on non-Event=1 rows: exactly zero, both numerically and as literal tokens
bad = [R['ln'] for R in rows if R['ev'] != 1 and np.any(R['T'] != 0.0)]
badtok = [R['ln'] for R in rows if R['ev'] != 1 and any(tok != '0' for tok in R['Ttok'])]
if not bad: ok('thrust columns numerically (0,0,0) on every Event=0/2/3/4 row' + ('' if not badtok else f' (but {len(badtok)} rows use a non-"0" literal e.g. -0 / 0.0)'))
else: finding('FAIL', f'non-zero thrust on non-Event=1 rows at text lines {bad[:10]}')
# Asteroid_ID column
bad = [R['ln'] for R in rows if (R['ev'] == 3) != (1 <= R['ast'] <= 300)]
bad2 = [R['ln'] for R in rows if R['ev'] != 3 and R['ast'] != 0]
if not bad and not bad2: ok('Asteroid_ID in 1..300 exactly on Event=3 rows and 0 elsewhere')
else: finding('FAIL', f'Asteroid_ID misuse at text lines {bad[:10]} {bad2[:10]}')
# mass column
if ms[0] <= M0MAX and ms[0] >= MDRY: ok(f'm0 (Event=0 mass column) = {float(ms[0])!r} kg <= 2000 (margin {M0MAX - ms[0]:.3g} kg; bound is inclusive, eq. 11)')
else: finding('FAIL', f'm0 = {ms[0]} outside [600, 2000]')
if np.all(ms >= MDRY): ok(f'mass >= 600 kg on every row; final mass {ms[-1]:.6f} kg, margin above dry mass {ms[-1] - MDRY:.4f} kg')
else: finding('FAIL', f'mass < 600 kg on {int((ms < MDRY).sum())} rows (min {ms.min()})')
if np.all(np.diff(ms) <= 0): ok('mass column monotonically non-increasing')
else: finding('FAIL', f'mass increases at {int((np.diff(ms) > 0).sum())} row transitions (max +{np.diff(ms).max():.3e} kg)')
print(f'  propellant used m0 - m_final = {ms[0] - ms[-1]:.4f} kg of {ms[0] - MDRY:.1f} kg available ({100*(ms[0]-ms[-1])/(ms[0]-MDRY):.2f} %)')
x = (ms[0] - MDRY) / 1400.0; print(f'  cost of this spacecraft J_i = 1 + x + x^2 = {1 + x + x*x:.6f} (x = {x:.4f})')

# ---------------------------------------------------------------- 3. launch
print('--- launch (Event=0) vs independent Table-3 Earth ephemeris')
L = rows[0]
rE, vE = elem_state(EARTH, EPH_EARTH_MJD, L['t'])
dpos = np.linalg.norm(L['r'] - rE); vinf_vec = L['v'] - vE; vinf = np.linalg.norm(vinf_vec)
if dpos <= 1.0: ok(f'launch position error {dpos:.3e} km (tolerance 1 km)')
else: finding('FAIL', f'launch position error {dpos:.3f} km > 1 km')
if vinf <= VINFMAX: ok(f'launch v_inf = {vinf:.6f} km/s <= 4 km/s (margin {VINFMAX - vinf:.4f} km/s; vector {np.array2string(vinf_vec, precision=6)})')
elif vinf <= VINFMAX + 0.01: finding('WARN', f'v_inf = {vinf:.6f} km/s exceeds 4 km/s but within the 0.01 km/s validator tolerance')
else: finding('FAIL', f'v_inf = {vinf:.6f} km/s > 4.01 km/s')
print(f'  launch t = {L["t"]:.0f} s = {L["t"]/DAY:.2f} d after t0; Earth r = {np.array2string(rE, precision=3)}')
try:
    sys.path.insert(0, str(ROOT)); from ctoc14.kepler import Ephemeris  # only as a cross-check of the independent ephemeris
    eph = Ephemeris(); rE2, vE2 = eph.earth_state(L['t'])
    print(f'  cross-check independent Earth ephemeris vs ctoc14.kepler: |dr| = {np.linalg.norm(rE - rE2):.3e} km, |dv| = {np.linalg.norm(vE - vE2):.3e} km/s')
except Exception as ex:
    eph = None; print(f'  (ctoc14.kepler cross-check unavailable: {ex})')

# ---------------------------------------------------------------- 4. flybys
print('--- Event=3 flyby rows vs independent asteroid ephemeris (MEA.txt, epoch MJD 61200)')
mea = load_mea(MEA)
fly = []
for R in rows:
    if R['ev'] != 3: continue
    ra, va = elem_state(mea[R['ast']], EPH_AST_MJD, R['t'])
    d = np.linalg.norm(R['r'] - ra); vrel = np.linalg.norm(R['v'] - va)
    d2 = np.linalg.norm(R['r'] - eph.ast_state(R['ast'] - 1, R['t'])[0]) if eph else float('nan')
    fly.append((R['ast'], R['t'], d, vrel, R['ln'], d2))
    if d > DFLY: finding('FAIL', f'text line {R["ln"]}: asteroid {R["ast"]} distance {d:.1f} km > 1000 km')
    if R['t'] > TMISSION: finding('FAIL', f'text line {R["ln"]}: flyby after mission end')
fly_ids = [f[0] for f in fly]
dup = sorted({a for a in fly_ids if fly_ids.count(a) > 1})
print(f'  {len(fly)} Event=3 rows, {len(set(fly_ids))} distinct asteroids; duplicates: {dup if dup else "none"}')
if fly:
    dmax = max(fly, key=lambda f: f[2])
    ok(f'all flyby distances <= 1000 km; worst {dmax[2]:.3f} km (asteroid {dmax[0]}, text line {dmax[4]}), margin {DFLY - dmax[2]:.1f} km')
    print(f'  max |own - ctoc14.kepler| asteroid position at flyby times: {max(abs(f[2]-f[5]) for f in fly):.3e} km')
    print('  five largest flyby distances (km):', ', '.join(f'ast {f[0]}: {f[2]:.3f}' for f in sorted(fly, key=lambda f: -f[2])[:5]))
    print(f'  flyby relative speeds: min {min(f[3] for f in fly):.3f}, max {max(f[3] for f in fly):.3f} km/s')
    # sanity: sensitivity of the distance to a time error (how fast the geometry changes)
    print(f'  distance growth per 1 s of time error (max |v_rel|): {max(f[3] for f in fly):.3f} km/s -> 1000 km reached after >= {DFLY/max(f[3] for f in fly):.0f} s of timing error')

# ---------------------------------------------------------------- 5. thrust arcs, spacing, samples
print('--- thrust arcs / sampling')
arcs = []   # list of (first_idx, last_idx) in rows (indices of Event=1 rows spanning the run)
k = 0
while k < N:
    if evs[k] == 1:
        j = k; last1 = k
        while j + 1 < N and evs[j + 1] in (1, 3):
            j += 1
            if evs[j] == 1: last1 = j
        arcs.append((k, last1)); k = last1 + 1
    else: k += 1
print(f'  {len(arcs)} thrust arc(s); samples per arc: {[int(sum(1 for i in range(a, b+1) if evs[i]==1)) for a, b in arcs]}')
arc_data = []
for (a, b) in arcs:
    idx = [i for i in range(a, b + 1) if evs[i] == 1]
    tk = ts[idx]; Tk = np.array([rows[i]['T'] for i in idx])
    arc_data.append((tk, Tk))
    if len(tk) < 4: finding('INFO', f'arc starting text line {rows[a]["ln"]} has only {len(tk)} samples (degree {len(tk)-1} interpolation)')
    if len(tk) > 1:
        gaps = np.diff(tk)
        if gaps.min() < DTMIN: finding('FAIL', f'arc starting text line {rows[a]["ln"]}: sample spacing {gaps.min():.3f} s < 8640 s')
        else: ok(f'arc starting text line {rows[a]["ln"]}: Event=1 spacing (ignoring Event=3 rows) min {gaps.min():.1f} s, max {gaps.max():.1f} s (>= 8640 s)')
    Tn = np.linalg.norm(Tk, axis=1)
    if Tn.max() <= TMAX: ok(f'sample |T| max {Tn.max():.6f} N <= 0.5 N (margin {TMAX - Tn.max():.6f} N); zero-thrust samples: {int((Tn == 0).sum())}')
    else: finding('FAIL', f'sample |T| = {Tn.max()} > 0.5 N')
n_e3_inside = sum(1 for i in range(N) if evs[i] == 3 and any(tk[0] < ts[i] < tk[-1] for tk, _ in arc_data))
print(f'  Event=3 rows inside an arc: {n_e3_inside} of {int((evs==3).sum())}')
gap_e3 = [min(abs(ts[i] - tk).min() for tk, _ in arc_data) for i in range(N) if evs[i] == 3]
if gap_e3: print(f'  closest Event=3 row to a thrust sample: {min(gap_e3):.3f} s (no rule; validator integrates a tiny segment)')

# ---------------------------------------------------------------- 6. per-segment interpolant, exact |T| max
def lagrange_coeffs(tau_nodes, Tw):
    """Monomial coefficients (deg+1, 3) of the componentwise Lagrange interpolant through (tau_nodes, Tw)."""
    kk = len(tau_nodes); C = np.zeros((kk, 3))
    for mm in range(kk):
        others = [tau_nodes[q] for q in range(kk) if q != mm]
        Lm = P.polyfromroots(others) / np.prod([tau_nodes[mm] - o for o in others]) if others else np.array([1.0])
        C[:len(Lm)] += Lm[:, None] * Tw[mm][None, :]
    return C

def segment_poly(tA, tB, reading='A', arcs_in=None):
    """Interpolant on [tA, tB] as (C, tauA, tauB, span): T(t) = sum_k C[k] tau^k, tau = (t - t_j)/(t_{j+1}-t_j).
    arcs_in: list of (tk, Tk) sample sets to use instead of arc_data (for checker-variant sensitivity tests)."""
    for tk, Tk in (arc_data if arcs_in is None else arcs_in):
        n = len(tk)
        if n < 2 or tA < tk[0] or tB > tk[-1] + 1e-9: continue
        j = int(np.clip(np.searchsorted(tk, tA, side='right') - 1, 0, n - 2))
        if tB > tk[j + 1] + 1e-9 or tA < tk[j] - 1e-9: raise RuntimeError('segment straddles an interval')
        if n < 4: win = list(range(n))
        else:
            l = int(np.clip(j - 1, 0, n - 4)) if reading == 'A' else int(np.clip(j, 0, n - 4))
            win = [l, l + 1, l + 2, l + 3]
        span = tk[j + 1] - tk[j]
        tau_nodes = [(tk[w] - tk[j]) / span for w in win]
        C = lagrange_coeffs(tau_nodes, Tk[win])
        return C, (tA - tk[j]) / span, (tB - tk[j]) / span, span
    return np.zeros((1, 3)), 0.0, 1.0, tB - tA

print('--- exact maximum of interpolated |T(t)| on every row-to-row segment (Reading A: 0-based centred window)')
seg = []
peak_all = 0.0; peak_where = None; grid_disagree = 0.0
for i in range(N - 1):
    C, ta, tb, span = segment_poly(ts[i], ts[i + 1])
    S = np.zeros(1)
    for c in range(3): S = P.polyadd(S, P.polymul(C[:, c], C[:, c]))
    cands = [ta, tb]
    if len(S) > 2:
        rts = P.polyroots(P.polyder(S))
        cands += [float(r.real) for r in rts if abs(r.imag) < 1e-9 and ta < r.real < tb]
    vals = np.sqrt(np.maximum(P.polyval(np.array(cands), S), 0.0))
    pk = float(vals.max()); tpk = ts[i] + (cands[int(vals.argmax())] - ta) * span
    # brute-force grid cross-check
    g = np.linspace(ta, tb, 2001); gv = np.sqrt(np.maximum(P.polyval(g, S), 0.0)).max()
    grid_disagree = max(grid_disagree, gv - pk)
    seg.append(dict(i=i, C=C, ta=ta, tb=tb, span=span, peak=pk, tpeak=tpk, coast=bool(np.all(C == 0.0))))
    if pk > peak_all: peak_all, peak_where = pk, (i, tpk)
i, tpk = peak_where
if peak_all <= TMAX: ok(f'exact max |T(t)| over the whole file = {peak_all:.9f} N at t = {tpk:.1f} s (segment from Line {rows[i]["line"]}); margin {TMAX - peak_all:.6f} N')
else: finding('FAIL', f'interpolated |T(t)| = {peak_all:.9f} N > 0.5 N at t = {tpk:.1f} s (text line {rows[i]["ln"]})')
print(f'  brute-force 2001-point grid never exceeds the exact peak (max grid - exact = {grid_disagree:.2e} N, must be <= 0)')
i, tpk = peak_where; Sg = np.zeros(1)
for c in range(3): Sg = P.polyadd(Sg, P.polymul(seg[i]['C'][:, c], seg[i]['C'][:, c]))
g400 = np.linspace(seg[i]['ta'], seg[i]['tb'], 401); print(f'  at the peak segment a 400-point grid (ctoc14.validator resolution) finds {np.sqrt(P.polyval(g400, Sg)).max():.9f} N vs exact {peak_all:.9f} N')
top = sorted(seg, key=lambda s: -s['peak'])[:5]
print('  five highest segment peaks (N):', ', '.join(f'{s["peak"]:.6f}@Line {rows[s["i"]]["line"]}' for s in top))
# characterise the worst overshoot: magnitude step and direction change across the samples bounding the peak segment
ii = peak_where[0]; ev1 = [q for q in range(ii, -1, -1) if evs[q] == 1][:3][::-1]; ev2 = [q for q in range(ii + 1, N) if evs[q] == 1][:2]
prof = [(rows[q]['line'], np.linalg.norm(rows[q]['T'])) for q in ev1 + ev2]
print('  |T| of the samples around the peak (Line, N):', ', '.join(f'({l}, {n:.4f})' for l, n in prof))
a_, b_ = ev1[-2], ev1[-1]
ua, ub = rows[a_]['T'] / max(np.linalg.norm(rows[a_]['T']), 1e-12), rows[b_]['T'] / max(np.linalg.norm(rows[b_]['T']), 1e-12)
print(f'  step into the peak interval: |T| {np.linalg.norm(rows[a_]["T"]):.4f} -> {np.linalg.norm(rows[b_]["T"]):.4f} N with a {np.degrees(np.arccos(np.clip(ua @ ub, -1, 1))):.1f} deg direction change between Lines {rows[a_]["line"]} and {rows[b_]["line"]}')
smax = max(np.linalg.norm(Tk, axis=1).max() for _, Tk in arc_data)
finding('INFO', f'sample cap is {smax:.4f} N but the cubic interpolant overshoots to {peak_all:.6f} N (+{100*(peak_all/smax-1):.2f} % above the cap) at a hard magnitude step with a large direction change; margin to the 0.5 N hard bound is only {TMAX-peak_all:.4f} N. A 180 deg turn at a step can overshoot +12.8 % (docs), i.e. 0.45 N samples are NOT provably safe in general - the exact per-interval check must stay in the pipeline (it passes here).')
over_sample = [(s['peak'] - max(np.linalg.norm(rows[s['i']]['T']), np.linalg.norm(rows[s['i']+1]['T']))) for s in seg]
print(f'  largest overshoot of |T(t)| above both bounding samples: {max(over_sample):.6f} N')
n_over = sum(1 for o in over_sample if o > 1e-12); print(f'  segments whose peak lies strictly between samples: {n_over} of {len(seg)}')
coast_segs = [s for s in seg if s['coast']]
if coast_segs:
    lc = max(coast_segs, key=lambda s: s['span'] if s['C'].shape[0] == 1 else 0)
    L_coast = max(ts[s['i'] + 1] - ts[s['i']] for s in coast_segs)
    print(f'  pure-coast segments (interpolant identically 0): {len(coast_segs)}; longest {L_coast:.0f} s = {L_coast/DAY:.3f} d')
Lmax = dts.max(); imax = int(dts.argmax())
vmax = max(np.linalg.norm(R['v']) for R in rows); rmax = max(np.linalg.norm(R['r']) for R in rows)
print(f'  longest row-to-row segment {Lmax:.0f} s = {Lmax/DAY:.3f} d (text line {rows[imax]["ln"]}); max |v| {vmax:.3f} km/s, max |r| {rmax/AU:.3f} AU')
print(f'  17-digit print rounding: |dv| <= {vmax*2**-53:.2e} km/s -> {vmax*2**-53*Lmax:.2e} km over the longest segment; |dr| rounding {rmax*2**-53:.2e} km  (vs 1 km tolerance)')
print(f'  10-digit rounding (PDF minimum) would give {vmax*5e-10*Lmax:.2e} km over the longest segment')
# Time-column rounding: 17 significant digits of a 4.5e8 s epoch is ~6e-8 s; at 36 km/s that is ~2e-6 km
tmax_tok = max(abs(R['t']) for R in rows)
print(f'  Time-token rounding: ulp({tmax_tok:.0f} s) = {np.spacing(tmax_tok):.1e} s -> {np.spacing(tmax_tok)*vmax:.1e} km position equivalent at max |v|')

def exact_peak(C, ta, tb):
    S = np.zeros(1)
    for c in range(3): S = P.polyadd(S, P.polymul(C[:, c], C[:, c]))
    cands = [ta, tb]
    if len(S) > 2:
        rts = P.polyroots(P.polyder(S))
        cands += [float(r.real) for r in rts if abs(r.imag) < 1e-9 and ta < r.real < tb]
    vals = np.sqrt(np.maximum(P.polyval(np.array(cands), S), 0.0))
    return float(vals.max()), cands[int(vals.argmax())]

# 6b. exact peak under the FORWARD window (Reading B) - what a checker with the other index convention would measure
pkB = 0.0; whereB = None
for i in range(N - 1):
    C, ta, tb, span = segment_poly(ts[i], ts[i + 1], reading='B')
    pk, tau = exact_peak(C, ta, tb)
    if pk > pkB: pkB, whereB = pk, (i, ts[i] + (tau - ta) * span)
print(f'  Reading-B (forward window) exact max |T(t)| = {pkB:.9f} N at t = {whereB[1]:.1f} s (segment from Line {rows[whereB[0]]["line"]}); margin {TMAX - pkB:+.6f} N')
if pkB > TMAX: finding('INFO', f'under the forward-window reading the interpolant would reach {pkB:.6f} N > 0.5 N (Reading A, the PDF-example-confirmed one, peaks at {peak_all:.6f} N)')

# 6c. thrust-profile shape statistics (what drives the cubic overshoot)
print('--- thrust profile statistics (all Event=1 samples)')
tk_all = np.concatenate([tk for tk, _ in arc_data]); Tk_all = np.concatenate([Tk for _, Tk in arc_data]); Tn_all = np.linalg.norm(Tk_all, axis=1)
hist, edges = np.histogram(Tn_all, bins=[0, 1e-9, 0.1, 0.2, 0.3, 0.4, 0.45 - 1e-9, 0.45 + 1e-9, 0.5])
print('  |T| histogram [N]: ' + ', '.join(f'[{edges[k]:.2f},{edges[k+1]:.2f}): {hist[k]}' for k in range(len(hist))))
print(f'  mean sample |T| {Tn_all.mean():.4f} N; samples at the 0.45 N cap: {int((np.abs(Tn_all - 0.45) < 1e-12).sum())} (within 1e-12 N), {int((Tn_all > 0.449).sum())} above 0.449 N; zero samples: {int((Tn_all == 0).sum())} at Lines {[rows[i]["line"] for i in range(N) if evs[i] == 1 and np.all(rows[i]["T"] == 0)]}')
steps = np.abs(np.diff(Tn_all))
un = Tk_all / np.maximum(Tn_all[:, None], 1e-15); cosang = np.einsum('ij,ij->i', un[:-1], un[1:]); ang = np.degrees(np.arccos(np.clip(cosang, -1, 1)))
both = (Tn_all[:-1] > 1e-6) & (Tn_all[1:] > 1e-6)
print(f'  magnitude steps between consecutive samples: > 0.1 N: {int((steps > 0.1).sum())}, > 0.2 N: {int((steps > 0.2).sum())}, > 0.3 N: {int((steps > 0.3).sum())}, max {steps.max():.4f} N')
print(f'  direction change between consecutive non-zero samples: max {ang[both].max():.1f} deg; > 30 deg: {int((ang[both] > 30).sum())}, > 60 deg: {int((ang[both] > 60).sum())}, > 90 deg: {int((ang[both] > 90).sum())}')
big_turn_hi = [(rows[[q for q in range(N) if evs[q]==1][k]]['line'], ang[k], Tn_all[k], Tn_all[k+1]) for k in range(len(ang)) if both[k] and ang[k] > 90 and max(Tn_all[k], Tn_all[k+1]) > 0.3]
print(f'  turns > 90 deg where either sample exceeds 0.3 N: {len(big_turn_hi)}' + (': ' + ', '.join(f'Line {l} ({a:.0f} deg, {n1:.3f}->{n2:.3f} N)' for l, a, n1, n2 in big_turn_hi[:6]) if big_turn_hi else ''))
# how a coast is encoded: interior zero samples inside the arc (instead of Event=2 nodes)
interior_zero = [rows[i]['line'] for i in range(N) if evs[i] == 1 and np.all(rows[i]['T'] == 0) and any(tk[0] < ts[i] < tk[-1] for tk, _ in arc_data)]
print(f'  Event=2 coast nodes: {int((evs == 2).sum())}; zero-thrust samples strictly inside an arc (coast encoded inside the arc): {len(interior_zero)} at Lines {interior_zero}')
# rule 6.1-4 misapplied to ALL consecutive rows inside the arc (a sloppy checker that does not skip Event=3 rows)
gaps_all = [(dts[i], rows[i]['line'], evs[i], evs[i + 1]) for i in range(N - 1) if any(tk[0] <= ts[i] and ts[i + 1] <= tk[-1] for tk, _ in arc_data)]
short = sorted(g for g in gaps_all if g[0] < DTMIN)
print(f'  row-to-row gaps inside the arc shorter than 8640 s (all involve an Event=3 row; rule 6.1-4 says to ignore them): {len(short)}; shortest {short[0][0]:.1f} s at Line {short[0][1]} (Event {short[0][2]}->{short[0][3]})' if short else '  no row-to-row gap inside the arc is shorter than 8640 s')
finding('INFO', f'every one of the {int((evs==3).sum())} flybys sits INSIDE the single {len(tk_all)}-sample thrust arc (no Event=2 node anywhere); {len(short)} row gaps are < 8640 s because of Event=3 rows. This is exactly the structure the PDF allows (6.1-2/6.1-4) and mirrors the PDF example, but the file has zero tolerance for a checker that (wrongly) breaks arcs at Event=3 rows, treats them as samples, or applies the 8640 s rule to all rows - quantified below.')

# ---------------------------------------------------------------- 7. own dynamics check (vectorised RK4, several step sizes)
print('--- own dynamics consistency: fixed-step RK4 from every row to the next (models a checker with another integrator)')
K = max(s['C'].shape[0] for s in seg)
Cs = np.zeros((N - 1, K, 3))
for s in seg: Cs[s['i'], :s['C'].shape[0]] = s['C']
tas = np.array([s['ta'] for s in seg]); tbs = np.array([s['tb'] for s in seg])
Y0 = np.array([np.concatenate([rows[i]['r'], rows[i]['v'], [rows[i]['m']]]) for i in range(N - 1)])
Y1 = np.array([np.concatenate([rows[i]['r'], rows[i]['v'], [rows[i]['m']]]) for i in range(1, N)])
Lseg = dts.copy()

def thrust_at(sfrac, Cs):
    tau = tas + sfrac * (tbs - tas)
    T = np.zeros((N - 1, 3))
    for k in range(Cs.shape[1] - 1, -1, -1): T = T * tau[:, None] + Cs[:, k]
    return T

def rhs(sfrac, Y, Cs):
    r = Y[:, 0:3]; v = Y[:, 3:6]; m = Y[:, 6:7]
    T = thrust_at(sfrac, Cs)
    rn = np.linalg.norm(r, axis=1, keepdims=True)
    acc = -MU * r / rn ** 3 + T / m * 1e-3
    return np.hstack([v, acc, -np.linalg.norm(T, axis=1, keepdims=True) / VEX])

def rk4_all(nsteps, Cs):
    Y = Y0.copy(); h = (Lseg / nsteps)[:, None]
    for k in range(nsteps):
        s0 = k / nsteps; s1 = (k + 0.5) / nsteps; s2 = (k + 1) / nsteps
        k1 = rhs(s0, Y, Cs); k2 = rhs(s1, Y + 0.5 * h * k1, Cs); k3 = rhs(s1, Y + 0.5 * h * k2, Cs); k4 = rhs(s2, Y + h * k3, Cs)
        Y = Y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return Y

for nsteps in (10, 96, 1440):
    Y = rk4_all(nsteps, Cs)
    ep = np.linalg.norm(Y[:, 0:3] - Y1[:, 0:3], axis=1); ev = np.linalg.norm(Y[:, 3:6] - Y1[:, 3:6], axis=1); em = np.abs(Y[:, 6] - Y1[:, 6])
    tag = f'RK4 {nsteps:5d} steps/segment (h = {Lmax/nsteps:7.1f} s on 1-day segments)'
    print(f'  {tag}: max |dr| {ep.max():.3e} km (line {rows[int(ep.argmax())]["ln"]}), max |dv| {ev.max()*1e3:.3e} m/s, max |dm| {em.max():.3e} kg')
    if nsteps == 1440:
        ep_A1440 = float(ep.max())
        if ep.max() <= 1.0 and ev.max() <= 1e-3 and em.max() <= 0.01: ok('all segments within 1 km / 1 m/s / 0.01 kg with an independent integrator')
        else: finding('FAIL', 'independent integrator disagrees with the file beyond the tolerances')
        print(f'  fraction of the 1 km / 1 m/s / 0.01 kg budgets used (worst segment): {100*ep.max():.2e} % / {100*ev.max()/1e-3:.2e} % / {100*em.max()/0.01:.2e} %')
# mass bookkeeping from the integrated |T|
print(f'  integrated propellant (sum of per-segment RK4 mass drops) = {(Y0[:, 6] - Y[:, 6]).sum():.4f} kg vs file {ms[0]-ms[-1]:.4f} kg')

# Reading-B sensitivity (forward window): what a checker with the other index convention would see
print('--- sensitivity: if the checker used the FORWARD window (Reading B, not the PDF-example-confirmed one)')
segB = [segment_poly(ts[i], ts[i + 1], reading='B') for i in range(N - 1)]
CsB = np.zeros_like(Cs)
for i, (C, ta, tb, span) in enumerate(segB): CsB[i, :C.shape[0]] = C
Y = rk4_all(96, CsB)
print(f'  interpolant differs between readings on {int(np.any(np.abs(CsB - Cs) > 1e-12, axis=(1, 2)).sum())} of {N-1} segments')
ep = np.linalg.norm(Y[:, 0:3] - Y1[:, 0:3], axis=1); ev = np.linalg.norm(Y[:, 3:6] - Y1[:, 3:6], axis=1); em = np.abs(Y[:, 6] - Y1[:, 6])
print(f'  Reading B: max |dr| {ep.max():.3e} km, max |dv| {ev.max()*1e3:.3e} m/s, max |dm| {em.max():.3e} kg; segments > 1 km: {int((ep > 1).sum())}, > 1 m/s: {int((ev > 1e-3).sum())}')
if (ep > 1).any() or (ev > 1e-3).any() or (em > 0.01).any():
    finding('INFO', f'validity depends on the centred (0-based) Lagrange window: with the forward window {int((ep > 1).sum())} segments exceed 1 km (max {ep.max():.0f} km), {int((ev > 1e-3).sum())} exceed 1 m/s. Reading A is the one that reproduces the PDF example to 2e-10 kg (docs/thrust_interpolation_rule.md), so this is a documented dependency, not a defect.')

# ---------------------------------------------------------------- 7b. naive-checker variants around Event=3 rows
print('--- sensitivity: checker variants that mishandle Event=3 rows inside an arc (both contradict PDF 6.1-2/6.1-3)')
def build_arcs(sample_events, breakers):
    """Arcs = maximal runs of rows whose Event is in sample_events; rows with Event in `breakers` end a run;
    every row of the run is a sample (Event=3 rows contribute their (0,0,0) thrust columns)."""
    out = []; k = 0
    while k < N:
        if evs[k] in sample_events and evs[k] not in breakers:
            j = k
            while j + 1 < N and evs[j + 1] in sample_events and evs[j + 1] not in breakers: j += 1
            idx = list(range(k, j + 1)); out.append((ts[idx], np.array([rows[i]['T'] for i in idx]))); k = j + 1
        else: k += 1
    return out
variants = {'C1: Event=3 row terminates the arc (like Event=2)': build_arcs({1}, set()),
            'C2: Event=3 row is a zero-thrust sample of the arc': build_arcs({1, 3}, set())}
naive_worst = {}
for name, arcs_v in variants.items():
    CsV = np.zeros_like(Cs)
    for i in range(N - 1):
        C, ta, tb, span = segment_poly(ts[i], ts[i + 1], arcs_in=arcs_v); CsV[i, :C.shape[0]] = C
    Y = rk4_all(96, CsV)
    ep = np.linalg.norm(Y[:, 0:3] - Y1[:, 0:3], axis=1); ev = np.linalg.norm(Y[:, 3:6] - Y1[:, 3:6], axis=1); em = np.abs(Y[:, 6] - Y1[:, 6])
    nd = int(np.any(np.abs(CsV - Cs) > 1e-12, axis=(1, 2)).sum())
    naive_worst[name] = (ep.max(), int((ep > 1).sum()), int((ev > 1e-3).sum()), int((em > 0.01).sum()))
    print(f'  {name}: interpolant differs on {nd} segments; max |dr| {ep.max():.3e} km, max |dv| {ev.max()*1e3:.3e} m/s, max |dm| {em.max():.3e} kg; segments failing 1 km: {int((ep > 1).sum())}, 1 m/s: {int((ev > 1e-3).sum())}, 0.01 kg: {int((em > 0.01).sum())}')
finding('INFO', 'a checker that breaks the arc at Event=3 rows or treats them as zero-thrust samples would fail ' +
        ', '.join(f'{v[1]} segments (max {v[0]:.0f} km) [{k.split(":")[0]}]' for k, v in naive_worst.items()) +
        '. The PDF text is explicit that neither is correct, and the PDF example itself is only consistent with Reading A, so this is a robustness note: moving flybys into short Event=2-bounded coasts would remove the dependency at some fuel cost.')

# ---------------------------------------------------------------- 7c. chained (non-restarting) integration drift
print('--- sensitivity: a checker that integrates continuously from the launch row WITHOUT restarting at each row')
from scipy.integrate import solve_ivp
def chained_drift():
    y = np.concatenate([rows[0]['r'], rows[0]['v'], [rows[0]['m']]]); worst = (0.0, 0, 0.0, 0.0); drift = []
    for i in range(N - 1):
        C = Cs[i]; ta, tb = tas[i], tbs[i]; L = Lseg[i]
        def f(t, y, C=C, ta=ta, tb=tb, L=L, t0=ts[i]):
            tau = ta + (t - t0) / L * (tb - ta)
            T = np.zeros(3)
            for k in range(C.shape[0] - 1, -1, -1): T = T * tau + C[k]
            r = y[0:3]; rn = np.linalg.norm(r)
            return np.concatenate([y[3:6], -MU * r / rn ** 3 + T / y[6] * 1e-3, [-np.linalg.norm(T) / VEX]])
        sol = solve_ivp(f, (ts[i], ts[i + 1]), y, method='DOP853', rtol=1e-13, atol=np.array([1e-9] * 3 + [1e-12] * 3 + [1e-10]))
        y = sol.y[:, -1]
        dr = np.linalg.norm(y[0:3] - rows[i + 1]['r']); dv = np.linalg.norm(y[3:6] - rows[i + 1]['v']); dm = abs(y[6] - rows[i + 1]['m'])
        drift.append((dr, dv, dm))
        if dr > worst[0]: worst = (dr, rows[i + 1]['line'], dv, dm)
    return np.array(drift), worst
drift, worst = chained_drift()
first_over = np.argmax(drift[:, 0] > 1.0) if (drift[:, 0] > 1.0).any() else None
print(f'  DOP853 rtol 1e-13, chained over all {N-1} segments: max |dr| {worst[0]:.3e} km at Line {worst[1]} (|dv| there {worst[2]*1e3:.3e} m/s, |dm| {worst[3]:.2e} kg); drift at the final row {drift[-1,0]:.3e} km / {drift[-1,1]*1e3:.3e} m/s / {drift[-1,2]:.2e} kg')
print(f'  rows where the accumulated drift exceeds 1 km: {int((drift[:, 0] > 1).sum())}' + (f' (first at Line {rows[first_over + 1]["line"]}, t = {ts[first_over + 1]/DAY:.1f} d)' if first_over is not None else ''))
print(f'  flyby rows: max accumulated drift {max(drift[i - 1, 0] for i in range(1, N) if evs[i] == 3):.3e} km (vs the 1000 km flyby ball)')
if (drift[:, 0] > 1).any():
    finding('INFO', f'a non-restarting checker would see up to {worst[0]:.1f} km of accumulated drift ({int((drift[:, 0] > 1).sum())} rows > 1 km); the PDF (Sec. 7 bullet 2, 从每一行状态出发) prescribes the per-row restart, so this is informational only.')
else:
    ok(f'even a non-restarting checker stays within 1 km over the whole 14.4-year file (max {worst[0]:.3e} km)')

# ---------------------------------------------------------------- 7d. fleet-level cross-check (read-only look at sibling files)
print('--- fleet context: sibling results/fleet_b300/sub_sc*.txt files (duplicate flybys across spacecraft are legal but wasted)')
import glob
sib = sorted(p for p in glob.glob(str(SUB.parent / 'sub_sc*.txt')) if Path(p) != SUB)
my_ids = set(fly_ids); seen = {}
for p in sib:
    ids_p = []
    for line in open(p, encoding='utf-8'):
        s = line.strip()
        if not s or s.startswith('#'): continue
        q = s.split()
        if q[2] == '3': ids_p.append(int(q[14]))
    ov = sorted(my_ids & set(ids_p))
    for a in ids_p: seen.setdefault(a, []).append(Path(p).name)
    print(f'  {Path(p).name}: {len(set(ids_p))} distinct flybys; overlap with this file: {ov if ov else "none"}')
ov_all = sorted(a for a in my_ids if a in seen)
if ov_all: finding('WARN', f'asteroids {ov_all} are flown by this spacecraft AND by {sorted({f for a in ov_all for f in seen[a]})} - the repeat earns nothing')
else: ok(f'none of the {len(my_ids)} asteroids of this file appears in any sibling file')
print(f'  NOTE for the merge: this file uses SC_ID=1 and Line 1..{N} (both must be renumbered), and its first text line is a "#" comment ("{lines[0][:60]}") which the PDF allows.')

# ---------------------------------------------------------------- 8. summary
print('=== SUMMARY ===')
fails = [f for f in FINDINGS if f[0] == 'FAIL']; warns = [f for f in FINDINGS if f[0] == 'WARN']; infos = [f for f in FINDINGS if f[0] == 'INFO']
print(f'FAIL: {len(fails)}  WARN: {len(warns)}  INFO: {len(infos)}')
for sev, txt in FINDINGS: print(f'  {sev}: {txt}')
print('margins:')
print(f'  v_inf            {vinf:.6f} km/s   (limit 4, margin {VINFMAX - vinf:.4f})')
print(f'  launch pos err   {dpos:.3e} km       (limit 1)')
print(f'  m0               {ms[0]:.6f} kg    (limit 2000 inclusive)')
print(f'  final mass       {ms[-1]:.6f} kg   (limit 600, margin {ms[-1]-MDRY:.4f})')
print(f'  max sample |T|   {max(np.linalg.norm(Tk, axis=1).max() for _, Tk in arc_data):.6f} N  (limit 0.5)')
print(f'  max interp |T|   {peak_all:.9f} N  (limit 0.5, margin {TMAX-peak_all:.6f}); forward-window reading would give {pkB:.6f} N')
if fly: print(f'  worst flyby dist {dmax[2]:.3f} km (ast {dmax[0]}; limit 1000, margin {DFLY-dmax[2]:.1f}); 2nd worst {sorted(fly, key=lambda f: -f[2])[1][2]:.3f} km')
print(f'  last row time    {ts[-1]:.3f} s ({(TMISSION-ts[-1])/DAY:.3f} d before window end)')
print(f'  longest segment  {Lmax:.0f} s; 17-digit rounding {vmax*2**-53*Lmax:.1e} km; own RK4-1440 worst segment {ep_A1440:.2e} km')
print(f'  chained drift    {worst[0]:.3e} km max (non-restarting checker)')
print(f'  naive Event=3    ' + '; '.join(f'{k.split(":")[0]} max {v[0]:.0f} km, {v[1]} segs > 1 km' for k, v in naive_worst.items()))
sys.exit(1 if fails else 0)
