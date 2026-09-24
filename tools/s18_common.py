"""Stage 18 (RHFA, rolling-horizon fleet auction): constants, class table, twin helpers, schemas, regrid law.

Every s18 module imports this file and nothing else of s18 (so the three work packages of docs/stage18/BUILD.md never
touch each other's files).  Existing tools (run_ialns, s14_twinbeam, s15b_regrid, ctoc14.*) are IMPORTED, never edited.

UNITS AND IDS (binding for every s18 file)
  * time            seconds of mission time internally (t = 0 is 2030-01-01); DAYS only in logs, JSON fields whose name
                    ends in `_d`, CLI arguments (--H 540 --L 360 --stop-T 1461) and results/s18/windows.json (`t_d`).
                    Convert with d2s / s2d.  Never store a bare "t" in days in a state or column file.
  * target ids      1-based ints 1..300 (MEA.txt order); 131 and 144 are UNREACHABLE; TARGETS = the 298 others.
                    Ephemeris calls take the 0-based index (id - 1) -- keep that conversion inside the function that
                    calls eph.ast_states_at, never in a schema.
  * craft index     0-based int 0..N-1 (log name c<i>); a ROOT column has craft = None.
  * twin file       run_ialns "ist" npz: tL (s), vinf (3,), ts (M,), Ts (M,3) [km/s], tf (n,), asts (n,) int.
  * tank            kg, = (600 + 1.5) exp(dv / ve) (ImpulsiveProblem.tank).  A tank is HONEST only after regrid + settle
                    (regrid_settle below).  Column tanks are planner-twin tanks (irregular nodes) and are so labelled.
  * J_i             cost(tank) = 1 + x + x^2, x = (tank - 600) / 1400.  Fleet J = sum J_i + (298 - covered) + 2.

THE REGRID LAW (law 8 of docs/stage17_brief.md)
  A twin born from the planner (s14_twinbeam.settle_child) carries irregular impulse nodes (junction impulses 1 d after a
  flyby, LinLeg profiles on a 10 d grid), so two impulses can sit in one 20 d bin, each at the cap: 1.1-1.6x the real
  0.43 N capability.  EVERY tank that is committed to a craft state, written to fleet.json, or reported as a result is
  taken from regrid_settle(ip) (s15b_regrid.regrid to regular 20 d bins + run_ialns.settle), and honest_check(ip) must
  report 0 irregular nodes and max_window_cap <= CAP_TOL.  Column tanks (pre-commit) are used only to RANK columns.
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, signal, pathlib, datetime
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import numpy as np
import run_ialns as RI
from s15b_regrid import regrid
from ctoc14.impulsive import ImpulsiveProblem
from ctoc14.search import UNREACHABLE
from ctoc14.constants import DAY, T_MISSION, M_DRY, FUEL_MAX, VE, AU, VINF_MAX

# ------------------------------------------------------------------ constants
PY = os.path.expanduser('~/.venvs/astro313/bin/python')
NICE = ['nice', '-n', '5']
S18 = ROOT / 'results/s18'
EVENTS_SRC = ROOT / 'results/s17/physics/events.json'     # list of {ast, t (DAYS), d (AU), ...}: torus-distance minima
WINDOWS_FILE = S18 / 'windows.json'
TARGETS = [t for t in range(1, 301) if t not in UNREACHABLE]   # the 298 reachable ids
N_TARGETS = len(TARGETS)
YEAR = 365.25 * DAY
DTAU = 20 * DAY               # twin node spacing (regular bins)
TOL_MISS = 150.0              # km: a twin is settled when every flyby misses by <= this
CAP_TOL = 1.02                # honest: max (impulse / tcap) per 20 d window (run_ialns.settle enforces 1.02)
M0P = 1600.0                  # planner launch mass used by s14_twinbeam scores (production convention)
TANK_DRY = M_DRY + 1.5        # 601.5 kg: tank of a zero-dv twin (ImpulsiveProblem margin 1.5)
D_WIN = 0.05                  # AU: a node event closer than this is a window (results/s17/physics classes)
D_WIN_FALLBACK = 0.08         # AU: 19 targets have no event < 0.05 AU; their windows are the events < 0.08 AU
PIN_MAX = 4                   # class pin: <= 4 windows (at D_WIN) in 15 yr
FILL_MIN = 10                 # class fil: >= 10 windows; med otherwise
CLASSES = ('pin', 'med', 'fil')
T_STOP_DEFAULT = T_MISSION - 60 * DAY     # (kept for reference; the loop now uses slice_H / T_END below)
T_END = T_MISSION - 2 * DAY               # last usable flyby epoch (epoch_candidates_h caps at the same value)
H_TAIL_MIN = 120 * DAY                    # a slice never leaves a tail shorter than this: the last slice absorbs it
FUEL_BAR_N8 = 3764.0                      # kg: the 09-26 bar for N = 8 covering 298 (mean m0 1070); docs/stage18_rhfa.md 0

# default parameters of a run (docs/stage18_rhfa.md section 2; days where the CLI says days)
# Fixer 09-24 (reviewer findings 3 + 4): m0max 1070 = 600 + 3764 / 8 is the bar itself (slack 0 for the auction; m0slack
# is now only the regrid tolerance of the commit), a per-craft budget RAMP (ramp_alpha / ramp_slack, ramp_cap below)
# paces the spend over the mission, a fleet knapsack row per slice (fuel_bar x H / T_MISSION x fleet_factor), lam
# starts at 1.0 and acts as the dual of that budget (s18_rhfa.update_lam); quotas 12 / 18 with the last 3 slices free.
DEF = dict(N=8, H=540.0, L=360.0, B=8, tries=16, beta=1.5, gamma=0.3, lam=1.0, qf=12, qm=18, m0max=1070.0,
           m0slack=10.0, nproc=8, seed=0, k_epochs=3, lin_cap=3.0, res_max=1.5e8, iters=100, timeout=240.0,
           grid_step=10.0, roots=60, max_levels=14, wall_job=2400.0, regrid_iters=150, regrid_timeout=600.0,
           beam_edf=1, proposer='segbeam', commit='regrid', stop_T=None, root_bin=60.0,
           # beam knobs added by the integrator (09-24): branch diversity cap (s18_segbeam.JOB_DEF) and the beam score
           # weights (None = s14 production w_t 1.0 / w_fuel 0.52; WP1 measured lighter columns at w_t 0.25)
           div_cap=3, w_t=None, w_fuel=None,
           # fixer 09-24: fuel pacing (finding 3), quota window (finding 4), beam stopping (finding 6), pool safety (1)
           fuel_bar=FUEL_BAR_N8, ramp_alpha=0.9, ramp_slack=60.0, fleet_factor=1.3, lam_max=6.0,
           quota_free_slices=3, look_stop=2, settle_cap=8, settle_ok_min=0.3, job_grace=600.0)


def slice_H(T, H):
    """Effective horizon [s] of the slice starting at T (finding 2: the last slice absorbs the mission tail instead
    of leaving it unused).  span = T_END - T; H_eff = span when span <= H + H_TAIL_MIN (this is the last slice, at most
    H + 120 d long), else H.  Returns 0.0 when nothing is left (T >= T_END)."""
    T = float(T); H = float(H); span = T_END - T
    if span <= 0.0:
        return 0.0
    return span if span <= H + H_TAIL_MIN else H


def ramp_cap(t_end, m0max, alpha=None, slack=None, prm=None):
    """Per-craft tank cap [kg] for a column whose committed window ends at t_end [s] (finding 3): TANK_DRY +
    (m0max - TANK_DRY) * min(1, t_end / T_MISSION)^alpha + slack, never above m0max (the bar is hard)."""
    alpha = float(DEF['ramp_alpha'] if alpha is None else alpha) if prm is None else float(prm.get('ramp_alpha', DEF['ramp_alpha']))
    slack = float(DEF['ramp_slack'] if slack is None else slack) if prm is None else float(prm.get('ramp_slack', DEF['ramp_slack']))
    m0max = float(m0max)
    f = min(1.0, max(0.0, float(t_end) / T_MISSION)) ** alpha
    return min(m0max, TANK_DRY + (m0max - TANK_DRY) * f + slack)


def slice_fuel_budget(H_eff, prm):
    """Fleet propellant budget [kg] of one slice (finding 3, the knapsack row): fuel_bar x H_eff / T_MISSION x
    fleet_factor."""
    return float(prm.get('fuel_bar', DEF['fuel_bar'])) * float(H_eff) / T_MISSION * float(prm.get('fleet_factor', DEF['fleet_factor']))


def quota_active(T, H_eff, prm):
    """False in the last `quota_free_slices` slices (T + H_eff > T_END - n H): the class quotas are lifted there
    (finding 4: with zero slack the quotas strand fillers at the end)."""
    n = int(prm.get('quota_free_slices', DEF['quota_free_slices'])); H = d2s(prm.get('H', DEF['H']))
    return float(T) + float(H_eff) <= T_END - n * H + 1e-6


def set_threads():
    """BLAS single-threaded in this process (the beam is parallel over craft, never inside a craft)."""
    for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[v] = '1'


def d2s(d):
    return float(d) * DAY


def s2d(s):
    return float(s) / DAY


def now():
    return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')


def cost(tank):
    """J_i of a craft with launch mass `tank` [kg]."""
    x = (float(tank) - M_DRY) / FUEL_MAX
    return 1.0 + x + x * x


def sum_Ji(tanks):
    return float(sum(cost(t) for t in tanks))


def fleet_J(tanks, covered):
    """Raw contest J of a fleet: sum J_i + misses (298 - covered) + 2 (the two unreachable targets)."""
    return sum_Ji(tanks) + (N_TARGETS - int(covered)) + 2.0


def fuel_of(tanks):
    """Total propellant [kg] = sum (tank - 600)."""
    return float(sum(float(t) - M_DRY for t in tanks))


# ------------------------------------------------------------------ files
def jdump(obj, path):
    """Atomic JSON write (tmp + rename); numpy scalars/arrays and sets are converted."""
    path = pathlib.Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f'.tmp{os.getpid()}')

    def cv(o):
        if isinstance(o, np.ndarray): return o.tolist()
        if isinstance(o, np.floating): return float(o)
        if isinstance(o, np.integer): return int(o)
        if isinstance(o, (set, frozenset)): return sorted(o)
        if isinstance(o, pathlib.Path): return str(o)
        raise TypeError(f'not JSON serialisable: {type(o).__name__}')
    with open(tmp, 'w') as f:
        json.dump(obj, f, indent=1, default=cv)
    os.replace(tmp, path)


def jload(path):
    with open(path) as f:
        return json.load(f)


def logger(path):
    """Line logger `say(s)` -> stdout + append to `path` with a HH:MM:SS stamp (run_ialns.logger)."""
    pathlib.Path(path).parent.mkdir(parents=True, exist_ok=True)
    return RI.logger(path)


def run_dir(name):
    """results/s18/<name> (a name may also be an absolute or ROOT-relative path)."""
    p = pathlib.Path(name)
    if p.is_absolute():
        return p
    if str(p).startswith('results/'):
        return ROOT / p
    return S18 / p


def slice_dir(rd, k):
    """results/s18/<run>/slices/s<kk>/ (per-slice jobs, columns.json, auction.json, commit.json)."""
    return pathlib.Path(rd) / 'slices' / f's{int(k):02d}'


def craft_npz(rd, i):
    """The committed, regridded, settled twin of craft i (run_ialns ist format)."""
    return pathlib.Path(rd) / f'craft_{int(i)}.npz'


class Timeout(Exception):
    pass


def _alarm(sig, frm):
    raise Timeout('timeout')


def with_timeout(fn, secs, *args, **kw):
    """Run fn(*args) under a SIGALRM wall limit (main thread of the calling process only)."""
    signal.signal(signal.SIGALRM, _alarm); signal.setitimer(signal.ITIMER_REAL, float(secs))
    try:
        return fn(*args, **kw)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


# ------------------------------------------------------------------ windows and classes
def build_windows(events=EVENTS_SRC, out=WINDOWS_FILE, d_win=D_WIN, d_fallback=D_WIN_FALLBACK):
    """results/s18/windows.json from results/s17/physics/events.json.

    Format: {"meta": {source, d_win, d_fallback, pin_max, fill_min, built, n_fallback},
             "targets": {"<id>": {"t_d": [event epochs, DAYS, ascending], "n": windows at d_win,
                                  "cls": "pin"|"med"|"fil", "fallback": bool}}}
    windows(t) = the target's node events with d < d_win; class from THAT count (pin <= 4, fil >= 10, med otherwise),
    exactly results/s17/physics/history_classes.py (95 pin / 123 med / 80 fil).  The 19 targets with no event under
    d_win (their minima sit at 0.050-0.060 AU) keep class pin but get their events under d_fallback as windows, so that
    the earliest-deadline bonus never degenerates to "now or never for ever".  Both unreachable ids are omitted."""
    ev = jload(events)
    by = {}
    for e in ev:
        by.setdefault(int(e['ast']), []).append((float(e['t']), float(e['d'])))
    T = {}; nfb = 0
    for X in TARGETS:
        lst = sorted(by.get(X, []))
        w = [t for t, d in lst if d < d_win]
        n = len(w); fb = False
        if not w:
            w = [t for t, d in lst if d < d_fallback]; fb = True; nfb += 1
        cls = 'pin' if n <= PIN_MAX else ('fil' if n >= FILL_MIN else 'med')
        T[str(X)] = dict(t_d=w, n=n, cls=cls, fallback=fb)
    W = dict(meta=dict(source=str(events), d_win=d_win, d_fallback=d_fallback, pin_max=PIN_MAX, fill_min=FILL_MIN,
                       built=now(), n_fallback=nfb), targets=T)
    if out is not None:
        jdump(W, out)
    return _index_windows(W)


def _index_windows(W):
    """Runtime table: {ast (int): {"t": np.array (SECONDS, ascending), "n": int, "cls": str, "fallback": bool}}."""
    out = {}
    for k, v in W['targets'].items():
        out[int(k)] = dict(t=np.asarray(v['t_d'], float) * DAY, n=int(v['n']), cls=str(v['cls']), fallback=bool(v['fallback']))
    out['_meta'] = W.get('meta', {})
    return out


def load_windows(path=WINDOWS_FILE):
    """Class table (see build_windows); built from EVENTS_SRC on first use."""
    path = pathlib.Path(path)
    if not path.exists():
        return build_windows(out=path)
    return _index_windows(jload(path))


def class_of(W, ast):
    return W[int(ast)]['cls']


def count_classes(W, targets):
    """{'pin': n, 'med': n, 'fil': n} of a target list."""
    c = dict(pin=0, med=0, fil=0)
    for t in targets:
        c[class_of(W, t)] += 1
    return c


def add_counts(a, b):
    return {k: int(a.get(k, 0)) + int(b.get(k, 0)) for k in CLASSES}


def windows_after(W, ast, t_s):
    """Number of windows of `ast` strictly after mission time t_s [SECONDS]."""
    return int(np.sum(W[int(ast)]['t'] > float(t_s)))


def edf(W, ast, t_after_s, beta):
    """Earliest-deadline bonus of taking `ast` in a slice that commits up to t_after_s (= T + H): w = windows after
    t_after_s; beta if w = 0 (now or never), beta / 2 if w = 1, beta / 4 if w = 2, else 0."""
    w = windows_after(W, ast, t_after_s)
    return float(beta) if w == 0 else (float(beta) / 2 if w == 1 else (float(beta) / 4 if w == 2 else 0.0))


# ------------------------------------------------------------------ twins
def load_twin(path):
    """npz (run_ialns ist format) -> state dict with tL as float and asts as ints."""
    z = np.load(path); st = {k: z[k] for k in z.files}
    st['tL'] = float(st['tL']); st['asts'] = np.asarray(st['asts'], int)
    return st


def save_twin(path, st):
    """Atomic npz write of an ist state dict."""
    path = pathlib.Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f'.tmp{os.getpid()}.npz')
    np.savez(tmp, **{k: np.asarray(v) for k, v in st.items()})
    os.replace(tmp, path)


def twin_of(st):
    """ImpulsiveProblem of an ist state dict (run_ialns.ipr)."""
    return RI.ipr(st)


def miss_of(ip):
    """Max flyby miss [km] of a twin as it stands (one integration, no optimisation)."""
    Yf, _ = ip.integrate(); dm, _ = ip.misses(Yf)
    return float(np.linalg.norm(dm, axis=1).max()) if len(dm) else 0.0


def order_of(st):
    """Flyby order (indices into st['tf'] / st['asts'] sorted by epoch)."""
    return np.argsort(np.asarray(st['tf'], float), kind='stable')


def twin_summary(st, ip=None):
    """Cheap per-route facts: n, targets (flyby order), epochs (s), t_launch, t_end, tank (as stored, NOT necessarily
    honest), miss (km), dv (km/s)."""
    ip = ip or twin_of(st); o = order_of(st)
    return dict(n=int(len(st['asts'])), targets=[int(st['asts'][i]) for i in o], epochs=[float(st['tf'][i]) for i in o],
                t_launch=float(st['tL']), t_end=float(np.max(st['tf'])) if len(st['tf']) else float(st['tL']),
                tank=float(ip.tank()), miss=miss_of(ip), dv=float(ip.dv()))


def nreg_of(ip, tol=1.0):
    """Number of regular node centres tL + (k + 1/2) dtau, k = 0.., present in ip.ts (s14_twinbeam's `nreg`: the
    prefix twin's regular nodes; child_problem appends centres from nreg on)."""
    ts = np.asarray(ip.ts, float)
    k = (ts - ip.tL) / ip.h - 0.5
    reg = np.abs(k - np.round(k)) * ip.h <= tol
    have = set(np.round(k[reg]).astype(int).tolist())
    n = 0
    while n in have:
        n += 1
    return int(n)


def pack_state(ip, miss=None, how='commit'):
    """s14_twinbeam twin STATE dict (st, nreg, tank, dv, miss, t_end, r_end, v_end, asts, tfs, how) of a settled
    ImpulsiveProblem: the object the segment beam expands (s14_twinbeam.search_state / child_problem / settle_child)."""
    import s14_twinbeam as TB
    if miss is None:
        miss = miss_of(ip)
    return TB.pack(ip, nreg_of(ip), float(miss), how)


def honest_check(ip, tol_s=1.0):
    """Honesty of a twin AS IT STANDS (no optimisation): {n, irregular, max_imp_cap, max_window_cap, miss_km, tank, dv,
    ok}.  irregular = nodes off the regular grid tL + (k + 1/2) dtau; max_window_cap = the largest sum of impulses inside
    one regular dtau window divided by tcap (the continuous-thrust capability at the tank mass); ok = irregular == 0 and
    max_window_cap <= CAP_TOL and miss <= TOL_MISS.  Same columns as results/s16/best/honest_check.txt."""
    ts = np.asarray(ip.ts, float); Ts = np.asarray(ip.Ts, float).reshape(-1, 3)
    k = (ts - ip.tL) / ip.h - 0.5
    irregular = int(np.sum(np.abs(k - np.round(k)) * ip.h > tol_s))
    imp = np.linalg.norm(Ts, axis=1) if len(Ts) else np.zeros(0)
    tcap = float(ip.tcap)
    if len(ts):
        edges = np.arange(ip.tL, max(ts.max(), float(np.max(ip.tf)) if len(ip.tf) else ts.max()) + 2 * ip.h, ip.h)
        b = np.clip(np.searchsorted(edges, ts, side='right') - 1, 0, len(edges) - 2)
        win = np.zeros(len(edges) - 1); np.add.at(win, b, imp)
        max_window = float(win.max() / tcap); max_imp = float(imp.max() / tcap)
    else:
        max_window = max_imp = 0.0
    miss = miss_of(ip)
    return dict(n=int(len(ip.asts)), irregular=irregular, max_imp_cap=round(max_imp, 4), max_window_cap=round(max_window, 4),
                miss_km=round(miss, 1), tank=float(ip.tank()), dv=float(ip.dv()),
                ok=bool(irregular == 0 and max_window <= CAP_TOL and miss <= TOL_MISS))


def regrid_settle(ip, dtau=DTAU, iters=None, timeout=None):
    """THE REGRID LAW.  s15b_regrid.regrid (every impulse summed into its regular dtau bin) + run_ialns.settle (restore,
    L1 optimise, restore, thrust-cap homotopy) on a COPY of the twin.  Returns
        dict(ip=<regridded ImpulsiveProblem>, st=<ist dict>, tank=<honest tank kg>, miss=<km>, ok=<miss <= 150 and
             honest ok>, tank_in=<tank before>, over_in=<max impulse / cap before>, over_out=<after>, honest=<honest_check>,
             secs=<wall>, error=<None | str>)
    ok False (miss > 150 km, timeout, exception) means the route is NOT usable at this tank; the caller keeps its previous
    state.  A twin that is already regular is still re-settled (cheap; idempotent to ~0.1 kg)."""
    iters = int(DEF['regrid_iters'] if iters is None else iters)
    timeout = float(DEF['regrid_timeout'] if timeout is None else timeout)
    tic = time.time(); tank_in = float(ip.tank())
    over_in = float((np.linalg.norm(ip.Ts, axis=1) / ip.tcap).max()) if len(ip.Ts) else 0.0
    g = regrid(ip, dtau); err = None; miss = np.inf
    try:
        miss = float(with_timeout(RI.settle, timeout, g, iters))
    except Timeout:
        err = 'timeout'
    except Exception as e:
        err = repr(e)[:120]
    over_out = float((np.linalg.norm(g.Ts, axis=1) / g.tcap).max()) if len(g.Ts) else 0.0
    hc = honest_check(g) if err is None else dict(ok=False, irregular=-1, max_window_cap=over_out, miss_km=miss)
    ok = bool(err is None and miss <= TOL_MISS and hc['ok'])
    return dict(ip=g, st=RI.ist(g), tank=float(g.tank()), miss=float(miss), ok=ok, tank_in=tank_in, over_in=over_in,
                over_out=over_out, honest=hc, secs=round(time.time() - tic, 1), error=err)


def truncate_twin(st, t_cut, iters=100, timeout=300.0):
    """Prefix of a route: the flybys with epoch <= t_cut [s] and the impulse nodes before the last of them, re-settled
    (s15b_prefix.make_prefix).  Returns (regrid_settle result dict, n_kept).  Used by the WP1 acceptance test and the WP3
    replay fixture to make realistic launched-craft prefixes out of existing routes (results/s16/best/fleet)."""
    ip = twin_of(st); o = order_of(st); keep = [i for i in o if float(st['tf'][i]) <= t_cut]
    if not keep:
        raise ValueError('no flyby before t_cut')
    tk = float(max(st['tf'][i] for i in keep)); m = ip.ts <= tk
    ip2 = ImpulsiveProblem(RI.eph(), ip.tL, ip.vinf, ip.ts[m], ip.Ts[m], np.asarray(st['tf'], float)[keep],
                           [int(st['asts'][i]) for i in keep])
    return regrid_settle(ip2, iters=iters, timeout=timeout), len(keep)


# ------------------------------------------------------------------ schemas
COLUMN_SCHEMA = """COLUMN (one proposal of one craft for one slice; JSON-serialisable dict; list in slices/s<k>/columns.json)
  id          str   unique in the run: "<job_id>_k<n>" e.g. "s03_c2_k5"
  craft       int | None   0-based craft index; None = ROOT column (a launch, interchangeable among unlaunched craft)
  kind        "craft" | "root"
  targets     [int] NEW targets of this column in flyby order (1-based; not in the craft's prefix)
  epochs      [float] their flyby epochs, SECONDS; every one in (T, T + H]  (the committed part)
  all_targets [int] the whole route after adoption (prefix + new), flyby order
  n_new       int   len(targets)
  npz         str   path (relative to the run dir) of the settled twin of the WHOLE route (ist format, planner-born:
                    irregular nodes allowed; commit regrids it).  A fixture "fake" column has npz = None.
  tank        float planner-twin tank [kg] of the whole route (NOT honest)
  tank_prev   float honest tank of the prefix before this column (root: 601.5)
  dtank       float tank - tank_prev   (ranking only)
  miss        float max flyby miss [km] of the twin (<= 150)
  t_launch    float s (root columns: the launch epoch of the proposal; craft columns: the prefix's)
  t_end       float s, last committed flyby
  potential   int   max number of extra flybys any beam descendant of this state reached in (T + H, T + H + L]
  counts      {"pin","med","fil"} classes of the NEW targets
  depth       int   beam level of the state (root = 1)
  score       float the beam's own score (for diagnostics)
  job_id      str
"""

STATE_SCHEMA = """STATE (results/s18/<run>/state.json; the resumable truth of a run; written atomically after every phase)
  run         str   run name;  created / updated  timestamps;  params  dict (DEF keys; H, L, grid_step, root_bin in DAYS as
                    given on the CLI, converted with d2s where used)
  N           int   number of craft
  T           float s, the current commit epoch (start of the current slice)
  slice       int   index of the current slice (0-based);  phase  "propose" | "auction" | "commit" | "done"
  craft       [ {i, launched (bool), route ("craft_<i>.npz" | None), targets [int] (flyby order), epochs [float s],
                 counts {"pin","med","fil"}, tank (HONEST, kg; 601.5 when unlaunched), t_launch (s | None), t_end (s | None),
                 history [ {slice, column_id, n_new, tank, dtank, secs} ] } ]  for i = 0..N-1
  remaining   [int] targets not yet taken (start: TARGETS; never contains 131/144)
  taken       {"<target>": craft index}
  log         [ per slice: {slice, T_d, n_jobs, n_columns, n_kept, chosen {"<craft>": column id}, rejected [..],
                 launched, covered, fuel, sumJi, wall_s} ]
  done        bool;  stop  {kind, T_d} when done
Per-slice files under slices/s<kk>/: jobs/<job_id>.json (JOB), <job>/result.json (RESULT), columns.json (all COLUMNs),
auction.json (the choice), commit.json (regrid_settle outcomes).  Resume = load state.json, then skip every phase of the
current slice whose file exists.
"""


def new_state(run, params, N):
    return dict(run=str(run), created=now(), updated=now(), params=dict(params), N=int(N), T=0.0, slice=0,
                phase='propose',
                craft=[dict(i=i, launched=False, route=None, targets=[], epochs=[], counts=dict(pin=0, med=0, fil=0),
                            tank=TANK_DRY, t_launch=None, t_end=None, history=[]) for i in range(int(N))],
                remaining=list(TARGETS), taken={}, log=[], done=False, stop=None)


def check_state(S):
    """Invariants of a STATE dict; raises AssertionError with the broken one."""
    assert S['N'] == len(S['craft'])
    rem = set(S['remaining']); tk = {int(k): v for k, v in S['taken'].items()}
    assert not (rem & set(UNREACHABLE)), 'unreachable in remaining'
    assert not (rem & set(tk)), 'taken target still remaining'
    assert rem | set(tk) == set(TARGETS), 'remaining + taken != TARGETS'
    for c in S['craft']:
        assert len(c['targets']) == len(set(c['targets'])), f'craft {c["i"]}: duplicate target'
        assert all(tk.get(t) == c['i'] for t in c['targets']), f'craft {c["i"]}: taken map disagrees'
        assert len(c['targets']) == len(c['epochs'])
        assert all(a < b for a, b in zip(c['epochs'], c['epochs'][1:])), f'craft {c["i"]}: epochs not increasing'
        assert c['launched'] == (c['route'] is not None)
        assert sum(c['counts'].values()) == len(c['targets'])
    assert sum(len(c['targets']) for c in S['craft']) == len(tk)
    return True


def save_state(S, rd):
    S['updated'] = now(); check_state(S)
    jdump(S, pathlib.Path(rd) / 'state.json')


def load_state(rd):
    S = jload(pathlib.Path(rd) / 'state.json'); check_state(S)
    return S


def new_column(**kw):
    """COLUMN dict with every key present (see COLUMN_SCHEMA); the caller fills what it knows."""
    c = dict(id=None, craft=None, kind='craft', targets=[], epochs=[], all_targets=[], n_new=0, npz=None, tank=TANK_DRY,
             tank_prev=TANK_DRY, dtank=0.0, miss=0.0, t_launch=None, t_end=None, potential=0,
             counts=dict(pin=0, med=0, fil=0), depth=1, score=0.0, job_id=None)
    c.update(kw); c['n_new'] = len(c['targets']); c['dtank'] = float(c['tank']) - float(c['tank_prev'])
    c['kind'] = 'root' if c['craft'] is None else 'craft'
    return c


def check_column(c, T=None, H=None):
    """Schema and horizon invariants of a COLUMN (epochs in (T, T + H] when given)."""
    for k in ('id', 'targets', 'epochs', 'all_targets', 'tank', 'tank_prev', 'miss', 'potential', 'counts'):
        assert k in c, f'column missing {k}'
    assert len(c['targets']) == len(c['epochs']) == c['n_new']
    assert len(set(c['targets'])) == len(c['targets']), 'duplicate target in column'
    assert not (set(c['targets']) & set(UNREACHABLE))
    assert set(c['targets']) <= set(c['all_targets'])
    assert all(a < b for a, b in zip(c['epochs'], c['epochs'][1:])), 'epochs not increasing'
    assert c['miss'] <= TOL_MISS, f'column miss {c["miss"]} > {TOL_MISS}'
    assert sum(c['counts'].values()) == c['n_new']
    if T is not None and H is not None and c['epochs']:
        assert min(c['epochs']) > T - 1e-6 and max(c['epochs']) <= T + H + 1e-6, 'column epoch outside (T, T+H]'
    # fixer 09-24 (reviewer minor finding): the launch precedes the first flyby; a root launch lies in [T, T + H]
    # (T itself only for the mission-start launch of slice 0)
    if c.get('t_launch') is not None and c['epochs']:
        assert float(c['t_launch']) < min(c['epochs']), 'launch after the first flyby'
        if c.get('craft') is None and T is not None and H is not None:
            assert T - 1e-6 <= float(c['t_launch']) <= T + H + 1e-6, 'root launch outside [T, T+H]'
    return True


if __name__ == '__main__':
    # self-test: build the class table and print the histogram
    W = load_windows()
    from collections import Counter
    print('windows:', Counter(class_of(W, t) for t in TARGETS), 'fallback', W['_meta'].get('n_fallback'))
    print('edf(1, T+H = 540 d, beta 1.5) =', edf(W, 1, d2s(540), 1.5), '; windows after 5000 d of target 2:', windows_after(W, 2, d2s(5000)))
