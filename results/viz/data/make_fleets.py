#!/usr/bin/env python
"""Write results/viz/data/fleets/<key>.json and results/viz/data/index.json (docs/viz_spec.md section 2) for every
fleet listed in results/catalog/viz_selection.json.

Run from the project root:  nice -n 5 ~/.venvs/astro313/bin/python results/viz/data/make_fleets.py [key ...]
Reads results/CTOC14_Result_*.txt (submission files) and results/s18/*/fleet/route_*.npz (impulsive twins);
writes ONLY under results/viz/data/fleets/ and results/viz/data/index.json.

Conversion rules (see docs/viz_fleets.md):
- frame: heliocentric ecliptic J2000; positions AU (6 significant digits), times days since MJD 62502.0
  (integer-day grid samples are written as integers, flyby epochs to 1e-6 d), mass kg (2 decimals).
- submission files: every row is a sample (Event 0 launch, 1 thrust sample at 1 d, 3 flyby, 4 end). Kept: launch, end,
  every flyby row (its OWN position, so the viewer's flyby sample is what the validator checked), the GRID_D-day grid and
  every 1 d row within NEAR_FB_D of a flyby. Thrust arcs = contiguous runs of Event=1 rows with |T| >= THRUST_ON_N
  (all rows are Event=1, coasts are rows with ~0 N, so a magnitude threshold is needed); arc = [t_first, t_next_after_last].
  m0, mass samples and fuel (m0 - m_end) come from the file's mass column.
- twin fleets (route_*.npz, run_ialns "ist" format): ImpulsiveProblem.integrate() gives the state at every impulse node
  (just before the impulse) and at every flyby; each Kepler arc from a node's post-impulse state is resampled with
  ctoc14.kepler.propagate_twobody on the grid; flyby samples are the integrate() states at ip.tf. Thrust arcs: impulse nodes
  with |dv| >= DV_MIN, half-width max(1 d, 10 d * |dv| / dv_cap) with dv_cap = 0.43 N * 20 d / m (the bin's continuous-thrust
  capability, the same quantity ImpulsiveProblem.tcap caps), so an impulse at the cap shows as the whole 20 d bin.
  m0 = ip.tank(), mass drops by the rocket equation at every node, fuel_kg = ip.tank() - 600.
- size: if a fleet file exceeds SIZE_LIMIT the grid step is raised (2 -> 3 -> 4 -> 5 d); flyby samples are never dropped.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]                 # .../ctoc14
sys.path.insert(0, str(ROOT))
from ctoc14.constants import AU, DAY, M_DRY, T0_MJD, VE, cost_sc          # noqa: E402
from ctoc14.kepler import Ephemeris, propagate_twobody                    # noqa: E402
from ctoc14.impulsive import ImpulsiveProblem                             # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / 'fleets'
SEL = ROOT / 'results' / 'catalog' / 'viz_selection.json'

GRID_STEPS_D = (2, 3, 4, 5)      # sample step candidates (days), first that fits SIZE_LIMIT wins
NEAR_FB_D = 4.0                  # 1 d samples within +- this of a flyby epoch
THRUST_ON_N = 0.05               # submission files: thrust "on" when |T| >= 10 % of Tmax (0.5 N)
DV_MIN_KMS = 1e-3                # twins: impulses below 1 m/s are numerical noise, not thrust arcs
BIN_D = 20.0                     # twins: impulse-node spacing (ImpulsiveProblem dtau)
T_CONT_N = 0.43                  # twins: continuous-thrust capability used by ImpulsiveProblem.tcap
SIZE_LIMIT = 1_500_000           # bytes per fleet file

_EPH = None
def eph():
    global _EPH
    if _EPH is None:
        _EPH = Ephemeris(ROOT / 'MEA.txt')
    return _EPH

# ------------------------------------------------------------------------------------------------ number formatting
def f6(x):
    """6 significant digits as a float (json prints the shortest repr, e.g. 0.987654 or 1.2e-05)."""
    return float(f'{x:.6g}')

def tday(t, flyby=False):
    """Time in days: flyby epochs to 1e-6 d; grid epochs to 1e-4 d, integers written as int."""
    if flyby:
        return round(float(t), 6)
    r = round(float(t), 4)
    return int(r) if r == int(r) else r

def merge_arcs(arcs):
    arcs = sorted((float(a), float(b)) for a, b in arcs if b > a)
    out = []
    for a, b in arcs:
        if out and a <= out[-1][1] + 1e-9:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [[tday(a), tday(b)] for a, b in out]

def thin_mask(t, keep, fb_t, step):
    """Greedy thinning: keep[k] rows always; other rows when >= step since the last kept row, or within NEAR_FB_D of a
    flyby and >= 1 d since the last kept row. t in days (sorted)."""
    fb_t = np.asarray(fb_t, float)
    n = len(t); mask = np.zeros(n, bool); last = -np.inf
    before_gap = np.zeros(n, bool); before_gap[:-1] = np.diff(t) > step + 1e-9      # last row before a sparse stretch
    for k in range(n):
        near = fb_t.size and np.min(np.abs(fb_t - t[k])) <= NEAR_FB_D
        if keep[k] or before_gap[k] or t[k] - last >= step - 1e-9 or (near and t[k] - last >= 1.0 - 1e-9):
            mask[k] = True; last = t[k]
    return mask

def strictly_increasing(t, flyby_flag):
    """Indices to keep so that t is strictly increasing; a flyby row wins over a non-flyby row at the same epoch."""
    keep = []
    for k in range(len(t)):
        if keep and t[k] <= t[keep[-1]] + 1e-12:
            if flyby_flag[k] and not flyby_flag[keep[-1]]:
                keep[-1] = k
            continue
        keep.append(k)
    return np.array(keep, int)

# ------------------------------------------------------------------------------------------------ submission files
def parse_submission(path):
    """Rows of a CTOC14 submission file grouped by spacecraft id (in order of first appearance)."""
    craft = {}
    with open(path) as fh:
        for line in fh:
            s = line.split()
            if len(s) < 15 or not s[0].lstrip('-').isdigit():
                continue
            sc = int(s[1]); ev = int(s[2])
            craft.setdefault(sc, []).append((float(s[3]) / DAY, ev, float(s[4]) / AU, float(s[5]) / AU, float(s[6]) / AU,
                                             float(s[10]), math.sqrt(float(s[11]) ** 2 + float(s[12]) ** 2 + float(s[13]) ** 2),
                                             int(s[14])))
    return craft

def convert_submission_craft(sc, rows, step):
    rows = sorted(rows, key=lambda r: r[0])                       # stable: file order within equal epochs
    t = np.array([r[0] for r in rows]); ev = np.array([r[1] for r in rows])
    xyz = np.array([[r[2], r[3], r[4]] for r in rows]); m = np.array([r[5] for r in rows])
    Tn = np.array([r[6] for r in rows]); ast = np.array([r[7] for r in rows])
    fb = ev == 3
    launch_rows = np.where(ev == 0)[0]; end_rows = np.where(ev == 4)[0]
    launch_d = t[launch_rows[0]] if len(launch_rows) else t[0]
    end_d = t[end_rows[-1]] if len(end_rows) else t[-1]
    m0 = m[launch_rows[0]] if len(launch_rows) else m[0]
    m_end = m[end_rows[-1]] if len(end_rows) else m[-1]
    # thrust arcs from the full-resolution rows: runs of consecutive Event=1 rows with |T| >= THRUST_ON_N, closed by the
    # next Event=1 row (flyby rows sit between two 1 d samples and do not break a run)
    e1 = np.where(ev == 1)[0]
    on = Tn[e1] >= THRUST_ON_N
    arcs = []; k = 0; n = len(rows); n1 = len(e1)
    while k < n1:
        if on[k]:
            j = k
            while j + 1 < n1 and on[j + 1]:
                j += 1
            t_off = t[e1[j + 1]] if j + 1 < n1 else end_d
            arcs.append((t[e1[k]], min(t_off, end_d)))
            k = j + 1
        else:
            k += 1
    # thinning
    keep = (ev == 0) | (ev == 3) | (ev == 4)
    mask = thin_mask(t, keep, t[fb], step)
    idx = np.where(mask)[0]
    idx = idx[strictly_increasing(t[idx], fb[idx])]
    tt = [tday(t[i], flyby=bool(fb[i])) for i in idx]
    flybys = [{'t': tday(t[i], flyby=True), 'ast': int(ast[i])} for i in np.where(fb)[0]]
    return {'id': int(sc), 'm0': round(float(m0), 6), 'launch_d': tday(launch_d), 'end_d': tday(end_d),
            'n_flybys': int(fb.sum()), 'fuel_kg': round(float(m0 - m_end), 3),
            't': tt, 'xyz': [[f6(v) for v in xyz[i]] for i in idx], 'm': [round(float(m[i]), 2) for i in idx],
            'thrust': merge_arcs(arcs), 'flybys': flybys,
            'n_rows_file': int(n), 'thrust_on_days': float(round(on.sum(), 1))}

def build_submission(sel, step):
    craft_rows = parse_submission(ROOT / sel['source'])
    craft = [convert_submission_craft(sc, rows, step) for sc, rows in craft_rows.items()]
    return craft

# ------------------------------------------------------------------------------------------------ twin fleets
def load_route(path):
    z = np.load(path); st = {k: z[k] for k in z.files}
    return ImpulsiveProblem(eph(), float(st['tL']), st['vinf'], st['ts'], st['Ts'], st['tf'], [int(a) for a in st['asts']])

def convert_twin_craft(cid, route_name, ip, step):
    Yf, Ys = ip.integrate()
    tank = float(ip.tank())
    t_end = float(ip.tf.max())
    # ordered node list on [tL, t_end]: (t [s], r [km], v_after [km/s], kind, dv, index)
    nodes = [(ip.tL, ip.rE.copy(), ip.vE + ip.vinf, 'L', 0.0, -1)]
    for k in range(len(ip.ts)):
        if ip.ts[k] <= t_end:
            nodes.append((float(ip.ts[k]), Ys[k, :3].copy(), Ys[k, 3:6] + ip.Ts[k], 'I', float(np.linalg.norm(ip.Ts[k])), k))
    for j in range(len(ip.tf)):
        nodes.append((float(ip.tf[j]), Yf[j, :3].copy(), Yf[j, 3:6].copy(), 'F', 0.0, j))
    nodes.sort(key=lambda q: (q[0], {'L': 0, 'I': 1, 'F': 2}[q[3]]))
    tL_d = ip.tL / DAY
    fb_d = ip.tf / DAY
    grid = tL_d + step * np.arange(0, int((t_end / DAY - tL_d) / step) + 2)
    fine = np.concatenate([np.arange(np.floor(x - NEAR_FB_D), np.ceil(x + NEAR_FB_D) + 1) for x in fb_d]) if len(fb_d) else np.array([])
    extra_d = np.unique(np.concatenate([grid, fine]))
    extra_d = extra_d[(extra_d > tL_d) & (extra_d < t_end / DAY)]
    node_d = np.array([q[0] / DAY for q in nodes])
    # drop grid epochs that coincide with a node epoch
    extra_d = extra_d[np.min(np.abs(extra_d[:, None] - node_d[None, :]), axis=1) > 1e-6]
    samples = []                                             # (t_d, xyz_AU, m, flyby_ast or None)
    m = tank; impulses = []; arcs = []
    max_arc_err_km = 0.0
    for q, (tn, r, v, kind, dv, idx) in enumerate(nodes):
        if kind == 'I' and dv > 0:
            m *= math.exp(-dv / VE)
            if dv >= DV_MIN_KMS:
                dv_cap = T_CONT_N * BIN_D * DAY / m / 1000.0            # km/s deliverable in one 20 d bin at mass m
                w = min(max(1.0, float(step), 0.5 * BIN_D * dv / dv_cap), 0.5 * BIN_D)   # >= grid step: the viewer needs a sample pair inside
                arcs.append((max(tn / DAY - w, tL_d), min(tn / DAY + w, t_end / DAY)))
                impulses.append([tday(tn / DAY), round(dv, 4)])
        samples.append((tn / DAY, r / AU, m, ip.asts[idx] if kind == 'F' else None, kind))
        if q + 1 < len(nodes):
            t_next = nodes[q + 1][0]
            sel = extra_d[(extra_d > tn / DAY) & (extra_d < t_next / DAY)]
            if len(sel):
                rr, _ = propagate_twobody(r, v, sel * DAY - tn)
                for td, rk in zip(sel, rr):
                    samples.append((float(td), rk / AU, m, None, 'G'))
            # consistency: the arc propagated to the next node must land on integrate()'s state there
            if t_next > tn:
                r_chk, _ = propagate_twobody(r, v, t_next - tn)
                max_arc_err_km = max(max_arc_err_km, float(np.linalg.norm(r_chk - nodes[q + 1][1])))
    samples.sort(key=lambda s: (s[0], 0 if s[3] is not None else 1))
    t = np.array([s[0] for s in samples]); fbf = np.array([s[3] is not None for s in samples])
    idx = strictly_increasing(t, fbf)
    samples = [samples[i] for i in idx]
    flybys = [{'t': tday(s[0], flyby=True), 'ast': int(s[3])} for s in samples if s[3] is not None]
    assert len(flybys) == len(ip.tf), (route_name, len(flybys), len(ip.tf))
    return {'id': int(cid), 'route': route_name, 'm0': round(tank, 6), 'launch_d': tday(tL_d), 'end_d': tday(t_end / DAY, flyby=True),
            'n_flybys': int(len(ip.tf)), 'fuel_kg': round(tank - M_DRY, 3), 'vinf_kms': round(float(np.linalg.norm(ip.vinf)), 4),
            'dv_kms': round(float(ip.dv()), 4),
            't': [tday(s[0], flyby=s[3] is not None) for s in samples],
            'xyz': [[f6(v) for v in s[1]] for s in samples], 'm': [round(float(s[2]), 2) for s in samples],
            'thrust': merge_arcs(arcs), 'impulses': impulses, 'flybys': flybys,
            'max_arc_err_km': round(max_arc_err_km, 6)}

def build_twin(sel, step):
    d = ROOT / sel['source']
    routes = sorted(d.glob('route_*.npz'), key=lambda p: (len(p.stem), p.stem))
    craft = []
    for k, p in enumerate(routes):
        craft.append(convert_twin_craft(k + 1, p.stem[6:], load_route(p), step))
    return craft

# ------------------------------------------------------------------------------------------------ fleet assembly
def assemble(sel, craft, step):
    targets = sorted({f['ast'] for c in craft for f in c['flybys']})
    covered = len(targets)
    missed = [k for k in range(1, 301) if k not in set(targets)]
    sumJi = float(sum(cost_sc(c['m0']) for c in craft))
    J_raw = sumJi + 300 - covered
    fleet = {'key': sel['key'], 'title': sel['title'], 'source': sel['source'], 'kind': sel['kind'],
             'n_craft': len(craft), 'covered': covered, 'missed': missed,
             'J_raw': round(J_raw, 6), 'sumJi': round(sumJi, 6), 'date': sel['date'],
             'n_flybys': int(sum(c['n_flybys'] for c in craft)),
             'fuel_kg': round(float(sum(c['fuel_kg'] for c in craft)), 1),
             'submitted': sel.get('submitted'), 'shown': sel.get('shown'), 'md5': sel.get('md5'),
             'note': sel.get('note', ''),
             'frame': 'heliocentric ecliptic J2000', 'units': {'length': 'AU', 'time': 'day since MJD 62502.0', 'mass': 'kg'},
             'grid_d': step, 'near_flyby_1d_window_d': NEAR_FB_D,
             'thrust_def': (f'contiguous Event=1 runs with |T| >= {THRUST_ON_N} N (10 % of Tmax); arc closed by the next row'
                            if sel['kind'] == 'submission' else
                            f'impulse nodes with |dv| >= {DV_MIN_KMS} km/s, half-width max(1 d, 10 d * |dv| / dv_cap), '
                            f'dv_cap = {T_CONT_N} N * {BIN_D:.0f} d / m (continuous-thrust capability of one bin)'),
             'mass_def': ('file mass column' if sel['kind'] == 'submission'
                          else 'm0 = ip.tank() (601.5 exp(dv/ve)); rocket equation at every impulse; fuel_kg = tank - 600'),
             'craft': craft}
    return fleet

def dumps(obj):
    return json.dumps(obj, separators=(',', ':'), allow_nan=False)

def write_fleet(sel):
    for step in GRID_STEPS_D:
        craft = build_submission(sel, step) if sel['kind'] == 'submission' else build_twin(sel, step)
        fleet = assemble(sel, craft, step)
        txt = dumps(fleet)
        if len(txt) <= SIZE_LIMIT:
            break
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"{sel['key']}.json"
    path.write_text(txt)
    n_s = sum(len(c['t']) for c in craft)
    d = {k: (fleet[k], sel.get(k)) for k in ('n_craft', 'covered', 'J_raw', 'sumJi', 'n_flybys')}
    flags = ' '.join(f'{k}:{a}!={b}' for k, (a, b) in d.items() if b is not None and (abs(a - b) > 1e-5 if isinstance(a, float) else a != b))
    print(f"{sel['key']:>5} {sel['kind']:>10} craft {fleet['n_craft']:>2} covered {fleet['covered']:>3} flybys {fleet['n_flybys']:>3} "
          f"J_raw {fleet['J_raw']:>10.6f} fuel {fleet['fuel_kg']:>8.1f} kg | step {step} d samples {n_s:>6} "
          f"({n_s / len(craft):.0f}/craft) size {len(txt) / 1e6:.3f} MB -> {path.relative_to(ROOT)}"
          + (f'  MISMATCH vs catalog: {flags}' if flags else ''))
    return fleet, path, len(txt), step

def main(argv):
    sel_all = json.loads(SEL.read_text())['fleets']
    keys = argv or [s['key'] for s in sel_all]
    sel_map = {s['key']: s for s in sel_all}
    entries = []
    for key in keys:
        fleet, path, size, step = write_fleet(sel_map[key])
        entries.append({'key': key, 'file': f'fleets/{key}.json', 'title': fleet['title'], 'n_craft': fleet['n_craft'],
                        'covered': fleet['covered'], 'n_flybys': fleet['n_flybys'], 'J_raw': fleet['J_raw'], 'sumJi': fleet['sumJi'],
                        'fuel_kg': fleet['fuel_kg'], 'shown': fleet['shown'], 'submitted': fleet['submitted'], 'date': fleet['date'],
                        'kind': fleet['kind'], 'source': fleet['source'], 'priority': sel_map[key].get('priority'),
                        'grid_d': step, 'size_bytes': size, 'note': fleet['note']})
    # index.json: every selected fleet, ordered by date (then selection order); a partial run keeps the other entries
    idx_path = HERE / 'index.json'
    old = json.loads(idx_path.read_text())['fleets'] if idx_path.exists() else []
    merged = {e['key']: e for e in old if e['key'] in sel_map}
    merged.update({e['key']: e for e in entries})
    order = {s['key']: k for k, s in enumerate(sel_all)}
    fleets = sorted(merged.values(), key=lambda e: (e['date'], order.get(e['key'], 99)))
    idx = {'generated_by': 'results/viz/data/make_fleets.py', 'frame': 'heliocentric ecliptic J2000, AU, days since MJD 62502.0',
           'selection': str(SEL.relative_to(ROOT)), 'fleets': fleets}
    idx_path.write_text(json.dumps(idx, indent=1))
    print(f'wrote {idx_path.relative_to(ROOT)} with {len(fleets)} fleets')

if __name__ == '__main__':
    main(sys.argv[1:])
