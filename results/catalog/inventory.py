"""CATALOG inventory: parse every submission file and every twin fleet directory under results/.
Writes results/catalog/inventory_raw.json (measured numbers only)."""
import os, sys, json, glob, hashlib, time, pathlib, re
import numpy as np
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
sys.path.insert(0, str(ROOT))
from ctoc14.constants import M_DRY, FUEL_MAX, VE, DAY

def cost(m0):
    x = (float(m0) - M_DRY) / FUEL_MAX
    return 1.0 + x + x * x

def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

def parse_validator(path):
    if not path or not os.path.exists(path):
        return None
    txt = open(path, errors='replace').read()
    out = dict(log=str(pathlib.Path(path).relative_to(ROOT)))
    m = re.search(r'VERDICT:\s*(\w+)', txt)
    out['verdict'] = m.group(1) if m else None
    m = re.search(r'N_covered\s*=\s*(\d+)\s+N_miss\s*=\s*(\d+)\s+sum J_i\s*=\s*([\d.]+)\s+J\s*=\s*([\d.]+)', txt)
    if m:
        out['covered'] = int(m.group(1)); out['miss'] = int(m.group(2))
        out['sumJi'] = float(m.group(3)); out['J'] = float(m.group(4))
    return out

# ---------------------------------------------------------------- submissions
def parse_submission(path):
    craft = {}
    order = []
    with open(path, errors='replace') as f:
        for line in f:
            s = line.strip()
            if not s or s[0] == '#':
                continue
            p = s.split()
            if len(p) < 15:
                continue
            sc = int(p[1]); ev = int(p[2])
            c = craft.get(sc)
            if c is None:
                c = craft[sc] = dict(id=sc, m0=None, t0=None, t_end=None, n_ev1=0, n_ev2=0, flybys=[], m_min=1e9, rows=0)
                order.append(sc)
            c['rows'] += 1
            t = float(p[3]); m = float(p[10])
            if m < c['m_min']:
                c['m_min'] = m
            if ev == 0:
                c['m0'] = m; c['t0'] = t
            elif ev == 1:
                c['n_ev1'] += 1
            elif ev == 2:
                c['n_ev2'] += 1
            elif ev == 3:
                c['flybys'].append((t, int(p[14])))
            elif ev == 4:
                c['t_end'] = t
    seen = set(); covered = []
    for sc in order:
        for t, a in craft[sc]['flybys']:
            if a not in seen:
                seen.add(a); covered.append(a)
    per = []
    for sc in order:
        c = craft[sc]
        per.append(dict(id=sc, m0=c['m0'], J_i=cost(c['m0']), n_flybys=len(c['flybys']),
                        n_unique=len(set(a for _, a in c['flybys'])),
                        launch_d=c['t0'] / DAY if c['t0'] is not None else None,
                        end_d=c['t_end'] / DAY if c['t_end'] is not None else None,
                        fuel_kg=c['m0'] - c['m_min'] if c['m0'] is not None else None,
                        n_thrust_rows=c['n_ev1'], n_coast_rows=c['n_ev2'], rows=c['rows']))
    sumJi = sum(x['J_i'] for x in per)
    missed = sorted(set(range(1, 301)) - seen)
    return dict(n_craft=len(order), n_flybys=sum(len(craft[s]['flybys']) for s in order), covered=len(seen),
                missed=missed, sumJi=sumJi, J_raw=sumJi + 300 - len(seen), craft=per,
                m0_total=sum(x['m0'] for x in per), fuel_total=sum(x['fuel_kg'] for x in per))

subs = []
files = sorted(glob.glob(str(ROOT / 'results/CTOC14_Result_*.txt')))
extra = pathlib.Path.home() / '.Trash/CTOC14_Result_Miki_s16E.txt'
if extra.exists():
    files.append(str(extra))
for f in files:
    t0 = time.time()
    tag = re.sub(r'^CTOC14_Result_', '', pathlib.Path(f).stem)
    st = os.stat(f)
    rec = dict(tag=tag, file=str(pathlib.Path(f).relative_to(ROOT)) if str(f).startswith(str(ROOT)) else f,
               bytes=st.st_size, mtime=time.strftime('%Y-%m-%d %H:%M', time.localtime(st.st_mtime)),
               md5=md5(f))
    rec.update(parse_submission(f))
    vlog = ROOT / f'results/validator2_{tag}.log'
    if not vlog.exists():
        vlog = ROOT / f'results/catalog/validator2_{tag}.log'
    rec['validator'] = parse_validator(str(vlog) if vlog.exists() else None)
    subs.append(rec)
    print(f"{tag:14s} craft={rec['n_craft']:2d} fb={rec['n_flybys']:4d} cov={rec['covered']:3d} sumJi={rec['sumJi']:.6f} "
          f"J={rec['J_raw']:.6f} md5={rec['md5'][:8]} {rec['mtime']} val={rec['validator'] and rec['validator'].get('verdict')} ({time.time()-t0:.1f}s)", flush=True)

# selection-only tags (never combined into a txt)
sel_only = []
for f in sorted(glob.glob(str(ROOT / 'results/CTOC14_Result_*.txt.selection'))):
    tag = re.sub(r'^CTOC14_Result_', '', pathlib.Path(f).name.replace('.txt.selection', ''))
    if not (ROOT / f'results/CTOC14_Result_{tag}.txt').exists():
        st = os.stat(f)
        sel_only.append(dict(tag=tag, file=str(pathlib.Path(f).relative_to(ROOT)), bytes=st.st_size,
                             mtime=time.strftime('%Y-%m-%d %H:%M', time.localtime(st.st_mtime)),
                             content=open(f).read().strip()[:400]))

# ---------------------------------------------------------------- twin fleets
def tank_of(st):
    Ts = np.asarray(st['Ts'], float).reshape(-1, 3)
    dv = float(np.linalg.norm(Ts, axis=1).sum()) if len(Ts) else 0.0
    return (M_DRY + 1.5) * np.exp(dv / VE), dv

def read_routes(d):
    routes = {}
    for f in sorted(pathlib.Path(d).glob('route_*.npz')):
        z = np.load(f, allow_pickle=True)
        if not all(k in z.files for k in ('tL', 'Ts', 'tf', 'asts')):
            continue
        st = {k: z[k] for k in z.files}
        tank, dv = tank_of(st)
        tf = np.asarray(st['tf'], float).ravel()
        routes[f.stem[6:]] = dict(tank=float(tank), J_i=cost(tank), dv_kms=dv, flybys=int(len(tf)),
                                  asts=[int(a) for a in np.asarray(st['asts']).ravel()],
                                  launch_d=float(st['tL']) / DAY,
                                  last_flyby_d=float(tf.max()) / DAY if len(tf) else None)
    return routes

twins = []
seen_dirs = set()
all_json = sorted(glob.glob(str(ROOT / 'results/**/fleet.json'), recursive=True))
for jf in all_json:
    d = pathlib.Path(jf).parent
    seen_dirs.add(str(d))
    try:
        js = json.load(open(jf))
    except Exception as e:
        twins.append(dict(dir=str(d.relative_to(ROOT)), error=repr(e))); continue
    st = os.stat(jf)
    rec = dict(dir=str(d.relative_to(ROOT)), stage=d.relative_to(ROOT).parts[1],
               mtime=time.strftime('%Y-%m-%d %H:%M', time.localtime(st.st_mtime)),
               has_fleet_json=True, note=js.get('note', ''))
    rr = js.get('routes') or js.get('craft') or {}
    fmt = 'routes' if 'routes' in js else ('craft' if 'craft' in js else 'none')
    rec['format'] = fmt
    rec['n'] = js.get('n', len(rr))
    rec['covered'] = js.get('covered')
    rec['J'] = js.get('J')
    rec['honest'] = js.get('honest')
    m0s = []
    for k, v in rr.items():
        m0s.append(v.get('tank', v.get('m0')))
    rec['m0'] = [float(x) for x in m0s if x is not None]
    rec['sumJi'] = js.get('sumJi', sum(cost(x) for x in rec['m0']) if rec['m0'] else None)
    rec['fuel'] = js.get('fuel', sum(rec['m0']) - M_DRY * len(rec['m0']) if rec['m0'] else None)
    rec['flybys'] = int(sum(int(v.get('flybys', 0)) for v in rr.values())) if rr else None
    if 'misses' in js:
        rec['misses'] = js['misses']
    npz_route = sorted(d.glob('route_*.npz')); npz_craft = sorted(d.glob('craft_*.npz'))
    rec['n_route_npz'] = len(npz_route); rec['n_craft_npz'] = len(npz_craft)
    if rec['covered'] is None and npz_route:
        routes = read_routes(d)
        cov = set(a for r in routes.values() for a in r['asts'])
        rec['covered'] = len(cov); rec['n'] = len(routes)
        rec['m0'] = [r['tank'] for r in routes.values()]; rec['sumJi'] = sum(r['J_i'] for r in routes.values())
        rec['J'] = rec['sumJi'] + 300 - len(cov); rec['flybys'] = sum(r['flybys'] for r in routes.values())
        rec['computed_from_npz'] = True
    twins.append(rec)

# route_*.npz directories without fleet.json (route libraries and unsaved fleets)
npz_dirs = sorted(set(str(pathlib.Path(f).parent) for f in glob.glob(str(ROOT / 'results/**/route_*.npz'), recursive=True)))
libs = []
for d in npz_dirs:
    if d in seen_dirs:
        continue
    dp = pathlib.Path(d)
    routes = read_routes(d)
    if not routes:
        continue
    cov = set(a for r in routes.values() for a in r['asts'])
    newest = max(os.stat(f).st_mtime for f in dp.glob('route_*.npz'))
    sumJi = sum(r['J_i'] for r in routes.values())
    libs.append(dict(dir=str(dp.relative_to(ROOT)), stage=dp.relative_to(ROOT).parts[1], has_fleet_json=False,
                     mtime=time.strftime('%Y-%m-%d %H:%M', time.localtime(newest)), n=len(routes), covered=len(cov),
                     flybys=sum(r['flybys'] for r in routes.values()), sumJi=sumJi, J=sumJi + 300 - len(cov),
                     m0=[r['tank'] for r in routes.values()], computed_from_npz=True,
                     kind='library' if dp.name == 'lib' or len(routes) > 12 else 'fleet_no_json'))

out = dict(generated=time.strftime('%Y-%m-%d %H:%M:%S'), submissions=subs, selection_only=sel_only, twin_fleets=twins,
           npz_dirs_without_fleet_json=libs)
json.dump(out, open(ROOT / 'results/catalog/inventory_raw.json', 'w'), indent=1)
print('submissions', len(subs), 'selection-only', len(sel_only), 'fleet.json', len(twins), 'npz dirs w/o json', len(libs))
