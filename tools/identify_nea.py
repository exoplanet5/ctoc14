#!/usr/bin/env python
"""Identify the 300 CTOC14 MEA.txt asteroids against JPL SBDB (catalog search by
orbital-element filters) and confirm with the JPL Horizons API at epoch MJD 61200.0.

Outputs:
  data/sbdb_neo_asteroids.json   raw SBDB NEO asteroid catalog (full precision)
  data/sbdb_neo_comets.json      raw SBDB NEO comet catalog
  data/nea_identification.csv    one row per MEA object: best SBDB match + Horizons check
  data/nea_identification.json   same, machine readable
"""
import json, math, re, sys, time, pathlib, urllib.parse, urllib.request, csv
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'; DATA.mkdir(exist_ok=True)
MEA = np.loadtxt(ROOT / 'MEA.txt')
EPOCH_JD = 2461200.5  # MJD 61200.0
DO_HORIZONS = '--no-horizons' not in sys.argv

def fetch(url, retries=6, timeout=180):
    last = None
    for k in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read().decode('utf-8')
        except Exception as e:
            last = e; time.sleep(3 * (k + 1))
    raise last

FIELDS = ("spkid,full_name,pdes,name,prefix,kind,neo,pha,H,diameter,albedo,orbit_id,epoch,"
          "e,a,q,i,om,w,ma,ad,per_y,moid,class,data_arc,first_obs,last_obs,n_obs_used,condition_code,rms,soln_date")

def sbdb_catalog(kind):
    cache = DATA / f'sbdb_neo_{"asteroids" if kind=="a" else "comets"}.json'
    if cache.exists():
        return json.loads(cache.read_text())
    q = {'fields': FIELDS, 'sb-kind': kind, 'sb-group': 'neo', 'full-prec': 'true'}
    url = 'https://ssd-api.jpl.nasa.gov/sbdb_query.api?' + urllib.parse.urlencode(q)
    txt = fetch(url)
    js = json.loads(txt)
    if 'data' not in js:
        print('SBDB error:', txt[:500]); sys.exit(1)
    cache.write_text(json.dumps(js))
    return js

def wrap(d):
    return (d + 180.0) % 360.0 - 180.0

def main():
    cats = []
    for kind in ('a', 'c'):
        js = sbdb_catalog(kind)
        f = js['fields']; idx = {n: k for k, n in enumerate(f)}
        rows = js['data']
        print(f'SBDB NEO {kind}: {len(rows)} rows', flush=True)
        cats.append((idx, rows))

    # build numeric arrays for matching
    def col(idx, rows, name):
        out = np.full(len(rows), np.nan)
        for k, r in enumerate(rows):
            v = r[idx[name]]
            try: out[k] = float(v)
            except (TypeError, ValueError): pass
        return out
    allrows = []; A=[];E=[];I=[];OM=[];W=[];MA=[];EP=[]
    for idx, rows in cats:
        allrows += [(idx, r) for r in rows]
        A.append(col(idx, rows, 'a')); E.append(col(idx, rows, 'e')); I.append(col(idx, rows, 'i'))
        OM.append(col(idx, rows, 'om')); W.append(col(idx, rows, 'w')); MA.append(col(idx, rows, 'ma')); EP.append(col(idx, rows, 'epoch'))
    A=np.concatenate(A);E=np.concatenate(E);I=np.concatenate(I);OM=np.concatenate(OM);W=np.concatenate(W);MA=np.concatenate(MA);EP=np.concatenate(EP)

    results = []
    for row in MEA:
        mid = int(row[0]); a, e, inc, om, w, ma = row[1:7]
        # normalized element distance (a: 1e-3 AU, e: 5e-4, angles: 0.02 deg)
        d = (np.abs(A - a) / 1e-3) ** 2 + (np.abs(E - e) / 5e-4) ** 2 + (np.abs(I - inc) / 0.02) ** 2 \
            + (wrap(OM - om) / 0.02) ** 2 + (wrap(W - w) / 0.02) ** 2
        d = np.where(np.isnan(d), np.inf, d)
        order = np.argsort(d)[:3]
        best = int(order[0]); idx, r = allrows[best]
        g = lambda n: r[idx[n]]
        rec = dict(ID=mid, a=a, e=e, i=inc, Omega=om, omega=w, M=ma,
                   q=a*(1-e), Q=a*(1+e), period_yr=a**1.5,
                   match_full_name=g('full_name').strip(), pdes=g('pdes'), spkid=g('spkid'), kind=g('kind'),
                   sbdb_class=g('class'), pha=g('pha'), H=g('H'), diameter_km=g('diameter'), moid_au=g('moid'),
                   data_arc_d=g('data_arc'), n_obs=g('n_obs_used'), condition_code=g('condition_code'),
                   first_obs=g('first_obs'), last_obs=g('last_obs'), soln_date=g('soln_date'), orbit_id=g('orbit_id'),
                   sbdb_epoch_jd=EP[best], sbdb_a=A[best], sbdb_e=E[best], sbdb_i=I[best], sbdb_om=OM[best], sbdb_w=W[best], sbdb_ma=MA[best],
                   match_dist=float(math.sqrt(d[best])),
                   second_best=allrows[int(order[1])][1][allrows[int(order[1])][0]['full_name']].strip(),
                   second_dist=float(math.sqrt(d[int(order[1])])),
                   dM_sbdb=float(wrap(MA[best]-ma)) if abs(EP[best]-EPOCH_JD) < 1e-6 else None)
        results.append(rec)
        print(f"{mid:4d} -> {rec['match_full_name']:32s} dist={rec['match_dist']:.3f} (2nd {rec['second_dist']:.1f}) "
              f"epoch={EP[best]:.1f} dM={rec['dM_sbdb']}", flush=True)

    if DO_HORIZONS:
        print('\n--- Horizons confirmation at JD 2461200.5 TDB (heliocentric ecliptic J2000 osculating elements) ---', flush=True)
        pat = {k: re.compile(r'(?<![A-Za-z])' + k + r'\s*=\s*([-+0-9.E]+)') for k in ['EC', 'IN', 'OM', 'W', 'MA', 'A', 'QR']}
        for rec in results:
            cmd = f"'{rec['spkid']}'" if rec['kind'] == 'an' or rec['kind'] == 'au' else f"'DES={rec['pdes']};'"
            # SPK-ID works for all small bodies; try it first, then fall back to designation
            for command in (f"'{rec['spkid']}'", f"'DES={rec['pdes']};'", f"'{rec['pdes']};'"):
                q = dict(format='text', COMMAND=command, OBJ_DATA='NO', MAKE_EPHEM='YES', EPHEM_TYPE='ELEMENTS',
                         CENTER="'500@10'", REF_PLANE='ECLIPTIC', TLIST="'2461200.5'", OUT_UNITS="'AU-D'", CSV_FORMAT='NO')
                url = 'https://ssd.jpl.nasa.gov/api/horizons.api?' + urllib.parse.urlencode(q)
                try:
                    txt = fetch(url, retries=3, timeout=90)
                except Exception as ex:
                    txt = f'ERROR {ex}'
                if '$$SOE' in txt: break
                time.sleep(0.5)
            hz = {}
            if '$$SOE' in txt:
                blk = txt.split('$$SOE')[1].split('$$EOE')[0]
                for k, p in pat.items():
                    m = p.search(blk); hz[k] = float(m.group(1)) if m else None
                rec.update(hz_a=hz['A'], hz_e=hz['EC'], hz_i=hz['IN'], hz_om=hz['OM'], hz_w=hz['W'], hz_ma=hz['MA'],
                           hz_da=hz['A']-rec['a'], hz_de=hz['EC']-rec['e'], hz_di=hz['IN']-rec['i'],
                           hz_dom=float(wrap(hz['OM']-rec['Omega'])), hz_dw=float(wrap(hz['W']-rec['omega'])), hz_dM=float(wrap(hz['MA']-rec['M'])))
                rec['hz_maxdev'] = max(abs(rec['hz_da'])/rec['a'], abs(rec['hz_de']), abs(rec['hz_di'])/max(1e-9,abs(rec['i'])), abs(rec['hz_dom']), abs(rec['hz_dw']), abs(rec['hz_dM']))
                print(f"{rec['ID']:4d} {rec['match_full_name']:30s} Horizons: da={rec['hz_da']:+.2e} de={rec['hz_de']:+.2e} di={rec['hz_di']:+.2e} "
                      f"dOm={rec['hz_dom']:+.2e} dw={rec['hz_dw']:+.2e} dM={rec['hz_dM']:+.2e}", flush=True)
            else:
                rec['hz_error'] = txt[:300].replace('\n', ' ')
                print(f"{rec['ID']:4d} Horizons FAILED: {rec['hz_error'][:150]}", flush=True)
            time.sleep(0.4)

    (DATA / 'nea_identification.json').write_text(json.dumps(results, indent=1, default=str))
    keys = list(results[0].keys())
    for r in results:
        for k in r:
            if k not in keys: keys.append(k)
    with open(DATA / 'nea_identification.csv', 'w', newline='') as f:
        wr = csv.DictWriter(f, fieldnames=keys); wr.writeheader()
        for r in results: wr.writerow(r)
    dists = np.array([r['match_dist'] for r in results])
    print(f'\nMatched {np.sum(dists < 1.0)} / {len(results)} within tight tolerance; max dist = {dists.max():.3f}')
    if DO_HORIZONS:
        devs = np.array([r.get('hz_maxdev', np.nan) for r in results])
        print(f'Horizons max relative/angle deviation over all objects: {np.nanmax(devs):.3e}; failures: {int(np.isnan(devs).sum())}')

if __name__ == '__main__':
    main()
