#!/usr/bin/env python
"""Match all 300 MEA.txt targets against the MPC NEA catalog (MPCORB NEA.txt, epoch K2669 = 2026-06-09 = MJD 61200)
and merge with the JPL SBDB match (data/nea_identification.json if present). Writes data/nea_catalog_ids.csv/json."""
import json, pathlib, csv, numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]; DATA = ROOT/'data'
rows = []
for line in open(DATA/'MPC_NEA.txt', errors='ignore'):
    if len(line) < 105: continue
    try:
        des=line[0:7].strip(); H=line[8:13].strip(); M=float(line[26:35]); w=float(line[37:46]); Om=float(line[48:57]); i=float(line[59:68]); e=float(line[70:79]); a=float(line[92:103]); ep=line[20:25]
        readable=line[166:194].strip(); nobs=line[117:122].strip(); arc=line[127:136].strip(); U=line[105:106].strip(); flags=line[161:165].strip()
    except ValueError: continue
    rows.append(dict(mpc_des=des, mpc_name=readable, a=a, e=e, i=i, Om=Om, w=w, M=M, epoch=ep, H=H, nobs=nobs, arc=arc, U=U, flags=flags))
A=np.array([r['a'] for r in rows]);E=np.array([r['e'] for r in rows]);I=np.array([r['i'] for r in rows]);OM=np.array([r['Om'] for r in rows]);W=np.array([r['w'] for r in rows]);MM=np.array([r['M'] for r in rows])
def wrap(d): return (d+180)%360-180
mea=np.loadtxt(ROOT/'MEA.txt')
sb = {}
p = DATA/'nea_identification.json'
if p.exists():
    for r in json.load(open(p)): sb[r['ID']] = r
out=[]
for row in mea:
    mid=int(row[0])
    d=(np.abs(A-row[1])/1e-4)**2+(np.abs(E-row[2])/1e-5)**2+(np.abs(I-row[3])/1e-3)**2+(wrap(OM-row[4])/1e-3)**2+(wrap(W-row[5])/1e-3)**2+(wrap(MM-row[6])/1e-3)**2
    k=int(np.argmin(d)); r=rows[k]
    rec=dict(ID=mid, a=row[1], e=row[2], i=row[3], Omega=row[4], omega=row[5], M=row[6], q=row[1]*(1-row[2]), Q=row[1]*(1+row[2]), P_yr=row[1]**1.5,
             mpc_des=r['mpc_des'], mpc_name=r['mpc_name'], mpc_H=r['H'], mpc_nobs=r['nobs'], mpc_arc=r['arc'], mpc_U=r['U'], mpc_flags=r['flags'], mpc_epoch=r['epoch'],
             mpc_dist=float(np.sqrt(d[k])), mpc_max_abs_diff=float(max(abs(A[k]-row[1]),abs(E[k]-row[2]),abs(I[k]-row[3]),abs(wrap(OM[k]-row[4])),abs(wrap(W[k]-row[5])),abs(wrap(MM[k]-row[6])))))
    s=sb.get(mid)
    if s:
        rec.update(sbdb_name=s['match_full_name'], sbdb_spkid=s['spkid'], sbdb_dist=s['match_dist'], sbdb_H=s['H'], sbdb_diameter_km=s['diameter_km'], sbdb_moid_au=s['moid_au'], sbdb_pha=s['pha'], sbdb_class=s['sbdb_class'],
                   sbdb_arc_d=s['data_arc_d'], sbdb_nobs=s['n_obs'], sbdb_cc=s['condition_code'], hz_maxdev=s.get('hz_maxdev'))
    out.append(rec)
keys=[]
for r in out:
    for k in r:
        if k not in keys: keys.append(k)
with open(DATA/'nea_catalog_ids.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=keys); w.writeheader(); [w.writerow(r) for r in out]
json.dump(out, open(DATA/'nea_catalog_ids.json','w'), indent=1)
md=np.array([r['mpc_max_abs_diff'] for r in out]); H=np.array([float(r['mpc_H'] or 'nan') for r in out])
print(f"MPC exact matches (max |diff| < 1e-6 in all 6 elements): {(md<1e-6).sum()}/300 ; worst {md.max():.2e}")
print(f"H: min {np.nanmin(H):.2f} max {np.nanmax(H):.2f} median {np.nanmedian(H):.2f}; H<=18: {(H<=18).sum()}, H<=17: {(H<=17).sum()}, H<=16: {(H<=16).sum()}, H<=15: {(H<=15).sum()}")
named=[r for r in out if r['mpc_name'].startswith('(')]; print(f"numbered: {len(named)}; unnumbered: {300-len(named)}")
flags=[r['mpc_flags'] for r in out]; print('MPC flags (hex, orbit type/PHA bits) top values:', {f:flags.count(f) for f in sorted(set(flags))})
if sb:
    pha=[r.get('sbdb_pha') for r in out]; print('SBDB PHA Y/N:', pha.count('Y'), pha.count('N'))
    D=np.array([float(r['sbdb_diameter_km'] or 'nan') if r.get('sbdb_diameter_km') not in (None,'') else np.nan for r in out]); print(f"SBDB diameters known for {np.isfinite(D).sum()}: min {np.nanmin(D):.2f} max {np.nanmax(D):.2f} km")
    moid=np.array([float(r['sbdb_moid_au']) if r.get('sbdb_moid_au') not in (None,'') else np.nan for r in out]); print(f"SBDB Earth MOID: min {np.nanmin(moid):.4f} max {np.nanmax(moid):.4f} AU; MOID<0.05: {(moid<0.05).sum()}")
print("\nbrightest 10:"); [print('  ', r['ID'], r['mpc_name'] or r['mpc_des'], r['mpc_H']) for r in sorted(out, key=lambda r: float(r['mpc_H'] or 99))[:10]]
