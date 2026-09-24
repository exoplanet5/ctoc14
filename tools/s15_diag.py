"""Stage 15 diagnostic (task step 3): are the ISOLATED targets' legs mispriced by the planner?

For every isolated target X on a t10d route (results/s12/cover_clean/fleet = the t10d partition twin) and for a control
sample of ordinary legs of the same routes, compare
  planner prices of the leg INTO X
    chain   : ctoc14.search.expand from the PLANNER's own state carried along the route (results/s14/A/chaincheck.json;
              production dv 1.2 / drmax 0.15, relaxed 2.0 / 0.25; 'reseed' = not admissible even relaxed, single-rev
              Lambert junction at the source epoch)
    twin-st : single-leg costs from the TWIN's state at the previous flyby (results/s13/linchpin/legcheck_t10d.json:
              lam1 = single-rev Lambert junction, lin = LinLeg)
  twin costs
    leg     : the twin's own impulses on that leg (legcheck dv_twin)
    removal : TRUE marginal = tank(route) - tank(route without X), re-settled (run_ialns.w_remove)
    lintwin : s14_twinbeam.lin_price of re-inserting X at its t10d epoch into the settled route-without-X (all impulses,
              epochs and v_inf free, linearised; one min-L1 solve)
Prices are converted to kg with the planner's leg formula at m = 1400 kg (dm = m (1 - exp(-kappa dv/ve)),
kappa = 1 + dv/(2 a tof)).
usage: s15_diag.py OUT.json [--control 24] [--nproc 8]"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, signal, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
for _p in (str(ROOT), str(ROOT / 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)
import numpy as np
import run_ialns as RI
from greedy_cover import _timeout
from ctoc14.constants import DAY, VE, TMAX
import s14_twinbeam as TB

ISO = [64, 119, 137, 138, 157, 158, 168, 216]
FLEET = ROOT / 'results/s12/cover_clean/fleet'


def kg(dv, tof_d, m=1400.0):
    if dv is None or not np.isfinite(dv):
        return None
    a = TMAX / m * 1e-3
    kap = 1.0 + dv / (2 * a * max(tof_d, 1.0) * DAY)
    return float(m * (1.0 - np.exp(-kap * dv / VE)))


def w_job(job):
    """Removal marginal + lintwin re-insertion price of (route, X)."""
    name, st, X = job
    tic = time.time(); out = dict(route=name, ast=X)
    try:
        signal.signal(signal.SIGALRM, _timeout); signal.setitimer(signal.ITIMER_REAL, 240)
        r = RI.w_remove((name, st, [X]))
    except Exception as e:
        r = dict(ok=False, error=type(e).__name__)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    tank0 = float(RI.ipr(st).tank())
    if not r.get('ok'):
        out.update(removal_ok=False, sec=round(time.time() - tic, 1)); return out
    out.update(removal_ok=True, tank0=tank0, tank_minus=r['tank'], removal_kg=tank0 - r['tank'],
               removal_dJ=RI.cost(tank0) - RI.cost(r['tank']))
    # lintwin price of re-inserting X at its t10d epoch (and +-3 d) into the settled route-without-X
    tX = float(st['tf'][list(int(a) for a in st['asts']).index(X)])
    ip = RI.ipr(r['st'])
    par = dict(st=r['st'], nreg=10 ** 6, dv=float(ip.dv()))
    best = None
    for dt in (0.0, -3.0, 3.0):
        try:
            dvm, lin, cm = TB.lin_price(par, X, tX + dt * DAY)
        except Exception:
            continue
        rec = (float(dvm), float(lin['res_km']), float(cm), dt)
        if best is None or rec[0] < best[0]:
            best = rec
    if best is not None:
        m = r['tank']
        out.update(lintwin_dv=best[0], lintwin_res_km=best[1], lintwin_coast_au=best[2], lintwin_dt_d=best[3],
                   lintwin_kg=float(m * (np.exp(max(best[0], 0.0) / VE) - 1.0)))
    out['sec'] = round(time.time() - tic, 1)
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--control', type=int, default=24)
    ap.add_argument('--nproc', type=int, default=8); ap.add_argument('--seed', type=int, default=15)
    a = ap.parse_args()
    F = RI.IFleet(FLEET)
    lc = json.load(open(ROOT / 'results/s13/linchpin/legcheck_t10d.json'))['routes']
    cc = json.load(open(ROOT / 'results/s14/A/chaincheck.json'))['routes']
    rows = {}
    for rn, r in lc.items():
        cmap = {d['ast']: d for d in cc[rn]['detail']}
        for l in r['legs']:
            if l['k'] == 0:
                continue
            c = cmap.get(l['ast'], {})
            chain_dv = c.get('dv_leg') if c.get('prod') or c.get('rel') else c.get('reseed_lam1_dv')
            rows[(rn, l['ast'])] = dict(route=rn, ast=l['ast'], k=l['k'], tof_d=l['tof_d'], iso=l['ast'] in ISO,
                                        chain_prod=c.get('prod'), chain_rel=c.get('rel'), chain_dv=chain_dv,
                                        twin_lam1=l['lam1_dv'], twin_lin=l['lin_dv'], coast_au=l['coast_miss_au'],
                                        ok_planner_twinstate=l['ok_planner'], ok_relaxed_twinstate=l['ok_relaxed'],
                                        leg_dv_twin=l['dv_twin'])
    iso_keys = [k for k, v in rows.items() if v['iso']]
    rng = np.random.default_rng(a.seed)
    others = [k for k, v in rows.items() if not v['iso'] and k[0] in {'04', '05', '07', '08'}]
    ctrl = [others[i] for i in rng.choice(len(others), size=min(a.control, len(others)), replace=False)]
    keys = iso_keys + ctrl
    RI.eph()
    tic = time.time()
    with mp.get_context('fork').Pool(a.nproc, maxtasksperchild=4) as pool:
        res = pool.map(w_job, [(rn, F.routes[rn]['st'], X) for rn, X in keys], chunksize=1)
    for (rn, X), r in zip(keys, res):
        v = rows[(rn, X)]; v.update(r)
        v['chain_kg'] = kg(v['chain_dv'], v['tof_d']); v['twin_lam1_kg'] = kg(v['twin_lam1'], v['tof_d'])
        v['twin_lin_kg'] = kg(v['twin_lin'], v['tof_d']); v['leg_twin_kg'] = kg(v['leg_dv_twin'], v['tof_d'])
    out = [rows[k] for k in keys]
    # summary
    def med(xs):
        xs = [x for x in xs if x is not None and np.isfinite(x)]
        return round(float(np.median(xs)), 2) if xs else None
    summ = {}
    for grp, sel in (('iso', [r for r in out if r['iso']]), ('control', [r for r in out if not r['iso']])):
        ratio = [r['chain_kg'] / max(r['removal_kg'], 1.0) for r in sel if r.get('removal_ok') and r['chain_kg'] is not None]
        summ[grp] = dict(n=len(sel), chain_admissible_prod=sum(1 for r in sel if r['chain_prod']),
                         chain_admissible_rel=sum(1 for r in sel if r['chain_rel']),
                         med_chain_kg=med([r['chain_kg'] for r in sel]), med_removal_kg=med([r.get('removal_kg') for r in sel]),
                         med_lintwin_kg=med([r.get('lintwin_kg') for r in sel]),
                         med_ratio_chain_over_removal=med(ratio),
                         frac_chain_over_2x=round(float(np.mean([x > 2 for x in ratio])), 3) if ratio else None)
    json.dump(dict(summary=summ, rows=out, wall_s=round(time.time() - tic, 1)), open(a.out, 'w'), indent=1, default=float)
    print(f'{"route":>5} {"ast":>4} {"iso":>3} {"tof":>6} | {"chain dv/kg":>14} {"P/R":>4} | {"twinst lam1/lin kg":>18} | '
          f'{"leg twin kg":>11} | {"removal kg":>10} | {"lintwin kg":>10} {"res km":>9}')
    for r in out:
        f = lambda x, d=1: '-' if x is None else f'{x:.{d}f}'
        print(f'{r["route"]:>5} {r["ast"]:>4} {"*" if r["iso"] else " ":>3} {r["tof_d"]:6.1f} | {f(r["chain_dv"], 3):>6}/{f(r["chain_kg"]):>6} '
              f'{("Y" if r["chain_prod"] else "n") + ("Y" if r["chain_rel"] else "n"):>4} | {f(r["twin_lam1_kg"]):>8}/{f(r["twin_lin_kg"]):>8} | '
              f'{f(r["leg_twin_kg"]):>11} | {f(r.get("removal_kg")):>10} | {f(r.get("lintwin_kg")):>10} {f(r.get("lintwin_res_km"), 0):>9}')
    print(json.dumps(summ, indent=1))


if __name__ == '__main__':
    main()
