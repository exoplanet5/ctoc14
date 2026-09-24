"""Fleet-level ANT COLONY OPTIMISATION over the WHOLE partition (NUDT's GTOC9 top level, on CTOC14-A).

WHY
---
Every fleet we have built so far is a sequence of independent decisions: build a route, then another, then force the
leftovers somewhere.  Measured consequence (stage 7-10): a target a route CHOOSES costs 0.032-0.042 J, a target FORCED
on a finished route costs 0.115-0.148 J.  NUDT placed 2nd in GTOC9 -- a partition problem with the same
fixed+variable per-mission cost -- with an ACO top level in which ONE ant builds ALL missions in ONE constructive
pass.  Nothing is ever forced onto a finished route (the ant assigns each target exactly once by construction) and
pheromone learns, across thousands of passes, which targets to SAVE for a later craft.  Greedy construction
structurally cannot do that: it gives craft 1 forty targets and craft 8 twenty.

WHAT AN ANT DOES
----------------
    R = all 298 reachable targets
    while R and n_craft < NMAX:
        carrier c ~ tau_c^alpha * eta_c(R)^beta          eta_c = pheromone-free reach of c into R
        B  = the carrier's DP base restricted to R       (real (time,lag) DP events, carrier_dp)
        repeat: t* = argmax over R of tau[c,t]^alpha * (1/dJ_hat(t))^beta
                accept while dJ_hat(t*) < BREAKEVEN (the J one flyby is worth in a J<13 fleet)
        R -= (B + accepted);  emit craft (c, S)
    J = sum_i J_i(tank_i) + |R| + 2

Pheromone lives on (carrier) and on (carrier, target); elite ants deposit 1/(J - JREF).

THE COST SURROGATE, AND WHY IT IS BUILT THIS WAY
------------------------------------------------
Regressing the tank of 96 real farm routes (results/s7/farmO) on route features gives: DEPTH ALONE R2 0.553, and
every target-identity feature together adds only 0.06 -- fuel ~ 17.9 kg per flyby - 238 kg, essentially independent
of WHICH targets the route flew, even though their carrier closest approaches span 0.025-0.53 AU.  That sample is
selection-truncated (the farm only ever accepted insertions under 45 kg), so a cost regression on it cannot price a
FORCED target and would let an ACO assign anything anywhere for nothing -- pure air.

The signal that survives the truncation is SELECTION, not price.  Over all 96 x 298 = 28 608 (route, target) pairs,
P(a free-choice route took the target) falls from 0.21 below 0.02 AU to 0.014 beyond 0.5 AU (non-base pairs, base
rate 0.061), and a DP-base target is taken 999 times in 1000.  So the surrogate prices a target by how RELUCTANTLY
routes pick it up: REL(d) = P(take at 0.01 AU) / P(take at d), scaled by KX and by the measured depth penalty
(depth/30)^QD.  F0, KB, KX, QD are then fitted to the 96 real tanks.

Independent check, on a fleet the fit never saw: the banked-equivalent 10-craft fleet (results/s7/cover1, real
sum J_i 12.3546) is priced at 12.3270 when the surrogate is given FREE carrier choice -- 0.028 J of optimism over
10 craft.  On the calibration routes themselves free carrier choice buys 0.056 J per craft of optimism, so read
any ACO J as "+0.0 to +0.4 optimistic" before the depth extrapolation is even considered.

MEASURED VERDICT (2026-09-21) -- THE PREMISE IS WRONG, AND HERE IS THE NUMBER
------------------------------------------------------------------------------
The ACO works as an algorithm.  Four campaigns (results/s10/aco/A-D) learn: elite-mean J falls 0.55-1.30 over
55-116 generations, best-ever 0.25-0.42, and every run reaches 298/298 coverage.  Best surrogate fleets:
8 craft x 38 flybys, sum J_i 10.259, J 12.259 (plan D, capped at depth 38 / 1100 kg to stay inside the measured
envelope); 7 x 45 gives 11.62 if the caps are lifted.

Then plan D was BUILT for real (results/s10/aco/buildD):

    surrogate     8 craft,  38 flybys each,  298/298 covered,  sum J_i 10.259,  J  12.259
    settled       7 craft,  14-28 flybys,    151/298 covered,  sum J_i  8.236,  J 157.236

One craft could not even start (its carrier had only 2 of its DP-base targets left -- hence --nbmin and --dp).
The other seven flew 151 of their 298 assigned targets: 51 %.  Five of the seven stopped because no remaining
ASSIGNED target could be inserted at all (three consecutive homotopy failures), not because of the time budget.

    free-choice routes, same 240-carrier library (results/s7/farmO, 96 routes)
        depth median 34 (max 41),  J/flyby median 0.0395, best 0.0338
    assigned-set routes (this run)
        depth median 22 (max 28),  J/flyby median 0.0541, best 0.0436

So the cost model was NOT the problem -- it is well calibrated (the banked-equivalent 10-craft fleet, never seen
by the fit, real sum J_i 12.3546, priced at 12.3270 with free carrier choice; at depth >= 36 the median optimism
over 101 real routes is 0.018 J per craft).  What failed is the assumption underneath the whole ACO idea: that
"assign each target exactly once by construction" avoids forcing.  IT IS FORCING.  Fixing a target's membership
before the trajectory exists is exactly the move that costs 3x, whether it happens first or last.  The premium
then shows up not as a price but as an INFEASIBILITY: the route simply runs out of insertable assigned targets at
depth 14-28, while the same carrier with free choice reaches 34.

The missing term is per-(carrier, target) FEASIBILITY, and it cannot be had from the ballistic carrier distance:
every target is within 0.055 AU of some carrier (median 0.013), so distance says everything is reachable, while
insert_model's P(settle) -- measured on 1225 real settles -- falls to 0.31 at depth 30-37 and 0.21 at 38-48.  An
ACO over raw targets is therefore searching a space whose feasibility it cannot see.  The repair is to run the
colony over PRE-BUILT FRAGMENTS (Tsinghua's GTOC11 architecture) rather than over raw targets, so that everything
an ant assembles is known to be flyable, or to give the ant a "hand the target back" move priced by predict_ok.

Usage
-----
    fleet_aco.py calibrate [--out results/s10/aco/route_model.json]
    fleet_aco.py diag                       structural limits (base packing, per-target reach)
    fleet_aco.py score <fleet dirs...>      surrogate vs real settled tanks, by depth band
    fleet_aco.py run  --gens 40 --ants 64 --nproc 8 --out results/s10/aco/A
    fleet_aco.py build results/s10/aco/A/best.json --out results/s10/aco/fleetA [--nproc 8]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, time, pickle, pathlib, argparse, multiprocessing as mp
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
from ctoc14.constants import DAY
from ctoc14.search import UNREACHABLE

ALL = np.array([t for t in range(1, 301) if t not in UNREACHABLE])
NT = len(ALL)
M_DRY, FUEL_MAX = 600.0, 1400.0
JREF = 10.0                      # deposit 1/(J - JREF); J < 13 is the target, ~11 is the theoretical floor


def cost_J(tank):
    x = (np.asarray(tank, float) - M_DRY) / FUEL_MAX
    return 1.0 + x + x * x


# ----------------------------------------------------------------------------- data
class Data:
    """Everything an ant needs: the reach table, the DP bases, and the index maps."""

    def __init__(self, lib=None, reach=None):
        lib = lib or str(ROOT / 'results/s7/libopt.pkl')
        reach = reach or str(ROOT / 'results/s10/reach.npz')
        Z = np.load(reach)
        self.D = Z['D'].astype(np.float64)          # (ncar, NT) closest approach [AU]
        self.E = Z['E'].astype(np.float64)          # (ncar, NT) epoch of it [s]
        self.ids = Z['ids'].astype(int)
        self.idx = {int(t): i for i, t in enumerate(self.ids)}
        self.lib = pickle.load(open(lib, 'rb'))
        self.ncar = len(self.lib)
        self.base = [np.array(sorted(self.idx[t] for t, _, _ in c['base'] if t in self.idx), int)
                     for c in self.lib]
        self.base_t = [np.array([c['base'][j][1] for j in range(len(c['base']))
                                 if c['base'][j][0] in self.idx], float) for c in self.lib]
        self.tL = np.array([c['tL'] for c in self.lib], float)
        self.vinf = np.array([c['vinf'] for c in self.lib], float)
        self.isbase = np.zeros((self.ncar, NT), bool)
        for c in range(self.ncar):
            self.isbase[c, self.base[c]] = True
        self.ev = None

    def load_events(self, eps=0.01):
        """Per-carrier (time, lag) event plane, so an ant can RE-RUN the carrier DP against whatever targets are
        still free instead of reusing a base that earlier craft have already eaten.  ~1000 events per carrier,
        essentially free to build; dp_route on one is 0.05-0.11 s."""
        if self.ev is not None: return self.ev
        import carrier_dp as CD
        from ctoc14.kepler import load_mea, Ephemeris
        from exp_carrier import rv2frame
        eph = Ephemeris(); ids, elem = load_mea()
        keep = np.array([int(x) in self.idx for x in ids]); ids = ids[keep]; elem = elem[keep]
        self.ev = []
        for c in range(self.ncar):
            rE, vE = eph.earth_state(self.tL[c])
            r = np.array(rE, float); v = np.array(vE, float) + self.vinf[c]
            t, lag, gap, ast, Pc = CD.events_for(dict(tL=float(self.tL[c]), r=r, v=v), ids, elem, eps)
            self.ev.append(dict(t=t, lag=lag, gap=gap, ast=ast,
                                ac=float(rv2frame(r[None], v[None])[3][0])))
        return self.ev


_D = None


def data():
    global _D
    if _D is None: _D = Data()
    return _D


# ----------------------------------------------------------------------------- route surrogate
DEFAULT_MODEL = dict(F0=-40.0, KB=13.0, KX=5.0, QD=0.954, SEL_A=-1.2, SEL_B=-0.55, D0=0.01,
                     RELMAX=25.0, MARGIN=90.0, DMAX=1.20, PFLOOR=0.20)


def load_model(path=None):
    p = pathlib.Path(path or (ROOT / 'results/s10/aco/route_model.json'))
    if p.exists(): return json.load(open(p))
    return dict(DEFAULT_MODEL)


class RouteModel:
    """tank(carrier, assigned set), in two measured pieces.

    1. SELECTION.  Over all 96 farm routes x 298 targets (28 608 pairs, 3 220 flown) the probability that a
       free-choice route TOOK a target falls with the closest approach of its ballistic carrier:
       P = 0.48 (d < 0.02 AU), 0.29 (0.02-0.04), 0.20 (0.04-0.06), 0.14 (0.06-0.15), 0.10 (0.2-0.5), 0.03 (> 0.5),
       against a base rate of 0.113; a DP-base target is taken 999 times in 1000.  Fit
       P(take) = sigma(SEL_A + SEL_B log10 d) and define the RELATIVE PRICE of a target as
       REL(d) = P(take at D0) / P(take at d).  This is the only honest distance signal in the data: the farm's
       kg costs are selection-truncated (it only ever accepted insertions under 45 kg), so a cost regression on
       them reports "distance does not matter" (R2 0.55 from depth alone, +0.06 for every identity feature).

    2. PRICE.  tank = 600 + F0 + KB*nb + sum over the extras, nearest first, of KX * REL(d) * (depth/30)^QD,
       with F0, KB, KX, QD fitted to the 96 REAL settled tanks.  QD carries the measured depth penalty
       (insert_model: kg ~ depth^0.95).
    """

    def __init__(self, m=None):
        m = m or load_model()
        self.F0, self.KB, self.KX, self.QD = m['F0'], m['KB'], m['KX'], m['QD']
        self.SEL_A, self.SEL_B, self.D0 = m['SEL_A'], m['SEL_B'], m['D0']
        self.RELMAX, self.MARGIN = m.get('RELMAX', 25.0), m['MARGIN']
        self.DMAX, self.PFLOOR = m.get('DMAX', 1.20), m.get('PFLOOR', 0.20)
        import insert_model as IM
        self.w = IM.model()['success']['w']
        self._lm = np.log10(max(self.MARGIN, 1.0))
        self._p0 = self.ptake(self.D0)

    def ptake(self, d):
        z = self.SEL_A + self.SEL_B * np.log10(np.maximum(d, 1e-4))
        return 1.0 / (1.0 + np.exp(-np.clip(z, -40, 40)))

    def rel(self, d):
        """Relative price of a target at carrier distance d (1.0 at D0, RELMAX at the far end)."""
        return np.minimum(self._p0 / np.maximum(self.ptake(d), 1e-6), self.RELMAX)

    def kg(self, d, depth, tank=0.0):
        return self.KX * self.rel(d) * (max(depth, 1.0) / 30.0) ** self.QD

    def kg_v(self, d, depth):
        return self.KX * self.rel(np.asarray(d)) * (max(depth, 1.0) / 30.0) ** self.QD

    def tank(self, c, S, dat=None):
        """S: iterable of TARGET INDICES (0..NT-1).  Returns (tank, n_flown, n_dropped)."""
        dat = dat or data()
        S = np.asarray(list(S), int)
        if len(S) == 0: return M_DRY, 0, 0
        isb = dat.isbase[c, S]
        nb = int(isb.sum())
        tank = M_DRY + self.F0 + self.KB * nb
        dd = np.sort(dat.D[c, S[~isb]])
        depth = max(nb, 1); drop = 0
        for d in dd:
            if d > self.DMAX: drop += 1; continue
            tank += self.kg(d, depth); depth += 1
        return max(M_DRY + 1.0, tank), len(S) - drop, drop

    def marginal(self, c, ti, depth, tank=0.0, dat=None):
        dat = dat or data()
        if dat.isbase[c, ti]: return self.KB
        return self.kg(dat.D[c, ti], depth)


# ----------------------------------------------------------------------------- calibration
def _farm_routes():
    """(carrier index, set of target indices, real settled tank) for every route of results/s7/farmO."""
    import run_ialns as RI
    dat = data(); out = []
    for f in sorted((ROOT / 'results/s7/farmO').glob('route_c*.npz')):
        k = int(f.stem[7:]); z = np.load(f); st = {kk: z[kk] for kk in z.files}
        S = [dat.idx[int(a)] for a in st['asts'] if int(a) in dat.idx]
        out.append((k, S, float(RI.ipr(st).tank())))
    return out


def cmd_calibrate(a):
    from scipy.optimize import minimize
    dat = data(); R = _farm_routes()
    print(f'{len(R)} real farm routes, depth {min(len(s) for _,s,_ in R)}-{max(len(s) for _,s,_ in R)}, '
          f'tank {min(t for _,_,t in R):.0f}-{max(t for _,_,t in R):.0f} kg')

    # --- stage 1: the SELECTION curve, P(a free-choice route takes t) vs carrier distance
    X = []; Y = []
    for k, S, _ in R:
        took = np.zeros(NT, bool); took[np.asarray(S, int)] = True
        m = ~dat.isbase[k]                     # DP-base pairs are taken 99.9 % of the time: excluded
        X.append(np.log10(np.maximum(dat.D[k][m], 1e-4))); Y.append(took[m])
    X = np.concatenate(X); Y = np.concatenate(Y).astype(float)
    def nll(w):
        z = np.clip(w[0] + w[1] * X, -40, 40)
        return float(np.mean(np.logaddexp(0, z) - Y * z))
    rs = minimize(nll, np.array([-1.2, -0.55]), method='Nelder-Mead', options=dict(maxiter=4000))
    SEL_A, SEL_B = float(rs.x[0]), float(rs.x[1])
    print(f'selection: P(take) = sigma({SEL_A:+.3f} {SEL_B:+.3f} log10 d) over {len(X)} non-base pairs '
          f'({Y.sum():.0f} taken, base rate {Y.mean():.3f})')
    for lo, hi in ((0, .02), (.02, .04), (.04, .06), (.06, .10), (.10, .20), (.20, .50), (.50, 99)):
        m = (10 ** X >= lo) & (10 ** X < hi)
        if m.sum(): print(f'   d {lo:.2f}-{hi:<5.2f}  observed {Y[m].mean():.3f}  '
                          f'fitted {1/(1+np.exp(-(SEL_A+SEL_B*X[m]))).mean():.3f}  n {m.sum():6d}')

    # --- stage 2: F0, KB, KX, QD against the 96 real tanks
    def unpack(x):
        m = dict(DEFAULT_MODEL)
        m['F0'], m['KB'], m['KX'] = float(x[0]), float(abs(x[1])), float(abs(x[2]))
        m['QD'] = float(min(4.5, abs(x[3])))
        m['SEL_A'], m['SEL_B'] = SEL_A, SEL_B
        return m

    def err(x):
        rm = RouteModel(unpack(x))
        e = np.array([rm.tank(k, S, dat)[0] - t for k, S, t in R])
        return float(np.mean(e ** 2))

    best = None
    for x0 in ([-240, 18, 6, 0.95], [-40, 13, 4, 0.5], [-150, 16, 8, 1.4], [0, 10, 3, 0.2], [-300, 20, 10, 1.8]):
        r = minimize(err, np.array(x0, float), method='Nelder-Mead',
                     options=dict(maxiter=8000, maxfev=8000, xatol=1e-3, fatol=1e-3))
        if best is None or r.fun < best.fun: best = r
    m = unpack(best.x); rm = RouteModel(m)
    pred = np.array([rm.tank(k, S, dat)[0] for k, S, _ in R]); true = np.array([t for _, _, t in R])
    e = pred - true
    m.update(n_calib=len(R), rmse_kg=float(np.sqrt(np.mean(e ** 2))), median_abs_kg=float(np.median(np.abs(e))),
             r2=float(1 - np.var(e) / np.var(true)), bias_kg=float(e.mean()),
             J_abs_err=float(np.mean(np.abs(cost_J(pred) - cost_J(true)))), n_pairs=int(len(X)),
             sel_base_rate=float(Y.mean()))
    print(f'\nF0 {m["F0"]:+.1f} kg  KB {m["KB"]:.2f} kg/base-fb  KX {m["KX"]:.2f} kg/extra at d={m["D0"]}  '
          f'QD {m["QD"]:.3f}')
    print(f'  RMSE {m["rmse_kg"]:.1f} kg, median |err| {m["median_abs_kg"]:.1f}, R2 {m["r2"]:.3f}, '
          f'bias {m["bias_kg"]:+.1f} kg, mean |dJ_i| {m["J_abs_err"]:.4f}')
    print('  price of one extra at depth 30:  ' + '  '.join(
        f'd={d:.2f}:{rm.kg(d,30):.0f}kg' for d in (0.01, 0.03, 0.06, 0.12, 0.25, 0.5, 1.0)))
    out = pathlib.Path(a.out or (ROOT / 'results/s10/aco/route_model.json'))
    out.parent.mkdir(parents=True, exist_ok=True); json.dump(m, open(out, 'w'), indent=1)
    print('->', out)


# ----------------------------------------------------------------------------- the ant
def ant(dat, rm, tau_c, tau_ct, rng, a):
    """One complete fleet.  Returns dict(crafts=[(carrier, [target idx])], J, tanks, missed)."""
    R = np.ones(NT, bool)
    used = np.zeros(dat.ncar, bool)
    crafts = []; tanks = []; allbases = []
    dJ_break = a.breakeven
    while R.any() and len(crafts) < a.nmax:
        # --- choose a carrier
        near = (dat.D <= a.dnear) & R[None, :]
        eta = near.sum(1).astype(float) + 2.0 * (dat.isbase & R[None, :]).sum(1)
        eta[used] = 0.0
        if eta.max() <= 0: break
        p = (tau_c ** a.alpha) * (eta ** a.beta)
        p[used] = 0.0
        if p.sum() <= 0: break
        c = int(rng.choice(dat.ncar, p=p / p.sum()))
        used[c] = True
        # --- the carrier's DP base.  --dp re-runs the exact (time,lag) DP against only the targets still free,
        #     so a late craft gets a REAL base of its own instead of the scraps of a base earlier craft ate.
        bases = None
        if a.dp:
            import carrier_dp as CD
            e = dat.ev[c]
            w = np.zeros(301)
            w[dat.ids[R]] = 1.0
            try:
                val, seq = CD.dp_route(e['t'], e['lag'], e['gap'], e['ast'], float(dat.tL[c]),
                                       a.lam, a.kgap, ac=e['ac'], w=w)
            except Exception:
                seq = []
            seen = set(); seq = [i for i in seq if not (e['ast'][i] in seen or seen.add(e['ast'][i]))]
            seq = [i for i in seq if int(e['ast'][i]) in dat.idx and R[dat.idx[int(e['ast'][i])]]]
            bases = [(int(e['ast'][i]), float(e['t'][i]), float(e['lag'][i])) for i in seq]
            S = [dat.idx[b[0]] for b in bases]
        else:
            S = [int(i) for i in dat.base[c] if R[i]]
        if len(S) < a.nbmin:                       # no viable base from this carrier: it cannot fly
            continue
        for i in S: R[i] = False
        tank = M_DRY + rm.F0 + rm.KB * len(S); depth = max(len(S), 1)
        # --- depth phase: append while the marginal target is worth less than a flyby
        cand = np.nonzero(R & (dat.D[c] <= rm.DMAX))[0]
        if len(cand):
            dd = dat.D[c, cand]
            order = np.argsort(dd)
            cand = cand[order]; dd = dd[order]
            rel = rm.rel(dd)
            tw = tau_ct[c, cand] ** a.alpha
            for _ in range(a.depth - len(S)):
                if not len(cand) or depth >= a.depth: break
                kg = rm.KX * rel * (max(depth, 1.0) / 30.0) ** rm.QD
                dJ = cost_J(tank + kg) - cost_J(tank)
                ok = (dJ < dJ_break) & (tank + kg <= a.maxtank)
                if not ok.any(): break
                sc = np.where(ok, tw * (dJ_break / np.maximum(dJ, 1e-6)) ** a.beta, 0.0)
                if a.greedy or rng.random() < a.q0:
                    j = int(sc.argmax())
                else:
                    j = int(rng.choice(len(sc), p=sc / sc.sum()))
                t = int(cand[j]); tank += kg[j]; depth += 1
                R[t] = False; S.append(t)
                keep = np.ones(len(cand), bool); keep[j] = False
                cand, dd, rel, tw = cand[keep], dd[keep], rel[keep], tw[keep]
        if len(S) < a.nmin:                      # too shallow to be worth a craft: give the targets back
            for i in S: R[i] = True
            continue
        crafts.append((c, S)); tanks.append(tank); allbases.append(bases)
    missed = int(R.sum())
    J = float(sum(cost_J(t) for t in tanks)) + missed + 2.0
    return dict(crafts=crafts, tanks=tanks, missed=missed, J=J, bases=allbases,
                sumJi=float(sum(cost_J(t) for t in tanks)), covered=NT - missed)


# ----------------------------------------------------------------------------- generation worker
_W = {}


def _init(a, mj):
    _W['dat'] = data(); _W['rm'] = RouteModel(mj); _W['a'] = a
    if a.dp: _W['dat'].load_events(a.eps)


def _gen(job):
    seed, tau_c, tau_ct = job
    dat, rm, a = _W['dat'], _W['rm'], _W['a']
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(a.chunk):
        out.append(ant(dat, rm, tau_c, tau_ct, rng, a))
    out.sort(key=lambda r: r['J'])
    return out[:max(1, a.chunk // 4)]


def cmd_run(a):
    dat = data(); mj = load_model(a.model); rm = RouteModel(mj)
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    log = open(out / 'log.txt', 'a')

    def say(s):
        print(s, flush=True); log.write(s + '\n'); log.flush()

    say(f'# fleet_aco {time.strftime("%Y-%m-%d %H:%M:%S")}  carriers {dat.ncar}, targets {NT}')
    say(f'# model F0 {rm.F0:+.1f} KB {rm.KB:.2f} KX {rm.KX:.2f} QD {rm.QD:.2f} '
        f'sel({rm.SEL_A:+.2f},{rm.SEL_B:+.2f})  (calib rmse {mj.get("rmse_kg",0):.0f} kg, R2 {mj.get("r2",0):.2f})')
    say(f'# alpha {a.alpha} beta {a.beta} rho {a.rho} ants {a.ants} gens {a.gens} '
        f'nmax {a.nmax} depth {a.depth} breakeven {a.breakeven}')
    tau_c = np.ones(dat.ncar); tau_ct = np.ones((dat.ncar, NT))
    best = None; hist = []
    rng = np.random.default_rng(a.seed)
    nchunk = max(1, a.ants // a.chunk)
    ctx = mp.get_context('fork')
    t0 = time.time()
    with ctx.Pool(a.nproc, initializer=_init, initargs=(a, mj)) as pl:
        for g in range(a.gens):
            jobs = [(int(rng.integers(1 << 30)), tau_c, tau_ct) for _ in range(nchunk)]
            elite = [r for part in pl.map(_gen, jobs) for r in part]
            elite.sort(key=lambda r: r['J'])
            if best is None or elite[0]['J'] < best['J']:
                best = elite[0]
                json.dump(dict(J=best['J'], sumJi=best['sumJi'], covered=best['covered'],
                               missed=best['missed'], gen=g, dp=bool(a.dp),
                               tanks=[float(x) for x in best['tanks']],
                               bases=best['bases'],
                               crafts=[[int(c), [int(dat.ids[i]) for i in S]] for c, S in best['crafts']]),
                          open(out / 'best.json', 'w'), indent=1)
            # --- pheromone
            tau_c *= (1 - a.rho); tau_ct *= (1 - a.rho)
            for r in elite[:a.nelite]:
                dep = 1.0 / max(r['J'] - JREF, 0.25)
                for c, S in r['crafts']:
                    tau_c[c] += dep
                    tau_ct[c, np.asarray(S, int)] += dep
            np.clip(tau_c, a.tmin, a.tmax, out=tau_c); np.clip(tau_ct, a.tmin, a.tmax, out=tau_ct)
            Js = np.array([r['J'] for r in elite])
            hist.append(dict(gen=g, best=float(Js.min()), mean=float(Js.mean()),
                             best_ever=float(best['J']), covered=int(elite[0]['covered']),
                             n=len(elite[0]['crafts'])))
            say(f'gen {g:3d}  best {Js.min():7.3f}  eliteMean {Js.mean():7.3f}  ever {best["J"]:7.3f} '
                f'({len(best["crafts"])} craft, {best["covered"]}/{NT} cov, sumJi {best["sumJi"]:.3f})  '
                f'[{time.time()-t0:.0f} s]')
    json.dump(hist, open(out / 'history.json', 'w'), indent=1)
    say(f'BEST J {best["J"]:.4f}  sumJi {best["sumJi"]:.4f}  covered {best["covered"]}/{NT}  '
        f'craft {len(best["crafts"])}  depths ' + ' '.join(str(len(S)) for _, S in best['crafts']))
    say(f'      tanks ' + ' '.join(f'{t:.0f}' for t in best['tanks']))


# ----------------------------------------------------------------------------- real build (validation)
def _build(job):
    """Settle carrier c's DP base restricted to S, then insert the rest of S nearest-first, for real."""
    import run_ialns as RI, carrier_dp as CD
    from ctoc14.globalopt import insert_block
    k, c, S, tmax, maxstep, maxtank, bases = job
    dat = data(); lib = dat.lib[c]
    want = set(int(dat.ids[i]) for i in S)
    t0 = time.time()
    full = [tuple(b) for b in bases] if bases else [b for b in lib['base'] if b[0] in want]
    full = [b for b in full if b[0] in want]
    ip = None
    for drop in [None] + list(range(len(full))):
        bs = [b for i, b in enumerate(full) if i != drop]
        if len(bs) < 4: break
        try:
            q = CD.seed_ip(RI.eph(), lib['tL'], np.array(lib['vinf']),
                           np.array([b[1] for b in bs]), np.array([b[2] for b in bs]))
            if insert_block(q, [(b[0], b[1]) for b in bs], log=RI.QUIET).max() > 1e4: continue
            if RI.settle(q, 100) <= 150: ip = q; break
        except Exception:
            continue
    if ip is None: return k, None, 'base failed'
    st = RI.ist(ip); tank = float(ip.tank())
    nfail = 0
    while time.time() - t0 < tmax:
        have = set(int(x) for x in st['asts'])
        left = [t for t in want if t not in have]
        if not left: break
        cands = RI.w_cands(('h', st, left, 0.35))
        if not cands: break
        tf = np.asarray(st['tf']); best = {}
        for x in sorted(cands, key=lambda x: x['dist']):
            if np.abs(tf - x['t']).min() > 12 * DAY: best.setdefault(x['ast'], x)
        trial = sorted(best.values(), key=lambda x: x['dist'])[:6]
        if not trial: break
        got = False
        for x in trial:
            r = RI.w_insert(('h', st, x['ast'], x['t']))
            if r.get('ok') and r['tank'] <= maxtank and r['tank'] - tank <= maxstep:
                st, tank = r['st'], r['tank']; got = True; break
        if not got:
            nfail += 1
            if nfail >= 3: break
        else:
            nfail = 0
    return k, (st, tank, time.time() - t0), None


def cmd_build(a):
    import run_ialns as RI
    dat = data()
    P = json.load(open(a.plan))
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    say = RI.logger(out / 'log.txt')
    say(f'plan {a.plan}: surrogate J {P["J"]:.4f} sumJi {P["sumJi"]:.4f} covered {P["covered"]} '
        f'({len(P["crafts"])} craft)')
    jobs = []
    BS = P.get('bases') or [None] * len(P['crafts'])
    for k, (c, ids) in enumerate(P['crafts']):
        S = [dat.idx[int(t)] for t in ids if int(t) in dat.idx]
        jobs.append((k, int(c), S, a.tmax, a.max_step, a.max_tank, BS[k]))
    FL = RI.IFleet(); t0 = time.time()
    with mp.get_context('fork').Pool(a.nproc) as pl:
        for k, r, why in pl.imap_unordered(_build, jobs):
            if r is None:
                say(f'  craft {k}: FAILED ({why})'); continue
            st, tank, dt = r
            np.savez(out / f'route_{k:02d}.npz', **st)
            FL.routes[f'{k:02d}'] = dict(st=st, tank=tank)
            want = len(P['crafts'][k][1])
            say(f'  craft {k}: {len(st["asts"]):2d}/{want} fb @ {tank:6.0f} kg  '
                f'{RI.cost(tank)/max(1,len(st["asts"])):.4f} J/fb  ({dt:.0f} s)')
    cov = FL.coverage()
    sumJi = sum(RI.cost(r['tank']) for r in FL.routes.values())
    J = sumJi + (NT - len(cov)) + 2.0
    say(f'SETTLED: {len(FL.routes)} craft, covered {len(cov)}/{NT}, sum J_i {sumJi:.4f}, J {J:.4f} '
        f'(surrogate said sumJi {P["sumJi"]:.4f}, J {P["J"]:.4f})   [{time.time()-t0:.0f} s]')
    FL.save(out / 'fleet', note=f'fleet_aco build of {a.plan}')
    json.dump(dict(settled_J=J, settled_sumJi=sumJi, covered=len(cov), n=len(FL.routes),
                   surrogate_J=P['J'], surrogate_sumJi=P['sumJi'], surrogate_covered=P['covered']),
              open(out / 'report.json', 'w'), indent=1)


# ----------------------------------------------------------------------------- diagnostics
def cmd_score(a):
    """Price REAL settled routes with the surrogate, giving it the same freedom the ACO has (it may pick the
    carrier).  This is the honesty test: the gap is exactly the optimism an ACO J carries."""
    import run_ialns as RI
    dat = data(); rm = RouteModel(load_model(a.model))
    rows = []; seen = set()
    for pat in a.dirs:
        for f in sorted(pathlib.Path(pat).rglob('route_*.npz')):
            try:
                z = np.load(f); st = {k: z[k] for k in z.files}
                if 'asts' not in st: continue
                S = [dat.idx[int(x)] for x in st['asts'] if int(x) in dat.idx]
                if len(S) < a.minfb: continue
                t = float(RI.ipr(st).tank()); key = (len(S), round(t, 1))
                if key in seen: continue
                seen.add(key)
                own = rm.tank(int(f.stem[7:]), S, dat)[0] if f.stem[7:].isdigit() else np.nan
                best = min(rm.tank(c, S, dat)[0] for c in range(dat.ncar))
                rows.append((len(S), t, best, own))
            except Exception:
                pass
    A = np.array([r[:3] for r in rows]); n, tr, pb = A[:, 0], A[:, 1], A[:, 2]
    print(f'{len(rows)} distinct real routes, depth {n.min():.0f}-{n.max():.0f}')
    print('depth band     n   pred-true [kg]   dJ_i median    dJ_i mean')
    for lo, hi in ((20, 26), (26, 30), (30, 34), (34, 38), (38, 42), (42, 60)):
        m = (n >= lo) & (n < hi)
        if m.sum() >= 3:
            dj = cost_J(pb[m]) - cost_J(tr[m])
            print(f'  {lo:2d}-{hi:<2d}      {m.sum():4d}      {np.median(pb[m]-tr[m]):+7.0f}      '
                  f'{np.median(dj):+.4f}      {np.mean(dj):+.4f}')


def cmd_diag(a):
    """How many craft can have a DP base of their own?  The library bases overlap heavily, so a disjoint packing
    runs out fast -- after that every further craft is built purely by insertion, which is the FORCED regime."""
    dat = data()
    R = np.ones(NT, bool); seq = []
    while True:
        n = (dat.isbase & R[None, :]).sum(1)
        for c, _ in seq: n[c] = 0
        c = int(n.argmax())
        if n[c] < 3: break
        seq.append((c, int(n[c]))); R[dat.isbase[c] & R] = False
    print('greedy disjoint packing of the 240 library DP bases:')
    print('  sizes      ', [s for _, s in seq])
    print('  cumulative ', np.cumsum([s for _, s in seq]).tolist())
    print(f'  union of all bases: {int(dat.isbase.any(0).sum())} of {NT} targets')
    mn = dat.D.min(0)
    print(f'\nper target, min over all {dat.ncar} carriers of the ballistic closest approach:')
    print('  percentiles(10,50,90,99) ' + ' '.join(f'{x:.3f}' for x in np.percentile(mn, [10, 50, 90, 99])) + ' AU')
    print(f'  every target is within {mn.max():.3f} AU of SOME carrier -> no target is structurally unreachable;'
          f'\n  the wall is combinatorial (median target has only {np.median((dat.D<0.05).sum(0)):.0f} carriers within 0.05 AU).')


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    c = sub.add_parser('calibrate'); c.add_argument('--out', default=''); c.set_defaults(fn=cmd_calibrate)
    r = sub.add_parser('run')
    r.add_argument('--out', required=True); r.add_argument('--model', default='')
    r.add_argument('--gens', type=int, default=40); r.add_argument('--ants', type=int, default=64)
    r.add_argument('--chunk', type=int, default=8); r.add_argument('--nproc', type=int, default=8)
    r.add_argument('--alpha', type=float, default=1.0); r.add_argument('--beta', type=float, default=2.0)
    r.add_argument('--rho', type=float, default=0.10); r.add_argument('--nelite', type=int, default=6)
    r.add_argument('--q0', type=float, default=0.5); r.add_argument('--greedy', action='store_true')
    r.add_argument('--tmin', type=float, default=0.05); r.add_argument('--tmax', type=float, default=20.0)
    r.add_argument('--nmax', type=int, default=12); r.add_argument('--nmin', type=int, default=8)
    r.add_argument('--depth', type=int, default=45, help='hard depth cap: the deepest route ever SETTLED is 48')
    r.add_argument('--maxtank', type=float, default=1200.0, help='the heaviest route ever settled is 1143 kg')
    r.add_argument('--dnear', type=float, default=0.08)
    r.add_argument('--breakeven', type=float, default=0.0369); r.add_argument('--seed', type=int, default=0)
    r.add_argument('--dp', action='store_true', help='re-run the exact carrier DP against the remaining targets '
                   'for every craft (0.05-0.11 s each) instead of reusing the library base')
    r.add_argument('--nbmin', type=int, default=5, help='a carrier whose base is smaller than this cannot fly')
    r.add_argument('--eps', type=float, default=0.01); r.add_argument('--lam', type=float, default=0.5)
    r.add_argument('--kgap', type=float, default=15.0)
    r.set_defaults(fn=cmd_run)
    b = sub.add_parser('build'); b.add_argument('plan'); b.add_argument('--out', required=True)
    b.add_argument('--nproc', type=int, default=8); b.add_argument('--tmax', type=float, default=1200.0)
    b.add_argument('--max-step', type=float, default=60.0); b.add_argument('--max-tank', type=float, default=1400.0)
    b.set_defaults(fn=cmd_build)
    s = sub.add_parser('score'); s.add_argument('dirs', nargs='+'); s.add_argument('--model', default='')
    s.add_argument('--minfb', type=int, default=20); s.set_defaults(fn=cmd_score)
    d = sub.add_parser('diag'); d.set_defaults(fn=cmd_diag)
    a = ap.parse_args(); a.fn(a)


if __name__ == '__main__':
    main()
