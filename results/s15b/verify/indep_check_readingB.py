#!/usr/bin/env python
"""Independent CTOC14 Problem A rules checker (written from CTOC14_problem.txt, no repo validator code).

Usage: indep_check.py RESULT_FILE [MEA.txt]
Only the Earth elements are read from ctoc14/constants.py (and cross-checked against Table 3 of the PDF).
"""
import os, sys, math
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
import numpy as np

# ---- constants (Table 2) ----
MU = 1.32712440018e11
AU = 149597870.7
G0 = 9.80665
ISP = 4000.0
TMAXN = 0.5
MDRY = 600.0
M0MAX = 2000.0
VINF = 4.0
DMAX = 1000.0
DTMIN = 8640.0
TWIN = 4.73364e8
T0_MJD = 62502.0
AST_EPH_MJD = 61200.0
EAR_EPH_MJD = 60676.0
VEX = ISP * G0            # m/s

# Earth elements: from ctoc14/constants.py, cross-checked against Table 3 literal values
sys.path.insert(0, "/Users/mickey/solarsystem/ctoc14")
from ctoc14.constants import EARTH_ELEMENTS as _EE  # noqa: E402
TABLE3 = np.array([1.0009175020, 0.017566762041, 0.002976847126, 189.953211282428, 273.196254000254, 357.4135031077])
assert np.array_equal(np.asarray(_EE, float), TABLE3), ("constants.py Earth elements differ from Table 3", _EE)
EARTH = TABLE3

violations = []
def viol(kind, msg):
    violations.append((kind, msg))


def kepE(M, e):
    M = np.mod(M, 2 * np.pi)
    E = np.where(e > 0.8, np.pi * np.ones_like(M), M + e * np.sin(M))
    for _ in range(60):
        f = E - e * np.sin(E) - M
        fp = 1 - e * np.cos(E)
        dE = -f / fp
        E = E + dE
        if np.all(np.abs(dE) < 1e-15):
            break
    return E


def elem_state(el, t_since_eph):
    """el: (...,6) a[AU] e i O w M0 [deg]; t_since_eph seconds. returns r (km), v (km/s)."""
    el = np.atleast_2d(el)
    a = el[:, 0] * AU
    e = el[:, 1]
    i, O, w, M0 = [np.radians(el[:, k]) for k in (2, 3, 4, 5)]
    n = np.sqrt(MU / a ** 3)
    M = M0 + n * t_since_eph
    E = kepE(M, e)
    cE, sE = np.cos(E), np.sin(E)
    b = a * np.sqrt(1 - e * e)
    xp = a * (cE - e)
    yp = b * sE
    rr = a * (1 - e * cE)
    Edot = n / (1 - e * cE)
    vxp = -a * sE * Edot
    vyp = b * cE * Edot
    cO, sO, cw, sw, ci, si = np.cos(O), np.sin(O), np.cos(w), np.sin(w), np.cos(i), np.sin(i)
    P = np.stack([cO * cw - sO * sw * ci, sO * cw + cO * sw * ci, sw * si], -1)
    Q = np.stack([-cO * sw - sO * cw * ci, -sO * sw + cO * cw * ci, cw * si], -1)
    r = xp[:, None] * P + yp[:, None] * Q
    v = vxp[:, None] * P + vyp[:, None] * Q
    return r, v


def kepler_prop(r0, v0, dt):
    """vectorised two-body propagation (elliptic), eccentric-anomaly-difference form."""
    r0n = np.linalg.norm(r0, axis=1)
    v02 = np.sum(v0 * v0, axis=1)
    a = 1.0 / (2.0 / r0n - v02 / MU)
    if np.any(a <= 0):
        raise RuntimeError("non-elliptic coast segment")
    sa = np.sqrt(a)
    n = np.sqrt(MU / a ** 3)
    sig0 = np.sum(r0 * v0, axis=1) / np.sqrt(MU)
    c1 = 1 - r0n / a
    c2 = sig0 / sa
    Mdt = n * dt
    dE = Mdt.copy()
    for _ in range(100):
        F = dE - c1 * np.sin(dE) + c2 * (1 - np.cos(dE)) - Mdt
        Fp = 1 - c1 * np.cos(dE) + c2 * np.sin(dE)
        step = -F / Fp
        step = np.clip(step, -1.0, 1.0)
        dE = dE + step
        if np.all(np.abs(step) < 1e-14):
            break
    cE, sE = np.cos(dE), np.sin(dE)
    f = 1 - a / r0n * (1 - cE)
    g = dt - (dE - sE) / n
    rv = f[:, None] * r0 + g[:, None] * v0
    rn = np.linalg.norm(rv, axis=1)
    fd = -np.sqrt(MU * a) / (rn * r0n) * sE
    gd = 1 - a / rn * (1 - cE)
    vv = fd[:, None] * r0 + gd[:, None] * v0
    return rv, vv


def main(path, mea_path):
    # ---------------- parse ----------------
    rows = []
    srcline = []
    with open(path, "r", encoding="utf-8") as fh:
        for ln, raw in enumerate(fh, 1):
            s = raw.strip()
            if not s or s.startswith("#"):
                continue
            tok = s.split()
            if len(tok) != 15:
                viol("format", f"file line {ln}: {len(tok)} columns (need 15)")
                continue
            try:
                ints = [int(tok[k]) for k in (0, 1, 2, 14)]
            except ValueError:
                viol("format", f"file line {ln}: non-integer in Line/SC_ID/Event/Asteroid_ID: {tok[0]} {tok[1]} {tok[2]} {tok[14]}")
                continue
            fl = [float(x) for x in tok[3:14]]
            if not all(math.isfinite(x) for x in fl):
                viol("format", f"file line {ln}: non-finite value")
            rows.append([ints[0], ints[1], ints[2]] + fl + [ints[3]])
            srcline.append(ln)
    A = np.array(rows, dtype=float)
    srcline = np.array(srcline)
    N = len(A)
    LINE = A[:, 0].astype(int)
    SC = A[:, 1].astype(int)
    EV = A[:, 2].astype(int)
    T = A[:, 3]
    R = A[:, 4:7]
    V = A[:, 7:10]
    M = A[:, 10]
    TH = A[:, 11:14]
    AID = A[:, 14].astype(int)
    print(f"parsed {N} data lines, all 15 columns: {len(rows)==N}")

    # line numbers
    if np.any(LINE <= 0):
        viol("format", "non-positive Line numbers")
    if not np.array_equal(LINE, np.arange(1, N + 1)):
        print("NOTE: Line column not exactly 1..N consecutive")
    else:
        print("Line column = 1..N consecutive")
    # events
    bad_ev = ~np.isin(EV, [0, 1, 2, 3, 4])
    for k in np.where(bad_ev)[0]:
        viol("format", f"file line {srcline[k]}: bad Event {EV[k]}")
    # SC ids
    if SC[0] != 1:
        viol("format", "first SC_ID != 1")
    dsc = np.diff(SC)
    if np.any(dsc < 0):
        viol("format", f"SC_ID not ascending at file lines {srcline[1:][dsc<0][:5]}")
    if np.any(dsc > 1):
        viol("format", f"SC_ID gaps at file lines {srcline[1:][dsc>1][:5]}")
    crafts = np.unique(SC)
    ncraft = len(crafts)
    if not np.array_equal(crafts, np.arange(1, ncraft + 1)):
        viol("format", f"SC_IDs not consecutive: {crafts}")

    # ---------------- per-craft structural checks ----------------
    # global
    for k in np.where((T < 0) | (T > TWIN))[0]:
        viol("window", f"file line {srcline[k]}: Time {T[k]!r} outside [0,{TWIN}]")
    for k in np.where(M < MDRY)[0]:
        viol("mass", f"file line {srcline[k]}: m={M[k]!r} < 600")
    nz = np.any(TH != 0.0, axis=1)
    for k in np.where(nz & (EV != 1))[0]:
        viol("thrust0", f"file line {srcline[k]}: Event {EV[k]} with nonzero thrust {TH[k]}")
    Tn = np.linalg.norm(TH, axis=1)
    for k in np.where((EV == 1) & (Tn > TMAXN))[0]:
        viol("thrust", f"file line {srcline[k]}: sample |T|={Tn[k]!r} > 0.5")
    for k in np.where((EV != 3) & (AID != 0))[0]:
        viol("format", f"file line {srcline[k]}: Event {EV[k]} with Asteroid_ID {AID[k]}")
    print(f"max sample |T| over Event=1 lines = {Tn[EV==1].max():.6f} N")

    starts = {}
    ends = {}
    arcs = []   # list of (craft, np.array of sample row indices)
    seg_thrust = np.zeros(N, dtype=int) - 1   # for row k (segment k->k+1): arc index if thrust-active
    min_gap_in_arc = np.inf
    min_gap_cross = np.inf
    cross_small = 0
    for c in crafts:
        idx = np.where(SC == c)[0]
        starts[c], ends[c] = idx[0], idx[-1]
        ev = EV[idx]
        if ev[0] != 0:
            viol("struct", f"SC {c}: first line Event={ev[0]} (file line {srcline[idx[0]]})")
        if ev[-1] != 4:
            viol("struct", f"SC {c}: last line Event={ev[-1]} (file line {srcline[idx[-1]]})")
        for k in idx[1:][ev[1:] == 0]:
            viol("struct", f"SC {c}: extra Event=0 at file line {srcline[k]}")
        for k in idx[:-1][ev[:-1] == 4]:
            viol("struct", f"SC {c}: Event=4 before end at file line {srcline[k]}")
        dt = np.diff(T[idx])
        for j in np.where(dt <= 0)[0]:
            viol("time", f"SC {c}: Time not strictly increasing at file line {srcline[idx[j+1]]} (dt={dt[j]!r})")
        dm = np.diff(M[idx])
        for j in np.where(dm > 0)[0]:
            viol("mass", f"SC {c}: mass increases by {dm[j]!r} kg at file line {srcline[idx[j+1]]}")
        # arcs
        cur = []
        last_sample_any = None
        for k in idx:
            e = EV[k]
            if e == 1:
                if last_sample_any is not None and not cur:
                    g = T[k] - T[last_sample_any]
                    min_gap_cross = min(min_gap_cross, g)
                    if g < DTMIN:
                        cross_small += 1
                cur.append(k)
                last_sample_any = k
            elif e == 3:
                pass
            else:
                if cur:
                    arcs.append((c, np.array(cur)))
                    cur = []
        if cur:
            arcs.append((c, np.array(cur)))
    for ai, (c, s) in enumerate(arcs):
        ts = T[s]
        if len(s) > 1:
            g = np.diff(ts)
            min_gap_in_arc = min(min_gap_in_arc, g.min())
            for j in np.where(g < DTMIN)[0]:
                viol("spacing", f"SC {c}: Event=1 gap {g[j]!r} s < 8640 at file line {srcline[s[j+1]]}")
        # mark thrust-active segments: rows k with T[k] >= t1 and T[k+1] <= tn, rows between s[0] and s[-1]
        for k in range(s[0], s[-1]):
            seg_thrust[k] = ai
    print(f"crafts={ncraft}, arcs={len(arcs)}, min Event=1 gap inside arcs={min_gap_in_arc:.3f} s, "
          f"min gap between last sample of an arc and first of next arc={min_gap_cross:.3f} s ({cross_small} < 8640)")

    # ---------------- launch ----------------
    print("\nLaunch checks:")
    m0 = {}
    for c in crafts:
        k = starts[c]
        te = T[k] + (T0_MJD - EAR_EPH_MJD) * 86400.0
        re, ve = elem_state(EARTH, np.array([te]))
        dr = np.linalg.norm(R[k] - re[0])
        vinf = np.linalg.norm(V[k] - ve[0])
        m0[c] = M[k]
        flag = ""
        if dr > 1.0:
            viol("launch", f"SC {c}: launch position off Earth by {dr:.6f} km")
            flag += " POS"
        if vinf > VINF + 0.01:
            viol("launch", f"SC {c}: v_inf {vinf:.6f} > 4.01")
            flag += " VINF"
        if not (MDRY <= M[k] <= M0MAX):
            viol("mass", f"SC {c}: m0 {M[k]} not in [600,2000]")
            flag += " M0"
        print(f"  SC {c}: t_launch={T[k]:.3f} s ({T[k]/86400:.4f} d)  |dr|={dr*1000:.3f} m  v_inf={vinf:.9f} km/s  m0={M[k]:.6f} kg{flag}")

    # ---------------- thrust interpolation: per-segment polynomial ----------------
    # For segment k (rows k,k+1) inside arc samples s (row idx), locate interval j (0-based) with
    # ts[j] <= T[k] and T[k+1] <= ts[j+1]; window l=clip(j-1,0,n-4) of 4 samples (n>=4) else all n.
    segk = np.where(seg_thrust >= 0)[0]
    # nodes for Gauss-Legendre / dense evaluation expressed through Lagrange basis on the window
    win_t = np.zeros((len(segk), 4))
    win_T = np.zeros((len(segk), 4, 3))
    win_n = np.zeros(len(segk), dtype=int)
    for q, k in enumerate(segk):
        c, s = arcs[seg_thrust[k]]
        ts = T[s]
        n = len(s)
        j = int(np.searchsorted(ts, T[k], side="right") - 1)
        j = min(max(j, 0), n - 2)
        assert ts[j] <= T[k] + 1e-9 and T[k + 1] <= ts[j + 1] + 1e-9, (k, j)
        if n >= 4:
            l = min(max(j, 0), n - 4)  # READING B
            w = s[l:l + 4]
        else:
            w = s
        win_n[q] = len(w)
        win_t[q, :len(w)] = T[w]
        win_T[q, :len(w)] = TH[w]

    def thrust_at(q_idx, t):
        """t: (len(q_idx), K) absolute times -> (len, K, 3) thrust in N."""
        out = np.zeros(t.shape + (3,))
        wt = win_t[q_idx]
        wT = win_T[q_idx]
        wn = win_n[q_idx]
        for nn in (1, 2, 3, 4):
            sel = wn == nn
            if not np.any(sel):
                continue
            tt = t[sel]
            ws = wt[sel, :nn]
            acc = np.zeros(tt.shape + (3,))
            for a in range(nn):
                L = np.ones_like(tt)
                for b in range(nn):
                    if b != a:
                        L *= (tt - ws[:, b:b + 1]) / (ws[:, a:a + 1] - ws[:, b:b + 1])
                acc += L[..., None] * wT[sel, a][:, None, :]
            out[sel] = acc
        return out

    # dense max |T| between samples
    K = 201
    fr = np.linspace(0.0, 1.0, K)
    maxT = 0.0
    argmax = None
    CH = 4000
    for st in range(0, len(segk), CH):
        qi = np.arange(st, min(st + CH, len(segk)))
        k = segk[qi]
        tt = T[k][:, None] + (T[k + 1] - T[k])[:, None] * fr[None, :]
        Tv = thrust_at(qi, tt)
        nrm = np.linalg.norm(Tv, axis=2)
        mx = nrm.max()
        if mx > maxT:
            maxT = mx
            qq, kk = np.unravel_index(np.argmax(nrm), nrm.shape)
            argmax = (srcline[k[qq]], tt[qq, kk])
        for q2 in np.where(nrm.max(axis=1) > TMAXN)[0]:
            viol("thrust", f"interpolated |T|={nrm[q2].max():.6f} > 0.5 between file lines {srcline[k[q2]]} and {srcline[k[q2]+1]}")
    print(f"\nmax interpolated |T| (201 pts/segment over {len(segk)} thrust segments) = {maxT:.6f} N at file line {argmax[0]} t={argmax[1]:.1f}")

    # ---------------- dynamics + mass integration (independent) ----------------
    print("\nDynamics consistency (own RK4 for thrust segments, analytic Kepler for coast segments):")
    allk = np.array([k for c in crafts for k in range(starts[c], ends[c])])
    coast = allk[seg_thrust[allk] < 0]
    thr = segk
    # coast: mass must be unchanged
    dmc = np.abs(M[coast + 1] - M[coast])
    rC, vC = kepler_prop(R[coast], V[coast], T[coast + 1] - T[coast])
    erC = np.linalg.norm(rC - R[coast + 1], axis=1)
    evC = np.linalg.norm(vC - V[coast + 1], axis=1)
    print(f"  coast segments {len(coast)}: max |dr|={erC.max():.4e} km  max |dv|={evC.max()*1000:.4e} m/s  max |dm|={dmc.max():.3e} kg"
          f"  longest={((T[coast+1]-T[coast]).max()/86400):.2f} d")
    for k in coast[erC > 1.0]:
        viol("dyn", f"coast file line {srcline[k]}->{srcline[k+1]}: |dr| {erC[coast==k][0]:.4f} km")
    for k in coast[evC > 1e-3]:
        viol("dyn", f"coast file line {srcline[k]}->{srcline[k+1]}: |dv| {evC[coast==k][0]*1000:.4f} m/s")
    for k in coast[dmc > 0.01]:
        viol("mass", f"coast file line {srcline[k]}->{srcline[k+1]}: |dm| {dmc[coast==k][0]:.4f} kg")

    # thrust: RK4 grouped by step count
    dts = T[thr + 1] - T[thr]
    nst = np.maximum(8, np.ceil(dts / 900.0)).astype(int)
    erT = np.zeros(len(thr)); evT = np.zeros(len(thr)); emT = np.zeros(len(thr))
    for ns in np.unique(nst):
        sel = np.where(nst == ns)[0]
        k = thr[sel]
        h = dts[sel] / ns
        r = R[k].copy(); v = V[k].copy(); m = M[k].copy(); t = T[k].copy()

        def f(t, r, v, m):
            Tv = thrust_at(sel, t[:, None])[:, 0, :]
            rn = np.linalg.norm(r, axis=1)
            acc = -MU * r / rn[:, None] ** 3 + Tv / m[:, None] / 1000.0
            md = -np.linalg.norm(Tv, axis=1) / VEX
            return v, acc, md

        for _ in range(ns):
            k1 = f(t, r, v, m)
            k2 = f(t + h / 2, r + h[:, None] / 2 * k1[0], v + h[:, None] / 2 * k1[1], m + h / 2 * k1[2])
            k3 = f(t + h / 2, r + h[:, None] / 2 * k2[0], v + h[:, None] / 2 * k2[1], m + h / 2 * k2[2])
            k4 = f(t + h, r + h[:, None] * k3[0], v + h[:, None] * k3[1], m + h * k3[2])
            r = r + h[:, None] / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            v = v + h[:, None] / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            m = m + h / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
            t = t + h
        erT[sel] = np.linalg.norm(r - R[k + 1], axis=1)
        evT[sel] = np.linalg.norm(v - V[k + 1], axis=1)
        emT[sel] = np.abs(m - M[k + 1])
    print(f"  thrust segments {len(thr)}: max |dr|={erT.max():.4e} km  max |dv|={evT.max()*1000:.4e} m/s  max |dm|={emT.max():.4e} kg"
          f"  longest={dts.max()/86400:.3f} d")
    for q in np.where(erT > 1.0)[0]:
        viol("dyn", f"thrust file line {srcline[thr[q]]}->{srcline[thr[q]+1]}: |dr| {erT[q]:.4f} km")
    for q in np.where(evT > 1e-3)[0]:
        viol("dyn", f"thrust file line {srcline[thr[q]]}->{srcline[thr[q]+1]}: |dv| {evT[q]*1000:.4f} m/s")
    for q in np.where(emT > 0.01)[0]:
        viol("mass", f"thrust file line {srcline[thr[q]]}->{srcline[thr[q]+1]}: |dm| {emT[q]:.5f} kg")

    if os.environ.get("IC_DEBUG"):
        for q in range(len(thr)):
            print(f"    DBG thrust seg file {srcline[thr[q]]}->{srcline[thr[q]+1]}: dr={erT[q]:.6e} km dv={evT[q]*1000:.6e} m/s dm={emT[q]:.6e} kg  m_end_file={M[thr[q]+1]!r}")
        for q in range(len(coast)):
            print(f"    DBG coast seg file {srcline[coast[q]]}->{srcline[coast[q]+1]}: dr={erC[q]:.6e} km dv={evC[q]*1000:.6e} m/s")
    for nm, arr, ks in (("dr km", erT, thr), ("dv km/s", evT, thr), ("dm kg", emT, thr), ("coast dr km", erC, coast), ("coast dv km/s", evC, coast)):
        top = np.argsort(arr)[-3:][::-1]
        print(f"  worst {nm}: " + ", ".join(f"{arr[q]:.4e} @file {srcline[ks[q]]}" for q in top))
    # min mass anywhere: mass is monotone so min is the last line per craft; also check final mass
    for c in crafts:
        print(f"  SC {c}: final mass {M[ends[c]]:.6f} kg, min over lines {M[starts[c]:ends[c]+1].min():.6f} kg")

    # ---------------- flybys ----------------
    print("\nFlyby checks:")
    mea = []
    with open(mea_path) as fh:
        for raw in fh:
            s = raw.strip()
            if not s or s.startswith("#"):
                continue
            tok = s.split()
            mea.append([float(x) for x in tok[:7]])
    mea = np.array(mea)
    assert len(mea) == 300 and np.array_equal(mea[:, 0], np.arange(1, 301)), "MEA ids"
    EL = mea[:, 1:7]
    fb = np.where(EV == 3)[0]
    ok_id = (AID[fb] >= 1) & (AID[fb] <= 300)
    for k in fb[~ok_id]:
        viol("flyby", f"file line {srcline[k]}: bad Asteroid_ID {AID[k]}")
    fbv = fb[ok_id]
    ra, va = elem_state(EL[AID[fbv] - 1], T[fbv] + (T0_MJD - AST_EPH_MJD) * 86400.0)
    d = np.linalg.norm(R[fbv] - ra, axis=1)
    vrel = np.linalg.norm(V[fbv] - va, axis=1)
    for q in np.where(d > DMAX)[0]:
        viol("flyby", f"file line {srcline[fbv[q]]}: SC {SC[fbv[q]]} ast {AID[fbv[q]]} distance {d[q]:.3f} km > 1000")
    print(f"  Event=3 lines: {len(fb)}, valid IDs: {ok_id.sum()}, max distance {d.max():.3f} km, "
          f"mean {d.mean():.3f} km, #d>1000: {(d>DMAX).sum()}, rel speed range {vrel.min():.3f}..{vrel.max():.3f} km/s")
    good = (d <= DMAX) & (T[fbv] >= 0) & (T[fbv] <= TWIN)
    order = np.argsort(T[fbv], kind="stable")
    covered = {}
    dup = 0
    for q in order:
        if not good[q]:
            continue
        a = AID[fbv[q]]
        if a in covered:
            dup += 1
        else:
            covered[a] = (SC[fbv[q]], T[fbv[q]])
    ncov = len(covered)
    missing = sorted(set(range(1, 301)) - set(covered))
    print(f"  distinct asteroids covered: {ncov}; repeat detections: {dup}; missing: {missing}")
    for c in crafts:
        sel = SC[fbv] == c
        firsts = sum(1 for a, (sc, _) in covered.items() if sc == c)
        print(f"  SC {c}: Event=3 lines {sel.sum()}, first-detections credited {firsts}, max d {d[sel].max():.3f} km")

    # ---------------- cost ----------------
    Js = {c: 1 + (m0[c] - MDRY) / 1400.0 + ((m0[c] - MDRY) / 1400.0) ** 2 for c in crafts}
    sumJ = sum(Js.values())
    J = sumJ + (300 - ncov)
    print("\nCost:")
    for c in crafts:
        print(f"  SC {c}: m0={m0[c]:.6f}  J_i={Js[c]:.6f}")
    print(f"  N={ncraft}  sum J_i = {sumJ:.6f}  N_miss = {300-ncov}  J = {J:.6f}")

    print(f"\nVIOLATIONS: {len(violations)}")
    from collections import Counter
    print("  by kind:", dict(Counter(k for k, _ in violations)))
    for k, m in violations[:200]:
        print(f"  [{k}] {m}")
    return 0


if __name__ == "__main__":
    p = sys.argv[1]
    mea = sys.argv[2] if len(sys.argv) > 2 else "/Users/mickey/solarsystem/ctoc14/MEA.txt"
    main(p, mea)
