#!/usr/bin/env python
"""Independent CTOC14 Problem A rules-compliance checker (lens: RULES), written for s16E.

Reads only the result file, MEA.txt and literal constants from CTOC14_problem.txt (Tables 2/3).
No repo code is imported.  Read-only on its inputs.

Usage: rules_indep_check.py RESULT_FILE [MEA.txt]
"""
import os, sys, math, re
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
from fractions import Fraction
from collections import Counter
import numpy as np

# ---------------- Table 2 / Table 3 literals ----------------
MU = 1.32712440018e11
AU = 149597870.7
G0 = 9.80665
ISP = 4000.0
VEX = ISP * G0                    # m/s  (T[N]/VEX[m/s] = kg/s)
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
EARTH = [1.0009175020, 0.017566762041, 0.002976847126, 189.953211282428, 273.196254000254, 357.4135031077]

V = []
def viol(kind, msg):
    V.append((kind, msg))

READING = os.environ.get("RC_READING", "A")
print(f"interpolation window reading: {READING}")
INT_RE = re.compile(r"^[+-]?\d+$")
FLT_RE = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$")


# ---------------- ephemerides: two independent Kepler solvers ----------------
def rot(i, O, w):
    cO, sO, cw, sw, ci, si = math.cos(O), math.sin(O), math.cos(w), math.sin(w), math.cos(i), math.sin(i)
    P = (cO * cw - sO * sw * ci, sO * cw + cO * sw * ci, sw * si)
    Q = (-cO * sw - sO * cw * ci, -sO * sw + cO * cw * ci, cw * si)
    return P, Q


def state_bisect(el, dt):
    """Scalar: Kepler equation by bisection (guaranteed bracket E in [M-e, M+e])."""
    a_au, e, i, O, w, M0 = el
    a = a_au * AU
    n = math.sqrt(MU / a ** 3)
    M = math.radians(M0) + n * dt
    M = math.fmod(M, 2 * math.pi)
    if M > math.pi:
        M -= 2 * math.pi
    if M < -math.pi:
        M += 2 * math.pi
    lo, hi = M - e - 1e-12, M + e + 1e-12
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid - e * math.sin(mid) - M > 0:
            hi = mid
        else:
            lo = mid
        if hi - lo < 1e-16:
            break
    E = 0.5 * (lo + hi)
    # one Newton polish
    E -= (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
    cE, sE = math.cos(E), math.sin(E)
    b = a * math.sqrt(1 - e * e)
    xp, yp = a * (cE - e), b * sE
    Ed = n / (1 - e * cE)
    vxp, vyp = -a * sE * Ed, b * cE * Ed
    P, Q = rot(math.radians(i), math.radians(O), math.radians(w))
    r = np.array([xp * P[k] + yp * Q[k] for k in range(3)])
    v = np.array([vxp * P[k] + vyp * Q[k] for k in range(3)])
    return r, v


def state_fg(el, dt):
    """Scalar: epoch state (M0) then universal-variable-free f&g propagation via Delta-E Newton (different path)."""
    a_au, e, i, O, w, M0 = el
    a = a_au * AU
    n = math.sqrt(MU / a ** 3)
    # epoch state from M0 (Newton with safe start)
    Mm = math.fmod(math.radians(M0), 2 * math.pi)
    E = math.pi if e > 0.8 else Mm
    for _ in range(100):
        d = (E - e * math.sin(E) - Mm) / (1 - e * math.cos(E))
        E -= d
        if abs(d) < 1e-15:
            break
    cE, sE = math.cos(E), math.sin(E)
    b = a * math.sqrt(1 - e * e)
    Ed = n / (1 - e * cE)
    P, Q = rot(math.radians(i), math.radians(O), math.radians(w))
    r0 = np.array([a * (cE - e) * P[k] + b * sE * Q[k] for k in range(3)])
    v0 = np.array([-a * sE * Ed * P[k] + b * cE * Ed * Q[k] for k in range(3)])
    # propagate by dt: solve dE - c1 sin dE + c2 (1-cos dE) = n dt  (mod 2pi)
    r0n = np.linalg.norm(r0)
    c1 = 1 - r0n / a
    c2 = float(np.dot(r0, v0)) / math.sqrt(MU * a)
    Mdt = math.fmod(n * dt, 2 * math.pi)
    dE = Mdt
    for _ in range(200):
        F = dE - c1 * math.sin(dE) + c2 * (1 - math.cos(dE)) - Mdt
        Fp = 1 - c1 * math.cos(dE) + c2 * math.sin(dE)
        st = max(-0.5, min(0.5, -F / Fp))
        dE += st
        if abs(st) < 1e-15:
            break
    f = 1 - a / r0n * (1 - math.cos(dE))
    g = Mdt / n - (dE - math.sin(dE)) / n
    r = f * r0 + g * v0
    rn = np.linalg.norm(r)
    fd = -math.sqrt(MU * a) / (rn * r0n) * math.sin(dE)
    gd = 1 - a / rn * (1 - math.cos(dE))
    v = fd * r0 + gd * v0
    return r, v


def kepler_prop_vec(r0, v0, dt):
    r0n = np.linalg.norm(r0, axis=1)
    a = 1.0 / (2.0 / r0n - np.sum(v0 * v0, axis=1) / MU)
    assert np.all(a > 0)
    n = np.sqrt(MU / a ** 3)
    c1 = 1 - r0n / a
    c2 = np.sum(r0 * v0, axis=1) / np.sqrt(MU * a)
    Mdt = n * dt
    dE = Mdt.copy()
    for _ in range(200):
        F = dE - c1 * np.sin(dE) + c2 * (1 - np.cos(dE)) - Mdt
        Fp = 1 - c1 * np.cos(dE) + c2 * np.sin(dE)
        st = np.clip(-F / Fp, -0.5, 0.5)
        dE += st
        if np.all(np.abs(st) < 1e-15):
            break
    f = 1 - a / r0n * (1 - np.cos(dE))
    g = dt - (dE - np.sin(dE)) / n
    r = f[:, None] * r0 + g[:, None] * v0
    rn = np.linalg.norm(r, axis=1)
    fd = -np.sqrt(MU * a) / (rn * r0n) * np.sin(dE)
    gd = 1 - a / rn * (1 - np.cos(dE))
    return r, fd[:, None] * r0 + gd[:, None] * v0


def main(path, mea_path):
    # ================= parse (strict) =================
    rows, src, raw_m0 = [], [], {}
    ntok = Counter()
    with open(path, "rb") as fh:
        data = fh.read()
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as ex:
        viol("format", f"not UTF-8: {ex}")
        text = data.decode("latin-1")
    for ln, raw in enumerate(text.splitlines(), 1):
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        tok = s.split()
        ntok[len(tok)] += 1
        if len(tok) != 15:
            viol("format", f"file line {ln}: {len(tok)} columns")
            continue
        bad = [k for k in (0, 1, 2, 14) if not INT_RE.match(tok[k])]
        bad += [k for k in range(3, 14) if not FLT_RE.match(tok[k])]
        if bad:
            viol("format", f"file line {ln}: malformed tokens at cols {[b+1 for b in bad]}")
            continue
        vals = [float(x) for x in tok[3:14]]
        if not all(math.isfinite(x) for x in vals):
            viol("format", f"file line {ln}: non-finite")
        rows.append([int(tok[0]), int(tok[1]), int(tok[2])] + vals + [int(tok[14])])
        src.append(ln)
        if int(tok[2]) == 0:
            raw_m0[int(tok[1])] = tok[10]
    print(f"token-count histogram over data lines: {dict(ntok)}")
    A = np.array(rows, dtype=float)
    src = np.array(src)
    N = len(A)
    LINE = A[:, 0].astype(np.int64); SC = A[:, 1].astype(int); EV = A[:, 2].astype(int)
    T = A[:, 3]; R = A[:, 4:7]; Vv = A[:, 7:10]; M = A[:, 10]; TH = A[:, 11:14]; AID = A[:, 14].astype(int)
    print(f"data lines {N}; Line col == 1..N: {np.array_equal(LINE, np.arange(1, N+1))}; Line>0: {bool(np.all(LINE>0))}")
    if np.any(LINE <= 0):
        viol("format", "Line <= 0")
    for k in np.where(~np.isin(EV, [0, 1, 2, 3, 4]))[0]:
        viol("format", f"file line {src[k]}: Event {EV[k]}")
    crafts = np.unique(SC)
    print(f"SC_IDs: {crafts.tolist()}  sorted ascending: {bool(np.all(np.diff(SC) >= 0))}")
    if not np.array_equal(crafts, np.arange(1, len(crafts) + 1)):
        viol("format", "SC_ID not consecutive from 1")
    if np.any(np.diff(SC) < 0):
        viol("format", "SC_ID not sorted ascending")
    print(f"Event histogram: {dict(sorted(Counter(EV.tolist()).items()))}")

    # ================= global line checks =================
    print(f"Time range: [{T.min()!r}, {T.max()!r}] (window [0, {TWIN!r}])")
    for k in np.where((T < 0) | (T > TWIN))[0]:
        viol("window", f"file line {src[k]}: Time {T[k]!r}")
    for k in np.where(M < MDRY)[0]:
        viol("mass", f"file line {src[k]}: m {M[k]!r} < 600")
    nz = np.any(TH != 0.0, axis=1)
    print(f"non-Event1 lines with any nonzero thrust component: {int(np.sum(nz & (EV != 1)))}")
    for k in np.where(nz & (EV != 1))[0]:
        viol("thrust0", f"file line {src[k]}: Event {EV[k]} thrust {TH[k]}")
    for k in np.where((EV != 3) & (AID != 0))[0]:
        viol("format", f"file line {src[k]}: Event {EV[k]} Asteroid_ID {AID[k]}")
    Tn = np.linalg.norm(TH, axis=1)
    print(f"max sample |T| on Event=1 lines: {Tn[EV==1].max()!r} N; count > 0.5: {int(np.sum((EV==1)&(Tn>TMAXN)))}")
    for k in np.where((EV == 1) & (Tn > TMAXN))[0]:
        viol("thrust", f"file line {src[k]}: sample |T| {Tn[k]!r}")

    # ================= per craft structure / arcs =================
    starts, ends, arcs = {}, {}, []
    segarc = -np.ones(N, dtype=int)
    min_in_arc, min_cross, n_cross_small = np.inf, np.inf, 0
    for c in crafts:
        idx = np.where(SC == c)[0]
        if not np.all(np.diff(idx) == 1):
            viol("format", f"SC {c}: lines not contiguous")
        starts[c], ends[c] = idx[0], idx[-1]
        ev = EV[idx]
        if ev[0] != 0:
            viol("struct", f"SC {c}: first Event {ev[0]}")
        if ev[-1] != 4:
            viol("struct", f"SC {c}: last Event {ev[-1]}")
        if np.sum(ev == 0) != 1:
            viol("struct", f"SC {c}: {np.sum(ev==0)} Event=0 lines")
        if np.sum(ev == 4) != 1:
            viol("struct", f"SC {c}: {np.sum(ev==4)} Event=4 lines")
        dt = np.diff(T[idx])
        for j in np.where(~(dt > 0))[0]:
            viol("time", f"SC {c}: non-increasing Time at file line {src[idx[j+1]]}")
        dm = np.diff(M[idx])
        for j in np.where(dm > 0)[0]:
            viol("mass", f"SC {c}: mass increases {dm[j]!r} at file line {src[idx[j+1]]}")
        cur, last1 = [], None
        for k in idx:
            e = EV[k]
            if e == 1:
                if not cur and last1 is not None:
                    g = T[k] - T[last1]
                    min_cross = min(min_cross, g)
                    n_cross_small += g < DTMIN
                cur.append(k); last1 = k
            elif e == 3:
                continue
            else:     # 0, 2, 4 separate arcs
                if cur:
                    arcs.append((c, np.array(cur))); cur = []
        if cur:
            arcs.append((c, np.array(cur)))
            viol("struct", f"SC {c}: arc not closed by 0/2/4")
    arc_n = Counter()
    for ai, (c, s) in enumerate(arcs):
        arc_n[c] += 1
        if len(s) > 1:
            g = np.diff(T[s])
            min_in_arc = min(min_in_arc, g.min())
            for j in np.where(g < DTMIN)[0]:
                viol("spacing", f"SC {c}: Event=1 spacing {g[j]!r} < 8640 at file line {src[s[j+1]]}")
        segarc[s[0]:s[-1]] = ai
    ns_hist = Counter(len(s) for _, s in arcs)
    print(f"arcs: {len(arcs)} (per craft {dict(arc_n)}); samples/arc min {min(len(s) for _,s in arcs)} max {max(len(s) for _,s in arcs)}")
    print(f"min Event=1 spacing inside arcs: {min_in_arc!r} s ; min gap last-sample->first-sample across an arc boundary: {min_cross!r} s ({n_cross_small} < 8640)")
    # all distinct Event=1 spacings inside arcs
    gaps = np.concatenate([np.diff(T[s]) for _, s in arcs if len(s) > 1])
    print(f"Event=1 spacing inside arcs: distinct values (rounded 1e-6) {sorted(set(np.round(gaps,6).tolist()))[:10]}")

    # ================= launch =================
    print("\nLaunch:")
    m0 = {}
    for c in crafts:
        k = starts[c]
        dte = T[k] + (T0_MJD - EAR_EPH_MJD) * 86400.0
        re1, ve1 = state_bisect(EARTH, dte)
        re2, ve2 = state_fg(EARTH, dte)
        dr = np.linalg.norm(R[k] - re1); vinf = np.linalg.norm(Vv[k] - ve1)
        dr2 = np.linalg.norm(R[k] - re2); vinf2 = np.linalg.norm(Vv[k] - ve2)
        m0[c] = M[k]
        fl = []
        if T[k] < 0: fl.append("T<0")
        if max(dr, dr2) > 1.0: fl.append("POS>1km"); viol("launch", f"SC {c}: |dr| {dr} km")
        if max(vinf, vinf2) > VINF: fl.append("VINF>4.000")
        if max(vinf, vinf2) > VINF + 0.01: viol("launch", f"SC {c}: v_inf {vinf}")
        if not (MDRY <= M[k] <= M0MAX): fl.append("M0"); viol("mass", f"SC {c}: m0 {M[k]}")
        print(f"  SC {c}: t={T[k]:.1f}s |dr|={dr*1e3:.4f} m (alt {dr2*1e3:.4f} m) v_inf={vinf:.12f} (alt {vinf2:.12f}) km/s m0={raw_m0[c]} {' '.join(fl) or 'ok'}")

    # ================= interpolated thrust: exact max =================
    segk = np.where(segarc >= 0)[0]
    W_t = np.zeros((len(segk), 4)); W_T = np.zeros((len(segk), 4, 3)); W_n = np.zeros(len(segk), dtype=int)
    for q, k in enumerate(segk):
        c, s = arcs[segarc[k]]
        ts = T[s]; n = len(s)
        j = int(np.searchsorted(ts, T[k], side="right") - 1)
        j = min(max(j, 0), n - 2)
        if not (ts[j] <= T[k] and T[k + 1] <= ts[j + 1]):
            viol("internal", f"segment {src[k]} not inside one sample interval")
        l0 = min(max(j - (0 if READING == 'B' else 1), 0), n - 4)
        w = s[l0:l0 + 4] if n >= 4 else s
        W_n[q] = len(w); W_t[q, :len(w)] = T[w]; W_T[q, :len(w)] = TH[w]

    def interp(qi, t):
        """qi: (Q,) segment indices; t: (Q,K) -> (Q,K,3)"""
        out = np.zeros(t.shape + (3,))
        for nn in (1, 2, 3, 4):
            sel = W_n[qi] == nn
            if not np.any(sel):
                continue
            tt = t[sel]; ws = W_t[qi][sel, :nn]; wv = W_T[qi][sel]
            acc = np.zeros(tt.shape + (3,))
            for a in range(nn):
                L = np.ones_like(tt)
                for b in range(nn):
                    if b != a:
                        L = L * (tt - ws[:, b:b+1]) / (ws[:, a:a+1] - ws[:, b:b+1])
                acc += L[..., None] * wv[:, a][:, None, :]
            out[sel] = acc
        return out

    # cubic coefficients in s = (t - t_k)/h on [0,1]: fit through 4 points (exact for deg <= 3)
    qi_all = np.arange(len(segk))
    h = T[segk + 1] - T[segk]
    sn = np.array([0.0, 1/3, 2/3, 1.0])
    vals = interp(qi_all, T[segk][:, None] + h[:, None] * sn[None, :])      # (Q,4,3)
    Vd = np.vander(sn, 4, increasing=True)                                    # (4,4)
    coef = np.einsum("ij,qjc->qic", np.linalg.inv(Vd), vals)                   # (Q,4,3) c0..c3
    # |T|^2 polynomial (deg 6) coefficients, then derivative (deg 5)
    P2 = np.zeros((len(segk), 7))
    for a in range(4):
        for b in range(4):
            P2[:, a + b] += np.sum(coef[:, a, :] * coef[:, b, :], axis=1)
    dP = P2[:, 1:] * np.arange(1, 7)[None, :]
    best = np.maximum(P2[:, 0], P2.sum(axis=1))   # endpoints s=0, s=1
    for q in range(len(segk)):
        d = dP[q]
        if not np.any(d):
            continue
        rts = np.roots(d[::-1])
        rr = rts.real[(np.abs(rts.imag) < 1e-7) & (rts.real > 0) & (rts.real < 1)]
        if rr.size:
            best[q] = max(best[q], np.polyval(P2[q][::-1], rr).max())
    Tmax_exact = np.sqrt(np.maximum(best, 0))
    # cross-check with dense sampling
    fr = np.linspace(0, 1, 401)
    dense = 0.0
    for st in range(0, len(segk), 3000):
        qi = np.arange(st, min(st + 3000, len(segk)))
        dense = max(dense, np.linalg.norm(interp(qi, T[segk[qi]][:, None] + h[qi][:, None] * fr[None, :]), axis=2).max())
    qm = int(np.argmax(Tmax_exact))
    print(f"\nthrust segments {len(segk)}; exact max interpolated |T| = {Tmax_exact.max():.9f} N at file line {src[segk[qm]]} "
          f"(window n={W_n[qm]}); dense-401 max = {dense:.9f} N; #segments > 0.5: {int(np.sum(Tmax_exact > TMAXN))}")
    for q in np.where(Tmax_exact > TMAXN)[0]:
        viol("thrust", f"interp |T| {Tmax_exact[q]:.6f} in segment file {src[segk[q]]}->{src[segk[q]+1]}")
    top = np.argsort(Tmax_exact)[-5:][::-1]
    print("  top-5 segments: " + ", ".join(f"{Tmax_exact[q]:.6f}@{src[segk[q]]}(SC{SC[segk[q]]})" for q in top))
    # overshoot vs sample: max ratio interpolated / max(adjacent samples)
    # ================= dynamics (own RK4 / Kepler) =================
    allk = np.concatenate([np.arange(starts[c], ends[c]) for c in crafts])
    coast = allk[segarc[allk] < 0]
    rC, vC = kepler_prop_vec(R[coast], Vv[coast], T[coast + 1] - T[coast])
    erC = np.linalg.norm(rC - R[coast + 1], axis=1); evC = np.linalg.norm(vC - Vv[coast + 1], axis=1)
    dmC = np.abs(M[coast + 1] - M[coast])
    print(f"\ncoast segments {len(coast)}: max|dr| {erC.max():.3e} km, max|dv| {evC.max()*1e3:.3e} m/s, max|dm| {dmC.max():.3e} kg")
    for arr, lim, nm in ((erC, 1.0, "dr"), (evC, 1e-3, "dv"), (dmC, 0.01, "dm")):
        for q in np.where(arr > lim)[0]:
            viol("dyn", f"coast {src[coast[q]]}: {nm} {arr[q]}")
    dts = h
    nst = np.maximum(8, np.ceil(dts / 600.0)).astype(int)
    erT = np.zeros(len(segk)); evT = np.zeros(len(segk)); emT = np.zeros(len(segk)); mEnd = np.zeros(len(segk))
    for nsv in np.unique(nst):
        sel = np.where(nst == nsv)[0]
        k = segk[sel]; hh = dts[sel] / nsv
        r = R[k].copy(); v = Vv[k].copy(); m = M[k].copy(); t = T[k].copy()
        def f(t, r, v, m):
            Tv = interp(sel, t[:, None])[:, 0, :]
            rn = np.linalg.norm(r, axis=1)
            return v, -MU * r / rn[:, None] ** 3 + Tv / m[:, None] / 1000.0, -np.linalg.norm(Tv, axis=1) / VEX
        for _ in range(nsv):
            a1 = f(t, r, v, m)
            a2 = f(t + hh/2, r + hh[:, None]/2 * a1[0], v + hh[:, None]/2 * a1[1], m + hh/2 * a1[2])
            a3 = f(t + hh/2, r + hh[:, None]/2 * a2[0], v + hh[:, None]/2 * a2[1], m + hh/2 * a2[2])
            a4 = f(t + hh, r + hh[:, None] * a3[0], v + hh[:, None] * a3[1], m + hh * a3[2])
            r = r + hh[:, None]/6 * (a1[0] + 2*a2[0] + 2*a3[0] + a4[0])
            v = v + hh[:, None]/6 * (a1[1] + 2*a2[1] + 2*a3[1] + a4[1])
            m = m + hh/6 * (a1[2] + 2*a2[2] + 2*a3[2] + a4[2])
            t = t + hh
        erT[sel] = np.linalg.norm(r - R[k + 1], axis=1); evT[sel] = np.linalg.norm(v - Vv[k + 1], axis=1)
        emT[sel] = np.abs(m - M[k + 1]); mEnd[sel] = m
    print(f"thrust segments {len(segk)}: max|dr| {erT.max():.3e} km, max|dv| {evT.max()*1e3:.3e} m/s, max|dm| {emT.max():.3e} kg; "
          f"longest {dts.max()/86400:.4f} d; min integrated end mass {mEnd.min():.6f} kg")
    for arr, lim, nm in ((erT, 1.0, "dr"), (evT, 1e-3, "dv"), (emT, 0.01, "dm")):
        for q in np.where(arr > lim)[0]:
            viol("dyn", f"thrust {src[segk[q]]}: {nm} {arr[q]}")
    if mEnd.min() < MDRY:
        viol("mass", f"integrated mass {mEnd.min()} < 600")
    for c in crafts:
        print(f"  SC {c}: lines {ends[c]-starts[c]+1}, t_end {T[ends[c]]:.1f} s, final mass {M[ends[c]]!r}")

    # ================= flybys =================
    mea = []
    for raw in open(mea_path):
        s = raw.strip()
        if s and not s.startswith("#"):
            mea.append([float(x) for x in s.split()[:7]])
    mea = np.array(mea)
    assert len(mea) == 300 and np.array_equal(mea[:, 0], np.arange(1, 301))
    fb = np.where(EV == 3)[0]
    print(f"\nEvent=3 lines: {len(fb)}")
    d1 = np.full(len(fb), np.inf); d2 = np.full(len(fb), np.inf); vr = np.zeros(len(fb))
    for q, k in enumerate(fb):
        a = AID[k]
        if not 1 <= a <= 300:
            viol("flyby", f"file line {src[k]}: Asteroid_ID {a}")
            continue
        dt = T[k] + (T0_MJD - AST_EPH_MJD) * 86400.0
        ra, va = state_bisect(mea[a - 1, 1:7], dt)
        rb, vb = state_fg(mea[a - 1, 1:7], dt)
        d1[q] = np.linalg.norm(R[k] - ra); d2[q] = np.linalg.norm(R[k] - rb); vr[q] = np.linalg.norm(Vv[k] - va)
    dd = np.maximum(d1, d2)
    print(f"  distance max {dd.max():.6f} km (bisection {d1.max():.6f}, f&g {d2.max():.6f}); max |d1-d2| {np.max(np.abs(d1-d2)):.3e} km; "
          f"#>1000: {int(np.sum(dd>DMAX))}; #>900: {int(np.sum(dd>900))}; rel-speed {vr.min():.3f}..{vr.max():.3f} km/s")
    for q in np.where(dd > DMAX)[0]:
        viol("flyby", f"file line {src[fb[q]]}: ast {AID[fb[q]]} d {dd[q]:.3f}")
    order = np.argsort(T[fb], kind="stable")
    cov = {}
    for q in order:
        k = fb[q]
        if dd[q] <= DMAX and 0 <= T[k] <= TWIN and AID[k] not in cov:
            cov[AID[k]] = (SC[k], T[k])
    rep = len(fb) - len(cov)
    missing = sorted(set(range(1, 301)) - set(cov))
    print(f"  distinct covered {len(cov)}; repeat/uncounted Event=3 lines {rep}; missing {missing}")
    for c in crafts:
        sel = SC[fb] == c
        print(f"  SC {c}: Event=3 {int(sel.sum())}, first-detections {sum(1 for x in cov.values() if x[0]==c)}, max d {dd[sel].max():.3f} km")

    # ================= cost =================
    sumJ = 0.0; sumF = Fraction(0)
    for c in crafts:
        x = (m0[c] - MDRY) / 1400.0
        sumJ += 1 + x + x * x
        xf = (Fraction(raw_m0[c]) - 600) / 1400
        sumF += 1 + xf + xf * xf
    nmiss = 300 - len(cov)
    print(f"\nN={len(crafts)} sum J_i = {sumJ:.9f} (exact-decimal {float(sumF):.12f}); N_miss={nmiss}; J = {sumJ+nmiss:.9f} (exact {float(sumF)+nmiss:.12f})")
    print(f"\nVIOLATIONS: {len(V)}  by kind {dict(Counter(k for k,_ in V))}")
    for k, m in V[:100]:
        print(f"  [{k}] {m}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "/Users/mickey/solarsystem/ctoc14/MEA.txt")
