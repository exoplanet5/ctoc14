"""Thrust model of the CTOC14 validator: 3rd-order sliding-window Lagrange interpolation of thrust samples.

An arc has n samples at strictly increasing times ts[0..n-1] with vectors Ts[0..n-1] (N). For t in [ts[j], ts[j+1]]
(0-based j in [0, n-2]) the window start is l = clip(j-1, 0, n-4) and the cubic Lagrange polynomial through samples
l..l+3 is evaluated component-wise. If n < 4 all n samples are used (degree n-1). Outside [ts[0], ts[-1]] thrust is 0.
The interpolant is continuous at the samples (each window passes through its samples) but its derivative is not.
Mass flow uses the NORM of the interpolated vector: mdot = -|T(t)| / (Isp g0).
"""
import numpy as np

def lagrange_basis(tw, t):
    """Lagrange basis weights of nodes tw (k,) at times t (...,). Returns (..., k)."""
    t = np.asarray(t, dtype=float)[..., None]
    k = len(tw)
    L = np.ones(t.shape[:-1] + (k,))
    for m in range(k):
        for q in range(k):
            if q != m:
                L[..., m] *= (t[..., 0] - tw[q]) / (tw[m] - tw[q])
    return L

def window_start(j, n):
    """Window start index for interval j (0-based) of an arc with n samples."""
    if n < 4:
        return 0
    return int(min(max(j - 1, 0), n - 4))

def thrust_interp(ts, Ts, t):
    """Evaluate the arc thrust vector at times t (scalar or array). ts: (n,), Ts: (n,3)."""
    ts = np.asarray(ts, dtype=float); Ts = np.asarray(Ts, dtype=float)
    t = np.asarray(t, dtype=float); scalar = t.ndim == 0; t = np.atleast_1d(t)
    n = len(ts)
    out = np.zeros(t.shape + (3,))
    inside = (t >= ts[0]) & (t <= ts[-1])
    if n == 1:
        out[inside] = Ts[0]
    elif n > 1 and np.any(inside):
        ti = t[inside]
        j = np.clip(np.searchsorted(ts, ti, side='right') - 1, 0, n - 2)
        k = n if n < 4 else 4
        l = np.zeros_like(j) if n < 4 else np.clip(j - 1, 0, n - 4)
        # gather window nodes (m, k)
        idx = l[:, None] + np.arange(k)[None, :]
        tw = ts[idx]; Tw = Ts[idx]                      # (m,k), (m,k,3)
        L = np.ones((len(ti), k))
        for a in range(k):
            for b in range(k):
                if a != b:
                    L[:, a] *= (ti - tw[:, b]) / (tw[:, a] - tw[:, b])
        out[inside] = np.einsum('mk,mkc->mc', L, Tw)
    return out[0] if scalar else out

def arc_max_thrust(ts, Ts, n_sub=400):
    """Max |T(t)| over the arc on a dense grid (n_sub points per interval) plus samples. Returns (max, t_at_max)."""
    ts = np.asarray(ts, dtype=float)
    if len(ts) == 1:
        return float(np.linalg.norm(Ts[0])), float(ts[0])
    grids = [np.linspace(ts[j], ts[j + 1], n_sub, endpoint=False) for j in range(len(ts) - 1)]
    tg = np.concatenate(grids + [ts[-1:]])
    Tn = np.linalg.norm(thrust_interp(ts, Ts, tg), axis=-1)
    k = int(np.argmax(Tn))
    return float(Tn[k]), float(tg[k])
