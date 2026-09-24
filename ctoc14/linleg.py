"""Linearised continuous-thrust leg model (position-only targeting).

For a spacecraft state (r, v) at time t coasting for TOF, the position at t+TOF responds linearly to small thrust
accelerations u_k applied over segments of length dtau centred at tau_k:  dr(TOF) = sum_k Phi_rv(TOF, tau_k) u_k dtau.
The minimum-fuel profile that moves the coast end-point onto a target position is approximated by IRLS on the L1 norm
(sum_k |u_k| dtau), with feasibility |u_k| <= ub * a (a = Tmax/m).  Arrival velocity follows from Phi_vv.
This replaces the Lambert-junction (impulsive) leg model of search.py, which is pessimistic for long legs
(> ~250 d: single-revolution Lambert) and optimistic for large junctions (leverage loss)."""
import numpy as np
from .kepler import propagate_twobody, propagate_batch
from .constants import DAY

class LinLeg:
    def __init__(self, r, v, tofs, dtau=10 * DAY, eps=1e-4):
        self.r = np.asarray(r, float); self.v = np.asarray(v, float); self.tofs = np.asarray(tofs, float); self.dtau = dtau
        nK = int(np.ceil(self.tofs.max() / dtau)); self.taus = (np.arange(nK) + 0.5) * dtau
        self.rc, self.vc = propagate_twobody(self.r, self.v, self.tofs)          # nominal coast (nJ,3)
        rk, vk = propagate_twobody(self.r, self.v, self.taus)
        nJ = len(self.tofs); self.Prv = np.zeros((nJ, nK, 3, 3)); self.Pvv = np.zeros((nJ, nK, 3, 3))
        # finite-difference sensitivities of every (segment k, axis i) impulse on every later end time j, in ONE batched
        # propagation (same result as the per-(k, i) loop with propagate_twobody, to rounding)
        J, Kk = np.nonzero(self.tofs[:, None] > self.taus[None, :])          # all (j, k) pairs with tau_k < tof_j
        if len(J):
            rows_j = np.repeat(J, 3); rows_k = np.repeat(Kk, 3); rows_i = np.tile(np.arange(3), len(J))
            V0 = vk[rows_k].copy(); V0[np.arange(len(rows_i)), rows_i] += eps
            rp, vp = propagate_batch(rk[rows_k], V0, self.tofs[rows_j] - self.taus[rows_k])
            self.Prv[rows_j, rows_k, :, rows_i] = (rp - self.rc[rows_j]) / eps * dtau
            self.Pvv[rows_j, rows_k, :, rows_i] = (vp - self.vc[rows_j]) / eps * dtau
        self.K = np.array([int(np.sum(self.taus < tof)) for tof in self.tofs])

    def solve(self, j, R, a, ub=0.8, iters=5, return_U=False):
        """Targets at positions R (T,3) at time t+tofs[j]. Returns cost dv [km/s] (inf if infeasible), arrival velocity (T,3),
        and the thrust profile norms (T,K).  a = Tmax/m in km/s^2."""
        K = self.K[j]; Mk = self.Prv[j, :K]; dR = R - self.rc[j]; T = len(R)
        if K < 1:
            return np.full(T, np.inf), np.broadcast_to(self.vc[j], (T, 3)), np.zeros((T, 0))
        MMk = np.einsum('kij,klj->kil', Mk, Mk)                    # M_k M_k^T  (K,3,3)
        w = np.ones((T, K)); best_cost = np.full(T, np.inf); best_U = np.zeros((T, K, 3))
        for it in range(iters):
            G = np.einsum('tk,kij->tij', 1 / w, MMk)
            try:
                y = np.linalg.solve(G, dR[..., None])[..., 0]
            except np.linalg.LinAlgError:
                y = np.einsum('tij,tj->ti', np.linalg.pinv(G), dR)
            U = np.einsum('tk,kji,tj->tki', 1 / w, Mk, y)
            un = np.linalg.norm(U, axis=-1)
            cost = un.sum(1) * self.dtau; feas = un.max(1) <= ub * a
            better = feas & (cost < best_cost)
            best_cost[better] = cost[better]; best_U[better] = U[better]
            w = np.maximum(un, 1e-9)
        dv_arr = np.einsum('kij,tkj->ti', self.Pvv[j, :K], best_U)      # arrival velocity change
        if return_U:
            return best_cost, self.vc[j] + dv_arr, np.linalg.norm(best_U, axis=-1), best_U
        return best_cost, self.vc[j] + dv_arr, np.linalg.norm(best_U, axis=-1)

    def check(self, j, U, R_target, mu=None):
        """Non-linear check: apply profile U (K,3) as impulses u_k*dtau at tau_k, propagate, return miss distance [km]."""
        K = self.K[j]; r, v = self.r.copy(), self.v.copy(); t = 0.0
        for k in range(K):
            r, v = propagate_twobody(r, v, np.array([self.taus[k] - t])); r, v = r[0], v[0]; t = self.taus[k]
            v = v + U[k] * self.dtau
        r, v = propagate_twobody(r, v, np.array([self.tofs[j] - t]))
        return float(np.linalg.norm(r[0] - R_target))
