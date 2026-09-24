"""Robustness check: validator2 with the alternative (forward-biased, reading B) Lagrange window {j..j+3} clipped to [0,n-4]."""
import sys, json
sys.path.insert(0, 'tools')
import numpy as np
import validator2 as v2

def window_B(self, t):
    n = self.n
    if n < 4:
        return list(range(n))
    j = int(np.searchsorted(self.t, t, side="right")) - 1
    j = min(max(j, 0), n - 2)
    l = min(max(j, 0), n - 4)
    return [l, l + 1, l + 2, l + 3]
v2.Arc.window = window_B
rep = v2.validate(sys.argv[1], 'MEA.txt', crosscheck=5, quiet=True)
print('reading B:', sys.argv[1])
print('verdict', rep['verdict'], 'N_cov', rep['N_covered'], 'J %.6f' % rep['J'])
print('max pos %.4e km at %s | vel %.4e m/s at %s | mass %.4e kg at %s | max interp |T| %.9f N'
      % (rep['max_pos_err_km'], rep['max_pos_err_at'], rep['max_vel_err_kms'] * 1e3, rep['max_vel_err_at'],
         rep['max_mass_err_kg'], rep['max_mass_err_at'], rep['max_T_interp']))
print('errors', len(rep['errors']), rep['errors'][:5], 'warnings', rep['warnings'])
