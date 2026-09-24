import os
os.environ['OMP_NUM_THREADS']='1'; os.environ['VECLIB_MAXIMUM_THREADS']='1'
import sys, time, numpy as np
sys.path.insert(0,'.')
from ctoc14.kepler import Ephemeris
from ctoc14.globalopt import problem_from_state, load_state, restore, tighten, write_craft, state_of, save_state
from ctoc14.constants import DAY
eph = Ephemeris(); path, out = sys.argv[1], sys.argv[2]; hs = float(sys.argv[3])
gp = problem_from_state(eph, load_state(path), hstep=hs*DAY)
Yf, _ = gp.integrate(dense_samples=True); d0 = np.linalg.norm(gp.misses(Yf)[0], axis=1)
print(f'hstep {hs} d: misses max {d0.max():.0f} km (m0 {gp.m0:.1f})', flush=True); tic = time.time()
d = restore(gp, tol_km=100, iters=30, log=lambda s: None); d = tighten(gp, 1.5, log=lambda s: None)
print(f'restored: max {d.max():.0f} km, m0 {gp.m0:.1f} ({time.time()-tic:.0f} s)', flush=True)
save_state(state_of(gp), path)
rep, fuel, pk = write_craft(gp, out, eph=eph)
print('valid', rep.ok, rep.errors[:3], 'flybys', len(rep.flybys), f'({time.time()-tic:.0f} s)')
