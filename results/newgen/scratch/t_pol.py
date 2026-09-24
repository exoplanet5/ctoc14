import os
os.environ['OMP_NUM_THREADS']='1'; os.environ['VECLIB_MAXIMUM_THREADS']='1'
import sys, time, numpy as np
sys.path.insert(0,'.')
from ctoc14.kepler import Ephemeris
from ctoc14.globalopt import problem_from_state, load_state, polish_ml, state_of, save_state
def cost(m0): x=(m0-600)/1400; return 1+x+x*x
eph = Ephemeris(); fleet, host, iters = sys.argv[1], sys.argv[2], int(sys.argv[3])
gp0 = problem_from_state(eph, load_state(f'{fleet}/craft_{host}.npz')); tic = time.time()
gt, d = polish_ml(gp0, iters, 1.5, verbose=True, log=lambda s: print(s, flush=True) if ' it ' in s and int(s.split()[1].rstrip(':')) % 10 == 0 else None)
print(f'{host}: m0 {gp0.m0:.1f} -> {gt.m0:.1f}, dJ {cost(gt.m0)-cost(gp0.m0):+.4f}, miss {d.max():.0f} ({time.time()-tic:.0f} s)')
save_state(state_of(gt), f'results/newgen/scratch/pol_{host}.npz')
