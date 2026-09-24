import os
os.environ['OMP_NUM_THREADS']='1'; os.environ['VECLIB_MAXIMUM_THREADS']='1'; os.environ['OPENBLAS_NUM_THREADS']='1'
import sys, time, numpy as np
sys.path.insert(0,'.')
from ctoc14.kepler import Ephemeris
from ctoc14.globalopt import problem_from_state, load_state, optimise, restore, add_flyby, extend_arc
from ctoc14.constants import DAY
eph = Ephemeris()
host, ast, t, iters = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
st = load_state(f'results/newgen/fleet_dd/craft_{host}.npz')
gp = problem_from_state(eph, st)
print('fuel before', gp.fuel(), 'm0', gp.m0, 'flybys', len(gp.asts))
extend_arc(gp, t*DAY + 5*DAY)
add_flyby(gp, ast, t*DAY)
d = restore(gp, tol_km=300, iters=40, rho=1e5); print("restored", gp.fuel(), d.max())
cur, hist = optimise(gp, iters=iters, rho0=float(sys.argv[5]) if len(sys.argv)>5 else 100.0)
d = restore(gp)
print('after insert: fuel', gp.fuel(), 'max miss', d.max())
