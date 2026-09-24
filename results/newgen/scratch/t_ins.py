import os
os.environ['OMP_NUM_THREADS']='1'; os.environ['VECLIB_MAXIMUM_THREADS']='1'; os.environ['OPENBLAS_NUM_THREADS']='1'
import sys, time, numpy as np
sys.path.insert(0,'.')
from ctoc14.kepler import Ephemeris
from ctoc14.globalopt import load_craft, GlobalProblem, optimise, restore, add_flyby, write_craft
from ctoc14.constants import DAY
eph = Ephemeris()
sc, ast, t, iters = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
c = load_craft('results/CTOC14_Result_TEAM.txt', sc)
gp = GlobalProblem(eph, c, hstep=0.25*DAY)
print('fuel before', gp.fuel(), 'flybys', len(gp.asts))
add_flyby(gp, ast, t*DAY)
cur, hist = optimise(gp, iters=iters, rho0=100.0, mu_km=2e4)
d = restore(gp, tol_km=50, iters=30)
print('after insert: fuel', gp.fuel(), 'max miss', d.max())
np.savez(f'results/newgen/scratch/ins_sc{sc}_{ast}.npz', Ts=gp.Ts, tf=gp.tf, vinf=gp.vinf, asts=np.array(gp.asts), m0=gp.m0)
