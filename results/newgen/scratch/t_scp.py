import sys, time, numpy as np
sys.path.insert(0,'.')
from ctoc14.kepler import Ephemeris
from ctoc14.globalopt import load_craft, GlobalProblem, optimise
from ctoc14.constants import DAY
eph = Ephemeris()
sc = int(sys.argv[1]); iters = int(sys.argv[2])
c = load_craft('results/CTOC14_Result_TEAM.txt', sc)
gp = GlobalProblem(eph, c, hstep=0.25*DAY)
cur, hist = optimise(gp, iters=iters, rho0=float(sys.argv[3]) if len(sys.argv)>3 else 1.0)
np.savez(f'results/newgen/scratch/scp_sc{sc}.npz', Ts=gp.Ts, tf=gp.tf, vinf=gp.vinf, ts=gp.ts, asts=np.array(gp.asts), tL=gp.tL, m0=gp.m0)
