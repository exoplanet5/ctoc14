import sys, time, numpy as np
sys.path.insert(0,'.')
from ctoc14.kepler import Ephemeris
from ctoc14.globalopt import load_craft, GlobalProblem, optimise, restore
from ctoc14.constants import DAY
eph = Ephemeris()
sc = int(sys.argv[1])
c = load_craft('results/CTOC14_Result_TEAM.txt', sc)
gp = GlobalProblem(eph, c, hstep=0.25*DAY)
z = np.load(f'results/newgen/scratch/scp_sc{sc}.npz')
gp.Ts, gp.tf, gp.vinf = z['Ts'], z['tf'], z['vinf']
d = restore(gp, iters=15)
np.savez(f'results/newgen/scratch/scp_sc{sc}_r.npz', Ts=gp.Ts, tf=gp.tf, vinf=gp.vinf, ts=gp.ts, asts=np.array(gp.asts), tL=gp.tL, m0=gp.m0)
