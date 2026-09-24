import sys, time, numpy as np
sys.path.insert(0,'.')
from ctoc14.kepler import Ephemeris
from ctoc14.globalopt import load_craft, GlobalProblem, write_craft
from ctoc14.constants import DAY
eph = Ephemeris()
sc = int(sys.argv[1]); src = sys.argv[2]
c = load_craft('results/CTOC14_Result_TEAM.txt', sc)
gp = GlobalProblem(eph, c, hstep=0.25*DAY)
z = np.load(src)
gp.Ts, gp.tf, gp.vinf = z['Ts'], z['tf'], z['vinf']
rep, fuel, pk = write_craft(gp, f'results/newgen/scratch/go_sc{sc}.txt')
print('ok', rep.ok, 'errors', rep.errors[:5], 'fuel', fuel, 'peak thrust', pk, 'covered', rep.n_covered, 'flybys', len(rep.flybys))
print('dists', sorted(round(v[2]) for v in rep.flybys.values())[-5:])
