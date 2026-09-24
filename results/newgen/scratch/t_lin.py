import sys, time, numpy as np
sys.path.insert(0,'.')
from ctoc14.kepler import Ephemeris
from ctoc14.globalopt import load_craft, GlobalProblem, irls_min_fuel
from ctoc14.constants import DAY, VE
eph = Ephemeris()
sc = int(sys.argv[1]) if len(sys.argv)>1 else 12
c = load_craft('results/CTOC14_Result_TEAM.txt', sc)
gp = GlobalProblem(eph, c, hstep=0.25*DAY)
tic=time.time(); Yf, Ys = gp.integrate(dense_samples=True); print('integrate', time.time()-tic)
dm, vr = gp.misses(Yf); print('ref misses km max', np.abs(np.linalg.norm(dm,axis=1)).max(), 'fuel model', gp.fuel(), 'file', c['m0']-c['m_end'])
tic=time.time(); B, Lv = gp.linearise(Yf, Ys); print('linearise', time.time()-tic, B.shape)
# FD check: perturb sample s by 0.02 N in x
for s in [100, 1500, 3000]:
    for ax in [0,1]:
        T2 = gp.Ts.copy(); T2[s, ax] += 0.02
        Yf2,_ = gp.integrate(Ts=T2)
        d_true = Yf2[:,0:3]-Yf[:,0:3]
        d_lin = B[:, s, :, ax]*0.02
        j = np.argmax(np.linalg.norm(d_true,axis=1))
        print(f's={s} ax={ax}: max true {np.linalg.norm(d_true,axis=1).max():.1f} km, lin err {np.linalg.norm(d_true-d_lin,axis=1).max():.2f} km')
dvv = np.array([0, 1e-3, 0]); Yf3,_ = gp.integrate(vinf=gp.vinf+dvv)
print('vinf check', np.linalg.norm(Yf3[:,0:3]-Yf[:,0:3],axis=1).max(), np.linalg.norm((Yf3[:,0:3]-Yf[:,0:3]) - Lv@dvv,axis=1).max())
np.savez('results/newgen/scratch/lin_sc%d.npz'%sc, B=B, Lv=Lv, Yf=Yf, Ys=Ys, dm=dm, vr=vr)
