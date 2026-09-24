"""Isolated test of closest_approach on a dumped window state."""
import sys, numpy as np; sys.path.insert(0, '/Users/mickey/solarsystem/ctoc14')
from ctoc14.kepler import Ephemeris
from ctoc14.lowthrust import integrate_rk4, closest_approach
from ctoc14.constants import DAY, T_MISSION
d = np.load(sys.argv[1]); eph = Ephemeris()
Y0, t0, ts, Ts, asts, t_noms, t_end = d['Y_start'], float(d['t_start']), d['ts'], d['Ts'], d['asts'], d['t_noms'], float(d['t_end'])
print('window start t=%.2f d, targets %s, t_noms %s d, t_end %.2f d' % (t0/DAY, asts, t_noms/DAY, t_end/DAY))
(Yend,), (ht, hY) = integrate_rk4(Y0[None], t0, [min(t_end, T_MISSION)], ts, Ts[None], 0.1*DAY, dense=True)
for a, tn in zip(asts, t_noms):
    t_star, Y, dvec, vrel, dist = closest_approach(eph, a-1, ht, hY, tn, ts=ts, Ts=Ts[None])
    print(f'ast {a}: closest_approach -> t*={t_star/DAY:.3f} d dist={dist:.1f} km |vrel|={np.linalg.norm(vrel):.2f} along-track={np.dot(dvec,vrel)/np.linalg.norm(vrel):.1f} km')
    # wide scan with a longer reference
    (Y2,), (ht2, hY2) = integrate_rk4(Y0[None], t0, [t_end + 150*DAY], ts, Ts[None], 0.1*DAY, dense=True)
    ra, va = eph.ast[a-1].state(ht2); dd = np.linalg.norm(hY2[:,0,:3]-ra, axis=1)
    k = int(np.argmin(dd)); print(f'   wide scan [{ht2[0]/DAY:.1f},{ht2[-1]/DAY:.1f}] d: min dist {dd[k]:.1f} km at t={ht2[k]/DAY:.2f} d; dist at t_end {dd[np.searchsorted(ht2,t_end)-1]:.1f} km; dist at t_nom {dd[np.searchsorted(ht2,tn)-1]:.1f}')
    # local minima list
    loc = [j for j in range(1,len(dd)-1) if dd[j]<dd[j-1] and dd[j]<dd[j+1]]
    print('   local minima (t d, dist km):', [(round(ht2[j]/DAY,1), round(dd[j])) for j in loc][:12])
