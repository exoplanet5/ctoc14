"""
Independent from-scratch cross-check of population_stats.csv (second implementation, no imports from
kepler.py / population_analysis.py / ctoc14). Own Newton Kepler solver, own perifocal->HEI rotation via the
argument of latitude, velocities from the radial/transverse components, 0.1-d brute-force grid over the window.
Checks: PDF sec. 6.2 sample lines; totals of perihelia / ecliptic crossings / in-band (0.8-1.6 AU) crossings and their
per-year histogram; inclination bins; r_min > 2.5 AU objects; closest Earth approaches; per-asteroid agreement with the CSV.
Run:  ~/.venvs/astro313/bin/python analysis/population/independent_check.py   (~4 s; output in independent_check.log)
"""
import numpy as np, pandas as pd
MU=1.32712440018e11; AU=149597870.7; DAY=86400.0
T0=62502.0; TEND=T0+4.73364e8/DAY
rows=[l.split() for l in open('/Users/mickey/solarsystem/ctoc14/MEA.txt') if l.strip() and not l.startswith('#')]
A=np.array([[float(x) for x in r] for r in rows]); ids=A[:,0].astype(int); el=A[:,1:7]
EARTH=np.array([1.0009175020,0.017566762041,0.002976847126,189.953211282428,273.196254000254,357.4135031077])
def kep(M,e):
    M=np.mod(M,2*np.pi); E=np.where(e<0.8,M,np.pi)*np.ones_like(M)
    for _ in range(100):
        d=(E-e*np.sin(E)-M)/(1-e*np.cos(E)); E=E-d
        if np.all(np.abs(d)<1e-14): break
    return E
def state(elm,teph,t):
    a,e,i,Om,w,M0=[np.asarray(x,float) for x in np.moveaxis(np.atleast_2d(elm),1,0)] if np.ndim(elm)==2 else [float(x) for x in elm]
    a=np.asarray(a)[...,None] if np.ndim(elm)==2 else a
    e=np.asarray(e)[...,None] if np.ndim(elm)==2 else e
    n=np.sqrt(MU/(a*AU)**3); M=np.radians(M0 if np.ndim(elm)<2 else np.asarray(M0)[...,None])+n*(np.asarray(t)-teph)*DAY
    E=kep(M,e); nu=2*np.arctan2(np.sqrt(1+e)*np.sin(E/2),np.sqrt(1-e)*np.cos(E/2))
    r=a*AU*(1-e*np.cos(E)); h=np.sqrt(MU*a*AU*(1-e*e))
    i,Om,w=[np.radians(np.asarray(x))[...,None] if np.ndim(elm)==2 else np.radians(x) for x in (i,Om,w)]
    u=w+nu
    x=r*(np.cos(Om)*np.cos(u)-np.sin(Om)*np.sin(u)*np.cos(i)); y=r*(np.sin(Om)*np.cos(u)+np.cos(Om)*np.sin(u)*np.cos(i)); z=r*np.sin(u)*np.sin(i)
    # velocity via vis-viva components
    vr=MU/h*e*np.sin(nu); vt=h/r
    vx=vr*x/r - vt*(np.cos(Om)*np.sin(u)+np.sin(Om)*np.cos(u)*np.cos(i))
    vy=vr*y/r - vt*(np.sin(Om)*np.sin(u)-np.cos(Om)*np.cos(u)*np.cos(i))
    vz=vr*z/r + vt*np.cos(u)*np.sin(i)
    return np.stack([x,y,z],-1), np.stack([vx,vy,vz],-1)
# 1. PDF sample lines
rE,vE=state(EARTH,60676.0,T0)
print('Earth t0 pos err km', np.linalg.norm(rE-np.array([-1.9500328647e+07,1.4581848347e+08,-7.6372048832e+03])))
v1=np.array([-2.9389477847e+01,-4.3166912970e+00,-3.9417153392e+00]); print('vinf km/s',np.linalg.norm(v1-vE))
ra,va=state(el[173],61200.0,T0+1.0008479152e+07/DAY)
print('ast174 err km',np.linalg.norm(ra-np.array([-1.1398611209e+08,-8.7213961221e+07,-1.6523817305e+07])))
# 2. fine-grid scan 0.1 d
t=np.arange(T0,TEND+1e-9,0.1); 
if t[-1]<TEND: t=np.append(t,TEND)
pe,_=state(EARTH,60676.0,t)
nper=np.zeros(300,int); ncross=np.zeros(300,int); inband_cross=[]; rmin=np.zeros(300); dmin=np.zeros(300); tdmin=np.zeros(300)
for c in range(0,300,50):
    sl=slice(c,c+50); p,_=state(el[sl],61200.0,t); r=np.linalg.norm(p,axis=-1)/AU; z=p[...,2]/AU
    nper[sl]=((r[:,1:-1]<r[:,:-2])&(r[:,1:-1]<=r[:,2:])).sum(1)
    sc=np.sign(z[:,:-1])*np.sign(z[:,1:])<0; ncross[sl]=sc.sum(1)
    rc=0.5*(r[:,:-1]+r[:,1:]); tc=0.5*(t[:-1]+t[1:])
    for k in range(sc.shape[0]):
        m=sc[k]&(rc[k]>=0.8)&(rc[k]<=1.6); inband_cross+= [(ids[c+k],tt,rr) for tt,rr in zip(tc[m],rc[k][m])]
    rmin[sl]=r.min(1); d=np.linalg.norm(p-pe[None],axis=-1)/AU; kk=d.argmin(1); dmin[sl]=d[np.arange(d.shape[0]),kk]; tdmin[sl]=t[kk]
ib=pd.DataFrame(inband_cross,columns=['id','t','r'])
edges=np.append(T0+np.arange(15)*365.25,TEND+1e-6)
print('perihelia total',nper.sum(),'crossings total',ncross.sum(),'in-band crossings',len(ib),'from',ib.id.nunique(),'asteroids')
print('in-band per year',np.histogram(ib.t,bins=edges)[0].tolist())
q=el[:,0]*(1-el[:,1]); Q=el[:,0]*(1+el[:,1]); inc=el[:,2]
print('q max',q.max(),'q<1',(q<1).sum(),'all cross 0.7-1.5',((q<=1.5)&(Q>=0.7)).all())
print('i<=10',(inc<=10).sum(),'10<i<=20',((inc>10)&(inc<=20)).sum(),'i>20',(inc>20).sum(),'retro',(inc>90).sum())
print('rmin>2.5 ids',ids[rmin>2.5].tolist(),'rmin values',rmin[rmin>2.5])
print('n_peri==0 ids',ids[nper==0].tolist())
o=np.argsort(dmin)[:6]; print('closest Earth approaches:',[(int(ids[j]),round(dmin[j],5),round(tdmin[j],2)) for j in o])
print('d_earth_min<0.1 asteroids',(dmin<0.1).sum(),'<0.05',(dmin<0.05).sum(),'<0.02',(dmin<0.02).sum())
# 3. 131/144 r at dates
for k in (131,144):
    tt=np.array([T0,T0+5*365.25,T0+10*365.25,TEND]); p,_=state(el[k-1],61200.0,tt); print(k,'r AU',np.round(np.linalg.norm(p,axis=1)/AU,3))
# 4. compare with stored CSV
S=pd.read_csv('/Users/mickey/solarsystem/ctoc14/analysis/population/population_stats.csv')
print('csv rows',len(S),'class counts',S.difficulty.value_counts().to_dict())
print('peri count diff ids',ids[S.n_perihelion_passages.values!=nper].tolist(),'cross count diff ids',ids[S.n_ecliptic_crossings.values!=ncross].tolist())
print('max |dmin - csv|',np.abs(S.d_earth_min_window_au.values-dmin).max(),'max |rmin-csv|',np.abs(S.r_min_window_au.values-rmin).max())
print('in-band crossing count per asteroid diff ids',ids[S.n_node_crossings_in_band.values!=ib.groupby('id').size().reindex(ids,fill_value=0).values].tolist())
