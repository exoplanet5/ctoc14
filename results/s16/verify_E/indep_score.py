# Independent J recomputation for a CTOC14 file (read-only). Different Kepler code from score_lens.py.
import sys, numpy as np
MU=1.32712440018e11; AU=149597870.7; T0=62502.0; TMAX=4.73364e8
fn=sys.argv[1]
D=np.loadtxt(fn, comments='#')
el=np.loadtxt('/Users/mickey/solarsystem/ctoc14/MEA.txt', comments='#')
elmap={int(r[0]):r[1:] for r in el}
def pos(elem, t, epoch):
    a,e,i,O,w,M0=elem; a*=AU; i,O,w,M0=np.radians([i,O,w,M0])
    n=np.sqrt(MU/a**3); M=M0+n*(t+(T0-epoch)*86400.0)
    M=np.mod(M,2*np.pi); E=np.where(e>0.8,np.pi*np.ones_like(M),M)
    for _ in range(60): E=E-(E-e*np.sin(E)-M)/(1-e*np.cos(E))
    nu=2*np.arctan2(np.sqrt(1+e)*np.sin(E/2),np.sqrt(1-e)*np.cos(E/2)); r=a*(1-e*np.cos(E))
    u=w+nu
    return r*np.array([np.cos(O)*np.cos(u)-np.sin(O)*np.sin(u)*np.cos(i),
                       np.sin(O)*np.cos(u)+np.cos(O)*np.sin(u)*np.cos(i),
                       np.sin(u)*np.sin(i)])
sc=D[:,1].astype(int); ev=D[:,2].astype(int); t=D[:,3]
Ji=[]; first={}
for s in np.unique(sc):
    R=D[sc==s]
    assert R[0,2]==0 and R[-1,2]==4
    m0=R[0,10]; x=(m0-600)/1400; Ji.append(1+x+x*x)
fb=D[ev==3]; dmax=0; nbad=0
for r in fb:
    aid=int(r[14]); d=np.linalg.norm(r[4:7]-pos(elmap[aid], r[3], 61200.0)); dmax=max(dmax,d)
    if d<=1000 and 0<=r[3]<=TMAX:
        if aid not in first or r[3]<first[aid]: first[aid]=r[3]
    else: nbad+=1
N=len(first); S=sum(Ji)
print('craft',len(Ji),'flyby rows',len(fb),'bad flybys',nbad,'max d %.3f km'%dmax,'covered',N,'missing',sorted(set(range(1,301))-set(first)))
print('J_i',' '.join('%.6f'%j for j in Ji))
print('sumJi %.9f  raw J %.9f'%(S,S+300-N))
