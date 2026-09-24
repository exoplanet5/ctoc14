# Independent score / format check of a CTOC14 submission file (read-only on the input).
import sys, math, re
import numpy as np
MU=1.32712440018e11; AU=149597870.7
T0_MJD=62502.0; TEND=473364000.0
fn=sys.argv[1]
# asteroid elements
ast={}
for L in open('/Users/mickey/solarsystem/ctoc14/MEA.txt'):
    if L.strip()=='' or L.lstrip().startswith('#'): continue
    p=L.split(); ast[int(p[0])]=[float(v) for v in p[1:7]]
def kep(el, t_s, eph_mjd):
    a,e,i,O,w,M0=el; a*=AU; i,O,w,M0=[math.radians(v) for v in (i,O,w,M0)]
    n=math.sqrt(MU/a**3); dt=t_s+(T0_MJD-eph_mjd)*86400.0
    M=math.fmod(M0+n*dt,2*math.pi)
    E=M if e<0.8 else math.pi
    for _ in range(100):
        f=E-e*math.sin(E)-M; d=f/(1-e*math.cos(E)); E-=d
        if abs(d)<1e-15: break
    cE,sE=math.cos(E),math.sin(E)
    x=a*(cE-e); y=a*math.sqrt(1-e*e)*sE
    r=a*(1-e*cE); vx=-math.sqrt(MU*a)/r*sE; vy=math.sqrt(MU*a)/r*math.sqrt(1-e*e)*cE
    cO,sO,ci,si,cw,sw=math.cos(O),math.sin(O),math.cos(i),math.sin(i),math.cos(w),math.sin(w)
    R=np.array([[cO*cw-sO*sw*ci,-cO*sw-sO*cw*ci],[sO*cw+cO*sw*ci,-sO*sw+cO*cw*ci],[sw*si,cw*si]])
    return R@np.array([x,y]), R@np.array([vx,vy])
EARTH=[1.0009175020,0.017566762041,0.002976847126,189.953211282428,273.196254000254,357.4135031077]
rows=[]; bad_tok=0; ncomment=0; ncols_bad=0; minsig=99; minsig_ex=None; sigcount={}
numre=re.compile(r'^[+-]?(\d*)\.?(\d*)(?:[eE][+-]?\d+)?$')
def sig(tok):
    m=re.match(r'^[+-]?([0-9.]+)',tok); mant=m.group(1).replace('.','').lstrip('0')
    return len(mant)
for ln,L in enumerate(open(fn),1):
    s=L.strip()
    if s=='' or s.startswith('#'): ncomment+=1; continue
    p=s.split()
    if len(p)!=15: ncols_bad+=1; continue
    for k,tok in enumerate(p):
        tl=tok.lower()
        if 'nan' in tl or 'inf' in tl or not numre.match(tok): bad_tok+=1
    v=[float(x) for x in p]
    if not all(math.isfinite(x) for x in v): bad_tok+=1
    # significant digits on state/mass fields that are not exactly-representable short values
    for k in range(4,11):
        tok=p[k]; sg=sig(tok)
        if sg<10:
            sigcount[k]=sigcount.get(k,0)+1
            if sg<minsig: minsig=sg; minsig_ex=(ln,k,tok)
    rows.append((ln,int(p[0]),int(p[1]),int(p[2]),v[3],np.array(v[4:7]),np.array(v[7:10]),v[10],np.array(v[11:14]),int(p[14])))
print('file',fn); print('data rows',len(rows),'comment/blank',ncomment,'bad-col rows',ncols_bad,'NaN/inf/garbage tokens',bad_tok)
print('Line column consecutive 1..N:', all(r[1]==i+1 for i,r in enumerate(rows)))
print('fields (cols 5-11) with <10 sig digits:',sigcount,'min',minsig,minsig_ex)
scs=sorted(set(r[2] for r in rows)); print('SC ids',scs, 'consecutive from 1:', scs==list(range(1,len(scs)+1)))
covered={}; Jsum=0; out=[]
for sc in scs:
    R=[r for r in rows if r[2]==sc]
    ev=[r[3] for r in R]; t=[r[4] for r in R]
    ok_struct = ev[0]==0 and ev[-1]==4 and all(b>a for a,b in zip(t,t[1:])) and ev.count(0)==1 and ev.count(4)==1
    tin = all(0<=x<=TEND for x in t)
    m0=R[0][7]; mmin=min(r[7] for r in R); x=(m0-600)/1400; Ji=1+x+x*x; Jsum+=Ji
    rE,vE=kep(EARTH,R[0][4],60676.0)
    dr=np.linalg.norm(R[0][5]-rE); vinf=np.linalg.norm(R[0][6]-vE)
    Tmaxs=max(np.linalg.norm(r[8]) for r in R)
    nf=0; nf_ok=0; dmax=0
    for r in R:
        if r[3]==3:
            nf+=1
            ra,_=kep(ast[r[9]],r[4],61200.0); d=np.linalg.norm(r[5]-ra); dmax=max(dmax,d)
            if d<=1000 and 0<=r[4]<=TEND:
                nf_ok+=1
                if r[9] not in covered or covered[r[9]][0]>r[4]: covered[r[9]]=(r[4],sc)
    out.append((sc,len(R),m0,Ji,mmin,R[-1][7],dr,vinf,Tmaxs,nf,nf_ok,dmax,ok_struct,tin,t[0]/86400,t[-1]/86400))
print('sc rows      m0         J_i      m_min     m_end   dr_launch[km] vinf[km/s] max|T_sample| nfb nfb<=1000 max_d[km] struct tin tL[d] tEnd[d]')
for o in out: print('%2d %6d %11.6f %9.6f %9.4f %9.4f %9.3e %9.6f %8.5f %4d %4d %8.2f %s %s %8.2f %8.2f'%o)
Nc=len(covered); Nmiss=300-Nc
print('distinct covered',Nc,'missing',sorted(set(range(1,301))-set(covered)))
print('sum J_i = %.6f  N_miss=%d  J = %.6f'%(Jsum,Nmiss,Jsum+Nmiss))
