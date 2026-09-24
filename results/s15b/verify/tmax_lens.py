import sys, numpy as np
fn=sys.argv[1]
d=np.loadtxt(fn,comments='#')
best=(0,None); narc=0; minsp=1e99
for sc in np.unique(d[:,1]).astype(int):
    R=d[d[:,1]==sc]; ev=R[:,2].astype(int)
    arcs=[]; cur=[]
    for r,e in zip(R,ev):
        if e==1: cur.append(r)
        elif e==3: continue
        else:
            if cur: arcs.append(np.array(cur)); cur=[]
    if cur: arcs.append(np.array(cur))
    for A in arcs:
        narc+=1; t=A[:,3]; T=A[:,11:14]; n=len(t)
        if n>1: minsp=min(minsp,np.min(np.diff(t)))
        for j in range(n-1):
            if n>=4: s=min(max(j-1,0),n-4); idx=range(s,s+4)
            else: idx=range(n)
            tt=np.linspace(t[j],t[j+1],41)
            val=np.zeros((41,3))
            for a in idx:
                L=np.ones(41)
                for b in idx:
                    if b!=a: L*= (tt-t[b])/(t[a]-t[b])
                val+=L[:,None]*T[a]
            m=np.linalg.norm(val,axis=1).max()
            if m>best[0]: best=(m,(sc,int(A[j,0]),t[j]))
print('arcs',narc,'min Event=1 spacing [s] %.1f'%minsp,'max interpolated |T| %.9f N at'%best[0],best[1])
