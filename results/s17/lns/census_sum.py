import json
an=json.load(open('results/s17/lns/anatomy.json'))
def cost(m): x=(m-600)/1400; return 1+x+x*x
out={}
for v in ('r9','r8','r7','r6','r5'):
    d=json.load(open(f'results/s17/lns/census_{v}.json'))
    by={}
    for c in d['cands']: by.setdefault(c['ast'],[]).append(c)
    best={X:min([c['lin_kg'] for c in by.get(X,[])]+[9999]) for X in d['targets']}
    band=lambda lo,hi: sum(1 for b in best.values() if lo<=b<hi)
    tot=sum(min(b,400) for b in best.values())
    cheap=sorted(X for X,b in best.items() if b<=40)
    hard=sorted(X for X,b in best.items() if b>150)
    tank=an[v]['tank']
    out[v]=dict(n=len(best),tank=tank,cJ=cost(tank),le40=band(0,40.0001),b40_100=band(40.0001,100),b100_150=band(100,150),gt150=band(150,1e9),sum_cap400=tot,cheap=cheap,hard=hard)
    print(f"{v}: n={len(best)} J_i={cost(tank):.3f} | best lin_kg <=40: {out[v]['le40']}, 40-100: {out[v]['b40_100']}, 100-150: {out[v]['b100_150']}, >150/none: {out[v]['gt150']} | sum(cap 400) {tot:.0f} kg | hard {hard}")
json.dump(out,open('results/s17/lns/census_summary.json','w'),indent=1)
