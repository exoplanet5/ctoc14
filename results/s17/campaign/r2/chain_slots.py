"""R2 CAMPAIGN zero-compute join: for each r9 hard-core target h and each surviving host with an approach < 0.2 AU,
count the host's own flybys within +-W days of the approach whose forced re-home price (best other host, not r9, not
self; redteam dissolve_*.jsonl lin screen) is <= 30 kg.  Those are the only droppable own targets that let a 1-for-1
re-plan link close under the 46 kg/target budget (F5).  Output r2/chain_slots.json."""
import json, collections
R = '/Users/mickey/solarsystem/ctoc14/results/s17/'
rows = [json.loads(l) for f in ('redteam/dissolve_r9_r8.jsonl', 'redteam/dissolve_r1_r2_r3_r4_r5_r6_r7.jsonl') for l in open(R + f)]
best = {}
for r in rows:
    if r['host'] in ('r9', r['V']): continue
    if r['ast'] not in best or r['lin_kg'] < best[r['ast']][0]: best[r['ast']] = (r['lin_kg'], r['host'])
an = json.load(open(R + 'lns/anatomy.json'))
c9 = json.load(open(R + 'lns/census_r9.json'))['cands']
hard = [17, 36, 75, 108, 118, 180, 212, 216, 219, 220, 243, 247, 288]
out = {}
for W in (150, 250):
    tab = {}
    for h in hard:
        slots = []
        for x in c9:
            if x['ast'] != h or x['dist'] >= 0.2: continue
            host = x['host']; t = x['t'] / 86400
            near = [(a, round(tf)) for a, tf in zip(an[host]['asts'], an[host]['tf']) if abs(tf - t) <= W and a in best and best[a][0] <= 30]
            slots.append(dict(host=host, t=round(t), dist=round(x['dist'], 3), forced_kg=round(x['lin_kg']), cheap_own_near=near))
        tab[h] = dict(n_hosts=len(slots), n_hosts_with_slot=len({s['host'] for s in slots if s['cheap_own_near']}), slots=slots)
    out[f'W{W}'] = tab
    ok = [h for h in hard if tab[h]['n_hosts_with_slot'] > 0]
    print(f'W={W} d: hard targets with >=1 host approach <0.2 AU: {sum(tab[h]["n_hosts"]>0 for h in hard)}/13; '
          f'with a cheap (<=30 kg) droppable own flyby near the approach: {len(ok)}/13 {ok}')
    for h in hard:
        print('  ', h, [(s['host'], s['t'], s['dist'], s['forced_kg'], s['cheap_own_near']) for s in tab[h]['slots']])
json.dump(out, open(R + 'campaign/r2/chain_slots.json', 'w'), indent=1)
