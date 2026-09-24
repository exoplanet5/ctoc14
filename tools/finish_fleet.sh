#!/bin/sh
# Finish a planned fleet: convert its new routes, tighten every route's best fragment twice, then a selection-only
# fragment MILP over ALL valid fragments.   Usage: tools/finish_fleet.sh fleetdir tag [nproc]
P=$HOME/.venvs/astro313/bin/python
cd /Users/mickey/solarsystem/ctoc14 || exit 1
N=${3:-3}
$P tools/convert_pool.py results/colgen/frags "$1" --nproc $N
$P tools/fleet_frags.py results/colgen/frags "$1" > results/colgen/tight_$2.txt
$P tools/tighten_frags.py results/colgen/tight_$2.txt --margin 6 --nproc $N --rounds 2 --min-left 8
ls results/*/frag_*_sc*.txt results/*/*/frag_*_sc*.txt results/*/*/*/frag_*_sc*.txt results/*/*_frag.txt results/colgen/frags/frag_*_p[12]*.txt 2>/dev/null \
  | sort -u | $P -c "
import sys, json, pathlib
for f in sys.stdin:
    f = f.strip(); i = pathlib.Path(f[:-4] + '_info.json')
    if i.exists():
        try:
            if not json.load(open(i)).get('ok', True): continue
        except Exception: pass
    print(f)" > results/colgen/frag_list_$2.txt
cat results/colgen/frag_list_$2.txt | xargs $P tools/select_frags_milp.py --dry results/CTOC14_Result_$2.txt
