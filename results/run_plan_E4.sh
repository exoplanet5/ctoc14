#!/bin/zsh
P=~/.venvs/astro313/bin/python
cd /Users/mickey/solarsystem/ctoc14
until grep -qE "FLEET|Traceback" results/plan_E3.log 2>/dev/null; do sleep 30; done
$P tools/run_fleet_search.py results/plan_E4 --beam 300 --wfuel 0.7 --mmargin 100 --dvmax 2.5 --improve-rounds 1 --eliminate --tail-from 6 --nproc 6 > results/plan_E4.log 2>&1
