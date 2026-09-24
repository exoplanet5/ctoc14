#!/bin/zsh
P=~/.venvs/astro313/bin/python
cd /Users/mickey/solarsystem/ctoc14
$P tools/run_fleet_search.py results/plan_E1 --beam 300 --wfuel 0.7 --mmargin 100 --dvmax 2.5 --improve-rounds 2 --eliminate --nproc 6 > results/plan_E1.log 2>&1
$P tools/run_fleet_search.py results/plan_E2 --beam 300 --wfuel 0.7 --mmargin 100 --dvmax 2.5 --wrare 0.5 --improve-rounds 2 --eliminate --nproc 6 > results/plan_E2.log 2>&1
$P tools/run_fleet_search.py results/plan_E3 --beam 800 --wfuel 0.7 --mmargin 100 --dvmax 2.5 --wrare 0.5 --improve-rounds 2 --eliminate --nproc 6 > results/plan_E3.log 2>&1
