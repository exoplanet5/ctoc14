#!/bin/zsh
PY=~/.venvs/astro313/bin/python
cd /Users/mickey/solarsystem/ctoc14/results/s13/planner_probe
while pgrep -f queue6.sh > /dev/null; do sleep 10; done
$PY probe_beam.py b100_p00_strict --pool routes:00 --mode strict --beam 100 --nproc 1 --settle 1 > b100_p00_strict.out 2>&1
$PY probe_beam.py b100_rand200_strict --pool pool_rand200.json --mode strict --beam 100 --nproc 1 --settle 1 > b100_rand200_strict.out 2>&1
$PY probe_beam.py b100_rand250_strict --pool pool_rand250.json --mode strict --beam 100 --nproc 1 --settle 1 > b100_rand250_strict.out 2>&1
