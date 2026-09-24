#!/bin/zsh
PY=~/.venvs/astro313/bin/python
cd /Users/mickey/solarsystem/ctoc14/results/s13/planner_probe
while pgrep -f queue7.sh > /dev/null; do sleep 10; done
$PY probe_beam.py b100_all_minusA46 --pool pool_all_minusA46.json --mode strict --beam 100 --nproc 1 --settle 2 > b100_all_minusA46.out 2>&1
