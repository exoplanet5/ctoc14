#!/bin/zsh
PY=~/.venvs/astro313/bin/python
cd /Users/mickey/solarsystem/ctoc14/results/s13/planner_probe
while pgrep -f queue5.sh > /dev/null; do sleep 10; done
$PY probe_beam.py b100_p0809_dr25 --pool routes:08,09 --mode strict --drmax 0.25 --beam 100 --nproc 1 --settle 2 > b100_p0809_dr25.out 2>&1
$PY probe_beam.py b100_p0809_dv16 --pool routes:08,09 --mode strict --dvmax 1.6 --beam 100 --nproc 1 --settle 2 > b100_p0809_dv16.out 2>&1
$PY probe_beam.py b300_p0809 --pool routes:08,09 --mode strict --beam 300 --nproc 1 --settle 2 > b300_p0809.out 2>&1
