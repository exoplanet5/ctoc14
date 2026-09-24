#!/bin/zsh
PY=~/.venvs/astro313/bin/python
cd /Users/mickey/solarsystem/ctoc14/results/s13/planner_probe
$PY probe_joint.py j2_b100_p0809_canon --pool routes:08,09 --K 2 --beam 100 --nproc 1 --canon --settle > j2_b100_p0809_canon.out 2>&1
$PY probe_joint.py j2_b100_p0809_canon_nowait --pool routes:08,09 --K 2 --beam 100 --nproc 1 --canon --wait 100000 --settle > j2_b100_p0809_canon_nowait.out 2>&1
