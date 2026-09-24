#!/bin/zsh
PY=~/.venvs/astro313/bin/python
cd /Users/mickey/solarsystem/ctoc14/results/s13/planner_probe
$PY probe_joint.py j2_b300_p0809_canon --pool routes:08,09 --K 2 --beam 300 --nproc 1 --canon > j2_b300_p0809_canon.out 2>&1
$PY settle_joint_tours.py j2_b300_p0809_canon > settle_j2_b300_50.out 2>&1
