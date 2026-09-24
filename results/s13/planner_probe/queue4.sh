#!/bin/zsh
PY=~/.venvs/astro313/bin/python
cd /Users/mickey/solarsystem/ctoc14/results/s13/planner_probe
$PY probe_beam.py b100_pool50_minusA --pool pool50_minusA.json --mode strict --beam 100 --nproc 1 --settle 1 > b100_pool50_minusA.out 2>&1
$PY probe_beam.py b100_pool77_minusA --pool pool77_minusA.json --mode strict --beam 100 --nproc 1 --settle 1 > b100_pool77_minusA.out 2>&1
$PY probe_joint.py j2_b100_p070809_canon --pool routes:07,08,09 --K 2 --beam 100 --nproc 1 --canon > j2_b100_p070809_canon.out 2>&1
$PY settle_joint_tours.py j2_b100_p070809_canon > settle_j2_77.out 2>&1
$PY probe_joint.py j3_b100_p070809_canon --pool routes:07,08,09 --K 3 --beam 100 --nproc 1 --canon > j3_b100_p070809_canon.out 2>&1
$PY probe_joint.py j3_b100_p070809_orig --pool routes:07,08,09 --K 3 --beam 100 --nproc 1 > j3_b100_p070809_orig.out 2>&1
