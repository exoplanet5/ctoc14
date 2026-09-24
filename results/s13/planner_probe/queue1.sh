#!/bin/zsh
PY=~/.venvs/astro313/bin/python
cd /Users/mickey/solarsystem/ctoc14/results/s13/planner_probe
$PY probe_beam.py b100_p0809_soft0 --pool routes:08,09 --mode soft --soft-prize 0.0 --beam 100 --nproc 1 --settle 2 > b100_p0809_soft0.out 2>&1
$PY probe_beam.py b100_p070809_strict --pool routes:07,08,09 --mode strict --beam 100 --nproc 1 --settle 2 > b100_p070809_strict.out 2>&1
$PY probe_beam.py b100_p070809_soft0 --pool routes:07,08,09 --mode soft --soft-prize 0.0 --beam 100 --nproc 1 --settle 2 > b100_p070809_soft0.out 2>&1
$PY probe_beam.py b100_p06_09_strict --pool routes:06,07,08,09 --mode strict --beam 100 --nproc 1 --settle 2 > b100_p06_09_strict.out 2>&1
$PY probe_beam.py b100_p04_09_strict --pool routes:04,05,06,07,08,09 --mode strict --beam 100 --nproc 1 --settle 2 > b100_p04_09_strict.out 2>&1
