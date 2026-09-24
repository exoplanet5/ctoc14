#!/bin/zsh
PY=~/.venvs/astro313/bin/python
cd /Users/mickey/solarsystem/ctoc14/results/s13/planner_probe
while pgrep -f queue1.sh > /dev/null; do sleep 10; done
$PY probe_beam.py b100_rand150_strict --pool pool_rand150.json --mode strict --beam 100 --nproc 1 --settle 2 > b100_rand150_strict.out 2>&1
$PY probe_beam.py b100_rand100_strict --pool pool_rand100.json --mode strict --beam 100 --nproc 1 --settle 2 > b100_rand100_strict.out 2>&1
cd /Users/mickey/solarsystem/ctoc14
/usr/bin/time -p $PY tools/skel_fill.py results/n8/skel8b results/s13/planner_probe/sf_free09 --only 09 --drop-wp all --mode windowed \
   --portfolio 100:0.15:45 --retry '' --pools results/s13/planner_probe/pools_0809_as09.json --pad-prize 0.0 \
   --m0 1600 --vinf 4.0 --dvmax 1.2 --try 2 --lam 0.4 --max-tank 1150 --nproc 1 > results/s13/planner_probe/sf_free09.out 2>&1
