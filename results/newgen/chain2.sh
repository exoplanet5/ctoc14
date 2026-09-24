#!/bin/zsh
cd /Users/mickey/solarsystem/ctoc14
P=~/.venvs/astro313/bin/python
$P tools/run_gfleet.py dedupe results/newgen/fleet_go1 results/newgen/fleet_dd2 --iters 60 --nproc 9 --keep-big
$P tools/run_gfleet.py dissolve results/newgen/fleet_dd2 results/newgen/fleet_ds2 --auto 4 --iters 60 --trials 6 --nproc 9
