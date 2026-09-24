#!/bin/zsh
# waits for go1, then: combined submission (go1), fleet init -> dedupe -> dissolve
cd /Users/mickey/solarsystem/ctoc14
P=~/.venvs/astro313/bin/python
while ps -eo args | grep "run_globalopt.py results/CTOC14_Result_TEAM.txt results/newgen/go1" | grep -v grep > /dev/null; do sleep 20; done
$P tools/run_gfleet.py init results/newgen/fleet_go1 --from-sub results/CTOC14_Result_TEAM.txt --from-ckpt results/newgen/go1
$P tools/run_gfleet.py dedupe results/newgen/fleet_go1 results/newgen/fleet_dd --iters 60 --nproc 9
$P tools/run_gfleet.py dissolve results/newgen/fleet_dd results/newgen/fleet_ds --auto 4 --iters 60 --trials 6 --nproc 9
