#!/bin/zsh
cd /Users/mickey/solarsystem/ctoc14
P=~/.venvs/astro313/bin/python
while ps -eo args | grep "export_d11.sh" | grep -v grep > /dev/null; do sleep 30; done
$P tools/run_ialns.py export results/newgen/ifleet_d08 results/newgen/efleet_d08 --nproc 2
$P tools/run_gfleet.py write results/newgen/efleet_d08 results/CTOC14_Result_d08.txt --nproc 2
