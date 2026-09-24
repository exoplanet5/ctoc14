#!/bin/zsh
cd /Users/mickey/solarsystem/ctoc14
P=~/.venvs/astro313/bin/python
$P tools/run_ialns.py export results/newgen/ifleet_d11 results/newgen/efleet_d11 --nproc 2
$P tools/run_gfleet.py write results/newgen/efleet_d11 results/CTOC14_Result_d11.txt --nproc 2
