#!/bin/zsh
cd /Users/mickey/solarsystem/ctoc14
P=~/.venvs/astro313/bin/python
$P tools/run_ialns.py resettle results/newgen/ifleet6 results/newgen/ifleet7 --nproc 8
$P tools/run_ialns.py search results/newgen/ifleet7 results/newgen/ifleet8 --nproc 8 --rounds 10 --relocate --retime --dissolve-k 0
