#!/bin/zsh
setopt nonomatch
cd /Users/mickey/solarsystem/ctoc14
P=~/.venvs/astro313/bin/python
$P tools/run_ialns.py export results/newgen/ifleet_t10d results/newgen/efleet_t10d --nproc 6
if [ -n "$(ls results/newgen/efleet_t10d/frags/ | grep INVALID)" ]; then echo "invalid craft present"; exit 1; fi
$P tools/combine_submission.py results/CTOC14_Result_t10d.txt $(ls results/newgen/efleet_t10d/frags/frag_*.txt | sort)
$P tools/validator2.py results/CTOC14_Result_t10d.txt > results/validator2_t10d.log 2>&1
tail -3 results/validator2_t10d.log
