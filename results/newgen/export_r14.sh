#!/bin/zsh
cd /Users/mickey/solarsystem/ctoc14
P=~/.venvs/astro313/bin/python
$P tools/run_ialns.py export results/newgen/ifleet_r14 results/newgen/efleet_r14 --nproc 2
if ls results/newgen/efleet_r14/frags/INVALID_* > /dev/null 2>&1; then echo "invalid craft present"; exit 1; fi
$P tools/combine_submission.py results/CTOC14_Result_r14.txt $(ls results/newgen/efleet_r14/frags/frag_*.txt | sort)
$P tools/validator2.py results/CTOC14_Result_r14.txt > results/validator2_r14.log 2>&1
tail -3 results/validator2_r14.log
