#!/bin/zsh
# N=9 fallback: fill the 9 rebalanced skeletons, insert leftovers, relocate, export, validate
setopt nonomatch
cd /Users/mickey/solarsystem/ctoc14
P=~/.venvs/astro313/bin/python
echo "=== fill (9 skeletons)"; date
$P tools/skel_fill.py results/n8/skel9b results/n9/fill1 --mode segmented --portfolio 300:0.15:45,300:0.15:90 --retry "" --nproc 8 --lam 0.1 || exit 1
grep -E "FINAL|sum J_i" results/n9/fill1/log.txt
echo "=== grow (insert leftovers)"; date
$P tools/grow_skeletons.py results/n9/fill1 results/n9/grow1 --max-dtank 80 --max-tank 1150 --nproc 8 --dmax 0.4 --dmax-wide 0.8 || exit 1
tail -2 results/n9/grow1.log
echo "=== export (pre-polish, bank it)"; date
$P tools/run_ialns.py export results/n9/grow1 results/n9/efleet1 --nproc 6 || exit 1
$P tools/combine_submission.py results/CTOC14_Result_n9a.txt $(ls results/n9/efleet1/frags/frag_*.txt | sort)
$P tools/validator2.py results/CTOC14_Result_n9a.txt > results/validator2_n9a.log 2>&1
tail -3 results/validator2_n9a.log; date
echo "=== relocate passes"; date
$P tools/run_ialns.py search results/n9/grow1 results/n9/ials1 --nproc 8 --relocate --rounds 2 --dissolve-k 0 || exit 1
tail -2 results/n9/ials1.log
echo "=== export (polished)"; date
$P tools/run_ialns.py export results/n9/ials1 results/n9/efleet2 --nproc 6 || exit 1
$P tools/combine_submission.py results/CTOC14_Result_n9b.txt $(ls results/n9/efleet2/frags/frag_*.txt | sort)
$P tools/validator2.py results/CTOC14_Result_n9b.txt > results/validator2_n9b.log 2>&1
tail -3 results/validator2_n9b.log; date
