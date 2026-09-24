#!/bin/zsh
# post-fill chain for results/n8/fill2: insert the leftovers, relocate passes, export, validate
setopt nonomatch
cd /Users/mickey/solarsystem/ctoc14
P=~/.venvs/astro313/bin/python
echo "=== refill (redistribute leftovers, skeletons fixed)"; date
$P tools/skel_fill.py results/n8/skel8b results/n8/refill2 --mode segmented --portfolio 300:0.15:45,300:0.15:90 --retry "" --nproc 8 --refill results/n8/fill2 --refill-rounds 2 || exit 1
grep -E "ACCEPTED|rejected|FINAL" results/n8/refill2/log.txt | tail -12
echo "=== grow (insert leftovers)"; date
$P tools/grow_skeletons.py results/n8/refill2 results/n8/grow2 --max-dtank 80 --max-tank 1150 --nproc 8 --dmax 0.4 --dmax-wide 0.8 || exit 1
tail -2 results/n8/grow2.log
echo "=== relocate passes"; date
$P tools/run_ialns.py search results/n8/grow2 results/n8/ials2 --nproc 8 --relocate --rounds 2 --dissolve-k 0 || exit 1
tail -2 results/n8/ials2.log
echo "=== export"; date
$P tools/run_ialns.py export results/n8/ials2 results/n8/efleet2 --nproc 6 || exit 1
if [ -n "$(ls results/n8/efleet2/frags/ 2>/dev/null | grep INVALID)" ]; then echo "invalid craft present"; fi
$P tools/combine_submission.py results/CTOC14_Result_n8a.txt $(ls results/n8/efleet2/frags/frag_*.txt | sort)
$P tools/validator2.py results/CTOC14_Result_n8a.txt > results/validator2_n8a.log 2>&1
tail -3 results/validator2_n8a.log; date
