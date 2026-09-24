#!/bin/zsh
# Combine + validate fleet E (stage 16 best, exact raw J 13.7466) from a STABLE copy of its frags.
setopt nonomatch
cd /Users/mickey/solarsystem/ctoc14
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
P=~/.venvs/astro313/bin/python
echo "[$(date +%H:%M:%S)] combine start"
$P tools/validator2.py /dev/null >/dev/null 2>&1
$P tools/combine_submission.py results/CTOC14_Result_s16E.txt $(ls results/s16/submit_E/frag_r*.txt | sort -V)
echo "[$(date +%H:%M:%S)] combined"
$P tools/validator2.py results/CTOC14_Result_s16E.txt > results/validator2_s16E.log 2>&1
echo "[$(date +%H:%M:%S)] validated"
tail -3 results/validator2_s16E.log
md5 results/CTOC14_Result_s16E.txt
