#!/bin/zsh
cd /Users/mickey/solarsystem/ctoc14
PY=~/.venvs/astro313/bin/python
D=results/s18/planbeam_test
nice -n 5 $PY tools/s18_planbeam.py test-slice $D/s3 $D/../g1_base 3 --crafts 0,1,4,7 --nproc 2 --regrid > $D/s3.out 2>&1
nice -n 5 $PY tools/s18_planbeam.py test-slice $D/var/wt3w300 $D/../g1_base 3 --crafts 4,7 --nproc 2 --width 300 > $D/var/wt3w300.out 2>&1
echo CHAIN2_DONE > $D/var/done2
