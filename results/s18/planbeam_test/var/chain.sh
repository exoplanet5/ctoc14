#!/bin/zsh
cd /Users/mickey/solarsystem/ctoc14
PY=~/.venvs/astro313/bin/python
D=results/s18/planbeam_test
nice -n 5 $PY tools/s18_planbeam.py test-slice $D/var/base $D/../g1_base 3 --crafts 4,7 --nproc 2 > $D/var/base.out 2>&1
nice -n 5 $PY tools/s18_planbeam.py test-slice $D/var/wt3 $D/../g1_base 3 --crafts 4,7 --nproc 2 --w-t 3.0 > $D/var/wt3.out 2>&1
nice -n 5 $PY tools/s18_planbeam.py test-slice $D/var/w300 $D/../g1_base 3 --crafts 4,7 --nproc 2 --width 300 > $D/var/w300.out 2>&1
echo CHAIN_DONE > $D/var/done
