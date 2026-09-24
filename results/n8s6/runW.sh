#!/bin/zsh
setopt nonomatch
cd /Users/mickey/solarsystem/ctoc14
PY=~/.venvs/astro313/bin/python
for h in 45 90 200 400 800; do
  for n in 12 06; do
    $PY tools/skel_fill.py results/n8/skel8b results/n8s6/W/${n}_h$h --mode segmented --only $n \
        --portfolio 300:0.15:$h --m0 1600 --vinf 4.0 --dvmax 1.2 --wp-dvmax 2.5 \
        --seg-k 30 --seg-depth 20 --max-tank 1200 --try 3 --nproc 8 > results/n8s6/W/${n}_h$h.out 2>&1
  done
done
