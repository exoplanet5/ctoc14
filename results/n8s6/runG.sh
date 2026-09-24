#!/bin/zsh
setopt nonomatch
cd /Users/mickey/solarsystem/ctoc14
PY=~/.venvs/astro313/bin/python
mkdir -p results/n8s6/G
for n in 10 06 09 05 03 13 12 01; do
  $PY tools/skel_fill.py results/n8/skel8b results/n8s6/G/$n --mode segmented --only $n \
      --portfolio 1000:0.15:45 --m0 1600 --vinf 4.0 --dvmax 1.2 --wp-dvmax 2.5 \
      --seg-k 50 --seg-depth 20 --max-tank 1150 --try 4 --lam 0.4 --nproc 8 \
      > results/n8s6/G/$n.out 2>&1
done
