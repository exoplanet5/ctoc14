#!/bin/zsh
setopt nonomatch
cd /Users/mickey/solarsystem/ctoc14
PY=~/.venvs/astro313/bin/python
mkdir -p results/n8s6/B
for n in 05 09 13 03 01 06 12 10; do
  $PY tools/skel_fill.py results/n8/skel8b results/n8s6/B/$n --mode segmented --only $n \
      --columns results/n8s6/B/cols_$n.jsonl --col-min 24 \
      --pool-fracs 1.0,0.7,0.5,0.35 --pool-seeds 5 \
      --portfolio 300:0.15:45,300:0.15:90 --m0 1600 --vinf 4.0 --dvmax 1.2 --wp-dvmax 2.5 \
      --seg-k 30 --seg-depth 14 --nproc 1 > results/n8s6/B/$n.out 2>&1 &
done
wait
