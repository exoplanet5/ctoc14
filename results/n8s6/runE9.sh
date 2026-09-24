#!/bin/zsh
setopt nonomatch
cd /Users/mickey/solarsystem/ctoc14
PY=~/.venvs/astro313/bin/python
mkdir -p results/n8s6/E9
for n in $(ls results/n8/skel9b/route_*.npz | sed 's/.*route_//;s/\.npz//'); do
  $PY tools/skel_fill.py results/n8/skel9b results/n8s6/E9/$n --mode segmented --only $n \
      --columns results/n8s6/E9/cols_$n.jsonl --col-min 20 --pool-fracs 1.0 --pool-seeds 1 \
      --portfolio 300:0.15:45 --m0 1600 --vinf 4.0 --dvmax 1.2 --wp-dvmax 2.5 \
      --seg-k 30 --seg-depth 20 --nproc 1 > results/n8s6/E9/$n.out 2>&1 &
done
wait
