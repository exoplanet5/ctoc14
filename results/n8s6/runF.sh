#!/bin/zsh
setopt nonomatch
cd /Users/mickey/solarsystem/ctoc14
PY=~/.venvs/astro313/bin/python
mkdir -p results/n8s6/F
for spec in 10:1000:50 06:1000:50 01:1000:50 10:600:40; do
  n=${spec%%:*}; rest=${spec#*:}; bw=${rest%%:*}; kk=${rest##*:}
  $PY tools/skel_fill.py results/n8/skel8b results/n8s6/F/${n}_b${bw}k${kk} --mode segmented --only $n \
      --columns results/n8s6/F/c_${n}_b${bw}k${kk}.jsonl --col-min 24 --pool-fracs 1.0 --pool-seeds 1 \
      --portfolio ${bw}:0.15:45 --m0 1600 --vinf 4.0 --dvmax 1.2 --wp-dvmax 2.5 \
      --seg-k $kk --seg-depth 20 --nproc 8 > results/n8s6/F/${n}_b${bw}k${kk}.out 2>&1
done
