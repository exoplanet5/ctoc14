#!/bin/zsh
setopt nonomatch
cd /Users/mickey/solarsystem/ctoc14
PY=~/.venvs/astro313/bin/python
run() {  # route, drop-spec, tag
  $PY tools/skel_fill.py results/n8/skel8b results/n8s6/D/$1_$3 --mode segmented --only $1 --drop-wp $2 \
      --columns results/n8s6/D/c_$1_$3.jsonl --col-min 20 --pool-fracs 1.0 --pool-seeds 1 \
      --portfolio 300:0.15:45 --m0 1600 --vinf 4.0 --dvmax 1.2 --wp-dvmax 2.5 \
      --seg-k 30 --seg-depth 20 --nproc 1 > results/n8s6/D/$1_$3.out 2>&1
}
# free-beam anchor (no waypoints at all) for two routes
run 10 all free &
run 01 all free &
# drop-one on the two weakest routes
for w in 18 25 62 116 166 169 199 282; do run 10 $w w$w & done
wait
for w in 16 34 86 108 145 224 240 254; do run 06 $w w$w & done
for w in 8 33 39 51 157 168 217 219 251; do run 09 $w w$w & done
wait
