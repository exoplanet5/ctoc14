#!/bin/zsh
# four parallel ACO campaigns with different exploration settings
cd /Users/mickey/solarsystem/ctoc14
PY=~/.venvs/astro313/bin/python
run() {  # name depth maxtank alpha beta rho q0 nmax
  nohup $PY tools/fleet_aco.py run --out results/s10/aco/$1 --gens 150 --ants 96 --chunk 12 --nproc 2 \
    --depth $2 --maxtank $3 --alpha $4 --beta $5 --rho $6 --q0 $7 --nmax $8 --seed $RANDOM \
    > results/s10/aco/$1.out 2>&1 &
}
run A 45 1200 1.0 2.0 0.10 0.5 12
run B 40 1150 1.0 3.0 0.15 0.7 12
run C 45 1200 2.0 1.0 0.05 0.3 12
run D 38 1100 1.0 2.0 0.10 0.5 14
wait
