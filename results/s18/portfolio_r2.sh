#!/bin/zsh
# Stage-18 portfolio round 2: waits for round 1, then launch cap/stagger and full-plan potential variants.
cd /Users/mickey/solarsystem/ctoc14
PY=~/.venvs/astro313/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
while ! grep -q "PORTFOLIO DONE" results/s18/portfolio_r1.log 2>/dev/null; do sleep 60; done
run() { name=$1; shift; mkdir -p results/s18/$name; echo "[$(date +%H:%M:%S)] START $name $@"; nice -n 5 $PY tools/s18_rhfa.py run results/s18/$name --N 8 --proposer planbeam --nproc 9 "$@" > results/s18/$name/run.log 2>&1; echo "[$(date +%H:%M:%S)] END $name: $(grep -E '^\[[0-9:]+\] \[slice [0-9]+ T=' results/s18/$name/run.log | tail -1 | cut -c1-110)"; nice -n 5 $PY tools/s18_eval.py results/s18/$name > results/s18/$name/eval.txt 2>&1; }
run p5_stagger --H 540 --L 360 --max-launch 3 --launch-stagger 60 --seed 5
S18_PLANBEAM='{"pot_mode":"full","pot_cap":8}' run p6_fullpot --H 540 --L 360 --gamma 0.15 --seed 6
S18_PLANBEAM='{"pot_mode":"full","pot_cap":8}' run p7_both --H 540 --L 360 --gamma 0.15 --max-launch 3 --launch-stagger 60 --lam 0.7 --seed 7
echo "[$(date +%H:%M:%S)] PORTFOLIO2 DONE"
