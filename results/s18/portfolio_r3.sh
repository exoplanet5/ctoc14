#!/bin/zsh
# Stage-18 portfolio round 3 (reordered 09:40): claims + miss-feedback iterations first, then long horizons.
cd /Users/mickey/solarsystem/ctoc14
PY=~/.venvs/astro313/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
while ! grep -q "PORTFOLIO2 DONE" results/s18/portfolio_r2.log 2>/dev/null; do sleep 60; done
run() { name=$1; shift; mkdir -p results/s18/$name; echo "[$(date +%H:%M:%S)] START $name $@"; nice -n 5 $PY tools/s18_rhfa.py run results/s18/$name --N 8 --proposer planbeam --nproc 9 "$@" > results/s18/$name/run.log 2>&1; echo "[$(date +%H:%M:%S)] END $name: $(grep -E '^\[[0-9:]+\] \[slice [0-9]+ T=' results/s18/$name/run.log | tail -1 | cut -c1-110)"; nice -n 5 $PY tools/s18_eval.py results/s18/$name > results/s18/$name/eval.txt 2>&1; }
run p10_claims --H 540 --L 360 --claims --claim-bonus 0.5 --claim-ttl 2 --seed 10
run p11_prizeG2 --H 540 --L 360 --prize-json results/s18/prizes_g2.json --seed 11
$PY tools/s18_rhfa.py misses-prize results/s18/p11_prizeG2 results/s18/prizes_p11.json --bonus 1.0 --prev results/s18/prizes_g2.json --decay 0.5
run p12_prizeIt2 --H 540 --L 360 --prize-json results/s18/prizes_p11.json --seed 12
S18_PLANBEAM='{"pot_mode":"full","pot_cap":10,"plan_wall":300,"job_wall":420}' run p8_H1080 --H 1080 --L 720 --gamma 0.15 --max-launch 4 --launch-stagger 90 --seed 8 --job-grace 900
echo "[$(date +%H:%M:%S)] PORTFOLIO3 DONE"
