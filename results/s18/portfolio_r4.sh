#!/bin/zsh
# Stage-18 portfolio round 4 (11:50): claims variants. Waits for p11_prizeG2 to finish.
cd /Users/mickey/solarsystem/ctoc14
PY=~/.venvs/astro313/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
while pgrep -f "s18_rhfa.py run results/s18/p11_prizeG2" >/dev/null; do sleep 30; done
nice -n 5 $PY tools/s18_eval.py results/s18/p11_prizeG2 > results/s18/p11_prizeG2/eval.txt 2>&1
echo "[$(date +%H:%M:%S)] END p11_prizeG2: $(grep -E '^\[[0-9:]+\] \[slice [0-9]+ T=' results/s18/p11_prizeG2/run.log | tail -1 | cut -c1-110)"
run() { name=$1; shift; mkdir -p results/s18/$name; echo "[$(date +%H:%M:%S)] START $name $@"; nice -n 5 $PY tools/s18_rhfa.py run results/s18/$name --N 8 --proposer planbeam --nproc 9 "$@" > results/s18/$name/run.log 2>&1; echo "[$(date +%H:%M:%S)] END $name: $(grep -E '^\[[0-9:]+\] \[slice [0-9]+ T=' results/s18/$name/run.log | tail -1 | cut -c1-110)"; nice -n 5 $PY tools/s18_eval.py results/s18/$name > results/s18/$name/eval.txt 2>&1; }
run p13_claims_ttl4 --H 540 --L 360 --claims --claim-bonus 1.0 --claim-ttl 4 --seed 13
S18_PLANBEAM='{"width": 300, "plan_wall": 300, "job_wall": 420}' run p14_claims_w300 --H 540 --L 360 --claims --claim-bonus 1.0 --claim-ttl 4 --seed 14 --job-grace 900
run p15_claims_prize --H 540 --L 360 --claims --claim-bonus 1.0 --claim-ttl 4 --prize-json results/s18/prizes_p10.json --seed 15
run p16_claims_ttl9 --H 540 --L 360 --claims --claim-bonus 1.5 --claim-ttl 9 --seed 16
echo "[$(date +%H:%M:%S)] PORTFOLIO4 DONE"
