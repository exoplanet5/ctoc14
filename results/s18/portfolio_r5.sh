#!/bin/zsh
# Stage-18 portfolio round 5: claims + wide + relaxed leg caps (spend the fuel slack on coverage). Waits for round 4.
cd /Users/mickey/solarsystem/ctoc14
PY=~/.venvs/astro313/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
while ! grep -q "PORTFOLIO4 DONE" results/s18/portfolio_r4.log 2>/dev/null; do sleep 60; done
run() { name=$1; shift; mkdir -p results/s18/$name; echo "[$(date +%H:%M:%S)] START $name $@"; nice -n 5 $PY tools/s18_rhfa.py run results/s18/$name --N 8 --proposer planbeam --nproc 9 "$@" > results/s18/$name/run.log 2>&1; echo "[$(date +%H:%M:%S)] END $name: $(grep -E '^\[[0-9:]+\] \[slice [0-9]+ T=' results/s18/$name/run.log | tail -1 | cut -c1-110)"; nice -n 5 $PY tools/s18_eval.py results/s18/$name > results/s18/$name/eval.txt 2>&1; }
S18_PLANBEAM='{"width": 300, "plan_wall": 300, "job_wall": 420, "dv": 1.6, "drmax": 0.2}' run p17_claims_relax --H 540 --L 360 --claims --claim-bonus 1.0 --claim-ttl 4 --seed 17 --job-grace 900
S18_PLANBEAM='{"width": 300, "plan_wall": 300, "job_wall": 420, "dv": 1.6, "drmax": 0.2}' run p18_claims_relax_lam05 --H 540 --L 360 --claims --claim-bonus 1.0 --claim-ttl 4 --lam 0.5 --seed 18 --job-grace 900
echo "[$(date +%H:%M:%S)] PORTFOLIO5 DONE"
