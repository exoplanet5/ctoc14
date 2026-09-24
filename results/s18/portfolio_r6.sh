#!/bin/zsh
# Stage-18 round 6: diagnostics. N=9 claims (does coverage scale with N?), and short slices with claims.
cd /Users/mickey/solarsystem/ctoc14
PY=~/.venvs/astro313/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
while ! grep -q "PORTFOLIO5 DONE" results/s18/portfolio_r5.log 2>/dev/null; do sleep 60; done
run() { name=$1; N=$2; shift 2; mkdir -p results/s18/$name; echo "[$(date +%H:%M:%S)] START $name N=$N $@"; nice -n 5 $PY tools/s18_rhfa.py run results/s18/$name --N $N --proposer planbeam --nproc 9 "$@" > results/s18/$name/run.log 2>&1; echo "[$(date +%H:%M:%S)] END $name: $(grep -E '^\[[0-9:]+\] \[slice [0-9]+ T=' results/s18/$name/run.log | tail -1 | cut -c1-110)"; nice -n 5 $PY tools/s18_eval.py results/s18/$name > results/s18/$name/eval.txt 2>&1; }
run p19_claims_N9 9 --H 540 --L 360 --claims --claim-bonus 1.0 --claim-ttl 4 --fuel-bar 4300 --seed 19
run p20_claims_H360 8 --H 360 --L 360 --claims --claim-bonus 1.0 --claim-ttl 6 --seed 20
echo "[$(date +%H:%M:%S)] PORTFOLIO6 DONE"
