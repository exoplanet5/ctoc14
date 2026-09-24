#!/bin/zsh
# Stage-18 portfolio round 1 (planbeam proposer): waits for g2_plan, then runs 4 variants sequentially on 9 workers.
cd /Users/mickey/solarsystem/ctoc14
PY=~/.venvs/astro313/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
while [ ! -f results/s18/g2_plan/fleet/fleet.json ] && kill -0 $(cat results/s18/g2_plan/pid) 2>/dev/null; do sleep 30; done
run() { name=$1; shift; mkdir -p results/s18/$name; echo "[$(date +%H:%M:%S)] START $name $@"; nice -n 5 $PY tools/s18_rhfa.py run results/s18/$name --N 8 --proposer planbeam --nproc 9 "$@" > results/s18/$name/run.log 2>&1; echo "[$(date +%H:%M:%S)] END $name: $(grep -E '^\[[0-9:]+\] \[slice [0-9]+ T=' results/s18/$name/run.log | tail -1 | cut -c1-110)"; nice -n 5 $PY tools/s18_eval.py results/s18/$name > results/s18/$name/eval.txt 2>&1; }
run p1_ration --H 540 --L 360 --qf 7 --qm 12 --beta 2.5 --lam 0.7 --gamma 0.5 --seed 1
run p2_count  --H 540 --L 360 --qf 12 --qm 18 --beta 1.5 --lam 0.3 --gamma 0.6 --seed 2
S18_PLANBEAM='{"width": 300, "plan_wall": 300, "job_wall": 420}' run p3_wide --H 540 --L 360 --qf 7 --qm 12 --beta 2.5 --lam 0.7 --gamma 0.5 --seed 3 --job-grace 900
run p4_long   --H 720 --L 540 --qf 7 --qm 12 --beta 2.5 --lam 0.7 --gamma 0.5 --seed 4
echo "[$(date +%H:%M:%S)] PORTFOLIO DONE"
