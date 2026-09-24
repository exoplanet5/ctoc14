#!/bin/zsh
# Track-A scheduler: start each job of a jobs file only while fewer than MAXP of track A's CPU-bound Python processes
# (s14_* tools and scratch calibration scripts) are running -- orphans of earlier runners count too.
# Jobs whose output dir already holds result.json are skipped (restart-safe: unfinished jobs resume from ckpt.pkl).
# Usage: sched.sh JOBS_FILE [MAXP=4]
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
cd /Users/mickey/solarsystem/ctoc14
MAXP=${2:-4}
count() { ps -eo command | grep -E "^/[^ ]*/(Python|python[0-9.]*) .*(s14_|calib_res)" | grep -v grep | wc -l | tr -d ' ' }
grep -v '^#' "$1" | grep -v '^$' | while IFS= read -r job; do
  outdir=$(echo "$job" | awk '{for(i=1;i<=NF;i++) if ($i ~ /^results\/s14\/A\/.*\/(prod|rel|lt|ltw|g)[A-Za-z0-9_]*$/) {print $i; exit}}')
  if [ -n "$outdir" ] && [ -f "$outdir/result.json" ]; then echo "[$(date +%H:%M:%S)] skip (done) $outdir"; continue; fi
  while [ "$(count)" -ge "$MAXP" ]; do sleep 10; done
  echo "[$(date +%H:%M:%S)] start $outdir"
  zsh -c "$job" &
  sleep 5
done
wait
echo "[$(date +%H:%M:%S)] all jobs done"
