#!/bin/zsh
# Track-A stage-2 scheduler: start each job line only while fewer than MAXP track-A CPU-bound Python processes
# (s14_* tools) are running.  A job whose output dir (the 3rd word after the tool path: tool MODE OUT) already holds
# result.json is skipped (restart-safe; unfinished free runs resume from ckpt.pkl).
# Usage: sched2.sh JOBS_FILE [MAXP=4]
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
cd /Users/mickey/solarsystem/ctoc14
MAXP=${2:-4}
count() { ps -eo command | grep -E "^/[^ ]*/(Python|python[0-9.]*) .*s14_" | grep -v grep | wc -l | tr -d ' ' }
grep -v '^#' "$1" | grep -v '^$' | while IFS= read -r job; do
  outdir=$(echo "$job" | awk '{for(i=1;i<=NF;i++) if ($i ~ /s14_[a-z_]*\.py$/) {print $(i+2); exit}}')
  if [ -n "$outdir" ] && [ -f "$outdir/result.json" ]; then echo "[$(date +%H:%M:%S)] skip (done) $outdir"; continue; fi
  while [ "$(count)" -ge "$MAXP" ]; do sleep 5; done
  echo "[$(date +%H:%M:%S)] start $outdir"
  zsh -c "$job" &
  sleep 2
done
wait
echo "[$(date +%H:%M:%S)] all jobs done"
