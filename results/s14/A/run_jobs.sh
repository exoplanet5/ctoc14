#!/bin/zsh
# run a jobs file (one shell command per line) with at most 4 concurrent jobs (track A CPU budget)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
cd /Users/mickey/solarsystem/ctoc14
grep -v '^#' "$1" | grep -v '^$' | xargs -P ${2:-4} -I{} zsh -c '{}'
