#!/bin/zsh
# Track B stage 2 (generations to G0.5).  Restart-safe: re-running continues from results/s13/gen/{gen.json,gen2.json}.
# Before re-running, check that no s13_gen2/s13_gen process is alive (ps -Ao pid,command | grep s13_gen); orphan
# s13_pass slots are waited for by s13_gen itself (slot .pid files).
# Slots run tools/s13_pass2.py (24-level front) with premium-set-only jitter; T7 dropped; value choice on the even slots.
cd /Users/mickey/solarsystem/ctoc14
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1 S13_JITTER_MODE=prem
nohup nice -n 5 ~/.venvs/astro313/bin/python tools/s13_gen2.py results/s13/gen --gens 3 \
  --alt full2,full4,full6,tail4,pilot2 --alt-args="--choice value" --deadline "2026-09-23 06:00" \
  --pass-tool tools/s13_pass2.py "--pass-args=--members T1,T2,T3,T5,T6" >> results/s13/gen/gen2.out 2>&1 &
echo $! > results/s13/gen/gen2.pid
