#!/bin/zsh
# SR (segment re-plan) closer probe: every N6a route re-planned segment by segment with ALL its own flybys pinned as
# waypoints (+-45 d) and the 35 pass leftovers as prize-1 targets; one process per host (independent pricing).
cd /Users/mickey/solarsystem/ctoc14
PY=~/.venvs/astro313/bin/python
SRC=/private/tmp/claude-502/-Users-mickey-solarsystem-ctoc14/5fde8c2f-167b-40ff-92f3-bf04233cf9a2/scratchpad/passes/N6a
for r in 01 02 03 04 05 06 07 08; do
  SR_TBIN=4 SR_NPT=6 SR_MC=150 nice -n 5 $PY results/s13/plan13/sr_fill.py $SRC results/s13/plan13/sr_N6a/$r --only $r --mode segmented \
     --portfolio 300:0.15:45 --retry '' --seg-k 30 --nproc 1 --try 6 --lam 0.04 --max-tank 1300 \
     > results/s13/plan13/sr_N6a/$r.out 2>&1 &
done
wait
echo SR_DONE
