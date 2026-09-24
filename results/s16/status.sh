#!/bin/zsh
# stage-16 status: pool, selections, pipelines (last log line each), CPU-bound process count
cd /Users/mickey/solarsystem/ctoc14
setopt nonomatch
echo "== $(date +%H:%M:%S)  CPU-bound python: $(ps -Ao pcpu,command | grep -iE 'Python.app' | awk '$1>50' | wc -l | tr -d ' ')"
tail -1 results/s16/honest_pool.out 2>/dev/null
for d in results/s16/select/*/; do echo "select $(basename $d): $(tail -1 $d/log.txt 2>/dev/null | cut -c1-220)"; done
for d in results/s16/pipe/*/; do echo "pipe $(basename $d): $(tail -1 $d/log.txt 2>/dev/null | cut -c1-220)"; done
