#!/bin/zsh
cd /Users/mickey/solarsystem/ctoc14
PY=~/.venvs/astro313/bin/python
for k in 1 2 3 4 5 6; do
  $PY tools/route_cover.py 'results/s7/farm*/route_*.npz' 'results/s7/cg[0-9]*/route_*.npz' 'results/**/route_*.npz' \
      --out results/s7/cgL$k --duals 2>&1 | tail -3
  $PY tools/carrier_farm.py results/s7/lib0.pkl results/s7/cgF$k --nroute 40 --frac 1.0 --ntrial 12 --nfail 3 \
      --tmax 1200 --top 150 --seed $((100+k)) --duals results/s7/cgL$k/duals.json --max-tank 1150 2>&1 | tail -2
done
$PY tools/route_cover.py 'results/s7/farm*/route_*.npz' 'results/s7/cgF*/route_*.npz' 'results/**/route_*.npz' \
    --out results/s7/cgFinal --duals 2>&1 | tail -3
