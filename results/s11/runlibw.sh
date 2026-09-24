#!/bin/zsh
cd /Users/mickey/solarsystem/ctoc14
P=~/.venvs/astro313/bin/python
$P tools/carrier_optw.py results/s11/libw_rare.pkl --n 40000 --starts 160 --iters 100 --nproc 3 --seed 811 --eps 0.012 --weights results/s11/prices_rare.json
$P tools/carrier_optw.py results/s11/libw_cmpl.pkl --n 40000 --starts 160 --iters 100 --nproc 3 --seed 907 --eps 0.012 --weights results/s11/prices_cmpl.json
$P tools/carrier_optw.py results/s11/libw_r1.pkl  --n 40000 --starts 130 --iters 100 --nproc 3 --seed 1013 --eps 0.014 --weights results/s11/prices_rnd1.json
$P tools/carrier_optw.py results/s11/libw_r2.pkl  --n 40000 --starts 130 --iters 100 --nproc 3 --seed 1117 --eps 0.010 --weights results/s11/prices_rnd2.json
$P tools/carrier_optw.py results/s11/libw_r3.pkl  --n 40000 --starts 130 --iters 100 --nproc 3 --seed 1223 --eps 0.016 --weights results/s11/prices_rnd3.json
