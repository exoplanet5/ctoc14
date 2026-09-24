#!/bin/zsh
cd /Users/mickey/solarsystem/ctoc14
P=~/.venvs/astro313/bin/python
# seed / eps / lam grid -> geometrically diverse optimised carriers
$P tools/carrier_opt.py results/s11/lib_a.pkl --n 60000 --starts 220 --iters 110 --nproc 4 --seed 211 --eps 0.012
$P tools/carrier_opt.py results/s11/lib_b.pkl --n 60000 --starts 220 --iters 110 --nproc 4 --seed 307 --eps 0.018
$P tools/carrier_opt.py results/s11/lib_c.pkl --n 60000 --starts 220 --iters 110 --nproc 4 --seed 419 --eps 0.008 --lam 0.7
$P tools/carrier_opt.py results/s11/lib_d.pkl --n 60000 --starts 220 --iters 110 --nproc 4 --seed 523 --eps 0.015 --lam 1.4
$P tools/carrier_opt.py results/s11/lib_e.pkl --n 60000 --starts 220 --iters 110 --nproc 4 --seed 631 --eps 0.012 --kgap 6
$P tools/carrier_opt.py results/s11/lib_f.pkl --n 60000 --starts 220 --iters 110 --nproc 4 --seed 733 --eps 0.020 --kphase 1.0
