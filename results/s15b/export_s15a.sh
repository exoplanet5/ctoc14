#!/bin/zsh
# Export of the stage-15b closed 9-craft fleet (twin sum J_i 11.9957, 298 covered). Mirrors results/newgen/export_t10d.sh
# with NEW paths only. Never touches results/CTOC14_Result_t10d.txt or results/validator2_t10d.log.
setopt nonomatch
cd /Users/mickey/solarsystem/ctoc14
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
P=~/.venvs/astro313/bin/python
echo "[$(date +%H:%M:%S)] export start"
nice -n 5 $P tools/run_ialns.py export results/s15b/best_closed/fleet results/s15b/efleet --nproc 9
echo "[$(date +%H:%M:%S)] export done"
ls results/s15b/efleet/frags/
if [ -n "$(ls results/s15b/efleet/frags/ | grep INVALID)" ]; then echo "INVALID craft present -- not combining"; exit 1; fi
$P tools/combine_submission.py results/CTOC14_Result_s15a.txt $(ls results/s15b/efleet/frags/frag_*.txt | sort)
echo "[$(date +%H:%M:%S)] combined"
$P tools/validator2.py results/CTOC14_Result_s15a.txt > results/validator2_s15a.log 2>&1
echo "[$(date +%H:%M:%S)] validated"
tail -8 results/validator2_s15a.log
