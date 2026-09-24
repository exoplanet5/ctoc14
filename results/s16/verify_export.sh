#!/bin/zsh
# Stage 16 VERIFICATION export (phase 1 only): per-craft exact conversion of results/s16/best/fleet into results/s16/efleet
# (craft_*.npz, frags/frag_*.txt, fleet.json with the exact J). Does NOT combine, does NOT write any results/CTOC14_Result_*.txt.
setopt nonomatch
cd /Users/mickey/solarsystem/ctoc14
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
P=~/.venvs/astro313/bin/python
SRC=results/s16/best/fleet; OUT=results/s16/efleet
MD5=$(cat $SRC/route_*.npz | md5 -q)
echo "[$(date +%H:%M:%S)] verify export of $SRC (md5 $MD5)"
rm -rf $OUT; mkdir -p $OUT
nice -n 5 $P tools/run_ialns.py export $SRC $OUT --nproc ${NPROC:-9}
echo $MD5 > $OUT/SOURCE_MD5
echo "[$(date +%H:%M:%S)] verify export done"; ls $OUT/frags/
