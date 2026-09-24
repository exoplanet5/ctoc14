#!/bin/zsh
# Stage 16 EXPORT (for the Export stage ONLY): results/s16/best/fleet -> results/CTOC14_Result_s16a.txt + validator2 log.
# Re-uses results/s16/efleet when its SOURCE_MD5 matches the current best fleet, else re-runs the per-craft export.
# Never touches results/CTOC14_Result_t10d.txt / _s15a.txt or their validator logs.
setopt nonomatch
cd /Users/mickey/solarsystem/ctoc14
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
P=~/.venvs/astro313/bin/python
SRC=results/s16/best/fleet; OUT=results/s16/efleet
MD5=$(cat $SRC/route_*.npz | md5 -q)
if [ "$(cat $OUT/SOURCE_MD5 2>/dev/null)" != "$MD5" ]; then
  echo "[$(date +%H:%M:%S)] efleet stale or missing -> export"; zsh results/s16/verify_export.sh
fi
if [ -n "$(ls $OUT/frags/ | grep INVALID)" ]; then echo "INVALID craft present -- not combining"; exit 1; fi
$P tools/combine_submission.py results/CTOC14_Result_s16a.txt $(ls $OUT/frags/frag_*.txt | sort)
echo "[$(date +%H:%M:%S)] combined"
$P tools/validator2.py results/CTOC14_Result_s16a.txt > results/validator2_s16a.log 2>&1
echo "[$(date +%H:%M:%S)] validated"; tail -8 results/validator2_s16a.log
