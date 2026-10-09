#!/bin/zsh
# Cadence sensitivity: loc_molusc with exptime 600 s and 200 s (actual cadences), seeds 1000-1002.
F=work/planet-hunt/deep/375939912/fpp
PYT=work/planet-hunt/tic240488148-companions/.venv/bin/python
cd $F/trc
for e in 600 200; do for k in 0 1 2; do
  EXPTIME_S=$e OMP_NUM_THREADS=2 $PYT -I -W ignore $F/scripts/run_fpp_exp.py loc_molusc $k >> $F/progress.log 2>&1
done; done
echo EXP_DONE >> $F/progress.log
