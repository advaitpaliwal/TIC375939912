#!/bin/zsh
# MOLUSC 4 x 50k (RUWE + Gaia contrast), M = 0.74 Msun (TIC 8.2), default age; then TRICERATOPS 10 runs x {base, loc}; loc_molusc after summary.
F=work/planet-hunt/deep/375939912/fpp
PYM=work/planet-hunt/tic240488148-companions/.venv-molusc/bin/python
PYT=work/planet-hunt/tic240488148-companions/.venv/bin/python
cp $F/scripts/molusc_run.py $F/molusc_run.py
cd $F
for s in 11 12 13 14; do
  SEED=$s $PYM -I -W ignore molusc_run.py -v --ruwe --gaia -a $F/molusc_out/m50k_s$s 20h54m44.552s -- +06d35m25.28s 50000 0.74 > molusc_out/m50k_s$s.log 2>&1 &
  [[ $s == 12 ]] && wait
done
wait; echo MOLUSC_DONE >> $F/progress.log
cd $F/trc
for m in base loc; do for k in 0 1 2 3 4 5 6 7 8 9; do
  OMP_NUM_THREADS=2 $PYT -I -W ignore $F/scripts/run_fpp.py $m $k >> $F/progress.log 2>&1
done; done
cd $F; $PYT -I -W ignore $F/scripts/molusc_summary.py > $F/molusc_summary.txt 2>&1
cd $F/trc
for k in 0 1 2 3 4 5 6 7 8 9; do OMP_NUM_THREADS=2 $PYT -I -W ignore $F/scripts/run_fpp.py loc_molusc $k >> $F/progress.log 2>&1; done
echo ALL_DONE >> $F/progress.log
