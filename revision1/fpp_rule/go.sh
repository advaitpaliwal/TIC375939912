#!/bin/sh
H=$(cd "$(dirname "$0")" && pwd); PY=~/Scratch/planet-hunt/tic240488148-companions/.venv/bin/python
seq 0 9 | xargs -P 10 -I '{}' sh -c 'cd "$1"/out && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "$2" -I -W ignore "$1"/run.py "$3" > "$1"/out/run_sel_"$3".log 2>&1' sh "$H" "$PY" '{}'
