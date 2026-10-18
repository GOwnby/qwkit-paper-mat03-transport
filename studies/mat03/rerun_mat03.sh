#!/bin/bash
# Reproduce every MAT-03 number and figure (about 1.8 h on 2 cores).
set -e
cd "$(dirname "$0")"
export PYTHONPATH=$(cd ../.. && pwd) OMP_NUM_THREADS=1 NPROC=${NPROC:-2}
python3 -W ignore subbands.py > subbands.log
python3 -W ignore hall_algebra.py > hall_algebra.log
python3 -W ignore consistency.py > consistency.log
python3 -W ignore mr.py > mr.log
python3 -W ignore coulomb.py > coulomb.log
python3 -W ignore coupled.py > coupled.log
python3 -W ignore roughness.py > roughness.log
python3 -W ignore conv3.py > conv3.log
python3 -W ignore trapping.py > trapping.log
python3 -W ignore localization.py > localization.log
python3 collect_numbers.py > /dev/null
python3 tables.py > /dev/null
python3 -W ignore paper_figures.py
