#!/bin/bash
# Analysis F, per pool: SweepFinder2 input (scripts/98_sweep_input.py), then the
# background folded spectrum per chromosome class with SweepFinder2 -f.
# Usage: 98_sweep_input.sh POOL
set -euo pipefail
unset LD_LIBRARY_PATH
P=/ru-auth/local/home/lzhao/vscode_projects/suzukii
O=/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC/sweeps/$1
/ru-auth/local/home/lzhao/anaconda3/envs/jupyter-ml/bin/python $P/scripts/98_sweep_input.py $1
for cls in autosome X; do
  $P/envs/popgen/bin/SweepFinder2 -f $O/$cls.all.sf2.txt $O/$cls.spect.txt
done
echo "$(date) done $1"
