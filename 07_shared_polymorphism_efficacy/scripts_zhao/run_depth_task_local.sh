#!/bin/bash
# Dispatcher for depth analyses that import scripts/depth_common.py (run from scripts/).
set -euo pipefail
unset LD_LIBRARY_PATH
set -- $*
P=/ru-auth/local/home/lzhao/vscode_projects/suzukii
cd $P/scripts
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 $P/envs/popgen/bin/python "$@"
