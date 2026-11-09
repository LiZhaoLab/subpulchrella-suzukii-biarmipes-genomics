#!/bin/bash
# Dispatcher for the depth analyses run as campaign inputs ("SCRIPT ARGS...").
set -euo pipefail
unset LD_LIBRARY_PATH
set -- $*
P=/ru-auth/local/home/lzhao/vscode_projects/suzukii
cd $P
script=$1; shift
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 $P/envs/popgen/bin/python scripts/$script "$@"
