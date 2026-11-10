#!/bin/bash
# Wrapper for campaign inputs "POOL CLS".
set -euo pipefail
unset LD_LIBRARY_PATH
set -- $*
/ru-auth/local/home/lzhao/vscode_projects/suzukii/envs/popgen/bin/python \
  /ru-auth/local/home/lzhao/vscode_projects/suzukii/scripts/100_sweep_null.py $1 $2
