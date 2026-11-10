#!/bin/bash
# One campaign step: polarized MK counts (and, for Dsub_L5_7, contamination evidence).
set -euo pipefail
unset LD_LIBRARY_PATH
exec /ru-auth/local/home/lzhao/anaconda3/envs/jupyter-ml/bin/python \
    /ru-auth/local/home/lzhao/vscode_projects/suzukii/scripts/60_polarized_mk.py "$1"
