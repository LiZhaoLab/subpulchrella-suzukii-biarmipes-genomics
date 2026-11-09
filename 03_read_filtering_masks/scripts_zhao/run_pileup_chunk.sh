#!/bin/bash
# One campaign step: per-site allele counts for SAMPLE:CHROM:START-END.
set -euo pipefail
unset LD_LIBRARY_PATH
exec /ru-auth/local/home/lzhao/anaconda3/envs/jupyter-ml/bin/python \
    /ru-auth/local/home/lzhao/vscode_projects/suzukii/scripts/01_pileup_counts.py "$1"
