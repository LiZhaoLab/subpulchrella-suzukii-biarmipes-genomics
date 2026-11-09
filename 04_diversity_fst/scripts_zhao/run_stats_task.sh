#!/bin/bash
# One campaign step: div:SAMPLE:CHROM:S-E | fst:CHROM:S-E | mk:ALIGN_DIR:INGROUP:NAME
set -euo pipefail
unset LD_LIBRARY_PATH
PY=/ru-auth/local/home/lzhao/anaconda3/envs/jupyter-ml/bin/python
S=/ru-auth/local/home/lzhao/vscode_projects/suzukii/scripts
case "$1" in
    div:*) exec $PY $S/32_chunk_stats.py div "${1#div:}" ;;
    fst:*) exec $PY $S/32_chunk_stats.py fst "${1#fst:}" ;;
    mk:*)  IFS=: read -r _ dir ing name <<< "$1"
           mkdir -p /ru-auth/local/home/lzhao/vscode_projects/suzukii/results/mk_recount
           exec $PY $S/20_mk_recount.py "$dir" "$ing" \
               /ru-auth/local/home/lzhao/vscode_projects/suzukii/results/mk_recount/$name.tsv ;;
    *) echo "unknown task $1" >&2; exit 2 ;;
esac
