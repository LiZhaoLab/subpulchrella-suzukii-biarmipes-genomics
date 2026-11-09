#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH -J dfam39_install
#SBATCH -o logs/dfam39_install.%j.out
#SBATCH -e logs/dfam39_install.%j.err

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/synteny_4species/reanalysis_20260824
DBDIR=${BASE}/databases/Dfam_3.9
RM_FAMDB=/ru-auth/local/home/asato/scratch/conda_envs/repeatmodeler/share/RepeatMasker/Libraries/famdb

GZ=${DBDIR}/dfam39_full.1.h5.gz
H5=${DBDIR}/dfam39_full.1.h5
LINK=${RM_FAMDB}/dfam39_full.1.h5

[[ -s "$GZ" ]] || {
    echo "ERROR: missing download: $GZ" >&2
    exit 1
}

echo "===== Decompress Brachycera partition ====="
date

gzip -dc "$GZ" > "${H5}.tmp.${SLURM_JOB_ID}"
mv "${H5}.tmp.${SLURM_JOB_ID}" "$H5"

[[ -s "$H5" ]] || {
    echo "ERROR: decompressed database is empty" >&2
    exit 1
}

ln -sfn "$H5" "$LINK"

ls -lh "$GZ" "$H5" "$LINK"

echo "FINAL: PASS"
date
