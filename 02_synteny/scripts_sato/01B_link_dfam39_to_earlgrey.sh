#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH -J dfam39_EG
#SBATCH -o logs/dfam39_EG.%j.out
#SBATCH -e logs/dfam39_EG.%j.err

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/synteny_4species/reanalysis_20260824

SOURCE=${BASE}/databases/Dfam_3.9/dfam39_full.1.h5
EG_FAMDB=/ru-auth/local/home/asato/scratch/conda_envs/earlgrey_env/share/RepeatMasker/Libraries/famdb
LINK=${EG_FAMDB}/dfam39_full.1.h5

FAMDB_PY=/ru-auth/local/home/asato/scratch/conda_envs/earlgrey_env/share/RepeatMasker/famdb.py
PYTHON=/ru-auth/local/home/asato/scratch/conda_envs/earlgrey_env/bin/python3

echo "===== Link Dfam 3.9 Brachycera partition to Earl Grey ====="
echo "source=${SOURCE}"
echo "link=${LINK}"

[[ -s "${SOURCE}" ]] || {
    echo "ERROR: source partition is missing or empty: ${SOURCE}" >&2
    exit 1
}

[[ -d "${EG_FAMDB}" ]] || {
    echo "ERROR: Earl Grey FamDB directory is missing: ${EG_FAMDB}" >&2
    exit 1
}

ln -sfn "${SOURCE}" "${LINK}"

[[ -L "${LINK}" ]] || {
    echo "ERROR: symbolic link was not created: ${LINK}" >&2
    exit 1
}

[[ "$(readlink -f "${LINK}")" == "$(readlink -f "${SOURCE}")" ]] || {
    echo "ERROR: symbolic link points to the wrong file" >&2
    exit 1
}

echo "===== Earl Grey FamDB information ====="

"${PYTHON}" "${FAMDB_PY}" \
    -i "${EG_FAMDB}" \
    info

echo "FINAL: PASS"
