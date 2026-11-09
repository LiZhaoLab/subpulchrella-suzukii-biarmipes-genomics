#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --array=0-1%2
#SBATCH -J Ne_fourfold
#SBATCH -o logs/Ne_fourfold.%A_%a.out
#SBATCH -e logs/Ne_fourfold.%A_%a.err

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/population_polymorphism/Ne_reanalysis_20260825
CONFIG=${BASE}/config/species.tsv
PYTHON=/ru-auth/local/home/asato/scratch/conda_envs/earlgrey_env/bin/python3
SCRIPT=${BASE}/scripts/03_extract_fourfold_sites.py

SPECIES_LIST=(Dsub Dsuz)
SPECIES=${SPECIES_LIST[$SLURM_ARRAY_TASK_ID]}

read -r FASTA GFF3 < <(
    awk -F $'\t' -v species="${SPECIES}" '
        NR > 1 && $1 == species {
            print $2, $3
            exit
        }
    ' "${CONFIG}"
)

[[ -n "${FASTA:-}" && -s "${FASTA}" ]] || {
    echo "ERROR: reference FASTA missing for ${SPECIES}" >&2
    exit 1
}

[[ -n "${GFF3:-}" && -s "${GFF3}" ]] || {
    echo "ERROR: GFF3 missing for ${SPECIES}" >&2
    exit 1
}

OUTDIR=${BASE}/sites/${SPECIES}
PREFIX=${OUTDIR}/${SPECIES}

mkdir -p "${OUTDIR}" "${BASE}/logs"

echo "===== Extract conservative four-fold sites ====="
echo "species=${SPECIES}"
echo "reference=${FASTA}"
echo "annotation=${GFF3}"
echo "definition=third positions of CCN, ACN, GTN, GCN and GGN codons"
echo "transcript_rule=longest valid complete CDS per gene"

"${PYTHON}" "${SCRIPT}" \
    --species "${SPECIES}" \
    --fasta "${FASTA}" \
    --gff3 "${GFF3}" \
    --output-prefix "${PREFIX}"

BED=${PREFIX}.fourfold_5families.strict.bed
SUMMARY=${PREFIX}.fourfold_5families.summary.tsv

[[ -s "${BED}" ]] || {
    echo "ERROR: output BED is missing or empty" >&2
    exit 1
}

[[ -s "${SUMMARY}" ]] || {
    echo "ERROR: output summary is missing or empty" >&2
    exit 1
}

awk '
    NF != 3 ||
    $2 !~ /^[0-9]+$/ ||
    $3 !~ /^[0-9]+$/ ||
    $2 < 0 ||
    $3 <= $2 {
        invalid++
    }
    END {
        if (invalid > 0) {
            print "ERROR: invalid BED rows=" invalid > "/dev/stderr"
            exit 1
        }
    }
' "${BED}"

echo "output_BED=${BED}"
echo "summary=${SUMMARY}"
echo "FINAL: PASS"
