#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --array=0-3%4
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH -J syn_RM4sp
#SBATCH -o logs/syn_RM4sp.%A_%a.out
#SBATCH -e logs/syn_RM4sp.%A_%a.err

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/synteny_4species/reanalysis_20260824

RM_ENV=/ru-auth/local/home/asato/scratch/conda_envs/repeatmodeler
PY_ENV=/ru-auth/local/home/asato/scratch/conda_envs/earlgrey_env
RM_DIR=${RM_ENV}/share/RepeatMasker
RM=${RM_DIR}/RepeatMasker

export PATH="${PY_ENV}/bin:${RM_ENV}/bin:/usr/bin:/bin"

SPECIES=(
  Dsuz
  Dsub
  Dbia
  Dmel
)

INPUTS=(
  "${BASE}/input_genomes/Dsuz.GCF_043229965.1.fna"
  "${BASE}/input_genomes/Dsub.GCF_014743375.2.fna"
  "${BASE}/input_genomes/Dbia.GCF_025231255.1.fna"
  "${BASE}/input_genomes/Dmel.GCF_000001215.4.fna"
)

SPECIES_NAME=${SPECIES[$SLURM_ARRAY_TASK_ID]}
INPUT=${INPUTS[$SLURM_ARRAY_TASK_ID]}
OUTDIR=${BASE}/repeatmasker/${SPECIES_NAME}

mkdir -p "$OUTDIR"

[[ -s "$INPUT" ]] || {
    echo "ERROR: missing input: $INPUT" >&2
    exit 1
}

[[ -s "${RM_DIR}/Libraries/famdb/dfam39_full.1.h5" ]] || {
    echo "ERROR: Dfam Brachycera partition is unavailable" >&2
    exit 1
}

echo "===== RepeatMasker common synteny mask ====="
echo "species=${SPECIES_NAME}"
echo "input=${INPUT}"
echo "output=${OUTDIR}"
echo "RepeatMasker_version=$("$RM" -v 2>&1 | head -1)"
echo "Dfam_version=3.9"
echo "Dfam_partition=1_Brachycera"
echo "masking=hard"
echo "sensitivity=default"
echo "search_engine=rmblast"
date

"$RM" \
  -e rmblast \
  -species 7215 \
  -pa 4 \
  -gff \
  -dir "$OUTDIR" \
  "$INPUT"

MASKED=${OUTDIR}/$(basename "$INPUT").masked

[[ -s "$MASKED" ]] || {
    echo "ERROR: masked FASTA missing: $MASKED" >&2
    exit 1
}

ls -lh "$MASKED"

echo "FINAL: PASS"
date
