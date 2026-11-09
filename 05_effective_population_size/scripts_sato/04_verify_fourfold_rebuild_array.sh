#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --job-name=Ne_4foldQC
#SBATCH --array=0-1%2
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/ru-auth/local/home/asato/scratch/population_polymorphism/Ne_reanalysis_20261023_repeatmasker_p98/logs/Ne_4foldQC.%A_%a.out
#SBATCH --error=/ru-auth/local/home/asato/scratch/population_polymorphism/Ne_reanalysis_20261023_repeatmasker_p98/logs/Ne_4foldQC.%A_%a.err

set -euo pipefail

NEW=/ru-auth/local/home/asato/scratch/population_polymorphism/Ne_reanalysis_20261023_repeatmasker_p98
POLY=/ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821
PYTHON=/ru-auth/local/home/asato/scratch/conda_envs/earlgrey_env/bin/python3
BEDTOOLS=/ru-auth/local/home/asato/scratch/conda_envs/earlgrey_env/bin/bedtools
EXTRACTOR="${NEW}/scripts/provenance_0825/03_extract_fourfold_sites.py"
SPECIES_LIST=(Dsub Dsuz)

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is not set}"
if (( TASK_ID < 0 || TASK_ID > 1 )); then
    echo "ERROR: array task ID must be 0 or 1" >&2
    exit 1
fi
SP="${SPECIES_LIST[$TASK_ID]}"

CONFIG="${NEW}/config/species.tsv"
read -r FASTA GFF3 < <(
    awk -F '\t' -v sp="$SP" 'NR>1 && $1==sp {print $2, $3; exit}' "$CONFIG"
)

for INPUT in "$CONFIG" "$FASTA" "$GFF3" "$EXTRACTOR" "$PYTHON" "$BEDTOOLS"; do
    if [[ ! -s "$INPUT" ]]; then
        echo "ERROR: missing or empty input: $INPUT" >&2
        exit 1
    fi
done

OUTDIR="${NEW}/checks/fourfold_rebuild/${SP}"
if [[ -e "$OUTDIR" ]]; then
    echo "ERROR: QC directory already exists; will not overwrite: $OUTDIR" >&2
    exit 1
fi
mkdir -p "$OUTDIR"
PREFIX="${OUTDIR}/${SP}"

echo "species=$SP"
echo "reference=$FASTA"
echo "annotation=$GFF3"
echo "extractor=$EXTRACTOR"
echo "QC_directory=$OUTDIR"
date

"$PYTHON" "$EXTRACTOR" \
    --species "$SP" \
    --fasta "$FASTA" \
    --gff3 "$GFF3" \
    --output-prefix "$PREFIX"

for SUFFIX in \
    fourfold_5families.strict.bed \
    fourfold_5families.strict.sites.tsv \
    fourfold_5families.summary.tsv \
    selected_longest_CDS.tsv; do
    REBUILT="${PREFIX}.${SUFFIX}"
    SAVED="${NEW}/sites/${SP}/${SP}.${SUFFIX}"
    if [[ ! -s "$REBUILT" || ! -s "$SAVED" ]] || ! cmp -s "$REBUILT" "$SAVED"; then
        echo "ERROR: regenerated file differs: $SUFFIX" >&2
        exit 1
    fi
    echo "MATCH: $SUFFIX"
done

for CLASS in autosome X; do
    REGION="${POLY}/config/regions/${SP}.${CLASS}.bed"
    REBUILT="${OUTDIR}/${SP}.${CLASS}.fourfold.bed"
    SAVED="${NEW}/sites/analysis_sets/${SP}.${CLASS}.fourfold.bed"

    if [[ ! -s "$REGION" || ! -s "$SAVED" ]]; then
        echo "ERROR: missing region BED or saved BED: $SP $CLASS" >&2
        exit 1
    fi
    "$BEDTOOLS" intersect \
        -a "${PREFIX}.fourfold_5families.strict.bed" \
        -b "$REGION" \
        -wa > "$REBUILT"
    if [[ ! -s "$REBUILT" ]] || ! cmp -s "$REBUILT" "$SAVED"; then
        echo "ERROR: regenerated class BED differs: $SP $CLASS" >&2
        exit 1
    fi
    echo "MATCH: ${SP}.${CLASS}.fourfold.bed"
done

echo "FINAL: PASS"
date
