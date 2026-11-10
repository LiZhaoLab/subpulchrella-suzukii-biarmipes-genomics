#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --job-name=TE_MerelIdx
#SBATCH --array=0-1%2
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only/logs/TE_MerelIdx.%A_%a.out
#SBATCH --error=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only/logs/TE_MerelIdx.%A_%a.err

set -euo pipefail

NEW=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only
BWA=/ru-auth/local/home/asato/scratch/conda_envs/gatk4/bin/bwa
SPECIES=(Dsub Dsuz)
INDEX=${SLURM_ARRAY_TASK_ID}

[[ "$INDEX" =~ ^[01]$ ]] || { echo "ERROR: invalid task ID" >&2; exit 1; }
SAMPLE=${SPECIES[$INDEX]}

REF="${NEW}/population_TE/reference/${SAMPLE}.Merel.TE_merged_reference.min500.fasta"
REF_QC="${NEW}/checks/reference_merge/${SAMPLE}.Merel.TE_merged_reference.min500.summary.tsv"
CHECKDIR="${NEW}/checks/reference_index"
SUMMARY="${CHECKDIR}/${SAMPLE}.Merel.index.summary.tsv"

for FILE in "$BWA" "$REF" "$REF_QC"; do
    [[ -s "$FILE" ]] || { echo "ERROR: missing input: $FILE" >&2; exit 1; }
done

awk -F '\t' -v ref="$REF" '
    $1=="merged_reference" && $2==ref {matching_ref=1}
    $1=="status" && $2=="PASS" {passed=1}
    END {exit !(matching_ref && passed)}
' "$REF_QC" || {
    echo "ERROR: merged-reference QC did not pass" >&2
    exit 1
}

mkdir -p "$CHECKDIR"
[[ ! -e "$SUMMARY" ]] || {
    echo "ERROR: summary already exists" >&2
    exit 1
}
for EXT in amb ann bwt pac sa; do
    [[ ! -e "${REF}.${EXT}" ]] || {
        echo "ERROR: index file already exists: ${REF}.${EXT}" >&2
        exit 1
    }
done

echo "sample=$SAMPLE"
echo "reference=$REF"
date

"$BWA" index -a bwtsw "$REF"

for EXT in amb ann bwt pac sa; do
    [[ -s "${REF}.${EXT}" ]] || {
        echo "ERROR: missing index file: ${REF}.${EXT}" >&2
        exit 1
    }
done

{
    printf 'metric\tvalue\n'
    printf 'sample\t%s\n' "$SAMPLE"
    printf 'reference\t%s\n' "$REF"
    printf 'index_algorithm\tbwtsw\n'
    printf 'status\tPASS\n'
} > "$SUMMARY"

cat "$SUMMARY"
echo "FINAL: PASS"
date
