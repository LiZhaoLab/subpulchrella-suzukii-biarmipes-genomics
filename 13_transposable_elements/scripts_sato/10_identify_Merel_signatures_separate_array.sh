#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --job-name=TE_MerelIdent
#SBATCH --array=0-1%2
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only/logs/TE_MerelIdent.%A_%a.out
#SBATCH --error=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only/logs/TE_MerelIdent.%A_%a.err

set -euo pipefail

NEW=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only
POPTE2=/ru-auth/local/home/asato/scratch/software/PoPoolationTE2/PoPoolationTE2.jar

ANALYSES=(Dsub_L5_7 Dsuz_4pop)
ANALYSIS="${ANALYSES[$SLURM_ARRAY_TASK_ID]}"

INPUT="${NEW}/population_TE/ppileup/${ANALYSIS}.Merel.raw.ppileup.gz"
OUTDIR="${NEW}/population_TE/signatures"
CHECKDIR="${NEW}/population_TE/checks/identify_signatures"
OUTPUT="${OUTDIR}/${ANALYSIS}.Merel.raw.separate.signatures"
SUMMARY="${CHECKDIR}/${ANALYSIS}.Merel.raw.separate.summary.tsv"
TMP="${OUTDIR}/.${ANALYSIS}.Merel.raw.separate.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.signatures"

[[ -s "$POPTE2" && -s "$INPUT" ]] || {
    echo "ERROR: missing PoPoolationTE2 or ppileup input" >&2
    exit 1
}

mkdir -p "$OUTDIR" "$CHECKDIR"

[[ ! -e "$OUTPUT" && ! -e "$SUMMARY" && ! -e "$TMP" ]] || {
    echo "ERROR: output or temporary file already exists; will not overwrite" >&2
    exit 1
}

trap 'rm -f -- "$TMP"' EXIT

echo "analysis=$ANALYSIS"
echo "input_ppileup=$INPUT"
echo "mode=separate"
echo "min_count=2"
date

java -Xmx48g -jar "$POPTE2" identifySignatures \
    --ppileup "$INPUT" \
    --mode separate \
    --min-count 2 \
    --output "$TMP"

[[ -s "$TMP" ]] || {
    echo "ERROR: empty signature output" >&2
    exit 1
}

LINES=$(wc -l < "$TMP")
mv "$TMP" "$OUTPUT"

{
    printf 'metric\tvalue\n'
    printf 'analysis\t%s\n' "$ANALYSIS"
    printf 'input_ppileup\t%s\n' "$INPUT"
    printf 'mode\tseparate\n'
    printf 'minimum_average_TE_support\t2\n'
    printf 'output_signatures\t%s\n' "$OUTPUT"
    printf 'signature_lines\t%s\n' "$LINES"
    printf 'status\tPASS\n'
} > "$SUMMARY"

cat "$SUMMARY"
echo "FINAL: PASS"
date
