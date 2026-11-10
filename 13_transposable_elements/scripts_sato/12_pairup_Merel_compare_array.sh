#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --job-name=TE_MerelPair
#SBATCH --array=0-1%2
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only/logs/TE_MerelPair.%A_%a.out
#SBATCH --error=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only/logs/TE_MerelPair.%A_%a.err

set -euo pipefail

NEW=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only
POPTE2=/ru-auth/local/home/asato/scratch/software/PoPoolationTE2/PoPoolationTE2.jar
HIER="${NEW}/population_TE/config/Merel.min500.hierarchy.order_resolved.tsv"

ANALYSES=(Dsub_L5_7 Dsuz_4pop)
SPECIES=(Dsub Dsuz)
ANALYSIS="${ANALYSES[$SLURM_ARRAY_TASK_ID]}"
SP="${SPECIES[$SLURM_ARRAY_TASK_ID]}"

REF="${NEW}/population_TE/reference/${SP}.Merel.TE_merged_reference.min500.fasta"
RAW="${NEW}/population_TE/frequency/${ANALYSIS}.Merel.raw.separate.freqsig"
FILTERED="${NEW}/population_TE/filtered_signatures/${ANALYSIS}.Merel.raw.separate.Merel_threshold.freqsig"
INPUT_SUMMARY="${NEW}/population_TE/checks/frequency_filter/${ANALYSIS}.Merel.raw.separate.frequency_filter.summary.tsv"

OUTDIR="${NEW}/population_TE/insertions"
CHECKDIR="${NEW}/population_TE/checks/pairup"
RAW_OUT="${OUTDIR}/${ANALYSIS}.Merel.raw.separate.unfiltered.teinsertions.txt"
FILTERED_OUT="${OUTDIR}/${ANALYSIS}.Merel.raw.separate.Merel_threshold.teinsertions.txt"
SUMMARY="${CHECKDIR}/${ANALYSIS}.Merel.pairup_comparison.tsv"

TMP_RAW="${OUTDIR}/.${ANALYSIS}.unfiltered.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt"
TMP_FILTERED="${OUTDIR}/.${ANALYSIS}.filtered.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt"
TMP_SUMMARY="${CHECKDIR}/.${ANALYSIS}.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.tsv"

for FILE in "$POPTE2" "$HIER" "$REF" "$RAW" "$FILTERED" "$INPUT_SUMMARY"; do
    [[ -s "$FILE" ]] || { echo "ERROR: missing input: $FILE" >&2; exit 1; }
done
grep -q $'^status\tPASS_PAIRUP_QC_PENDING$' "$INPUT_SUMMARY" || {
    echo "ERROR: frequency/filter summary is not PASS" >&2
    exit 1
}

mkdir -p "$OUTDIR" "$CHECKDIR"
for FILE in "$RAW_OUT" "$FILTERED_OUT" "$SUMMARY" \
            "$TMP_RAW" "$TMP_FILTERED" "$TMP_SUMMARY"; do
    [[ ! -e "$FILE" ]] || {
        echo "ERROR: output already exists; will not overwrite: $FILE" >&2
        exit 1
    }
done
trap 'rm -f -- "$TMP_RAW" "$TMP_FILTERED" "$TMP_SUMMARY"' EXIT

echo "analysis=$ANALYSIS"
echo "reference=$REF"
echo "hierarchy=$HIER"
echo "pair_distance_bp=-200_to_300"
date

java -Xmx48g -jar "$POPTE2" pairupSignatures \
    --signature "$RAW" \
    --ref-genome "$REF" \
    --hier "$HIER" \
    --min-distance -200 \
    --max-distance 300 \
    --output "$TMP_RAW"

java -Xmx48g -jar "$POPTE2" pairupSignatures \
    --signature "$FILTERED" \
    --ref-genome "$REF" \
    --hier "$HIER" \
    --min-distance -200 \
    --max-distance 300 \
    --output "$TMP_FILTERED"

[[ -s "$TMP_RAW" && -s "$TMP_FILTERED" ]] || {
    echo "ERROR: empty pairup output" >&2
    exit 1
}

{
    printf 'analysis\tcondition\ttotal\tFR\tF_only\tR_only\tFR_percent\tstatus\n'
    awk -F '\t' -v analysis="$ANALYSIS" -v condition="unfiltered" '
        $7=="FR" {fr++}
        $7=="F"  {f++}
        $7=="R"  {r++}
        END {
            total=fr+f+r
            printf "%s\t%s\t%d\t%d\t%d\t%d\t%.2f\tQC_PENDING\n",
                analysis, condition, total, fr, f, r,
                total ? 100*fr/total : 0
        }
    ' "$TMP_RAW"
    awk -F '\t' -v analysis="$ANALYSIS" -v condition="Merel_threshold" '
        $7=="FR" {fr++}
        $7=="F"  {f++}
        $7=="R"  {r++}
        END {
            total=fr+f+r
            printf "%s\t%s\t%d\t%d\t%d\t%d\t%.2f\tQC_PENDING\n",
                analysis, condition, total, fr, f, r,
                total ? 100*fr/total : 0
        }
    ' "$TMP_FILTERED"
} > "$TMP_SUMMARY"

mv "$TMP_RAW" "$RAW_OUT"
mv "$TMP_FILTERED" "$FILTERED_OUT"
mv "$TMP_SUMMARY" "$SUMMARY"

cat "$SUMMARY"
echo "FINAL: PAIRUP COMPARISON CREATED; QC PENDING"
date
