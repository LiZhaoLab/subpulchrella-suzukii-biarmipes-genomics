#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --job-name=TE_MerelFreqDsuz
#SBATCH --array=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=192G
#SBATCH --time=2-00:00:00
#SBATCH --output=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only/logs/TE_MerelFreqDsuz.%A_%a.out
#SBATCH --error=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only/logs/TE_MerelFreqDsuz.%A_%a.err

set -euo pipefail

NEW=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only
POPTE2=/ru-auth/local/home/asato/scratch/software/PoPoolationTE2/PoPoolationTE2.jar
ANALYSIS=Dsuz_4pop

PPILEUP="${NEW}/population_TE/ppileup/${ANALYSIS}.Merel.raw.ppileup.gz"
SIGNATURES="${NEW}/population_TE/signatures/${ANALYSIS}.Merel.raw.separate.signatures"
IDENT_SUMMARY="${NEW}/population_TE/checks/identify_signatures/${ANALYSIS}.Merel.raw.separate.summary.tsv"

FREQDIR="${NEW}/population_TE/frequency"
FILTERDIR="${NEW}/population_TE/filtered_signatures"
CHECKDIR="${NEW}/population_TE/checks/frequency_filter"
FREQ="${FREQDIR}/${ANALYSIS}.Merel.raw.separate.freqsig"
FILTERED="${FILTERDIR}/${ANALYSIS}.Merel.raw.separate.Merel_threshold.freqsig"
SUMMARY="${CHECKDIR}/${ANALYSIS}.Merel.raw.separate.frequency_filter.summary.tsv"

TMP_FREQ="${FREQDIR}/.${ANALYSIS}.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.freqsig"
TMP_FILTER="${FILTERDIR}/.${ANALYSIS}.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.freqsig"

for FILE in "$POPTE2" "$PPILEUP" "$SIGNATURES" "$IDENT_SUMMARY"; do
    [[ -s "$FILE" ]] || {
        echo "ERROR: missing input: $FILE" >&2
        exit 1
    }
done

grep -q $'^status\tPASS$' "$IDENT_SUMMARY" || {
    echo "ERROR: identifySignatures summary is not PASS" >&2
    exit 1
}

mkdir -p "$FREQDIR" "$FILTERDIR" "$CHECKDIR"

for FILE in "$FREQ" "$FILTERED" "$SUMMARY" "$TMP_FREQ" "$TMP_FILTER"; do
    [[ ! -e "$FILE" ]] || {
        echo "ERROR: result already exists; will not overwrite: $FILE" >&2
        exit 1
    }
done

trap 'rm -f -- "$TMP_FREQ" "$TMP_FILTER"' EXIT

echo "analysis=$ANALYSIS"
echo "java_heap=160g"
date

java -Xmx160g -jar "$POPTE2" frequency \
    --ppileup "$PPILEUP" \
    --signature "$SIGNATURES" \
    --output "$TMP_FREQ"

[[ -s "$TMP_FREQ" ]] || {
    echo "ERROR: empty frequency output" >&2
    exit 1
}

java -Xmx160g -jar "$POPTE2" filterSignatures \
    --input "$TMP_FREQ" \
    --min-coverage 10 \
    --max-otherte-count 2 \
    --max-structvar-count 2 \
    --output "$TMP_FILTER"

[[ -s "$TMP_FILTER" ]] || {
    echo "ERROR: empty filtered output" >&2
    exit 1
}

RAW_LINES=$(wc -l < "$TMP_FREQ")
FILTERED_LINES=$(wc -l < "$TMP_FILTER")

mv "$TMP_FREQ" "$FREQ"
mv "$TMP_FILTER" "$FILTERED"

{
    printf 'metric\tvalue\n'
    printf 'analysis\t%s\n' "$ANALYSIS"
    printf 'input_ppileup\t%s\n' "$PPILEUP"
    printf 'input_signatures\t%s\n' "$SIGNATURES"
    printf 'frequency_output\t%s\n' "$FREQ"
    printf 'frequency_lines\t%s\n' "$RAW_LINES"
    printf 'filtered_output\t%s\n' "$FILTERED"
    printf 'filtered_lines\t%s\n' "$FILTERED_LINES"
    printf 'min_coverage\t10\n'
    printf 'max_otherte_count\t2\n'
    printf 'max_structvar_count\t2\n'
    printf 'status\tPASS_PAIRUP_QC_PENDING\n'
} > "$SUMMARY"

cat "$SUMMARY"
echo "FINAL: FREQUENCY AND FILTER COMPLETE; PAIRUP QC PENDING"
date
