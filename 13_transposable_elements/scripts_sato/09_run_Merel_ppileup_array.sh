#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --job-name=TE_ppileup_Merel
#SBATCH --array=0-1%2
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00
#SBATCH --output=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only/logs/TE_ppileup_Merel.%A_%a.out
#SBATCH --error=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only/logs/TE_ppileup_Merel.%A_%a.err

set -euo pipefail

NEW=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only
POPTE2=/ru-auth/local/home/asato/scratch/software/PoPoolationTE2/PoPoolationTE2.jar
SAMTOOLS=/ru-auth/local/home/asato/scratch/conda_envs/gatk4/bin/samtools
HIER="${NEW}/population_TE/config/Merel.min500.hierarchy.order_resolved.tsv"
OUTDIR="${NEW}/population_TE/ppileup"
CHECKDIR="${NEW}/population_TE/checks/ppileup"

for FILE in "$POPTE2" "$SAMTOOLS" "$HIER"; do
    [[ -s "$FILE" ]] || { echo "ERROR: missing $FILE" >&2; exit 1; }
done

case "$SLURM_ARRAY_TASK_ID" in
    0)
        ANALYSIS=Dsub_L5_7
        SAMPLES=(Dsub_L5_7)
        ;;
    1)
        ANALYSIS=Dsuz_4pop
        SAMPLES=(JP-Kan JP-Sap JP-Sap2 JP-Tok)
        ;;
    *)
        echo "ERROR: unexpected array task ID" >&2
        exit 1
        ;;
esac

mkdir -p "$OUTDIR" "$CHECKDIR"
OUTPUT="${OUTDIR}/${ANALYSIS}.Merel.raw.ppileup.gz"
SUMMARY="${CHECKDIR}/${ANALYSIS}.Merel.ppileup.summary.tsv"
TMP="${OUTDIR}/.${ANALYSIS}.Merel.raw.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.ppileup.gz"

[[ ! -e "$OUTPUT" && ! -e "$SUMMARY" && ! -e "$TMP" ]] || {
    echo "ERROR: output already exists; will not overwrite" >&2
    exit 1
}
trap 'rm -f -- "$TMP"' EXIT

BAM_ARGS=()
for SAMPLE in "${SAMPLES[@]}"; do
    BAM="${NEW}/population_TE/bam_bwasw/${SAMPLE}/${SAMPLE}.Merel.bwasw.se2pe.sorted.bam"
    [[ -s "$BAM" && -s "${BAM}.bai" ]] || {
        echo "ERROR: missing BAM or index: $BAM" >&2
        exit 1
    }
    "$SAMTOOLS" quickcheck -v "$BAM"
    BAM_ARGS+=(--bam "$BAM")
done

echo "analysis=$ANALYSIS"
echo "sample_order=${SAMPLES[*]}"
echo "hierarchy=$HIER"
echo "map_qual=15"
echo "output=$OUTPUT"
date

java -Xmx112g -jar "$POPTE2" ppileup \
    "${BAM_ARGS[@]}" \
    --map-qual 15 \
    --hier "$HIER" \
    --output "$TMP"

[[ -s "$TMP" ]] || { echo "ERROR: empty ppileup" >&2; exit 1; }
gzip -t "$TMP"
mv "$TMP" "$OUTPUT"
trap - EXIT

SAMPLE_ORDER=$(IFS=,; echo "${SAMPLES[*]}")
{
    printf 'metric\tvalue\n'
    printf 'analysis\t%s\n' "$ANALYSIS"
    printf 'sample_order\t%s\n' "$SAMPLE_ORDER"
    printf 'number_of_samples\t%s\n' "${#SAMPLES[@]}"
    printf 'minimum_chromosome_mapping_quality\t15\n'
    printf 'hierarchy\t%s\n' "$HIER"
    printf 'output_ppileup\t%s\n' "$OUTPUT"
    printf 'status\tPASS\n'
} > "$SUMMARY"

cat "$SUMMARY"
echo "FINAL: PASS"
date
