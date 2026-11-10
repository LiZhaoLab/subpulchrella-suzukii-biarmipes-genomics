#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --job-name=TE_MerelMap
#SBATCH --array=0-9%6
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only/logs/TE_MerelMap.%A_%a.out
#SBATCH --error=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only/logs/TE_MerelMap.%A_%a.err

set -euo pipefail

NEW=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only
CONFIG="${NEW}/population_TE/config/samples.tsv"
BWA=/ru-auth/local/home/asato/scratch/conda_envs/gatk4/bin/bwa
SAMTOOLS=/ru-auth/local/home/asato/scratch/conda_envs/gatk4/bin/samtools

TASK=${SLURM_ARRAY_TASK_ID}
[[ "$TASK" =~ ^[0-9]$ ]] || {
    echo "ERROR: invalid task ID: $TASK" >&2
    exit 1
}

ROW=$((TASK / 2 + 2))
MATE_NUMBER=$((TASK % 2 + 1))
IFS=$'\t' read -r SAMPLE SPECIES FASTQ1 FASTQ2 < <(
    awk -F '\t' -v row="$ROW" \
        'NR==row {print $1 "\t" $2 "\t" $3 "\t" $4}' "$CONFIG"
)

case "$SAMPLE:$SPECIES" in
    Dsub_L5_7:Dsub|JP-Kan:Dsuz|JP-Sap:Dsuz|JP-Sap2:Dsuz|JP-Tok:Dsuz) ;;
    *) echo "ERROR: unexpected sample or species at row $ROW" >&2; exit 1 ;;
esac

MATE="R${MATE_NUMBER}"
if [[ "$MATE_NUMBER" -eq 1 ]]; then
    FASTQ=$FASTQ1
else
    FASTQ=$FASTQ2
fi

REF="${NEW}/population_TE/reference/${SPECIES}.Merel.TE_merged_reference.min500.fasta"
INDEX_QC="${NEW}/checks/reference_index/${SPECIES}.Merel.index.summary.tsv"
OUTDIR="${NEW}/population_TE/mapping_single_end_bwasw/${SAMPLE}"
CHECKDIR="${NEW}/checks/mapping_single_end_bwasw/${SAMPLE}"
OUTPUT="${OUTDIR}/${SAMPLE}.${MATE}.bwasw.single_end.bam"
TEMP="${OUTPUT}.tmp.${SLURM_JOB_ID}_${TASK}"
FLAGSTAT="${CHECKDIR}/${SAMPLE}.${MATE}.bwasw.flagstat.txt"
SUMMARY="${CHECKDIR}/${SAMPLE}.${MATE}.bwasw.mapping_summary.tsv"

for FILE in "$CONFIG" "$FASTQ" "$REF" "${REF}.bwt" "$INDEX_QC" "$BWA" "$SAMTOOLS"; do
    [[ -s "$FILE" ]] || { echo "ERROR: missing input: $FILE" >&2; exit 1; }
done

awk -F '\t' -v ref="$REF" '
    $1=="reference" && $2==ref {matching_ref=1}
    $1=="status" && $2=="PASS" {passed=1}
    END {exit !(matching_ref && passed)}
' "$INDEX_QC" || {
    echo "ERROR: reference index QC did not pass" >&2
    exit 1
}

[[ ! -e "$OUTPUT" && ! -e "$SUMMARY" ]] || {
    echo "ERROR: output already exists; will not overwrite" >&2
    exit 1
}

mkdir -p "$OUTDIR" "$CHECKDIR"
trap 'rm -f -- "$TEMP"' EXIT

echo "sample=$SAMPLE"
echo "mate=$MATE"
echo "FASTQ=$FASTQ"
echo "reference=$REF"
echo "output=$OUTPUT"
date

"$BWA" bwasw -t 12 "$REF" "$FASTQ" |
    "$SAMTOOLS" view -@ 4 -b -o "$TEMP" -

"$SAMTOOLS" quickcheck -v "$TEMP"
ln -- "$TEMP" "$OUTPUT"
rm -- "$TEMP"

"$SAMTOOLS" flagstat -@ 4 "$OUTPUT" > "$FLAGSTAT"

{
    printf 'metric\tvalue\n'
    printf 'sample\t%s\n' "$SAMPLE"
    printf 'mate\t%s\n' "$MATE"
    printf 'FASTQ\t%s\n' "$FASTQ"
    printf 'reference\t%s\n' "$REF"
    printf 'output_BAM\t%s\n' "$OUTPUT"
    printf 'aligner\tbwa bwasw -t 12\n'
    printf 'status\tPASS\n'
} > "$SUMMARY"

cat "$SUMMARY"
echo "FINAL: PASS"
date
