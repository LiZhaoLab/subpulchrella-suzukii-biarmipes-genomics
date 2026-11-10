#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --job-name=TE_se2pe_Merel
#SBATCH --array=0-4%2
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00
#SBATCH --output=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only/logs/TE_se2pe_Merel.%A_%a.out
#SBATCH --error=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only/logs/TE_se2pe_Merel.%A_%a.err

set -euo pipefail

NEW=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only
CONFIG="${NEW}/population_TE/config/samples.tsv"
POPTE2=/ru-auth/local/home/asato/scratch/software/PoPoolationTE2/PoPoolationTE2.jar
SAMTOOLS=/ru-auth/local/home/asato/scratch/conda_envs/gatk4/bin/samtools

ROW=$((SLURM_ARRAY_TASK_ID + 2))
IFS=$'\t' read -r SAMPLE SPECIES FASTQ1 FASTQ2 < <(
    awk -F '\t' -v row="$ROW" \
        'NR == row {print $1 "\t" $2 "\t" $3 "\t" $4}' "$CONFIG"
)

case "$SAMPLE:$SPECIES" in
    Dsub_L5_7:Dsub|JP-Kan:Dsuz|JP-Sap:Dsuz|JP-Sap2:Dsuz|JP-Tok:Dsuz) ;;
    *) echo "ERROR: unexpected sample or species at config row $ROW" >&2; exit 1 ;;
esac

MAPDIR="${NEW}/population_TE/mapping_single_end_bwasw/${SAMPLE}"
BAM1="${MAPDIR}/${SAMPLE}.R1.bwasw.single_end.bam"
BAM2="${MAPDIR}/${SAMPLE}.R2.bwasw.single_end.bam"

OUTDIR="${NEW}/population_TE/bam_bwasw/${SAMPLE}"
CHECKDIR="${NEW}/checks/se2pe_bwasw/${SAMPLE}"
OUTPUT="${OUTDIR}/${SAMPLE}.Merel.bwasw.se2pe.sorted.bam"
SUMMARY="${CHECKDIR}/${SAMPLE}.Merel.bwasw.se2pe.summary.tsv"
TMP="${OUTDIR}/.${SAMPLE}.Merel.se2pe.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.bam"

STREAM="${OUTDIR}/stream.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
FIFO1="${STREAM}/R1.fastq"
FIFO2="${STREAM}/R2.fastq"
P1=""
P2=""

for FILE in "$CONFIG" "$FASTQ1" "$FASTQ2" "$BAM1" "$BAM2" "$POPTE2" "$SAMTOOLS"; do
    [[ -s "$FILE" ]] || {
        echo "ERROR: missing input: $FILE" >&2
        exit 1
    }
done

[[ ! -e "$OUTPUT" && ! -e "${OUTPUT}.bai" && ! -e "$SUMMARY" ]] || {
    echo "ERROR: output already exists; will not overwrite" >&2
    exit 1
}

mkdir -p "$OUTDIR" "$CHECKDIR" "$STREAM"

cleanup() {
    [[ -z "$P1" ]] || kill "$P1" 2>/dev/null || true
    [[ -z "$P2" ]] || kill "$P2" 2>/dev/null || true
    rm -f -- "$TMP" "${TMP}.bai" "$FIFO1" "$FIFO2"
    rmdir "$STREAM" 2>/dev/null || true
}
trap cleanup EXIT

"$SAMTOOLS" quickcheck -v "$BAM1" "$BAM2"
mkfifo "$FIFO1" "$FIFO2"

# FASTQヘッダーの空白以降を除き、BAMのリード名と合わせる。
gzip -cd "$FASTQ1" |
    awk 'NR % 4 == 1 {sub(/[[:space:]].*$/, "")} {print}' > "$FIFO1" &
P1=$!

gzip -cd "$FASTQ2" |
    awk 'NR % 4 == 1 {sub(/[[:space:]].*$/, "")} {print}' > "$FIFO2" &
P2=$!

echo "sample=$SAMPLE"
echo "species=$SPECIES"
echo "BAM1=$BAM1"
echo "BAM2=$BAM2"
echo "output=$OUTPUT"
date

java -Xmx80g -jar "$POPTE2" se2pe \
    --fastq1 "$FIFO1" \
    --fastq2 "$FIFO2" \
    --bam1 "$BAM1" \
    --bam2 "$BAM2" \
    --sort --index \
    --output "$TMP"

wait "$P1"; P1=""
wait "$P2"; P2=""

"$SAMTOOLS" quickcheck -v "$TMP"

if [[ ! -s "${TMP}.bai" ]]; then
    "$SAMTOOLS" index -@ 3 "$TMP"
fi
[[ -s "${TMP}.bai" ]] || {
    echo "ERROR: BAM index missing" >&2
    exit 1
}
"$SAMTOOLS" idxstats "$TMP" > /dev/null

mv "$TMP" "$OUTPUT"
mv "${TMP}.bai" "${OUTPUT}.bai"

{
    printf 'metric\tvalue\n'
    printf 'sample\t%s\n' "$SAMPLE"
    printf 'input_R1_BAM\t%s\n' "$BAM1"
    printf 'input_R2_BAM\t%s\n' "$BAM2"
    printf 'output_BAM\t%s\n' "$OUTPUT"
    printf 'status\tBAM_CREATED_PAIR_QC_PENDING\n'
} > "$SUMMARY"

cat "$SUMMARY"
echo "FINAL: BAM CREATED; PAIR QC PENDING"
date
