#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --job-name=Ne_P98
#SBATCH --array=0-19%2
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --output=/ru-auth/local/home/asato/scratch/population_polymorphism/Ne_reanalysis_20261023_repeatmasker_p98/logs/Ne_P98.%A_%a.out
#SBATCH --error=/ru-auth/local/home/asato/scratch/population_polymorphism/Ne_reanalysis_20261023_repeatmasker_p98/logs/Ne_P98.%A_%a.err

set -euo pipefail

NEW=/ru-auth/local/home/asato/scratch/population_polymorphism/Ne_reanalysis_20261023_repeatmasker_p98
POLY=/ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821
P98=/ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20261012_repeatmasker_p98

SAMPLETAB="${POLY}/config/samples.tsv"
DEPTHTAB="${P98}/config/depth_thresholds.tsv"

SAMPLES=(Dsub_L5_7 JP-Kan JP-Sap JP-Sap2 JP-Tok)
CLASSES=(autosome X)
SITE_TYPES=(short_intron fourfold)

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is not set}"
if (( TASK_ID < 0 || TASK_ID > 19 )); then
    echo "ERROR: array task ID must be 0–19" >&2
    exit 1
fi

sample_index=$((TASK_ID / 4))
remainder=$((TASK_ID % 4))
class_index=$((remainder / 2))
site_index=$((remainder % 2))

SAMPLE="${SAMPLES[$sample_index]}"
CLASS="${CLASSES[$class_index]}"
SITE_TYPE="${SITE_TYPES[$site_index]}"

sample_row=$(
    awk -F '\t' -v sample="$SAMPLE" '
        NR > 1 && $1 == sample {
            print $2 "\t" $4 "\t" $5 "\t" $6
            exit
        }
    ' "$SAMPLETAB"
)
[[ -n "$sample_row" ]] || {
    echo "ERROR: sample not found: $SAMPLE" >&2
    exit 1
}

IFS=$'\t' read -r SPECIES REFERENCE_FASTA AUTOSOMAL_POOL_SIZE X_POOL_SIZE \
    <<< "$sample_row"

if [[ "$CLASS" == autosome ]]; then
    POOL_SIZE="$AUTOSOMAL_POOL_SIZE"
else
    POOL_SIZE="$X_POOL_SIZE"
fi

depth_row=$(
    awk -F '\t' -v sample="$SAMPLE" -v class="$CLASS" '
        NR > 1 && $1 == sample && $2 == class {
            print $3 "\t" $4
            exit
        }
    ' "$DEPTHTAB"
)
[[ -n "$depth_row" ]] || {
    echo "ERROR: P98 depth thresholds not found: $SAMPLE $CLASS" >&2
    exit 1
}

IFS=$'\t' read -r MIN_DEPTH MAX_DEPTH <<< "$depth_row"

BAM="${POLY}/bam_filtered/${SAMPLE}.${CLASS}.filtered.bam"
MASK_BED="${P98}/masks/diversity/${SAMPLE}.${CLASS}.indel_repeatmasker.bed"

if [[ "$SITE_TYPE" == short_intron ]]; then
    NEUTRAL_BED="${NEW}/sites/analysis_sets/${SPECIES}.${CLASS}.short_intron.len40to65.bed"
    OUTDIR="${NEW}/grenedalf_neutral_len40to65/${CLASS}/short_intron"
    SITE_LABEL=short_intron40to65
else
    NEUTRAL_BED="${NEW}/sites/analysis_sets/${SPECIES}.${CLASS}.fourfold.bed"
    OUTDIR="${NEW}/grenedalf_neutral/${CLASS}/fourfold"
    SITE_LABEL=fourfold
fi

LOGDIR="${OUTDIR}/logs"
SUFFIX="${SAMPLE}_${CLASS}_${SITE_LABEL}_20kb_DP${MIN_DEPTH}-${MAX_DEPTH}_MAPQ20_BQ20_indelS3F0p05flank5_repeatmasker_P98"
OUTPUT="${OUTDIR}/diversity-${SUFFIX}.csv"

for FILE in \
    "$SAMPLETAB" \
    "$DEPTHTAB" \
    "$BAM" \
    "${BAM}.bai" \
    "$REFERENCE_FASTA" \
    "$MASK_BED" \
    "$NEUTRAL_BED"
do
    [[ -s "$FILE" ]] || {
        echo "ERROR: missing or empty input: $FILE" >&2
        exit 1
    }
done

[[ ! -e "$OUTPUT" ]] || {
    echo "ERROR: output already exists; will not overwrite: $OUTPUT" >&2
    exit 1
}

mkdir -p "$OUTDIR" "$LOGDIR"

source /ru-auth/local/home/asato/miniconda3/etc/profile.d/conda.sh
conda activate /ru-auth/local/home/asato/scratch/conda_envs/grenedalf_env
GRENEDALF=$(command -v grenedalf)

echo "sample=$SAMPLE"
echo "species=$SPECIES"
echo "class=$CLASS"
echo "neutral_site_type=$SITE_TYPE"
echo "neutral_site_definition=$SITE_LABEL"
echo "bam=$BAM"
echo "reference=$REFERENCE_FASTA"
echo "neutral_site_BED=$NEUTRAL_BED"
echo "indel_repeatmasker_mask=$MASK_BED"
echo "P98_depth=${MIN_DEPTH}-${MAX_DEPTH}"
echo "pool_size=$POOL_SIZE"
echo "output=$OUTPUT"
date

"$GRENEDALF" diversity \
    --sam-path "$BAM" \
    --sam-min-map-qual 20 \
    --sam-min-base-qual 20 \
    --sam-flags-include-all 2 \
    --sam-flags-exclude-any 3844 \
    --reference-genome-fasta "$REFERENCE_FASTA" \
    --filter-region-bed "$NEUTRAL_BED" \
    --filter-mask-total-bed "$MASK_BED" \
    --filter-sample-min-count 2 \
    --filter-sample-min-read-depth "$MIN_DEPTH" \
    --filter-sample-max-read-depth "$MAX_DEPTH" \
    --pool-sizes "$POOL_SIZE" \
    --window-type interval \
    --window-interval-width 20000 \
    --window-interval-stride 20000 \
    --window-average-policy valid-loci \
    --tajima-d-denominator-policy empirical-min-read-depth \
    --na-entry nan \
    --out-dir "$OUTDIR" \
    --file-suffix "-${SUFFIX}" \
    --threads "$SLURM_CPUS_PER_TASK" \
    > "${LOGDIR}/${SUFFIX}.grenedalf.log" 2>&1

[[ -s "$OUTPUT" ]] || {
    echo "ERROR: missing or empty grenedalf output: $OUTPUT" >&2
    exit 1
}

ROWS=$(awk 'END {print NR-1}' "$OUTPUT")
[[ "$ROWS" -gt 0 ]] || {
    echo "ERROR: no data rows: $OUTPUT" >&2
    exit 1
}

echo "data_rows=$ROWS"
echo "FINAL: PASS"
date
