#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --array=0-9%4
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=0-08:00:00
#SBATCH -J poly_indelMaskFinal
#SBATCH -o /ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821/logs/poly_indelMaskFinal.%A_%a.out
#SBATCH -e /ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821/logs/poly_indelMaskFinal.%A_%a.err

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821

SAMPLETAB=${BASE}/config/samples.tsv
REGIONDIR=${BASE}/config/regions
BAMDIR=${BASE}/bam_rmdup
OUTDIR=${BASE}/indel_masks
CHECKDIR=${BASE}/checks/indel_masks_final
PYTHON_SCRIPT=${BASE}/scripts/build_fraction_aware_indel_mask.py

SAMTOOLS=/ru-auth/local/home/asato/scratch/conda_envs/augustus_env/bin/samtools

MIN_MAPQ=20
MIN_BASEQ=20
REQUIRE_FLAGS=2
EXCLUDE_FLAGS=3844

MIN_SUPPORT=3
MIN_FRACTION=0.05
FLANK=5

SAMPLES=(
    Dsub_L5_7
    JP-Kan
    JP-Sap
    JP-Sap2
    JP-Tok
)

CLASSES=(
    autosome
    X
)

sample_index=$((SLURM_ARRAY_TASK_ID / 2))
class_index=$((SLURM_ARRAY_TASK_ID % 2))

SAMPLE=${SAMPLES[$sample_index]}
CLASS=${CLASSES[$class_index]}

sample_row=$(
    awk -F'\t' -v sample="${SAMPLE}" '
        NR > 1 && $1 == sample {
            print $2 "\t" $4
            exit
        }
    ' "${SAMPLETAB}"
)

[[ -n "${sample_row}" ]] || {
    echo "ERROR: sample not found: ${SAMPLE}" >&2
    exit 1
}

IFS=$'\t' read -r SPECIES REFERENCE_FASTA <<< "${sample_row}"

BAM=${BAMDIR}/${SAMPLE}.rmdup.bam
BED=${REGIONDIR}/${SPECIES}.${CLASS}.bed

PREFIX=${SAMPLE}.${CLASS}
MASK_BED=${OUTDIR}/${PREFIX}.indel_mask.support3.fraction0p05.flank5.bed
EVENT_TABLE=${CHECKDIR}/${PREFIX}.passing_indel_events.tsv
SUMMARY=${CHECKDIR}/${PREFIX}.indel_mask.summary.tsv

mkdir -p \
    "${OUTDIR}" \
    "${CHECKDIR}"

for required in \
    "${BAM}" \
    "${BAM}.bai" \
    "${BED}" \
    "${REFERENCE_FASTA}" \
    "${PYTHON_SCRIPT}"
do
    [[ -s "${required}" ]] || {
        echo "ERROR: missing input: ${required}" >&2
        exit 1
    }
done

echo "sample=${SAMPLE}"
echo "species=${SPECIES}"
echo "analysis_class=${CLASS}"
echo "bam=${BAM}"
echo "reference=${REFERENCE_FASTA}"
echo "region_bed=${BED}"
echo "minimum_MAPQ=${MIN_MAPQ}"
echo "minimum_base_quality=${MIN_BASEQ}"
echo "required_flags=${REQUIRE_FLAGS}"
echo "excluded_flags=${EXCLUDE_FLAGS}"
echo "minimum_supporting_reads=${MIN_SUPPORT}"
echo "minimum_supporting_fraction=${MIN_FRACTION}"
echo "flank_bp=${FLANK}"

"${SAMTOOLS}" mpileup \
    --fasta-ref "${REFERENCE_FASTA}" \
    --positions "${BED}" \
    --min-MQ "${MIN_MAPQ}" \
    --min-BQ "${MIN_BASEQ}" \
    --incl-flags "${REQUIRE_FLAGS}" \
    --excl-flags "${EXCLUDE_FLAGS}" \
    --max-depth 1000000 \
    --no-output-ends \
    "${BAM}" |
python3 "${PYTHON_SCRIPT}" \
    --sample "${SAMPLE}" \
    --species "${SPECIES}" \
    --analysis-class "${CLASS}" \
    --region-bed "${BED}" \
    --minimum-support "${MIN_SUPPORT}" \
    --minimum-fraction "${MIN_FRACTION}" \
    --flank "${FLANK}" \
    --mask-bed "${MASK_BED}" \
    --event-output "${EVENT_TABLE}" \
    --summary-output "${SUMMARY}"

[[ -s "${MASK_BED}" ]] || {
    echo "ERROR: final mask BED is empty: ${MASK_BED}" >&2
    exit 1
}

[[ -s "${SUMMARY}" ]] || {
    echo "ERROR: summary missing: ${SUMMARY}" >&2
    exit 1
}

cat "${SUMMARY}"

echo "mask_bed=${MASK_BED}"
echo "event_table=${EVENT_TABLE}"
echo "FINAL: PASS"
