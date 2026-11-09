#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --array=0-9%2
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=0-08:00:00
#SBATCH -J poly_depthQC
#SBATCH -o /ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821/logs/poly_depthQC.%A_%a.out
#SBATCH -e /ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821/logs/poly_depthQC.%A_%a.err

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821

CONFIG=${BASE}/config/samples.tsv
BAMDIR=${BASE}/bam_rmdup
BEDDIR=${BASE}/config/regions
OUTDIR=${BASE}/checks/depth_qc
LOGDIR=${BASE}/logs

SAMTOOLS=/ru-auth/local/home/asato/scratch/conda_envs/augustus_env/bin/samtools
PYTHON=python3

mkdir -p \
    "${OUTDIR}" \
    "${LOGDIR}"

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

SAMPLE_INDEX=$((SLURM_ARRAY_TASK_ID / 2))
CLASS_INDEX=$((SLURM_ARRAY_TASK_ID % 2))

SAMPLE=${SAMPLES[$SAMPLE_INDEX]}
ANALYSIS_CLASS=${CLASSES[$CLASS_INDEX]}

SPECIES=$(
    awk -F'\t' -v sample="${SAMPLE}" \
      '$1==sample {print $2; exit}' \
      "${CONFIG}"
)

if [[ -z "${SPECIES}" ]]; then
    echo "ERROR: species not found for ${SAMPLE}" >&2
    exit 1
fi

BAM=${BAMDIR}/${SAMPLE}.rmdup.bam
BED=${BEDDIR}/${SPECIES}.${ANALYSIS_CLASS}.bed

SUMMARY=${OUTDIR}/${SAMPLE}.${ANALYSIS_CLASS}.depth_summary.tsv
HISTOGRAM=${OUTDIR}/${SAMPLE}.${ANALYSIS_CLASS}.depth_histogram.tsv

SUMMARY_TMP=${SUMMARY}.tmp.${SLURM_JOB_ID}
HISTOGRAM_TMP=${HISTOGRAM}.tmp.${SLURM_JOB_ID}

for input in \
    "${CONFIG}" \
    "${BAM}" \
    "${BAM}.bai" \
    "${BED}"
do
    if [[ ! -s "${input}" ]]; then
        echo "ERROR: missing or empty input: ${input}" >&2
        exit 1
    fi
done

echo "===== Depth QC ====="
echo "sample=${SAMPLE}"
echo "species=${SPECIES}"
echo "analysis_class=${ANALYSIS_CLASS}"
echo "bam=${BAM}"
echo "bed=${BED}"
echo "minimum_MAPQ=20"
echo "minimum_base_quality=20"
echo "require_flags=2"
echo "exclude_flags=3844"
echo "threads=${SLURM_CPUS_PER_TASK}"

"${SAMTOOLS}" quickcheck -v "${BAM}"

rm -f \
    "${SUMMARY_TMP}" \
    "${HISTOGRAM_TMP}"

"${SAMTOOLS}" depth \
    -@ "${SLURM_CPUS_PER_TASK}" \
    -a \
    -b "${BED}" \
    -q 20 \
    -Q 20 \
    --require-flags 2 \
    --excl-flags 3844 \
    "${BAM}" |
"${PYTHON}" -c '
import sys
from collections import Counter

summary_path = sys.argv[1]
histogram_path = sys.argv[2]
sample = sys.argv[3]
analysis_class = sys.argv[4]

histogram = Counter()

total_sites = 0
covered_sites = 0
depth_sum = 0
maximum_depth = 0

for line in sys.stdin:
    fields = line.rstrip("\n").split("\t")

    if len(fields) < 3:
        continue

    depth = int(fields[2])

    histogram[depth] += 1
    total_sites += 1
    depth_sum += depth

    if depth > 0:
        covered_sites += 1

    if depth > maximum_depth:
        maximum_depth = depth

zero_depth_sites = total_sites - covered_sites

def calculate_quantile(counts, number_of_sites, probability, minimum_depth):
    if number_of_sites == 0:
        return "NA"

    rank = max(
        1,
        int(probability * number_of_sites + 0.999999999)
    )

    cumulative = 0

    for depth in sorted(counts):
        if depth < minimum_depth:
            continue

        cumulative += counts[depth]

        if cumulative >= rank:
            return str(depth)

    return str(max(counts))

if total_sites > 0:
    mean_all = depth_sum / total_sites
    zero_percent = 100.0 * zero_depth_sites / total_sites
else:
    mean_all = 0.0
    zero_percent = 0.0

if covered_sites > 0:
    mean_covered = depth_sum / covered_sites
else:
    mean_covered = 0.0

quantiles = [
    (0.01, "P01"),
    (0.05, "P05"),
    (0.10, "P10"),
    (0.25, "P25"),
    (0.50, "median"),
    (0.75, "P75"),
    (0.90, "P90"),
    (0.95, "P95"),
    (0.99, "P99")
]

with open(summary_path, "w") as output:
    output.write("metric\tvalue\n")
    output.write(f"sample\t{sample}\n")
    output.write(f"analysis_class\t{analysis_class}\n")
    output.write(f"total_reference_sites\t{total_sites}\n")
    output.write(f"covered_sites\t{covered_sites}\n")
    output.write(f"zero_depth_sites\t{zero_depth_sites}\n")
    output.write(f"zero_depth_percent\t{zero_percent:.6f}\n")
    output.write(f"mean_depth_all_sites\t{mean_all:.6f}\n")
    output.write(
        f"mean_depth_covered_sites\t{mean_covered:.6f}\n"
    )
    output.write(f"maximum_depth\t{maximum_depth}\n")

    for probability, label in quantiles:
        all_sites_value = calculate_quantile(
            histogram,
            total_sites,
            probability,
            0
        )

        covered_sites_value = calculate_quantile(
            histogram,
            covered_sites,
            probability,
            1
        )

        output.write(
            f"{label}_depth_all_sites\t"
            f"{all_sites_value}\n"
        )

        output.write(
            f"{label}_depth_covered_sites\t"
            f"{covered_sites_value}\n"
        )

with open(histogram_path, "w") as output:
    output.write(
        "depth\tsite_count\t"
        "fraction_all_sites\tfraction_covered_sites\n"
    )

    for depth in sorted(histogram):
        count = histogram[depth]

        if total_sites > 0:
            fraction_all = count / total_sites
        else:
            fraction_all = 0.0

        if depth > 0 and covered_sites > 0:
            fraction_covered = count / covered_sites
        else:
            fraction_covered = 0.0

        output.write(
            f"{depth}\t{count}\t"
            f"{fraction_all:.12f}\t"
            f"{fraction_covered:.12f}\n"
        )
' \
"${SUMMARY_TMP}" \
"${HISTOGRAM_TMP}" \
"${SAMPLE}" \
"${ANALYSIS_CLASS}"

if [[ ! -s "${SUMMARY_TMP}" ]]; then
    echo "ERROR: empty depth summary" >&2
    exit 1
fi

if [[ ! -s "${HISTOGRAM_TMP}" ]]; then
    echo "ERROR: empty depth histogram" >&2
    exit 1
fi

EXPECTED_SITES=$(
    awk '
      {
          total += $3-$2
      }
      END {
          print total+0
      }
    ' "${BED}"
)

OBSERVED_SITES=$(
    awk -F'\t' \
      '$1=="total_reference_sites" {print $2}' \
      "${SUMMARY_TMP}"
)

if [[ "${OBSERVED_SITES}" -ne "${EXPECTED_SITES}" ]]; then
    echo \
      "ERROR: expected ${EXPECTED_SITES} sites," \
      "observed ${OBSERVED_SITES}" >&2

    exit 1
fi

mv "${SUMMARY_TMP}" "${SUMMARY}"
mv "${HISTOGRAM_TMP}" "${HISTOGRAM}"

echo "===== Depth summary ====="
cat "${SUMMARY}"

echo "FINAL: PASS"
