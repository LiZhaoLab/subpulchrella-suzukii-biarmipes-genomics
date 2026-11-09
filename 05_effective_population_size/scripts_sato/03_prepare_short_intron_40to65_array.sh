#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --job-name=Ne_SI40to65_BED
#SBATCH --array=0-1%2
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/ru-auth/local/home/asato/scratch/population_polymorphism/Ne_reanalysis_20261023_repeatmasker_p98/logs/Ne_SI40to65_BED.%A_%a.out
#SBATCH --error=/ru-auth/local/home/asato/scratch/population_polymorphism/Ne_reanalysis_20261023_repeatmasker_p98/logs/Ne_SI40to65_BED.%A_%a.err

set -euo pipefail

NEW=/ru-auth/local/home/asato/scratch/population_polymorphism/Ne_reanalysis_20261023_repeatmasker_p98
POLY=/ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821

PYTHON=/ru-auth/local/home/asato/scratch/conda_envs/earlgrey_env/bin/python3
BEDTOOLS=/ru-auth/local/home/asato/scratch/conda_envs/earlgrey_env/bin/bedtools
EXTRACT="${NEW}/scripts/02_extract_short_intron_sites_40to65.py"

SPECIES=(Dsub Dsuz)
SP="${SPECIES[$SLURM_ARRAY_TASK_ID]}"

GFF=$(
    awk -F '\t' -v sp="$SP" '
        NR > 1 && $1 == sp {
            print $3
            exit
        }
    ' "${NEW}/config/species.tsv"
)

PREFIX="${NEW}/sites/${SP}/${SP}.len40to65"
INTRONS="${PREFIX}.short_introns.tsv"
RAW_BED="${PREFIX}.short_intron_8to30.strict.bed"
SUMMARY="${PREFIX}.short_intron.summary.tsv"

OUT_AUTOSOME="${NEW}/sites/analysis_sets/${SP}.autosome.short_intron.len40to65.bed"
OUT_X="${NEW}/sites/analysis_sets/${SP}.X.short_intron.len40to65.bed"

for INPUT in "$PYTHON" "$BEDTOOLS" "$EXTRACT" "$GFF" \
    "${POLY}/config/regions/${SP}.autosome.bed" \
    "${POLY}/config/regions/${SP}.X.bed"
do
    [[ -s "$INPUT" ]] || {
        echo "ERROR: missing input: $INPUT" >&2
        exit 1
    }
done

for OUTPUT in "$INTRONS" "$RAW_BED" "$SUMMARY" "$OUT_AUTOSOME" "$OUT_X"
do
    [[ ! -e "$OUTPUT" ]] || {
        echo "ERROR: output already exists; will not overwrite: $OUTPUT" >&2
        exit 1
    }
done

mkdir -p "${NEW}/sites/${SP}" "${NEW}/sites/analysis_sets"

echo "species=$SP"
echo "gff=$GFF"
echo "definition=positions_8_to_30_of_introns_40_to_65bp"
date

"$PYTHON" "$EXTRACT" \
    --gff "$GFF" \
    --species "$SP" \
    --output-prefix "$PREFIX"

[[ -s "$INTRONS" && -s "$RAW_BED" && -s "$SUMMARY" ]] || {
    echo "ERROR: missing extraction output" >&2
    exit 1
}

awk -F '\t' '
    NR == 1 {next}
    {
        n++
        if ($5 < 40 || $5 > 65) bad++
        if ($4 == "+" && $3 - $8 < 10) bad++
        if ($4 == "-" && $7 - $2 < 10) bad++
    }
    END {
        print "checked_introns=" n+0
        print "invalid_introns=" bad+0
        if (n == 0 || bad > 0) exit 1
    }
' "$INTRONS"

for CLASS in autosome X
do
    REGION="${POLY}/config/regions/${SP}.${CLASS}.bed"
    OUT="${NEW}/sites/analysis_sets/${SP}.${CLASS}.short_intron.len40to65.bed"
    TMP="${OUT}.tmp.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

    "$BEDTOOLS" intersect \
        -a "$RAW_BED" \
        -b "$REGION" \
        -wa > "$TMP"

    [[ -s "$TMP" ]] || {
        echo "ERROR: empty analysis BED: $CLASS" >&2
        exit 1
    }

    mv "$TMP" "$OUT"

    SITES=$(awk '{n += $3 - $2} END {print n+0}' "$OUT")
    echo "${SP}_${CLASS}_sites=$SITES"
done

cat "$SUMMARY"
echo "FINAL: PASS"
date
