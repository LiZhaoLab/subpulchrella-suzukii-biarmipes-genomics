#!/bin/bash

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/population_polymorphism/Ne_reanalysis_20260825
POLYBASE=/ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821
BEDTOOLS=/ru-auth/local/home/asato/scratch/conda_envs/earlgrey_env/bin/bedtools

OUTDIR=${BASE}/sites/analysis_sets
CHECKDIR=${BASE}/checks/site_sets

mkdir -p "${OUTDIR}" "${CHECKDIR}"

SUMMARY=${CHECKDIR}/neutral_site_sets_by_analysis_class.tsv
TMP_SUMMARY=${SUMMARY}.tmp.$$

printf \
  "species\tanalysis_class\tsite_type\tsites\tBED_intervals\toutput_BED\n" \
  > "${TMP_SUMMARY}"

for species in Dsub Dsuz
do
    for analysis_class in autosome X
    do
        REGION_BED=${POLYBASE}/config/regions/${species}.${analysis_class}.bed

        [[ -s "${REGION_BED}" ]] || {
            echo "ERROR: missing region BED: ${REGION_BED}" >&2
            exit 1
        }

        for site_type in short_intron fourfold
        do
            case "${site_type}" in
                short_intron)
                    INPUT_BED=${BASE}/sites/${species}/${species}.short_intron_8to30.strict.bed
                    ;;
                fourfold)
                    INPUT_BED=${BASE}/sites/${species}/${species}.fourfold_5families.strict.bed
                    ;;
            esac

            [[ -s "${INPUT_BED}" ]] || {
                echo "ERROR: missing neutral-site BED: ${INPUT_BED}" >&2
                exit 1
            }

            OUTPUT_BED=${OUTDIR}/${species}.${analysis_class}.${site_type}.bed
            TMP_BED=${OUTPUT_BED}.tmp.$$

"${BEDTOOLS}" intersect \
  -a "${INPUT_BED}" \
  -b "${REGION_BED}" \
  -wa \
  > "${TMP_BED}"

            [[ -s "${TMP_BED}" ]] || {
                echo "ERROR: empty output: ${OUTPUT_BED}" >&2
                exit 1
            }

            INVALID=$(
              awk '
                NF != 3 ||
                $2 !~ /^[0-9]+$/ ||
                $3 !~ /^[0-9]+$/ ||
                $2 < 0 ||
                $3 <= $2 {
                    invalid++
                }
                END {
                    print invalid + 0
                }
              ' "${TMP_BED}"
            )

            [[ "${INVALID}" -eq 0 ]] || {
                echo "ERROR: invalid BED rows in ${OUTPUT_BED}: ${INVALID}" >&2
                exit 1
            }

            SITES=$(
              awk '{
                  total += $3 - $2
              }
              END {
                  print total + 0
              }' "${TMP_BED}"
            )

            INTERVALS=$(awk 'END {print NR + 0}' "${TMP_BED}")

            mv "${TMP_BED}" "${OUTPUT_BED}"

            printf \
              "%s\t%s\t%s\t%s\t%s\t%s\n" \
              "${species}" \
              "${analysis_class}" \
              "${site_type}" \
              "${SITES}" \
              "${INTERVALS}" \
              "${OUTPUT_BED}" \
              >> "${TMP_SUMMARY}"
        done
    done
done

mv "${TMP_SUMMARY}" "${SUMMARY}"

echo "===== Final neutral-site sets ====="
column -t -s $'\t' "${SUMMARY}"

echo "FINAL: PASS"
