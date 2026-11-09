#!/bin/bash

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821

CHROM_CONFIG=${BASE}/config/chromosomes.tsv
OUTDIR=${BASE}/config/regions
CHECKDIR=${BASE}/checks
LOGDIR=${BASE}/logs

DSUB_FAI=/ru-auth/local/home/asato/scratch/genome/Dsub_v1.1/GCF_014743375.2_RU_Dsub_v1.1_genomic.fna.fai
DSUZ_FAI=/ru-auth/local/home/asato/scratch/genome/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.fna.fai

SUMMARY=${CHECKDIR}/analysis_region_beds.summary.tsv

mkdir -p \
    "${OUTDIR}" \
    "${CHECKDIR}" \
    "${LOGDIR}"

for input in \
    "${CHROM_CONFIG}" \
    "${DSUB_FAI}" \
    "${DSUZ_FAI}"
do
    if [[ ! -s "${input}" ]]; then
        echo "ERROR: missing or empty input: ${input}" >&2
        exit 1
    fi
done

make_region_bed()
{
    local species="$1"
    local analysis_class="$2"
    local fai="$3"
    local output="$4"
    local expected_regions="$5"
    local tmp="${output}.tmp.$$"

    rm -f "${tmp}"

    awk -F'\t' \
      -v species="${species}" \
      -v analysis_class="${analysis_class}" '
        NR==FNR {
            if (FNR > 1 && $1 == species && $5 == analysis_class && $8 == 1) {
                keep[$2] = 1
            }
            next
        }

        $1 in keep {
            print $1 "\t0\t" $2
            found[$1] = 1
        }

        END {
            missing = 0

            for (contig in keep) {
                if (!(contig in found)) {
                    print "ERROR: contig absent from FAI: " contig > "/dev/stderr"
                    missing = 1
                }
            }

            if (missing == 1) {
                exit 1
            }
        }
      ' \
      "${CHROM_CONFIG}" \
      "${fai}" \
      > "${tmp}"

    if [[ ! -s "${tmp}" ]]; then
        echo "ERROR: failed to create ${output}" >&2
        rm -f "${tmp}"
        exit 1
    fi

    local observed_regions

    observed_regions=$(
        awk 'NF {n++} END {print n+0}' "${tmp}"
    )

    if [[ "${observed_regions}" -ne "${expected_regions}" ]]; then
        echo \
          "ERROR: ${species} ${analysis_class}: expected" \
          "${expected_regions} regions, observed" \
          "${observed_regions}" >&2

        rm -f "${tmp}"
        exit 1
    fi

    if awk '
        $2 < 0 || $3 <= $2 {
            exit 1
        }
      ' "${tmp}"
    then
        :
    else
        echo "ERROR: invalid BED coordinates in ${tmp}" >&2
        rm -f "${tmp}"
        exit 1
    fi

    if cut -f1 "${tmp}" | sort | uniq -d | grep -q .
    then
        echo "ERROR: duplicated contig in ${tmp}" >&2
        rm -f "${tmp}"
        exit 1
    fi

    mv "${tmp}" "${output}"
}

echo "===== Build analysis-region BED files ====="

make_region_bed \
    Dsub \
    autosome \
    "${DSUB_FAI}" \
    "${OUTDIR}/Dsub.autosome.bed" \
    6

make_region_bed \
    Dsub \
    X \
    "${DSUB_FAI}" \
    "${OUTDIR}/Dsub.X.bed" \
    1

make_region_bed \
    Dsuz \
    autosome \
    "${DSUZ_FAI}" \
    "${OUTDIR}/Dsuz.autosome.bed" \
    4

make_region_bed \
    Dsuz \
    X \
    "${DSUZ_FAI}" \
    "${OUTDIR}/Dsuz.X.bed" \
    1

{
    echo -e \
      "species\tanalysis_class\tregions\ttotal_length_bp\tBED"

    for species in Dsub Dsuz
    do
        for analysis_class in autosome X
        do
            bed=${OUTDIR}/${species}.${analysis_class}.bed

            regions=$(
                awk 'NF {n++} END {print n+0}' "${bed}"
            )

            total_length=$(
                awk '
                  {
                      total += $3-$2
                  }
                  END {
                      print total+0
                  }
                ' "${bed}"
            )

            echo -e \
              "${species}\t${analysis_class}\t${regions}\t${total_length}\t${bed}"
        done
    done
} > "${SUMMARY}"

echo "===== BED contents ====="

for bed in \
    "${OUTDIR}/Dsub.autosome.bed" \
    "${OUTDIR}/Dsub.X.bed" \
    "${OUTDIR}/Dsuz.autosome.bed" \
    "${OUTDIR}/Dsuz.X.bed"
do
    echo "----- ${bed} -----"
    column -t -s $'\t' "${bed}"
done

echo "===== Summary ====="

column -t -s $'\t' "${SUMMARY}"

echo "FINAL: PASS"
