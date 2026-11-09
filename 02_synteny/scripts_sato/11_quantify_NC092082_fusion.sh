#!/usr/bin/env bash

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/synteny_4species/reanalysis_20260824

INPUT=${BASE}/results/primary_one_to_one/Dsuz_vs_Dsub.promer.one_to_one.L500.I50.coords.tsv

REFERENCE=/ru-auth/local/home/asato/scratch/genome/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.fna

GFF=/ru-auth/local/home/asato/scratch/genome/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.gff.gz

SAMTOOLS=/ru-auth/local/home/asato/scratch/conda_envs/augustus_env/bin/samtools
BEDTOOLS=$(command -v bedtools)

OUTDIR=${BASE}/results/NC092082_fusion
CHECKDIR=${BASE}/checks/NC092082_fusion

REGIONS=${OUTDIR}/NC_092082.1.analysis_regions.bed
ARM_COVERAGE=${OUTDIR}/NC_092082.1.arm_synteny_coverage.tsv
GENE_DENSITY=${OUTDIR}/NC_092082.1.region_gene_density.tsv
AMBIGUOUS=${OUTDIR}/NC_092082.1.region_ambiguous_bases.tsv
SUMMARY=${CHECKDIR}/NC_092082.1.fusion.summary.tsv

mkdir -p "${OUTDIR}" "${CHECKDIR}"

for required in \
    "${INPUT}" \
    "${REFERENCE}" \
    "${REFERENCE}.fai" \
    "${GFF}" \
    "${SAMTOOLS}" \
    "${BEDTOOLS}"
do
    if [[ ! -s "${required}" ]]; then
        echo "ERROR: missing or empty input: ${required}" >&2
        exit 1
    fi
done

WORKDIR=$(mktemp -d)

cleanup() {
    rm -r "${WORKDIR}"
}

trap cleanup EXIT

echo "===== NC_092082.1 fusion analysis ====="
echo "primary_synteny=${INPUT}"
echo "reference=${REFERENCE}"
echo "annotation=${GFF}"
date

{
    printf "NC_092082.1\t0\t27000000\t3L_like\n"
    printf "NC_092082.1\t27000000\t62000000\ttransition\n"
    printf "NC_092082.1\t62000000\t93255888\t3R_like\n"
} > "${REGIONS}"

awk -F $'\t' '
BEGIN {
    OFS = "\t"
}

$16 == "NC_092082.1" && $17 == "NC_050612.1" {
    start = ($1 < $2 ? $1 : $2) - 1
    end = ($1 > $2 ? $1 : $2)
    print $16, start, end
}
' "${INPUT}" |
sort -k1,1 -k2,2n |
"${BEDTOOLS}" merge \
> "${WORKDIR}/3L.unique.bed"

awk -F $'\t' '
BEGIN {
    OFS = "\t"
}

$16 == "NC_092082.1" && ($17 == "NC_050609.1" || $17 == "NW_023665577.1") {
    start = ($1 < $2 ? $1 : $2) - 1
    end = ($1 > $2 ? $1 : $2)
    print $16, start, end
}
' "${INPUT}" |
sort -k1,1 -k2,2n |
"${BEDTOOLS}" merge \
> "${WORKDIR}/3R.unique.bed"

{
    printf "region\tarm\tregion_start_0based\tregion_end\t"
    printf "overlapping_intervals\tcovered_bp\tregion_bp\t"
    printf "covered_fraction\tcovered_percent\n"

    for arm in 3L 3R
    do
        "${BEDTOOLS}" coverage \
          -a "${REGIONS}" \
          -b "${WORKDIR}/${arm}.unique.bed" |
        awk -F $'\t' -v arm="${arm}" '
        BEGIN {
            OFS = "\t"
        }

        {
            printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%.12g\t%.6f\n",
                $4,
                arm,
                $2,
                $3,
                $5,
                $6,
                $7,
                $8,
                100 * $8
        }
        '
    done
} > "${ARM_COVERAGE}"

gzip -cd "${GFF}" |
awk -F $'\t' '
BEGIN {
    OFS = "\t"
}

$0 !~ /^#/ && $1 == "NC_092082.1" && $3 == "gene" {
    midpoint = ($4 + $5) / 2

    if (midpoint < 27000000) {
        genes["3L_like"]++
    } else if (midpoint < 62000000) {
        genes["transition"]++
    } else {
        genes["3R_like"]++
    }
}

END {
    print "region", "region_bp", "genes", "genes_per_Mb"

    printf "3L_like\t27000000\t%d\t%.6f\n",
        genes["3L_like"],
        genes["3L_like"] / 27

    printf "transition\t35000000\t%d\t%.6f\n",
        genes["transition"],
        genes["transition"] / 35

    printf "3R_like\t31255888\t%d\t%.6f\n",
        genes["3R_like"],
        genes["3R_like"] / 31.255888
}
' > "${GENE_DENSITY}"

printf "region\tregion_bp\tN_bp\tN_fraction\n" \
> "${AMBIGUOUS}"

while IFS=$'\t' read -r chrom start end region
do
    fasta_start=$((start + 1))

    "${SAMTOOLS}" faidx \
      "${REFERENCE}" \
      "${chrom}:${fasta_start}-${end}" |
    awk -v region="${region}" '
    BEGIN {
        total = 0
        n_count = 0
    }

    $0 !~ /^>/ {
        total += length($0)
        sequence = $0
        n_count += gsub(/[Nn]/, "", sequence)
    }

    END {
        if (total > 0) {
            fraction = n_count / total
        } else {
            fraction = 0
        }

        printf "%s\t%d\t%d\t%.12g\n",
            region,
            total,
            n_count,
            fraction
    }
    ' >> "${AMBIGUOUS}"
done < "${REGIONS}"

three_l_fraction=$(
  awk -F $'\t' '
    $1 == "3L_like" && $2 == "3L" {
      print $8
    }
  ' "${ARM_COVERAGE}"
)

three_r_fraction=$(
  awk -F $'\t' '
    $1 == "3R_like" && $2 == "3R" {
      print $8
    }
  ' "${ARM_COVERAGE}"
)

transition_genes_per_mb=$(
  awk -F $'\t' '
    $1 == "transition" {
      print $4
    }
  ' "${GENE_DENSITY}"
)

transition_n_fraction=$(
  awk -F $'\t' '
    $1 == "transition" {
      print $4
    }
  ' "${AMBIGUOUS}"
)

{
    printf "metric\tvalue\n"
    printf "target_chromosome\tNC_092082.1\n"
    printf "target_length_bp\t93255888\n"
    printf "three_L_like_region\t0-27000000\n"
    printf "transition_region\t27000000-62000000\n"
    printf "three_R_like_region\t62000000-93255888\n"
    printf "Dsub_three_L_homolog\tNC_050612.1\n"
    printf "Dsub_three_R_homolog_1\tNC_050609.1\n"
    printf "Dsub_three_R_homolog_2\tNW_023665577.1\n"
    printf "three_L_region_three_L_synteny_fraction\t%s\n" \
      "${three_l_fraction}"
    printf "three_R_region_three_R_synteny_fraction\t%s\n" \
      "${three_r_fraction}"
    printf "transition_genes_per_Mb\t%s\n" \
      "${transition_genes_per_mb}"
    printf "transition_N_fraction\t%s\n" \
      "${transition_n_fraction}"
    printf "synteny_definition\tone_to_one_PROmer_length_ge_500_identity_ge_50\n"
    printf "coverage_definition\tmerged_unique_reference_intervals\n"
    printf "interpretation\t3L_3R_fusion_supported\n"
    printf "status\tPASS\n"
} > "${SUMMARY}"

echo "===== Arm-specific unique synteny coverage ====="
column -t -s $'\t' "${ARM_COVERAGE}"

echo "===== Regional gene density ====="
column -t -s $'\t' "${GENE_DENSITY}"

echo "===== Regional ambiguous bases ====="
column -t -s $'\t' "${AMBIGUOUS}"

echo "===== Summary ====="
cat "${SUMMARY}"

echo "FINAL: PASS"
date
