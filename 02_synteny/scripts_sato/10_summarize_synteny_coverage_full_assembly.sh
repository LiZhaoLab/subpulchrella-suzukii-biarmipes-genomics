#!/usr/bin/env bash

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/synteny_4species/reanalysis_20260824
SAMTOOLS=/ru-auth/local/home/asato/scratch/conda_envs/gatk4/bin/samtools

comparisons=(
  Dsuz_vs_Dsub
  Dsub_vs_Dbia
  Dbia_vs_Dmel
)

references=(
  "${BASE}/repeatmasker/Dsuz/Dsuz.GCF_043229965.1.fna.masked"
  "${BASE}/repeatmasker/Dsub/Dsub.GCF_014743375.2.fna.masked"
  "${BASE}/repeatmasker/Dbia/Dbia.GCF_025231255.1.fna.masked"
)

queries=(
  "${BASE}/repeatmasker/Dsub/Dsub.GCF_014743375.2.fna.masked"
  "${BASE}/repeatmasker/Dbia/Dbia.GCF_025231255.1.fna.masked"
  "${BASE}/repeatmasker/Dmel/Dmel.GCF_000001215.4.fna.masked"
)

OUTROOT=${BASE}/results/primary_one_to_one_coverage
CHECKROOT=${BASE}/checks/primary_one_to_one_coverage
SUMMARY=${CHECKROOT}/primary_one_to_one.L500.I50.full_assembly_coverage.tsv

mkdir -p "${OUTROOT}" "${CHECKROOT}"

printf \
  "comparison\tgenome_side\tassembly_records\tassembly_bp\tunique_aligned_bp\tcoverage_percent\toutput\n" \
  > "${SUMMARY}"

summarize_side() {
  local comparison=$1
  local side=$2
  local fasta=$3

  local outdir=${OUTROOT}/${comparison}
  local merged=${outdir}/${comparison}.${side}.merged.bed
  local output=${outdir}/${comparison}.${side}.full_assembly.coverage.tsv

  [[ -s "${fasta}" ]] || {
    echo "ERROR: missing FASTA: ${fasta}" >&2
    exit 1
  }

  [[ -s "${merged}" ]] || {
    echo "ERROR: missing merged BED: ${merged}" >&2
    exit 1
  }

  if [[ ! -s "${fasta}.fai" ]]; then
    "${SAMTOOLS}" faidx "${fasta}"
  fi

  awk -F $'\t' -v OFS=$'\t' '
    NR == FNR {
      seq_length[$1] = $2
      next
    }

    {
      if (!($1 in seq_length)) {
        print "ERROR: aligned sequence absent from FASTA index: " $1 \
          > "/dev/stderr"
        exit 1
      }

      covered[$1] += $3 - $2
    }

    END {
      print "sequence", "length_bp", "unique_aligned_bp", "coverage_percent"

      for (sequence in seq_length) {
        value = covered[sequence] + 0

        if (value > seq_length[sequence]) {
          print "ERROR: coverage exceeds length: " sequence \
            > "/dev/stderr"
          exit 1
        }

        printf "%s\t%d\t%d\t%.6f\n",
          sequence,
          seq_length[sequence],
          value,
          100 * value / seq_length[sequence]
      }
    }
  ' "${fasta}.fai" "${merged}" |
  {
    read -r header
    printf "%s\n" "${header}"
    sort -t $'\t' -k3,3nr
  } > "${output}"

  read -r records assembly_bp aligned_bp coverage_percent < <(
    awk -F $'\t' '
      NR > 1 {
        records++
        assembly += $2
        aligned += $3
      }

      END {
        printf "%d %.0f %.0f %.6f\n",
          records,
          assembly,
          aligned,
          100 * aligned / assembly
      }
    ' "${output}"
  )

  printf \
    "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${comparison}" \
    "${side}" \
    "${records}" \
    "${assembly_bp}" \
    "${aligned_bp}" \
    "${coverage_percent}" \
    "${output}" \
    >> "${SUMMARY}"
}

for index in 0 1 2
do
  summarize_side \
    "${comparisons[$index]}" \
    reference \
    "${references[$index]}"

  summarize_side \
    "${comparisons[$index]}" \
    query \
    "${queries[$index]}"
done

cat "${SUMMARY}"

echo "FINAL: PASS"
