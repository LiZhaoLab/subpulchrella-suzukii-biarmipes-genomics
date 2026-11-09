#!/usr/bin/env bash

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/synteny_4species/reanalysis_20260824
OUTROOT=${BASE}/results/primary_one_to_one_coverage
CHECKROOT=${BASE}/checks/primary_one_to_one_coverage

mkdir -p "${OUTROOT}" "${CHECKROOT}"

SUMMARY=${CHECKROOT}/primary_one_to_one.L500.I50.coverage_summary.tsv

printf \
  "comparison\tgenome_side\tsequences\tgenome_bp\tunique_aligned_bp\tcoverage_percent\toutput\n" \
  > "${SUMMARY}"

for comparison in \
  Dsuz_vs_Dsub \
  Dsub_vs_Dbia \
  Dbia_vs_Dmel
do
  INPUT=${BASE}/results/primary_one_to_one/${comparison}.promer.one_to_one.L500.I50.coords.tsv
  OUTDIR=${OUTROOT}/${comparison}

  mkdir -p "${OUTDIR}"

  [[ -s "${INPUT}" ]] || {
    echo "ERROR: missing input: ${INPUT}" >&2
    exit 1
  }

  REF_BED=${OUTDIR}/${comparison}.reference.raw.bed
  QUERY_BED=${OUTDIR}/${comparison}.query.raw.bed

  REF_MERGED=${OUTDIR}/${comparison}.reference.merged.bed
  QUERY_MERGED=${OUTDIR}/${comparison}.query.merged.bed

  REF_LENGTHS=${OUTDIR}/${comparison}.reference.lengths.tsv
  QUERY_LENGTHS=${OUTDIR}/${comparison}.query.lengths.tsv

  REF_COVERAGE=${OUTDIR}/${comparison}.reference.coverage.tsv
  QUERY_COVERAGE=${OUTDIR}/${comparison}.query.coverage.tsv

  awk -F $'\t' -v OFS=$'\t' '
    {
      start = ($1 < $2 ? $1 : $2) - 1
      end   = ($1 > $2 ? $1 : $2)
      print $16, start, end
    }
  ' "${INPUT}" > "${REF_BED}"

  awk -F $'\t' -v OFS=$'\t' '
    {
      start = ($3 < $4 ? $3 : $4) - 1
      end   = ($3 > $4 ? $3 : $4)
      print $17, start, end
    }
  ' "${INPUT}" > "${QUERY_BED}"

  awk -F $'\t' -v OFS=$'\t' '
    !seen[$16]++ {
      print $16, $10
    }
  ' "${INPUT}" |
  sort -k1,1 > "${REF_LENGTHS}"

  awk -F $'\t' -v OFS=$'\t' '
    !seen[$17]++ {
      print $17, $11
    }
  ' "${INPUT}" |
  sort -k1,1 > "${QUERY_LENGTHS}"

  sort -k1,1 -k2,2n "${REF_BED}" |
  bedtools merge -i - \
    > "${REF_MERGED}"

  sort -k1,1 -k2,2n "${QUERY_BED}" |
  bedtools merge -i - \
    > "${QUERY_MERGED}"

  awk -F $'\t' -v OFS=$'\t' '
    NR == FNR {
      seq_length[$1] = $2
      next
    }

    {
      covered[$1] += $3 - $2
    }

    END {
      print "sequence", "length_bp", "unique_aligned_bp", "coverage_percent"

      for (sequence in seq_length) {
        value = covered[sequence] + 0

        if (value > seq_length[sequence]) {
          print "ERROR: coverage exceeds sequence length: " sequence \
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
  ' "${REF_LENGTHS}" "${REF_MERGED}" |
  {
    read -r header
    printf "%s\n" "${header}"
    sort -t $'\t' -k3,3nr
  } > "${REF_COVERAGE}"

  awk -F $'\t' -v OFS=$'\t' '
    NR == FNR {
      seq_length[$1] = $2
      next
    }

    {
      covered[$1] += $3 - $2
    }

    END {
      print "sequence", "length_bp", "unique_aligned_bp", "coverage_percent"

      for (sequence in seq_length) {
        value = covered[sequence] + 0

        if (value > seq_length[sequence]) {
          print "ERROR: coverage exceeds sequence length: " sequence \
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
  ' "${QUERY_LENGTHS}" "${QUERY_MERGED}" |
  {
    read -r header
    printf "%s\n" "${header}"
    sort -t $'\t' -k3,3nr
  } > "${QUERY_COVERAGE}"

  for side in reference query
  do
    COVERAGE=${OUTDIR}/${comparison}.${side}.coverage.tsv

    read -r sequence_count genome_bp aligned_bp coverage_percent < <(
      awk -F $'\t' '
        NR > 1 {
          count++
          genome += $2
          aligned += $3
        }

        END {
          printf "%d %.0f %.0f %.6f\n",
            count,
            genome,
            aligned,
            100 * aligned / genome
        }
      ' "${COVERAGE}"
    )

    printf \
      "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "${comparison}" \
      "${side}" \
      "${sequence_count}" \
      "${genome_bp}" \
      "${aligned_bp}" \
      "${coverage_percent}" \
      "${COVERAGE}" \
      >> "${SUMMARY}"
  done
done

cat "${SUMMARY}"

echo "FINAL: PASS"
