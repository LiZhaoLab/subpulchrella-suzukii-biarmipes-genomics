#!/usr/bin/env bash

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/synteny_4species/reanalysis_20260824

OUTDIR=${BASE}/results/primary_one_to_one
CHECKDIR=${BASE}/checks/primary_one_to_one

mkdir -p "${OUTDIR}" "${CHECKDIR}"

SUMMARY=${CHECKDIR}/primary_one_to_one.L500.I50.summary.tsv

printf \
  "comparison\tinput_blocks\tretained_blocks\tretained_percent\tsame_orientation_blocks\topposite_orientation_blocks\toutput\n" \
  > "${SUMMARY}"

for comparison in \
  Dsuz_vs_Dsub \
  Dsub_vs_Dbia \
  Dbia_vs_Dmel
do
  INPUT=${BASE}/alignments/${comparison}/${comparison}.promer.one_to_one.L100.coords.tsv
  OUTPUT=${OUTDIR}/${comparison}.promer.one_to_one.L500.I50.coords.tsv

  [[ -s "${INPUT}" ]] || {
    echo "ERROR: missing input: ${INPUT}" >&2
    exit 1
  }

  awk -F $'\t' -v OFS=$'\t' '
    NF != 17 {
      invalid++
      next
    }

    ($5 + 0) >= 500 &&
    ($6 + 0) >= 500 &&
    ($7 + 0) >= 50 {
      print
    }

    END {
      if (invalid > 0) {
        print "ERROR: invalid rows=" invalid > "/dev/stderr"
        exit 1
      }
    }
  ' "${INPUT}" > "${OUTPUT}"

  input_blocks=$(wc -l < "${INPUT}")
  retained_blocks=$(wc -l < "${OUTPUT}")

  retained_percent=$(
    awk -v retained="${retained_blocks}" -v total="${input_blocks}" \
      'BEGIN {printf "%.6f", 100 * retained / total}'
  )

  same_orientation=$(
    awk -F $'\t' \
      '($14 + 0) * ($15 + 0) > 0 {count++}
       END {print count + 0}' \
      "${OUTPUT}"
  )

  opposite_orientation=$(
    awk -F $'\t' \
      '($14 + 0) * ($15 + 0) < 0 {count++}
       END {print count + 0}' \
      "${OUTPUT}"
  )

  printf \
    "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${comparison}" \
    "${input_blocks}" \
    "${retained_blocks}" \
    "${retained_percent}" \
    "${same_orientation}" \
    "${opposite_orientation}" \
    "${OUTPUT}" \
    >> "${SUMMARY}"
done

cat "${SUMMARY}"

echo "FINAL: PASS"
