#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-00:30:00
#SBATCH -J filter_Dsub_MK
#SBATCH -o /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/filter_Dsub_MK.%j.out
#SBATCH -e /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/filter_Dsub_MK.%j.err

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817
INDIR=${BASE}/cmktest_results/Dsub_pairwise
CHECKDIR=${BASE}/checks
SUMMARY=${CHECKDIR}/Dsub_cmktest_filter.summary.tsv

mkdir -p "${CHECKDIR}" "${BASE}/logs"

python3 - <<'PY'
from pathlib import Path
import csv
import math

base = Path(
    "/ru-auth/local/home/asato/scratch/MKtest/"
    "reanalysis_longest_fix_20260817"
)

indir = base / "cmktest_results" / "Dsub_pairwise"
checkdir = base / "checks"
summary_file = checkdir / "Dsub_cmktest_filter.summary.tsv"

comparisons = [
    "Dsub_vs_Dsuz",
    "Dsub_vs_Dbia",
    "Dsub_vs_Dmel",
]

required_columns = [
    "Name",
    "Ncodon",
    "Npopulation",
    "Pn",
    "Dn",
    "Ps",
    "Ds",
    "Alpha",
    "Pvalue",
]


def finite_number(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None

    if not math.isfinite(number):
        return None

    return number


summary_rows = []

for comparison in comparisons:
    input_file = indir / f"{comparison}_mktest.csv"
    output_file = indir / f"{comparison}_mktest.filtered.csv"
    excluded_file = checkdir / f"{comparison}.cmktest_filter_excluded.tsv"

    if not input_file.is_file():
        raise FileNotFoundError(input_file)

    total = 0
    retained = 0

    excluded_nonfinite_count = 0
    excluded_empty_nonsynonymous = 0
    excluded_empty_synonymous = 0
    excluded_no_polymorphism = 0
    excluded_no_divergence = 0

    with (
        input_file.open(newline="") as source,
        output_file.open("w", newline="") as destination,
        excluded_file.open("w", newline="") as excluded,
    ):
        reader = csv.DictReader(source)

        if reader.fieldnames != required_columns:
            raise RuntimeError(
                f"{input_file}: unexpected columns: {reader.fieldnames}"
            )

        writer = csv.DictWriter(
            destination,
            fieldnames=required_columns,
            lineterminator="\n",
        )
        writer.writeheader()

        excluded_writer = csv.writer(
            excluded,
            delimiter="\t",
            lineterminator="\n",
        )
        excluded_writer.writerow(
            [
                "Name",
                "reason",
                "Pn",
                "Dn",
                "Ps",
                "Ds",
                "Alpha",
                "Pvalue",
            ]
        )

        for row in reader:
            total += 1

            pn = finite_number(row["Pn"])
            dn = finite_number(row["Dn"])
            ps = finite_number(row["Ps"])
            ds = finite_number(row["Ds"])
            alpha = finite_number(row["Alpha"])
            pvalue = finite_number(row["Pvalue"])

            reasons = []

            if None in (pn, dn, ps, ds, alpha, pvalue):
                reasons.append("nonfinite_or_nonnumeric")
                excluded_nonfinite_count += 1
            else:
                if pn + dn <= 0:
                    reasons.append("no_nonsynonymous_change")
                    excluded_empty_nonsynonymous += 1

                if ps + ds <= 0:
                    reasons.append("no_synonymous_change")
                    excluded_empty_synonymous += 1

                if pn + ps <= 0:
                    reasons.append("no_polymorphism")
                    excluded_no_polymorphism += 1

                if dn + ds <= 0:
                    reasons.append("no_fixed_difference")
                    excluded_no_divergence += 1

            if reasons:
                excluded_writer.writerow(
                    [
                        row["Name"],
                        ";".join(reasons),
                        row["Pn"],
                        row["Dn"],
                        row["Ps"],
                        row["Ds"],
                        row["Alpha"],
                        row["Pvalue"],
                    ]
                )
                continue

            writer.writerow(row)
            retained += 1

    excluded_count = total - retained

    summary_rows.append(
        [
            comparison,
            total,
            retained,
            excluded_count,
            excluded_nonfinite_count,
            excluded_empty_nonsynonymous,
            excluded_empty_synonymous,
            excluded_no_polymorphism,
            excluded_no_divergence,
            str(input_file),
            str(output_file),
            str(excluded_file),
        ]
    )

    print(f"===== {comparison} =====")
    print(f"input_rows\t{total}")
    print(f"retained_rows\t{retained}")
    print(f"excluded_rows\t{excluded_count}")
    print(f"output\t{output_file}")
    print(f"excluded_table\t{excluded_file}")

with summary_file.open("w", newline="") as out:
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")

    writer.writerow(
        [
            "comparison",
            "input_rows",
            "retained_rows",
            "excluded_rows",
            "nonfinite_or_nonnumeric",
            "no_nonsynonymous_change",
            "no_synonymous_change",
            "no_polymorphism",
            "no_fixed_difference",
            "input_CSV",
            "filtered_CSV",
            "excluded_table",
        ]
    )

    writer.writerows(summary_rows)

print(f"summary\t{summary_file}")
print("FINAL\tPASS")
PY
