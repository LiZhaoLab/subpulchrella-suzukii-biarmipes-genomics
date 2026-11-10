#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-01:00:00
#SBATCH -J filter_Dsuz_MK
#SBATCH -o /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/filter_Dsuz_MK.%j.out
#SBATCH -e /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/filter_Dsuz_MK.%j.err

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817
INDIR=${BASE}/cmktest_results/Dsuz_pairwise
CHECKDIR=${BASE}/checks

mkdir -p "${CHECKDIR}" "${BASE}/logs"

SUMMARY=${CHECKDIR}/Dsuz_cmktest_filter.summary.tsv

python3 - <<'PY'
from pathlib import Path
import csv
import math

base = Path("/ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817")
indir = base / "cmktest_results" / "Dsuz_pairwise"
checkdir = base / "checks"

summary_file = checkdir / "Dsuz_cmktest_filter.summary.tsv"

comparisons = [
    ("JP-Kan",  "Dsuz_vs_Dsub"),
    ("JP-Kan",  "Dsuz_vs_Dbia"),
    ("JP-Kan",  "Dsuz_vs_Dmel"),
    ("JP-Sap",  "Dsuz_vs_Dsub"),
    ("JP-Sap",  "Dsuz_vs_Dbia"),
    ("JP-Sap",  "Dsuz_vs_Dmel"),
    ("JP-Sap2", "Dsuz_vs_Dsub"),
    ("JP-Sap2", "Dsuz_vs_Dbia"),
    ("JP-Sap2", "Dsuz_vs_Dmel"),
    ("JP-Tok",  "Dsuz_vs_Dsub"),
    ("JP-Tok",  "Dsuz_vs_Dbia"),
    ("JP-Tok",  "Dsuz_vs_Dmel"),
]

required_cols = [
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

def as_float(row, key):
    value = row.get(key, "")
    try:
        return float(value)
    except ValueError:
        return float("nan")

def finite_number(row, key):
    value = as_float(row, key)
    return math.isfinite(value)

summary_rows = []

for pop, comp in comparisons:
    input_csv = indir / pop / f"{comp}_mktest.csv"
    output_csv = indir / pop / f"{comp}_mktest.filtered.csv"
    excluded_tsv = checkdir / f"{pop}.{comp}.cmktest_filter_excluded.tsv"

    if not input_csv.exists():
        raise FileNotFoundError(f"Missing input CSV: {input_csv}")

    retained = []
    excluded = []

    counters = {
        "input_rows": 0,
        "retained_rows": 0,
        "excluded_rows": 0,
        "nonfinite_or_nonnumeric": 0,
        "no_nonsynonymous_change": 0,
        "no_synonymous_change": 0,
        "no_polymorphism": 0,
        "no_fixed_difference": 0,
    }

    with input_csv.open(newline="") as fh:
        reader = csv.DictReader(fh)

        if reader.fieldnames != required_cols:
            raise RuntimeError(
                f"Unexpected columns in {input_csv}: {reader.fieldnames}"
            )

        for row in reader:
            counters["input_rows"] += 1

            pn = as_float(row, "Pn")
            dn = as_float(row, "Dn")
            ps = as_float(row, "Ps")
            ds = as_float(row, "Ds")

            reasons = []

            numeric_ok = (
                math.isfinite(pn)
                and math.isfinite(dn)
                and math.isfinite(ps)
                and math.isfinite(ds)
                and finite_number(row, "Alpha")
                and finite_number(row, "Pvalue")
            )

            if not numeric_ok:
                reasons.append("nonfinite_or_nonnumeric")
                counters["nonfinite_or_nonnumeric"] += 1
            else:
                if pn + dn <= 0:
                    reasons.append("no_nonsynonymous_change")
                    counters["no_nonsynonymous_change"] += 1

                if ps + ds <= 0:
                    reasons.append("no_synonymous_change")
                    counters["no_synonymous_change"] += 1

                if pn + ps <= 0:
                    reasons.append("no_polymorphism")
                    counters["no_polymorphism"] += 1

                if dn + ds <= 0:
                    reasons.append("no_fixed_difference")
                    counters["no_fixed_difference"] += 1

            if reasons:
                excluded.append((row, ";".join(reasons)))
            else:
                retained.append(row)

    counters["retained_rows"] = len(retained)
    counters["excluded_rows"] = len(excluded)

    with output_csv.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=required_cols)
        writer.writeheader()
        writer.writerows(retained)

    with excluded_tsv.open("w", newline="") as out:
        out.write(
            "\t".join(
                [
                    "population",
                    "comparison",
                    *required_cols,
                    "exclusion_reason",
                ]
            )
            + "\n"
        )

        for row, reason in excluded:
            out.write(
                "\t".join(
                    [
                        pop,
                        comp,
                        *[row[col] for col in required_cols],
                        reason,
                    ]
                )
                + "\n"
            )

    print(f"===== {pop} {comp} =====")
    print(f"input_rows\t{counters['input_rows']}")
    print(f"retained_rows\t{counters['retained_rows']}")
    print(f"excluded_rows\t{counters['excluded_rows']}")
    print(f"output\t{output_csv}")
    print(f"excluded_table\t{excluded_tsv}")

    summary_rows.append(
        [
            pop,
            comp,
            counters["input_rows"],
            counters["retained_rows"],
            counters["excluded_rows"],
            counters["nonfinite_or_nonnumeric"],
            counters["no_nonsynonymous_change"],
            counters["no_synonymous_change"],
            counters["no_polymorphism"],
            counters["no_fixed_difference"],
            str(input_csv),
            str(output_csv),
            str(excluded_tsv),
        ]
    )

with summary_file.open("w") as out:
    out.write(
        "\t".join(
            [
                "population",
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
        + "\n"
    )

    for row in summary_rows:
        out.write("\t".join(map(str, row)) + "\n")

print(f"summary\t{summary_file}")
print("FINAL\tPASS")
PY
