#!/usr/bin/env python3

from pathlib import Path
import csv
import math


BASE = Path(
    "/ru-auth/local/home/asato/scratch/"
    "synteny_4species/reanalysis_20260824"
)

ALIGNMENT_ROOT = BASE / "alignments"
OUTPUT_DIR = BASE / "checks" / "alignment_distributions"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

COMPARISONS = [
    "Dsuz_vs_Dsub",
    "Dsub_vs_Dbia",
    "Dbia_vs_Dmel",
]

FILTER_TYPES = {
    "many_to_many": "many_to_many",
    "one_to_one": "one_to_one",
}

LENGTH_THRESHOLDS = [100, 500, 1000, 5000, 10000]
IDENTITY_THRESHOLDS = [0, 40, 50, 60, 70, 80]


def percentile(values, probability):
    if not values:
        return math.nan

    ordered = sorted(values)

    if len(ordered) == 1:
        return ordered[0]

    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)

    if lower == upper:
        return ordered[lower]

    weight = position - lower

    return (
        ordered[lower] * (1 - weight)
        + ordered[upper] * weight
    )


def read_coords(path):
    records = []

    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n\r")

            if not line:
                continue

            fields = line.split("\t")

            if len(fields) < 7:
                raise RuntimeError(
                    f"Unexpected coordinate row in "
                    f"{path}:{line_number}: {fields}"
                )

            try:
                ref_start = int(fields[0])
                ref_end = int(fields[1])
                qry_start = int(fields[2])
                qry_end = int(fields[3])
                ref_length = int(fields[4])
                qry_length = int(fields[5])
                identity = float(fields[6])
            except ValueError as error:
                raise RuntimeError(
                    f"Non-numeric coordinate row in "
                    f"{path}:{line_number}: {fields[:7]}"
                ) from error

            alignment_length = min(
                abs(ref_length),
                abs(qry_length),
            )

            records.append(
                {
                    "ref_start": ref_start,
                    "ref_end": ref_end,
                    "qry_start": qry_start,
                    "qry_end": qry_end,
                    "alignment_length": alignment_length,
                    "identity": identity,
                }
            )

    return records


summary_rows = []
threshold_rows = []

for comparison in COMPARISONS:
    for filter_name, filename_label in FILTER_TYPES.items():
        path = (
            ALIGNMENT_ROOT
            / comparison
            / (
                f"{comparison}.promer."
                f"{filename_label}.L100.coords.tsv"
            )
        )

        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(
                f"Missing or empty coordinate file: {path}"
            )

        records = read_coords(path)

        lengths = [
            record["alignment_length"]
            for record in records
        ]

        identities = [
            record["identity"]
            for record in records
        ]

        summary_rows.append(
            {
                "comparison": comparison,
                "filter": filter_name,
                "alignments": len(records),
                "length_min": min(lengths),
                "length_P05": percentile(lengths, 0.05),
                "length_median": percentile(lengths, 0.50),
                "length_P95": percentile(lengths, 0.95),
                "length_P99": percentile(lengths, 0.99),
                "length_max": max(lengths),
                "identity_min": min(identities),
                "identity_P05": percentile(identities, 0.05),
                "identity_median": percentile(
                    identities, 0.50
                ),
                "identity_P95": percentile(
                    identities, 0.95
                ),
                "identity_P99": percentile(
                    identities, 0.99
                ),
                "identity_max": max(identities),
                "source_file": str(path),
            }
        )

        total = len(records)

        for minimum_length in LENGTH_THRESHOLDS:
            for minimum_identity in IDENTITY_THRESHOLDS:
                retained = sum(
                    record["alignment_length"]
                    >= minimum_length
                    and record["identity"]
                    >= minimum_identity
                    for record in records
                )

                threshold_rows.append(
                    {
                        "comparison": comparison,
                        "filter": filter_name,
                        "minimum_length": minimum_length,
                        "minimum_identity": minimum_identity,
                        "input_alignments": total,
                        "retained_alignments": retained,
                        "retained_percent": (
                            100 * retained / total
                            if total > 0
                            else math.nan
                        ),
                    }
                )

summary_path = (
    OUTPUT_DIR / "promer_alignment_distribution_summary.tsv"
)

threshold_path = (
    OUTPUT_DIR / "promer_filter_threshold_candidates.tsv"
)

with summary_path.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(summary_rows[0].keys()),
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(summary_rows)

with threshold_path.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(threshold_rows[0].keys()),
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(threshold_rows)

print(f"comparisons\t{len(COMPARISONS)}")
print(f"summary_rows\t{len(summary_rows)}")
print(f"threshold_rows\t{len(threshold_rows)}")
print(f"summary\t{summary_path}")
print(f"thresholds\t{threshold_path}")
print("FINAL: PASS")
