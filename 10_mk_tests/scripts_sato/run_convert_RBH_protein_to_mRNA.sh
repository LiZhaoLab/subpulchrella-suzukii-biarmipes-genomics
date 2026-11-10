#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-00:15:00
#SBATCH -J RBH_to_mRNA
#SBATCH -o /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/RBH_to_mRNA.%j.out
#SBATCH -e /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/RBH_to_mRNA.%j.err

set -euo pipefail

NEW_BASE=/ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817
OLD_BASE=/ru-auth/local/home/asato/scratch/MKtest

Dsub_MAP=${NEW_BASE}/pep_longest/Dsub.longest_per_gene.map.tsv
Dsuz_MAP=${OLD_BASE}/orth_longest/Dsuz.IsoJpt1.0.protein_to_mrna.tsv
Dbia_MAP=${OLD_BASE}/orth_longest/Dbia.protein_to_mrna.tsv
Dmel_MAP=${OLD_BASE}/orth_longest/Dmel.FBpp_to_mrna_gene_product.tsv

RBH_Dsuz=${NEW_BASE}/orth_longest/Dsub_Dsuz_IsoJpt1.0.RBH.longest.protein.tsv
RBH_Dbia=${NEW_BASE}/orth_longest/Dsub_Dbia.RBH.longest.protein.tsv
RBH_Dmel=${NEW_BASE}/orth_longest/Dsub_Dmel.RBH.longest.protein.tsv

export NEW_BASE
export Dsub_MAP
export Dsuz_MAP
export Dbia_MAP
export Dmel_MAP
export RBH_Dsuz
export RBH_Dbia
export RBH_Dmel

for f in \
  "${Dsub_MAP}" \
  "${Dsuz_MAP}" \
  "${Dbia_MAP}" \
  "${Dmel_MAP}" \
  "${RBH_Dsuz}" \
  "${RBH_Dbia}" \
  "${RBH_Dmel}"
do
    if [[ ! -s "${f}" ]]; then
        echo "ERROR: missing or empty input: ${f}" >&2
        exit 1
    fi
done

python3 - <<'PY'
import csv
import os
import sys
from pathlib import Path

base = Path(os.environ["NEW_BASE"])

dsub_map_file = Path(os.environ["Dsub_MAP"])
dsuz_map_file = Path(os.environ["Dsuz_MAP"])
dbia_map_file = Path(os.environ["Dbia_MAP"])
dmel_map_file = Path(os.environ["Dmel_MAP"])

comparisons = {
    "Dsub_Dsuz_IsoJpt1.0": {
        "rbh": Path(os.environ["RBH_Dsuz"]),
        "outgroup_map": dsuz_map_file,
    },
    "Dsub_Dbia": {
        "rbh": Path(os.environ["RBH_Dbia"]),
        "outgroup_map": dbia_map_file,
    },
    "Dsub_Dmel": {
        "rbh": Path(os.environ["RBH_Dmel"]),
        "outgroup_map": dmel_map_file,
    },
}


def normalize_protein(value):
    value = value.strip().split()[0]

    for prefix in ("cds-", "protein-"):
        if value.startswith(prefix):
            value = value[len(prefix):]

    return value


def normalize_mrna(value):
    value = value.strip().split()[0]

    if value.startswith("rna-"):
        value = value[4:]

    return value


def output_mrna(value):
    value = normalize_mrna(value)
    return f"rna-{value}"


# ============================================================
# Dsub protein-to-mRNA map
#
# Columns:
# gene_key, GeneID, gene_symbol, mRNA_accession,
# protein_accession, ...
# ============================================================

dsub_protein_to_mrna = {}

with dsub_map_file.open() as fh:
    reader = csv.DictReader(fh, delimiter="\t")

    required = {
        "mRNA_accession",
        "protein_accession",
    }

    if reader.fieldnames is None or not required.issubset(reader.fieldnames):
        raise RuntimeError(
            f"Unexpected Dsub mapping header: {reader.fieldnames}"
        )

    for row in reader:
        protein = normalize_protein(row["protein_accession"])
        mrna = output_mrna(row["mRNA_accession"])

        old = dsub_protein_to_mrna.get(protein)

        if old is not None and old != mrna:
            raise RuntimeError(
                f"Conflicting Dsub mapping for {protein}: "
                f"{old} versus {mrna}"
            )

        dsub_protein_to_mrna[protein] = mrna


# ============================================================
# Generic outgroup protein-to-mRNA map reader
#
# The first two columns are protein ID and mRNA ID.
# Header rows, when present, are skipped.
# ============================================================

def read_outgroup_map(path):
    mapping = {}

    with path.open() as fh:
        reader = csv.reader(fh, delimiter="\t")

        for row in reader:
            if len(row) < 2:
                continue

            first = row[0].strip().lower()
            second = row[1].strip().lower()

            if (
                first in {
                    "protein_id",
                    "protein_accession",
                    "fbpp",
                }
                or second in {
                    "mrna_id",
                    "mrna_accession",
                    "transcript_id",
                }
            ):
                continue

            protein = normalize_protein(row[0])
            mrna = output_mrna(row[1])

            old = mapping.get(protein)

            if old is not None and old != mrna:
                raise RuntimeError(
                    f"Conflicting mapping in {path}: "
                    f"{protein}: {old} versus {mrna}"
                )

            mapping[protein] = mrna

    return mapping


# ============================================================
# Convert each RBH table
# ============================================================

summary_rows = []

for comparison, info in comparisons.items():
    rbh_file = info["rbh"]
    outgroup_map = read_outgroup_map(info["outgroup_map"])

    output_file = (
        base
        / "orth_longest"
        / f"{comparison}.RBH.longest.mRNA.tsv"
    )

    unmapped_file = (
        base
        / "checks"
        / f"{comparison}.RBH.protein_to_mRNA.unmapped.tsv"
    )

    converted = []
    unmapped = []

    with rbh_file.open() as fh:
        reader = csv.reader(fh, delimiter="\t")

        for line_number, row in enumerate(reader, start=1):
            if len(row) < 2:
                continue

            outgroup_protein = normalize_protein(row[0])
            dsub_protein = normalize_protein(row[1])

            outgroup_mrna = outgroup_map.get(outgroup_protein)
            dsub_mrna = dsub_protein_to_mrna.get(dsub_protein)

            reasons = []

            if outgroup_mrna is None:
                reasons.append("outgroup_protein_unmapped")

            if dsub_mrna is None:
                reasons.append("Dsub_protein_unmapped")

            if reasons:
                unmapped.append(
                    (
                        str(line_number),
                        outgroup_protein,
                        dsub_protein,
                        ";".join(reasons),
                    )
                )
                continue

            converted.append(
                (
                    outgroup_mrna,
                    dsub_mrna,
                    outgroup_protein,
                    dsub_protein,
                )
            )

    # Ensure that the converted table remains one-to-one.
    outgroup_mrnas = [row[0] for row in converted]
    dsub_mrnas = [row[1] for row in converted]

    duplicate_outgroup = (
        len(outgroup_mrnas) - len(set(outgroup_mrnas))
    )
    duplicate_dsub = (
        len(dsub_mrnas) - len(set(dsub_mrnas))
    )

    with output_file.open("w") as out:
        for (
            outgroup_mrna,
            dsub_mrna,
            outgroup_protein,
            dsub_protein,
        ) in converted:
            out.write(
                f"{outgroup_mrna}\t{dsub_mrna}\n"
            )

    with unmapped_file.open("w") as out:
        out.write(
            "RBH_line\toutgroup_protein\tDsub_protein\treason\n"
        )

        for row in unmapped:
            out.write("\t".join(row) + "\n")

    summary_rows.append(
        (
            comparison,
            len(converted) + len(unmapped),
            len(converted),
            len(unmapped),
            duplicate_outgroup,
            duplicate_dsub,
            len(outgroup_map),
        )
    )

    print(f"===== {comparison} =====")
    print(f"input_RBH_pairs\t{len(converted) + len(unmapped)}")
    print(f"converted_pairs\t{len(converted)}")
    print(f"unmapped_pairs\t{len(unmapped)}")
    print(f"duplicated_outgroup_mRNA\t{duplicate_outgroup}")
    print(f"duplicated_Dsub_mRNA\t{duplicate_dsub}")
    print(f"output\t{output_file}")

    if unmapped:
        print(
            f"ERROR: {comparison} has {len(unmapped)} unmapped pairs.",
            file=sys.stderr,
        )


summary_file = (
    base
    / "checks"
    / "RBH_protein_to_mRNA_conversion.summary.tsv"
)

with summary_file.open("w") as out:
    out.write(
        "comparison\tinput_RBH_pairs\tconverted_pairs\t"
        "unmapped_pairs\tduplicated_outgroup_mRNA\t"
        "duplicated_Dsub_mRNA\toutgroup_mapping_entries\n"
    )

    for row in summary_rows:
        out.write("\t".join(map(str, row)) + "\n")


if any(row[3] > 0 for row in summary_rows):
    raise RuntimeError(
        "One or more RBH pairs could not be converted. "
        "See the unmapped tables in checks/."
    )


# ============================================================
# Orco check
# ============================================================

orco_mrna = "rna-XM_037873416.1"

print("===== Orco mRNA check =====")

for comparison in comparisons:
    output_file = (
        base
        / "orth_longest"
        / f"{comparison}.RBH.longest.mRNA.tsv"
    )

    found = []

    with output_file.open() as fh:
        for line in fh:
            row = line.rstrip("\n").split("\t")

            if len(row) >= 2 and row[1] == orco_mrna:
                found.append(line.rstrip("\n"))

    print(comparison)

    if found:
        for line in found:
            print(line)
    else:
        print("Orco not present")


print(f"summary\t{summary_file}")
print("Done.")
PY
