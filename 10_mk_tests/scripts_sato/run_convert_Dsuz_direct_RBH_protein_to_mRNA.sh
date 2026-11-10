#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0-02:00:00
#SBATCH -J Dsuz_RBH_to_mRNA
#SBATCH -o /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/Dsuz_RBH_to_mRNA.%j.out
#SBATCH -e /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/Dsuz_RBH_to_mRNA.%j.err

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817
MKBASE=/ru-auth/local/home/asato/scratch/MKtest

mkdir -p \
  "${BASE}/orth_longest" \
  "${BASE}/checks" \
  "${BASE}/logs"

python3 - <<'PY'
from pathlib import Path
import csv
from collections import Counter

base = Path(
    "/ru-auth/local/home/asato/scratch/MKtest/"
    "reanalysis_longest_fix_20260817"
)
mkbase = Path(
    "/ru-auth/local/home/asato/scratch/MKtest"
)

orth = base / "orth_longest"
checks = base / "checks"

old_dsub_dsuz = (
    orth
    / "Dsub_Dsuz_IsoJpt1.0.RBH.longest.protein.tsv"
)

dsuz_dsub_protein = (
    orth
    / "Dsuz_Dsub.direct.RBH.longest.protein.tsv"
)

dsuz_dbia_protein = (
    orth
    / "Dsuz_Dbia.direct.RBH.longest.protein.tsv"
)

dsuz_dmel_protein = (
    orth
    / "Dsuz_Dmel.direct.RBH.longest.protein.tsv"
)

dsub_map_file = (
    base
    / "pep_longest"
    / "Dsub.longest_per_gene.map.tsv"
)

dsuz_map_file = (
    mkbase
    / "orth_longest"
    / "Dsuz.IsoJpt1.0.protein_to_mrna.tsv"
)

dbia_map_file = (
    mkbase
    / "orth_longest"
    / "Dbia.protein_to_mrna.tsv"
)

dmel_map_file = (
    mkbase
    / "orth_longest"
    / "Dmel.FBpp_to_mrna_gene_product.tsv"
)

required = [
    old_dsub_dsuz,
    dsuz_dbia_protein,
    dsuz_dmel_protein,
    dsub_map_file,
    dsuz_map_file,
    dbia_map_file,
    dmel_map_file,
]

for path in required:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Required input missing or empty: {path}"
        )


def normalize_mrna(x):
    x = x.strip()

    if not x:
        return ""

    if x.startswith("rna-"):
        return x

    if x.startswith(("NM_", "XM_", "NR_", "XR_")):
        return "rna-" + x

    return x


def read_two_column_map(path, has_header=False):
    mapping = {}
    conflicts = []

    with path.open() as fh:
        reader = csv.reader(fh, delimiter="\t")

        if has_header:
            next(reader, None)

        for row in reader:
            if len(row) < 2:
                continue

            protein = row[0].strip()
            mrna = normalize_mrna(row[1])

            if not protein or not mrna:
                continue

            previous = mapping.get(protein)

            if previous is not None and previous != mrna:
                conflicts.append(
                    (protein, previous, mrna)
                )
                continue

            mapping[protein] = mrna

    return mapping, conflicts


def read_dsub_map(path):
    mapping = {}
    conflicts = []

    with path.open() as fh:
        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        for row in reader:
            protein = row["protein_accession"].strip()
            mrna = normalize_mrna(
                row["mRNA_accession"]
            )

            if not protein or not mrna:
                continue

            previous = mapping.get(protein)

            if previous is not None and previous != mrna:
                conflicts.append(
                    (protein, previous, mrna)
                )
                continue

            mapping[protein] = mrna

    return mapping, conflicts


dsub_map, dsub_conflicts = read_dsub_map(
    dsub_map_file
)

dsuz_map, dsuz_conflicts = read_two_column_map(
    dsuz_map_file,
    has_header=True,
)

dbia_map, dbia_conflicts = read_two_column_map(
    dbia_map_file,
    has_header=False,
)

dmel_map, dmel_conflicts = read_two_column_map(
    dmel_map_file,
    has_header=False,
)

all_conflicts = (
    [("Dsub",) + x for x in dsub_conflicts]
    + [("Dsuz",) + x for x in dsuz_conflicts]
    + [("Dbia",) + x for x in dbia_conflicts]
    + [("Dmel",) + x for x in dmel_conflicts]
)

conflict_file = (
    checks
    / "Dsuz_direct_RBH_protein_to_mRNA."
      "mapping_conflicts.tsv"
)

with conflict_file.open("w") as out:
    out.write(
        "species\tprotein_id\tmRNA_1\tmRNA_2\n"
    )

    for row in all_conflicts:
        out.write("\t".join(row) + "\n")

if all_conflicts:
    raise RuntimeError(
        f"Protein-to-mRNA mapping conflicts detected: "
        f"{len(all_conflicts)}"
    )


# ----------------------------------------------------------
# Reverse the new Dsub-Dsuz RBH table for Dsuz focal use
#
# Original:
#   Dsuz protein, Dsub protein
#
# New:
#   Dsub outgroup protein, Dsuz focal protein
# ----------------------------------------------------------

reversed_pairs = []

with old_dsub_dsuz.open() as fh:
    for line in fh:
        fields = line.rstrip("\n").split("\t")

        if len(fields) < 2:
            continue

        dsuz_protein = fields[0].strip()
        dsub_protein = fields[1].strip()

        reversed_pairs.append(
            (dsub_protein, dsuz_protein)
        )

reversed_pairs = sorted(set(reversed_pairs))

with dsuz_dsub_protein.open("w") as out:
    for dsub_protein, dsuz_protein in reversed_pairs:
        out.write(
            f"{dsub_protein}\t{dsuz_protein}\n"
        )


def read_rbh(path):
    pairs = []

    with path.open() as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 2:
                continue

            outgroup_protein = fields[0].strip()
            dsuz_protein = fields[1].strip()

            if outgroup_protein and dsuz_protein:
                pairs.append(
                    (outgroup_protein, dsuz_protein)
                )

    return pairs


def duplicated_values(values):
    counts = Counter(values)

    return sorted(
        value
        for value, count in counts.items()
        if count > 1
    )


def convert_comparison(
    comparison,
    protein_rbh,
    outgroup_map,
):
    protein_pairs = read_rbh(protein_rbh)

    output_mrna = (
        orth
        / f"{comparison}.direct.RBH."
          "longest.mRNA.tsv"
    )

    missing_file = (
        checks
        / f"{comparison}.direct_RBH."
          "protein_to_mRNA_missing.tsv"
    )

    converted = []
    missing = []

    for outgroup_protein, dsuz_protein in protein_pairs:
        outgroup_mrna = outgroup_map.get(
            outgroup_protein
        )
        dsuz_mrna = dsuz_map.get(
            dsuz_protein
        )

        reasons = []

        if outgroup_mrna is None:
            reasons.append(
                "outgroup_protein_not_mapped"
            )

        if dsuz_mrna is None:
            reasons.append(
                "Dsuz_protein_not_mapped"
            )

        if reasons:
            missing.append(
                (
                    outgroup_protein,
                    dsuz_protein,
                    ";".join(reasons),
                )
            )
            continue

        converted.append(
            (outgroup_mrna, dsuz_mrna)
        )

    with output_mrna.open("w") as out:
        for outgroup_mrna, dsuz_mrna in converted:
            out.write(
                f"{outgroup_mrna}\t{dsuz_mrna}\n"
            )

    with missing_file.open("w") as out:
        out.write(
            "outgroup_protein\tDsuz_protein\treason\n"
        )

        for row in missing:
            out.write("\t".join(row) + "\n")

    duplicated_outgroup = duplicated_values(
        x[0] for x in converted
    )
    duplicated_dsuz = duplicated_values(
        x[1] for x in converted
    )

    summary = {
        "comparison": comparison,
        "input_protein_RBH_pairs": len(protein_pairs),
        "converted_mRNA_pairs": len(converted),
        "unmapped_pairs": len(missing),
        "duplicated_outgroup_mRNA": len(
            duplicated_outgroup
        ),
        "duplicated_Dsuz_mRNA": len(
            duplicated_dsuz
        ),
        "protein_RBH": str(protein_rbh),
        "mRNA_RBH": str(output_mrna),
        "missing_table": str(missing_file),
    }

    if missing:
        raise RuntimeError(
            f"{comparison}: {len(missing)} "
            "protein RBH pairs could not be converted"
        )

    if duplicated_outgroup or duplicated_dsuz:
        raise RuntimeError(
            f"{comparison}: duplicated mRNA IDs detected"
        )

    return summary


summaries = []

summaries.append(
    convert_comparison(
        "Dsuz_Dsub",
        dsuz_dsub_protein,
        dsub_map,
    )
)

summaries.append(
    convert_comparison(
        "Dsuz_Dbia",
        dsuz_dbia_protein,
        dbia_map,
    )
)

summaries.append(
    convert_comparison(
        "Dsuz_Dmel",
        dsuz_dmel_protein,
        dmel_map,
    )
)

summary_file = (
    checks
    / "Dsuz_direct_RBH_protein_to_mRNA."
      "summary.tsv"
)

fields = [
    "comparison",
    "input_protein_RBH_pairs",
    "converted_mRNA_pairs",
    "unmapped_pairs",
    "duplicated_outgroup_mRNA",
    "duplicated_Dsuz_mRNA",
    "protein_RBH",
    "mRNA_RBH",
    "missing_table",
]

with summary_file.open("w", newline="") as out:
    writer = csv.DictWriter(
        out,
        fieldnames=fields,
        delimiter="\t",
    )

    writer.writeheader()

    for row in summaries:
        writer.writerow(row)

print(summary_file.read_text())


# ----------------------------------------------------------
# Orco check
# ----------------------------------------------------------

orco_dsub = "rna-XM_037873416.1"
orco_dsuz = "rna-NM_001328601.1"

dsub_mrna_file = (
    orth
    / "Dsuz_Dsub.direct.RBH.longest.mRNA.tsv"
)

orco_found = False

with dsub_mrna_file.open() as fh:
    for line in fh:
        fields = line.rstrip("\n").split("\t")

        if (
            len(fields) >= 2
            and fields[0] == orco_dsub
            and fields[1] == orco_dsuz
        ):
            orco_found = True
            break

print("===== Orco check =====")
print(
    f"Dsub={orco_dsub}\t"
    f"Dsuz={orco_dsuz}\t"
    f"present={orco_found}"
)

if not orco_found:
    raise RuntimeError(
        "Expected Dsub-Dsuz Orco RBH pair was not found"
    )

print("Done.")
PY
