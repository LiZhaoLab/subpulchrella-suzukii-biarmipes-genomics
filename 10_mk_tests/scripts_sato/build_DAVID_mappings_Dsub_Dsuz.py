#!/usr/bin/env python3

from pathlib import Path
from urllib.parse import unquote
from collections import Counter
import csv
import re

BASE = Path(
    "/ru-auth/local/home/asato/scratch/MKtest/"
    "reanalysis_longest_fix_20260817"
)

GFF = Path(
    "/ru-auth/local/home/asato/scratch/genome/Dmel/"
    "GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff"
)

RBH_FILES = {
    "Dsub": BASE / "orth_longest/Dsub_Dmel.RBH.longest.mRNA.tsv",
    "Dsuz": BASE / "orth_longest/Dsuz_Dmel.direct.RBH.longest.mRNA.tsv",
}

OUTDIR = BASE / "david_mapping"
OUTDIR.mkdir(parents=True, exist_ok=True)


def parse_attributes(text):
    result = {}

    for item in text.rstrip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = unquote(value)

    return result


def extract_fbgn(attributes):
    for key in ("Dbxref", "Parent", "ID"):
        match = re.search(
            r"FBgn[0-9]+",
            attributes.get(key, "")
        )

        if match:
            return match.group(0)

    return ""


# Dmel mRNA -> FBgn and official symbol
dmel_map = {}

with GFF.open() as handle:
    for line in handle:
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) != 9:
            continue

        if fields[2] not in {"mRNA", "transcript"}:
            continue

        attributes = parse_attributes(fields[8])

        dmel_mrna = re.sub(
            r"^rna-",
            "",
            attributes.get("ID", "")
        )

        fbgn = extract_fbgn(attributes)

        symbol = attributes.get(
            "gene",
            attributes.get("Name", "")
        )

        if dmel_mrna and fbgn:
            dmel_map[dmel_mrna] = (fbgn, symbol)

print(
    "Dmel_mRNA_to_FBgn_entries\t"
    f"{len(dmel_map)}"
)


for focal_species, rbh_file in RBH_FILES.items():
    pairs = []

    with rbh_file.open() as handle:
        for line in handle:
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 2:
                raise RuntimeError(
                    f"Invalid RBH row: {line.rstrip()}"
                )

            # column 1 = Dmel mRNA
            # column 2 = focal-species mRNA
            dmel_mrna = re.sub(r"^rna-", "", fields[0])
            focal_mrna = re.sub(r"^rna-", "", fields[1])

            pairs.append((focal_mrna, dmel_mrna))

    duplicated = [
        focal_mrna
        for focal_mrna, count in Counter(
            x[0] for x in pairs
        ).items()
        if count > 1
    ]

    if duplicated:
        raise RuntimeError(
            f"{focal_species}: duplicated focal IDs: "
            f"{len(duplicated)}"
        )

    output_file = (
        OUTDIR /
        f"{focal_species}_to_Dmel_FBgn.tsv"
    )

    missing_file = (
        OUTDIR /
        f"{focal_species}_to_Dmel_FBgn.missing.tsv"
    )

    mapped = 0
    missing = 0

    with (
        output_file.open("w", newline="") as output_handle,
        missing_file.open("w", newline="") as missing_handle
    ):
        writer = csv.writer(
            output_handle,
            delimiter="\t",
            lineterminator="\n"
        )

        missing_writer = csv.writer(
            missing_handle,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writerow([
            "focal_mrna",
            "dmel_mrna",
            "Dmel_FBgn",
            "Dmel_symbol_official",
        ])

        missing_writer.writerow([
            "focal_mrna",
            "dmel_mrna",
            "reason",
        ])

        for focal_mrna, dmel_mrna in pairs:
            annotation = dmel_map.get(dmel_mrna)

            if annotation is None:
                missing_writer.writerow([
                    focal_mrna,
                    dmel_mrna,
                    "Dmel_mRNA_not_mapped_to_FBgn_in_GFF",
                ])

                missing += 1
                continue

            fbgn, symbol = annotation

            writer.writerow([
                focal_mrna,
                dmel_mrna,
                fbgn,
                symbol,
            ])

            mapped += 1

    print(f"\n===== {focal_species} =====")
    print(f"input_RBH_pairs\t{len(pairs)}")
    print(f"mapped_rows\t{mapped}")
    print(f"missing_rows\t{missing}")
    print(f"output\t{output_file}")
    print(f"missing\t{missing_file}")

print("\nDone.")
