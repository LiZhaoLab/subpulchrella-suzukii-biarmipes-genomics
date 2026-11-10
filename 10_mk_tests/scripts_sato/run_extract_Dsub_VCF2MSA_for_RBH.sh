#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-01:00:00
#SBATCH -J extract_Dsub_MSA
#SBATCH -o /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/extract_Dsub_MSA.%j.out
#SBATCH -e /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/extract_Dsub_MSA.%j.err

set -euo pipefail

NEW_BASE=/ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817

CDS_MAP=/ru-auth/local/home/jpeng/scratch/Collaborations/Dsuz/Airi/vcf2msa/Dsub.cds_map.csv
MSA_INDEX=/ru-auth/local/home/asato/scratch/MKtest/Dsubpop/vcf2msa_2026/Dsub_msa/cds_snp_index.csv
RAW_MSA_DIR=/ru-auth/local/home/asato/scratch/MKtest/Dsubpop/vcf2msa_2026/Dsub_msa

export NEW_BASE
export CDS_MAP
export MSA_INDEX
export RAW_MSA_DIR

for f in "${CDS_MAP}" "${MSA_INDEX}"; do
    if [[ ! -s "${f}" ]]; then
        echo "ERROR: missing or empty input: ${f}" >&2
        exit 1
    fi
done

if [[ ! -d "${RAW_MSA_DIR}" ]]; then
    echo "ERROR: missing MSA directory: ${RAW_MSA_DIR}" >&2
    exit 1
fi

python3 - <<'PY'
import csv
import os
import shutil
from collections import defaultdict
from pathlib import Path

base = Path(os.environ["NEW_BASE"])
cds_map_file = Path(os.environ["CDS_MAP"])
index_file = Path(os.environ["MSA_INDEX"])
raw_msa_dir = Path(os.environ["RAW_MSA_DIR"])

comparisons = {
    "Dsub_vs_Dsuz": (
        base / "orth_longest"
        / "Dsub_Dsuz_IsoJpt1.0.RBH.longest.mRNA.tsv"
    ),
    "Dsub_vs_Dbia": (
        base / "orth_longest"
        / "Dsub_Dbia.RBH.longest.mRNA.tsv"
    ),
    "Dsub_vs_Dmel": (
        base / "orth_longest"
        / "Dsub_Dmel.RBH.longest.mRNA.tsv"
    ),
}


def normalize_mrna(value):
    value = value.strip().split()[0]

    if value.startswith("rna-"):
        value = value[4:]

    return value


def output_mrna(value):
    return f"rna-{normalize_mrna(value)}"


# ============================================================
# Read transcript -> CDS-coordinate mapping
#
# Important: this remains a forward mapping.
# We do not create a coordinate -> single-transcript mapping.
# ============================================================

transcript_to_coordinate = {}
mapping_conflicts = []

with cds_map_file.open() as fh:
    reader = csv.DictReader(fh, delimiter="\t")

    if reader.fieldnames is None:
        raise RuntimeError("CDS mapping file has no header")

    for row in reader:
        transcript = output_mrna(row["name"])
        coordinate = row["coor"].strip()

        old = transcript_to_coordinate.get(transcript)

        if old is not None and old != coordinate:
            mapping_conflicts.append(
                (transcript, old, coordinate)
            )
        else:
            transcript_to_coordinate[transcript] = coordinate

if mapping_conflicts:
    raise RuntimeError(
        f"Found {len(mapping_conflicts)} transcript-coordinate conflicts"
    )


# ============================================================
# Read CDS coordinate -> all VCF2MSA indices
#
# Multiple indices are intentionally retained for shared
# CDS coordinates.
# ============================================================

coordinate_to_indices = defaultdict(list)
index_snp = {}

with index_file.open() as fh:
    reader = csv.DictReader(fh)

    required = {"index", "header", "SNP"}

    if reader.fieldnames is None or not required.issubset(reader.fieldnames):
        raise RuntimeError(
            f"Unexpected index header: {reader.fieldnames}"
        )

    for row in reader:
        index = row["index"].strip()
        coordinate = row["header"].strip()
        snp = row["SNP"].strip()

        coordinate_to_indices[coordinate].append(index)
        index_snp[index] = snp

for coordinate in coordinate_to_indices:
    coordinate_to_indices[coordinate].sort(
        key=lambda value: int(value)
    )


# ============================================================
# Extract each RBH set
# ============================================================

summary_rows = []
fatal_problems = 0

for comparison, rbh_file in comparisons.items():
    if not rbh_file.is_file() or rbh_file.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty RBH mRNA table: {rbh_file}"
        )

    outdir = (
        base
        / "vcf2msa_longest_RBH"
        / comparison
    )
    outdir.mkdir(parents=True, exist_ok=True)

    # Remove only previous FASTA outputs from this new-analysis folder.
    for old_fasta in outdir.glob("rna-*.fasta"):
        old_fasta.unlink()

    manifest_file = (
        base
        / "checks"
        / f"{comparison}.VCF2MSA_extraction_manifest.tsv"
    )

    missing_file = (
        base
        / "checks"
        / f"{comparison}.VCF2MSA_not_extracted.tsv"
    )

    dsub_transcripts = []

    with rbh_file.open() as fh:
        reader = csv.reader(fh, delimiter="\t")

        for row in reader:
            if len(row) < 2:
                continue

            dsub_transcripts.append(
                output_mrna(row[1])
            )

    if len(dsub_transcripts) != len(set(dsub_transcripts)):
        raise RuntimeError(
            f"Duplicated Dsub mRNA IDs in {rbh_file}"
        )

    manifest = []
    missing = []

    n_multiple_index_coordinates = 0

    for transcript in dsub_transcripts:
        coordinate = transcript_to_coordinate.get(transcript)

        if coordinate is None:
            missing.append(
                (
                    transcript,
                    "transcript_not_in_CDS_map",
                    "",
                    "",
                )
            )
            fatal_problems += 1
            continue

        candidate_indices = coordinate_to_indices.get(
            coordinate,
            [],
        )

        if not candidate_indices:
            # vcf2msa_2026.py only emitted CDS records meeting its
            # retained-SNP criterion.
            missing.append(
                (
                    transcript,
                    "no_VCF2MSA_index_or_no_retained_SNP",
                    coordinate,
                    "",
                )
            )
            continue

        if len(candidate_indices) > 1:
            n_multiple_index_coordinates += 1

        existing = []

        for index in candidate_indices:
            fasta = raw_msa_dir / f"msa_{index}.fasta"

            if fasta.is_file() and fasta.stat().st_size > 0:
                existing.append((index, fasta))

        if not existing:
            candidate_snp_values = [
                str(index_snp.get(idx, ""))
                for idx in candidate_indices
            ]

            if (
                candidate_snp_values
                and all(snp == "0" for snp in candidate_snp_values)
            ):
                missing.append(
                    (
                        transcript,
                        "no_VCF2MSA_index_or_no_retained_SNP",
                        coordinate,
                        ";".join(candidate_indices),
                    )
                )
            else:
                missing.append(
                    (
                        transcript,
                        "index_present_but_MSA_file_missing",
                        coordinate,
                        ";".join(candidate_indices),
                    )
                )
                fatal_problems += 1

            continue

        # For shared coordinates, all candidates correspond to the
        # same CDS coordinates. Select the lowest available index
        # deterministically while preserving the RBH transcript ID.
        selected_index, source_fasta = existing[0]

        destination = outdir / f"{transcript}.fasta"
        shutil.copy2(source_fasta, destination)

        manifest.append(
            (
                transcript,
                coordinate,
                ";".join(candidate_indices),
                selected_index,
                index_snp.get(selected_index, ""),
                str(source_fasta),
                str(destination),
            )
        )

    with manifest_file.open("w") as out:
        out.write(
            "Dsub_mRNA\tCDS_coordinate\tcandidate_indices\t"
            "selected_index\tSNP\tsource_MSA\toutput_MSA\n"
        )

        for row in manifest:
            out.write("\t".join(row) + "\n")

    with missing_file.open("w") as out:
        out.write(
            "Dsub_mRNA\treason\tCDS_coordinate\t"
            "candidate_indices\n"
        )

        for row in missing:
            out.write("\t".join(row) + "\n")

    summary_rows.append(
        (
            comparison,
            len(dsub_transcripts),
            len(manifest),
            len(missing),
            sum(
                row[1] == "transcript_not_in_CDS_map"
                for row in missing
            ),
            sum(
                row[1] == "no_VCF2MSA_index_or_no_retained_SNP"
                for row in missing
            ),
            sum(
                row[1] == "index_present_but_MSA_file_missing"
                for row in missing
            ),
            n_multiple_index_coordinates,
        )
    )

    print(f"===== {comparison} =====")
    print(f"RBH_Dsub_transcripts\t{len(dsub_transcripts)}")
    print(f"extracted_MSA_files\t{len(manifest)}")
    print(f"not_extracted\t{len(missing)}")
    print(
        "shared_coordinate_with_multiple_indices\t"
        f"{n_multiple_index_coordinates}"
    )
    print(f"output_directory\t{outdir}")


summary_file = (
    base
    / "checks"
    / "Dsub_VCF2MSA_RBH_extraction.summary.tsv"
)

with summary_file.open("w") as out:
    out.write(
        "comparison\tRBH_Dsub_transcripts\textracted_MSA_files\t"
        "not_extracted\ttranscript_not_in_CDS_map\t"
        "no_index_or_no_retained_SNP\t"
        "index_present_but_MSA_missing\t"
        "shared_coordinate_multiple_indices\n"
    )

    for row in summary_rows:
        out.write("\t".join(map(str, row)) + "\n")


# ============================================================
# Orco check
# ============================================================

orco = "rna-XM_037873416.1"

print("===== Orco extraction check =====")

for comparison in comparisons:
    orco_file = (
        base
        / "vcf2msa_longest_RBH"
        / comparison
        / f"{orco}.fasta"
    )

    print(
        f"{comparison}\t"
        f"{'present' if orco_file.is_file() else 'absent'}\t"
        f"{orco_file}"
    )


print(f"summary\t{summary_file}")

if fatal_problems:
    raise RuntimeError(
        f"Found {fatal_problems} unexpected extraction failures. "
        "See checks/*.VCF2MSA_not_extracted.tsv"
    )

print("Done.")
PY
