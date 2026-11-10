#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --job-name=TE_MerelRef
#SBATCH --array=0-1%2
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only/logs/TE_MerelRef.%A_%a.out
#SBATCH --error=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only/logs/TE_MerelRef.%A_%a.err

set -euo pipefail

NEW=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only
SPECIES=(Dsub Dsuz)
GENOMES=(
    /ru-auth/local/home/asato/scratch/genome/Dsub_v1.1/GCF_014743375.2_RU_Dsub_v1.1_genomic.fna
    /ru-auth/local/home/asato/scratch/genome/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.fna
)

INDEX=${SLURM_ARRAY_TASK_ID}
[[ "$INDEX" =~ ^[01]$ ]] || {
    echo "ERROR: unexpected array task ID: $INDEX" >&2
    exit 1
}

SAMPLE=${SPECIES[$INDEX]}
GENOME=${GENOMES[$INDEX]}
GENOME_NAME=$(basename "$GENOME")
MASKED="${NEW}/shared_library_annotation_fullmask/${SAMPLE}/${GENOME_NAME}.masked"
LIBRARY="${NEW}/shared_library/Merel.min500.fa"
HIERARCHY="${NEW}/population_TE/config/Merel.min500.hierarchy.tsv"
MASK_QC="${NEW}/checks/repeatmasker_fullmask/${SAMPLE}.summary.tsv"
HIER_QC="${NEW}/checks/hierarchy/Merel.min500.hierarchy.summary.tsv"

OUTDIR="${NEW}/population_TE/reference"
CHECKDIR="${NEW}/checks/reference_merge"
OUTPUT="${OUTDIR}/${SAMPLE}.Merel.TE_merged_reference.min500.fasta"
SUMMARY="${CHECKDIR}/${SAMPLE}.Merel.TE_merged_reference.min500.summary.tsv"
TEMP="${OUTPUT}.tmp.${SLURM_JOB_ID}_${INDEX}"

for FILE in "$GENOME" "$MASKED" "$LIBRARY" "$HIERARCHY" "$MASK_QC" "$HIER_QC"; do
    [[ -s "$FILE" ]] || {
        echo "ERROR: missing or empty input: $FILE" >&2
        exit 1
    }
done

for QC in "$MASK_QC" "$HIER_QC"; do
    awk -F '\t' '$1=="status" && $2=="PASS" {ok=1} END {exit !ok}' "$QC" || {
        echo "ERROR: input QC did not pass: $QC" >&2
        exit 1
    }
done

[[ ! -e "$OUTPUT" && ! -e "$SUMMARY" ]] || {
    echo "ERROR: output already exists; will not overwrite" >&2
    exit 1
}

mkdir -p "$OUTDIR" "$CHECKDIR"
trap 'rm -f -- "$TEMP"' EXIT

python3 - "$SAMPLE" "$GENOME" "$MASKED" "$LIBRARY" \
    "$HIERARCHY" "$MASK_QC" "$TEMP" "$OUTPUT" "$SUMMARY" <<'PY'
import csv
import hashlib
import itertools
import os
import sys

(sample, genome_path, masked_path, library_path, hierarchy_path,
 mask_qc_path, temp_path, output_path, summary_path) = sys.argv[1:]

expected_genomes = {
    "Dsub": (299, 265153789),
    "Dsuz": (55, 282311783),
}
expected_records, expected_bases = expected_genomes[sample]

def records(path):
    name = None
    parts = []
    with open(path, "rb") as handle:
        for raw in handle:
            if raw.startswith(b">"):
                if name is not None:
                    yield name, b"".join(parts).upper()
                name = raw[1:].strip().split()[0].decode("ascii")
                parts = []
            elif raw.strip():
                if name is None:
                    raise RuntimeError(f"Sequence before FASTA header: {path}")
                parts.append(raw.strip())
    if name is not None:
        yield name, b"".join(parts).upper()

with open(mask_qc_path) as handle:
    mask_qc = dict(
        line.rstrip("\n").split("\t", 1)
        for line in handle if "\t" in line
    )
if mask_qc.get("status") != "PASS" or mask_qc.get("masked_genome") != masked_path:
    raise RuntimeError("Mask QC does not match the intended masked genome")

genome_ids = set()
genome_records = genome_bases = original_n = newly_masked = 0

for original, masked in itertools.zip_longest(
    records(genome_path), records(masked_path)
):
    if original is None or masked is None:
        raise RuntimeError("Genome and masked FASTA have different record counts")

    original_id, original_seq = original
    masked_id, masked_seq = masked
    if original_id != masked_id or len(original_seq) != len(masked_seq):
        raise RuntimeError(f"Genome ID or length changed: {original_id}")
    if original_id in genome_ids:
        raise RuntimeError(f"Duplicate genome ID: {original_id}")

    genome_ids.add(original_id)
    genome_records += 1
    genome_bases += len(original_seq)
    original_n += original_seq.count(b"N")

    for before, after in zip(original_seq, masked_seq):
        if before != after:
            if after != ord("N"):
                raise RuntimeError(f"Unexpected non-N change: {original_id}")
            newly_masked += 1

if (genome_records, genome_bases) != (expected_records, expected_bases):
    raise RuntimeError(
        f"Unexpected genome dimensions: {genome_records}, {genome_bases}"
    )

te_ids = set()
te_records = 0
for te_id, sequence in records(library_path):
    if te_id in te_ids or te_id in genome_ids:
        raise RuntimeError(f"Duplicate or conflicting TE ID: {te_id}")
    if len(sequence) < 500:
        raise RuntimeError(f"TE under 500 bp: {te_id}")
    te_ids.add(te_id)
    te_records += 1

if te_records != 2960:
    raise RuntimeError(f"Expected 2960 TE records; found {te_records}")

with open(hierarchy_path, newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    if reader.fieldnames != ["id", "family", "order"]:
        raise RuntimeError("Unexpected hierarchy columns")
    hierarchy_rows = list(reader)

hierarchy_ids = [row["id"] for row in hierarchy_rows]
if len(hierarchy_ids) != len(set(hierarchy_ids)):
    raise RuntimeError("Duplicate hierarchy ID")
if set(hierarchy_ids) != te_ids:
    raise RuntimeError("TE FASTA IDs and hierarchy IDs differ")
if any(not row["family"] or not row["order"] for row in hierarchy_rows):
    raise RuntimeError("Empty family or order")

with open(masked_path, "rb") as handle:
    handle.seek(-1, os.SEEK_END)
    if handle.read(1) != b"\n":
        raise RuntimeError("Masked FASTA does not end with a newline")

digest = hashlib.sha256()
with open(temp_path, "wb") as output:
    for path in (masked_path, library_path):
        with open(path, "rb") as source:
            while chunk := source.read(1024 * 1024):
                output.write(chunk)
                digest.update(chunk)

# 既存ファイルを上書きせず、検証済みの一時ファイルを公開する。
os.link(temp_path, output_path)
os.unlink(temp_path)

with open(summary_path, "x") as handle:
    handle.write("metric\tvalue\n")
    values = [
        ("sample", sample),
        ("original_genome", genome_path),
        ("masked_genome", masked_path),
        ("genome_records", genome_records),
        ("genome_bases", genome_bases),
        ("original_N_bases", original_n),
        ("newly_masked_bases", newly_masked),
        ("TE_library", library_path),
        ("TE_records", te_records),
        ("hierarchy", hierarchy_path),
        ("merged_reference", output_path),
        ("merged_reference_records", genome_records + te_records),
        ("merged_reference_sha256", digest.hexdigest()),
        ("status", "PASS"),
    ]
    for key, value in values:
        handle.write(f"{key}\t{value}\n")

print(open(summary_path).read(), end="")
print("FINAL: PASS")
PY

date
