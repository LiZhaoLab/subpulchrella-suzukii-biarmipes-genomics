#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --job-name=TE_MerelHier
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only/logs/TE_MerelHier.%j.out
#SBATCH --error=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only/logs/TE_MerelHier.%j.err

set -euo pipefail

NEW=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only
ANNOTATION="${NEW}/shared_library/Merel.all_TEs.original.txt"
FULL="${NEW}/shared_library/Merel.all_TEs.original.fa"
MIN500="${NEW}/shared_library/Merel.min500.fa"
OUTDIR="${NEW}/population_TE/config"
CHECKDIR="${NEW}/checks/hierarchy"
HIER="${OUTDIR}/Merel.min500.hierarchy.tsv"
SUMMARY="${CHECKDIR}/Merel.min500.hierarchy.summary.tsv"

for FILE in "$ANNOTATION" "$FULL" "$MIN500"; do
    [[ -s "$FILE" ]] || {
        echo "ERROR: missing input: $FILE" >&2
        exit 1
    }
done

[[ ! -e "$HIER" && ! -e "$SUMMARY" ]] || {
    echo "ERROR: output already exists; will not overwrite" >&2
    exit 1
}

mkdir -p "$OUTDIR" "$CHECKDIR"

python3 - "$ANNOTATION" "$FULL" "$MIN500" "$HIER" "$SUMMARY" <<'PY'
import csv
import sys

annotation_path, full_path, min_path, hier_path, summary_path = sys.argv[1:]

annotation = {}
with open(annotation_path, newline="") as handle:
    for line_number, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
        if len(row) != 7:
            raise RuntimeError(f"Annotation line {line_number}: expected 7 columns, found {len(row)}")
        original_id = row[0]
        order = row[3]
        family = row[4]
        if not original_id or not order or not family:
            raise RuntimeError(f"Annotation line {line_number}: empty ID, order, or family")
        if original_id in annotation:
            raise RuntimeError(f"Duplicate annotation ID: {original_id}")
        annotation[original_id] = (family, order)

def read_fasta(path):
    records = {}
    name = None
    parts = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    if name in records:
                        raise RuntimeError(f"Duplicate FASTA ID: {name}")
                    records[name] = "".join(parts).upper()
                name = line[1:].split()[0]
                parts = []
            elif line.strip():
                if name is None:
                    raise RuntimeError(f"Sequence before header in {path}")
                parts.append(line.strip())
        if name is not None:
            if name in records:
                raise RuntimeError(f"Duplicate FASTA ID: {name}")
            records[name] = "".join(parts).upper()
    return records

full = read_fasta(full_path)
min500 = read_fasta(min_path)

if len(annotation) != 3975 or len(full) != 3975 or len(min500) != 2960:
    raise RuntimeError(
        f"Unexpected counts: annotation={len(annotation)}, "
        f"full={len(full)}, min500={len(min500)}"
    )

full_original_ids = {name.split("#", 1)[0] for name in full}
if full_original_ids != set(annotation):
    raise RuntimeError("Original FASTA IDs do not match annotation IDs")

rows = []
for full_id, sequence in min500.items():
    original_id = full_id.split("#", 1)[0]
    if full_id not in full or full[full_id] != sequence:
        raise RuntimeError(f"min500 FASTA differs from full library: {full_id}")
    if len(sequence) < 500:
        raise RuntimeError(f"Sequence shorter than 500 bp: {full_id}")
    family, order = annotation[original_id]
    rows.append((full_id, family, order))

with open(hier_path, "x") as handle:
    handle.write("id\tfamily\torder\n")
    for row in rows:
        handle.write("\t".join(row) + "\n")

with open(summary_path, "x") as handle:
    handle.write("metric\tvalue\n")
    handle.write(f"annotation_records\t{len(annotation)}\n")
    handle.write(f"full_library_records\t{len(full)}\n")
    handle.write(f"min500_hierarchy_records\t{len(rows)}\n")
    handle.write(f"distinct_families\t{len({row[1] for row in rows})}\n")
    handle.write(f"output_hierarchy\t{hier_path}\n")
    handle.write("family_source\tall_TEs.txt_column_5\n")
    handle.write("order_source\tall_TEs.txt_column_4\n")
    handle.write("status\tPASS\n")

print(open(summary_path).read(), end="")
print("FINAL: PASS")
PY
