#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0-01:00:00
#SBATCH -J Dsub_longest_GFF
#SBATCH -o /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/Dsub_longest_GFF.%j.out
#SBATCH -e /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/Dsub_longest_GFF.%j.err

set -euo pipefail

GENOME_DIR=/ru-auth/local/home/asato/scratch/genome/Dsub_v1.1
NEW_BASE=/ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817

GFF=${GENOME_DIR}/GCF_014743375.2_RU_Dsub_v1.1_genomic.gff
PROTEIN_FASTA=${GENOME_DIR}/GCF_014743375.2_RU_Dsub_v1.1_protein.faa

OUT_FASTA=${NEW_BASE}/pep_longest/Dsub.longest_per_gene.pep.fa
OUT_MAP=${NEW_BASE}/pep_longest/Dsub.longest_per_gene.map.tsv
AUDIT=${NEW_BASE}/checks/Dsub.longest_per_gene.all_proteins.tsv
SUMMARY=${NEW_BASE}/checks/Dsub.longest_per_gene.summary.tsv
UNMAPPED=${NEW_BASE}/checks/Dsub.longest_per_gene.unmapped_proteins.txt
CONFLICTS=${NEW_BASE}/checks/Dsub.longest_per_gene.mapping_conflicts.tsv

mkdir -p \
  "${NEW_BASE}/pep_longest" \
  "${NEW_BASE}/checks" \
  "${NEW_BASE}/logs"

for f in "${GFF}" "${PROTEIN_FASTA}"; do
    if [[ ! -s "${f}" ]]; then
        echo "ERROR: missing or empty input file: ${f}" >&2
        exit 1
    fi
done

export GFF
export PROTEIN_FASTA
export OUT_FASTA
export OUT_MAP
export AUDIT
export SUMMARY
export UNMAPPED
export CONFLICTS

python3 - <<'PY'
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

gff_path = Path(os.environ["GFF"])
protein_fasta = Path(os.environ["PROTEIN_FASTA"])

out_fasta = Path(os.environ["OUT_FASTA"])
out_map = Path(os.environ["OUT_MAP"])
audit_file = Path(os.environ["AUDIT"])
summary_file = Path(os.environ["SUMMARY"])
unmapped_file = Path(os.environ["UNMAPPED"])
conflicts_file = Path(os.environ["CONFLICTS"])


def parse_attributes(text):
    attributes = {}

    for field in text.rstrip().split(";"):
        if "=" not in field:
            continue

        key, value = field.split("=", 1)
        attributes[key] = unquote(value)

    return attributes


def strip_prefix(value, prefix):
    if value.startswith(prefix):
        return value[len(prefix):]
    return value


def extract_geneid(dbxref):
    if not dbxref:
        return ""

    match = re.search(r"(?:^|,)GeneID:([0-9]+)(?:,|$)", dbxref)
    return match.group(1) if match else ""


# ============================================================
# Parse GFF
# ============================================================

mrna_to_gene = {}
mrna_to_geneid = {}
mrna_to_symbol = {}
mrna_to_product = {}

protein_to_mrna = {}
mapping_conflicts = []

with gff_path.open() as fh:
    for line_number, line in enumerate(fh, start=1):
        if not line or line.startswith("#"):
            continue

        columns = line.rstrip("\n").split("\t")

        if len(columns) != 9:
            continue

        feature = columns[2]
        attrs = parse_attributes(columns[8])

        if feature in {"mRNA", "transcript"}:
            mrna_feature_id = attrs.get("ID", "")
            parent_gene = attrs.get("Parent", "")
            transcript_id = attrs.get("transcript_id", "")

            if not mrna_feature_id or not parent_gene:
                continue

            mrna_feature_id = strip_prefix(mrna_feature_id, "rna-")
            parent_gene = strip_prefix(parent_gene, "gene-")

            # Use the accession as the mRNA identifier when available.
            mrna_id = transcript_id or mrna_feature_id

            gene_id = extract_geneid(attrs.get("Dbxref", ""))
            gene_symbol = attrs.get("gene", "")
            product = attrs.get("product", "")

            mrna_to_gene[mrna_id] = parent_gene
            mrna_to_geneid[mrna_id] = gene_id
            mrna_to_symbol[mrna_id] = gene_symbol
            mrna_to_product[mrna_id] = product

        elif feature == "CDS":
            protein_id = attrs.get("protein_id", "")
            parent_mrna = attrs.get("Parent", "")

            if not protein_id or not parent_mrna:
                continue

            parent_mrna = strip_prefix(parent_mrna, "rna-")

            old = protein_to_mrna.get(protein_id)

            if old is not None and old != parent_mrna:
                mapping_conflicts.append(
                    (
                        protein_id,
                        old,
                        parent_mrna,
                        line_number,
                    )
                )
            else:
                protein_to_mrna[protein_id] = parent_mrna


with conflicts_file.open("w") as out:
    out.write(
        "protein_id\tprevious_mRNA\tnew_mRNA\tGFF_line\n"
    )

    for row in mapping_conflicts:
        out.write("\t".join(map(str, row)) + "\n")


if mapping_conflicts:
    raise RuntimeError(
        f"Found {len(mapping_conflicts)} conflicting protein-to-mRNA "
        f"mappings. See {conflicts_file}"
    )


# ============================================================
# Read protein FASTA
# ============================================================

proteins = {}

header = None
sequence_parts = []


def save_record(current_header, parts):
    if current_header is None:
        return

    protein_id = current_header.split()[0]
    sequence = "".join(parts).replace(" ", "").replace("\t", "")
    length = len(sequence.rstrip("*"))

    if protein_id in proteins:
        raise RuntimeError(
            f"Duplicated protein FASTA identifier: {protein_id}"
        )

    proteins[protein_id] = {
        "protein_id": protein_id,
        "header": current_header,
        "sequence": sequence,
        "length": length,
    }


with protein_fasta.open() as fh:
    for line in fh:
        line = line.rstrip("\n")

        if line.startswith(">"):
            save_record(header, sequence_parts)
            header = line[1:]
            sequence_parts = []
        else:
            sequence_parts.append(line.strip())

save_record(header, sequence_parts)


# ============================================================
# Connect proteins to genes
# ============================================================

mapped_records = []
unmapped_records = []

for protein_id, protein in proteins.items():
    mrna_id = protein_to_mrna.get(protein_id)

    if mrna_id is None:
        unmapped_records.append((protein_id, "no_CDS_mapping"))
        continue

    gene_key = mrna_to_gene.get(mrna_id)

    if gene_key is None:
        unmapped_records.append(
            (protein_id, f"mRNA_not_mapped_to_gene:{mrna_id}")
        )
        continue

    record = dict(protein)
    record.update(
        {
            "mrna_id": mrna_id,
            "gene_key": gene_key,
            "GeneID": mrna_to_geneid.get(mrna_id, ""),
            "gene_symbol": mrna_to_symbol.get(mrna_id, ""),
            "product": mrna_to_product.get(mrna_id, ""),
        }
    )

    mapped_records.append(record)


with unmapped_file.open("w") as out:
    out.write("protein_id\treason\n")

    for row in unmapped_records:
        out.write("\t".join(row) + "\n")


if unmapped_records:
    print(
        f"ERROR: {len(unmapped_records)} protein records could not be "
        f"mapped to a gene.",
        file=sys.stderr,
    )
    print(f"See: {unmapped_file}", file=sys.stderr)
    sys.exit(1)


# ============================================================
# Select longest protein per gene
#
# Selection:
#   1. greatest amino-acid length
#   2. lexicographically smallest mRNA accession
#   3. lexicographically smallest protein accession
# ============================================================

records_by_gene = defaultdict(list)

for record in mapped_records:
    records_by_gene[record["gene_key"]].append(record)


selected_by_gene = {}

for gene_key, records in records_by_gene.items():
    ranked = sorted(
        records,
        key=lambda record: (
            -record["length"],
            record["mrna_id"],
            record["protein_id"],
        ),
    )

    selected_by_gene[gene_key] = ranked[0]


# ============================================================
# Write selected FASTA and mapping table
# ============================================================

with out_fasta.open("w") as fasta_out, out_map.open("w") as map_out:
    map_out.write(
        "gene_key\tGeneID\tgene_symbol\tmRNA_accession\t"
        "protein_accession\tprotein_length\tn_protein_isoforms\t"
        "product\n"
    )

    for gene_key in sorted(selected_by_gene):
        record = selected_by_gene[gene_key]

        fasta_out.write(
            f">{record['protein_id']} "
            f"gene={record['gene_key']} "
            f"GeneID={record['GeneID']} "
            f"mRNA={record['mrna_id']}\n"
        )

        sequence = record["sequence"]

        for start in range(0, len(sequence), 80):
            fasta_out.write(sequence[start:start + 80] + "\n")

        map_out.write(
            "\t".join(
                [
                    record["gene_key"],
                    record["GeneID"],
                    record["gene_symbol"],
                    record["mrna_id"],
                    record["protein_id"],
                    str(record["length"]),
                    str(len(records_by_gene[gene_key])),
                    record["product"],
                ]
            )
            + "\n"
        )


# ============================================================
# Write audit table containing all proteins
# ============================================================

with audit_file.open("w") as out:
    out.write(
        "gene_key\tGeneID\tgene_symbol\tmRNA_accession\t"
        "protein_accession\tprotein_length\tn_protein_isoforms\t"
        "selected\tproduct\n"
    )

    for gene_key in sorted(records_by_gene):
        selected = selected_by_gene[gene_key]

        ranked = sorted(
            records_by_gene[gene_key],
            key=lambda record: (
                -record["length"],
                record["mrna_id"],
                record["protein_id"],
            ),
        )

        for record in ranked:
            out.write(
                "\t".join(
                    [
                        record["gene_key"],
                        record["GeneID"],
                        record["gene_symbol"],
                        record["mrna_id"],
                        record["protein_id"],
                        str(record["length"]),
                        str(len(records_by_gene[gene_key])),
                        "1" if record["protein_id"]
                        == selected["protein_id"] else "0",
                        record["product"],
                    ]
                )
                + "\n"
            )


# ============================================================
# Summary and Orco check
# ============================================================

n_multisoform_genes = sum(
    len(records) > 1
    for records in records_by_gene.values()
)

n_removed_isoforms = sum(
    len(records) - 1
    for records in records_by_gene.values()
)

orco_candidates = [
    record
    for record in mapped_records
    if record["mrna_id"] in {
        "XM_037873416.1",
        "XM_037873417.1",
    }
]

orco_selected = [
    record
    for record in selected_by_gene.values()
    if record["mrna_id"] in {
        "XM_037873416.1",
        "XM_037873417.1",
    }
]

summary = [
    ("input_protein_records", len(proteins)),
    ("GFF_protein_to_mRNA_mappings", len(protein_to_mrna)),
    ("GFF_mRNA_to_gene_mappings", len(mrna_to_gene)),
    ("mapped_protein_records", len(mapped_records)),
    ("unmapped_protein_records", len(unmapped_records)),
    ("unique_gene_keys", len(records_by_gene)),
    ("selected_longest_proteins", len(selected_by_gene)),
    ("genes_with_multiple_protein_isoforms", n_multisoform_genes),
    ("extra_protein_isoforms_removed", n_removed_isoforms),
    ("mapping_conflicts", len(mapping_conflicts)),
    ("Orco_candidate_isoforms", len(orco_candidates)),
    ("Orco_selected_isoforms", len(orco_selected)),
]

with summary_file.open("w") as out:
    out.write("category\tcount_or_value\n")

    for category, value in summary:
        out.write(f"{category}\t{value}\n")

    for record in orco_candidates:
        out.write(
            "Orco_candidate_detail\t"
            f"{record['mrna_id']}|{record['protein_id']}|"
            f"length={record['length']}|gene={record['gene_key']}\n"
        )

    for record in orco_selected:
        out.write(
            "Orco_selected_detail\t"
            f"{record['mrna_id']}|{record['protein_id']}|"
            f"length={record['length']}|gene={record['gene_key']}\n"
        )


print("Dsub GFF-based longest-protein selection completed.")

for category, value in summary:
    print(f"{category}\t{value}")

for record in orco_selected:
    print(
        "Orco_selected\t"
        f"{record['mrna_id']}\t"
        f"{record['protein_id']}\t"
        f"length={record['length']}\t"
        f"gene={record['gene_key']}"
    )

print(f"output_fasta\t{out_fasta}")
print(f"output_map\t{out_map}")
print(f"audit_file\t{audit_file}")
print(f"summary_file\t{summary_file}")
PY

echo
echo "Output FASTA records:"
grep -c '^>' "${OUT_FASTA}"

echo
echo "Output files:"
ls -lh \
  "${OUT_FASTA}" \
  "${OUT_MAP}" \
  "${AUDIT}" \
  "${SUMMARY}" \
  "${UNMAPPED}" \
  "${CONFLICTS}"

echo
echo "Done."
