#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=0-04:00:00
#SBATCH -J make_combined_new
#SBATCH -o /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/make_combined_new.%j.out
#SBATCH -e /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/make_combined_new.%j.err

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817

mkdir -p \
  "${BASE}/logs" \
  "${BASE}/checks" \
  "${BASE}/cds_longest_RBH" \
  "${BASE}/combined_fasta"

python3 - <<'PY'
from pathlib import Path
import csv
import re
import shutil

base = Path(
    "/ru-auth/local/home/asato/scratch/MKtest/"
    "reanalysis_longest_fix_20260817"
)

configs = {
    "Dsub_vs_Dsuz": {
        "label": "dsuz",
        "protein_rbh": (
            base / "orth_longest"
            / "Dsub_Dsuz_IsoJpt1.0.RBH.longest.protein.tsv"
        ),
        "mrna_rbh": (
            base / "orth_longest"
            / "Dsub_Dsuz_IsoJpt1.0.RBH.longest.mRNA.tsv"
        ),
        "cds_source": Path(
            "/ru-auth/local/home/asato/scratch/genome/"
            "Dsuz_IsoJpt1.0/"
            "GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_"
            "cds_from_genomic.fna"
        ),
        "cds_key_type": "protein_id",
    },

    "Dsub_vs_Dbia": {
        "label": "dbia",
        "protein_rbh": (
            base / "orth_longest"
            / "Dsub_Dbia.RBH.longest.protein.tsv"
        ),
        "mrna_rbh": (
            base / "orth_longest"
            / "Dsub_Dbia.RBH.longest.mRNA.tsv"
        ),
        "cds_source": Path(
            "/ru-auth/local/home/asato/scratch/genome/"
            "Dbia/"
            "GCF_025231255.1_RU_DBia_V1.1_"
            "cds_from_genomic.fna"
        ),
        "cds_key_type": "protein_id",
    },

    "Dsub_vs_Dmel": {
        "label": "dmel",
        "protein_rbh": (
            base / "orth_longest"
            / "Dsub_Dmel.RBH.longest.protein.tsv"
        ),
        "mrna_rbh": (
            base / "orth_longest"
            / "Dsub_Dmel.RBH.longest.mRNA.tsv"
        ),
        "cds_source": Path(
            "/ru-auth/local/home/asato/scratch/genome/"
            "Dmel/"
            "GCF_000001215.4_Release_6_plus_ISO1_MT_"
            "cds_from_genomic.fna"
        ),
        "cds_key_type": "flybase_fbpp",
    },
}


def read_two_column_table(path):
    rows = []

    with path.open() as fh:
        for line_number, line in enumerate(fh, start=1):
            line = line.rstrip("\n")

            if not line:
                continue

            fields = line.split("\t")

            if len(fields) < 2:
                raise RuntimeError(
                    f"Fewer than two columns in {path}, "
                    f"line {line_number}"
                )

            rows.append((fields[0], fields[1]))

    return rows


def read_fasta(path, key_type):
    sequences = {}
    duplicate_keys = set()
    conflicting_keys = set()

    header = None
    sequence_parts = []

    def save_record(record_header, sequence):
        if record_header is None:
            return

        if key_type == "protein_id":
            match = re.search(
                r"\[protein_id=([^\]]+)\]",
                record_header,
            )

            if match is None:
                return

            key = match.group(1)

        elif key_type == "flybase_fbpp":
            match = re.search(
                r"FLYBASE:(FBpp[0-9]+)",
                record_header,
            )

            if match is None:
                return

            key = match.group(1)

        else:
            raise RuntimeError(
                f"Unknown CDS key type: {key_type}"
            )

        sequence = sequence.upper().replace(" ", "")

        if key in sequences:
            duplicate_keys.add(key)

            if sequences[key] != sequence:
                conflicting_keys.add(key)
        else:
            sequences[key] = sequence

    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")

            if line.startswith(">"):
                if header is not None:
                    save_record(
                        header,
                        "".join(sequence_parts),
                    )

                header = line[1:]
                sequence_parts = []
            else:
                sequence_parts.append(line.strip())

        if header is not None:
            save_record(
                header,
                "".join(sequence_parts),
            )

    if conflicting_keys:
        examples = sorted(conflicting_keys)[:10]

        raise RuntimeError(
            f"Conflicting CDS sequences for duplicate keys "
            f"in {path}: {examples}"
        )

    return sequences, duplicate_keys


def read_fasta_records(path):
    records = []
    header = None
    sequence_parts = []

    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")

            if line.startswith(">"):
                if header is not None:
                    records.append(
                        (header, "".join(sequence_parts))
                    )

                header = line[1:]
                sequence_parts = []
            else:
                sequence_parts.append(line.strip())

        if header is not None:
            records.append(
                (header, "".join(sequence_parts))
            )

    return records


def write_sequence(handle, header, sequence):
    handle.write(f">{header}\n")

    for start in range(0, len(sequence), 80):
        handle.write(
            sequence[start:start + 80] + "\n"
        )


summary_rows = []
fatal_problems = 0

for comparison, config in configs.items():
    print(f"\n===== {comparison} =====")

    required = [
        config["protein_rbh"],
        config["mrna_rbh"],
        config["cds_source"],
    ]

    for path in required:
        if not path.is_file():
            raise FileNotFoundError(path)

    protein_rows = read_two_column_table(
        config["protein_rbh"]
    )
    mrna_rows = read_two_column_table(
        config["mrna_rbh"]
    )

    if len(protein_rows) != len(mrna_rows):
        raise RuntimeError(
            f"{comparison}: protein and mRNA RBH tables "
            f"have different row counts: "
            f"{len(protein_rows)} versus {len(mrna_rows)}"
        )

    # The mRNA RBH table was generated directly from the
    # corresponding protein RBH table while preserving row order.
    rbh_records = []
    seen_dsub_mrna = set()
    duplicated_dsub_mrna = set()

    for protein_row, mrna_row in zip(
        protein_rows,
        mrna_rows,
    ):
        outgroup_protein, dsub_protein = protein_row
        outgroup_mrna, dsub_mrna = mrna_row

        if dsub_mrna in seen_dsub_mrna:
            duplicated_dsub_mrna.add(dsub_mrna)

        seen_dsub_mrna.add(dsub_mrna)

        rbh_records.append(
            {
                "outgroup_protein": outgroup_protein,
                "dsub_protein": dsub_protein,
                "outgroup_mrna": outgroup_mrna,
                "dsub_mrna": dsub_mrna,
            }
        )

    if duplicated_dsub_mrna:
        raise RuntimeError(
            f"{comparison}: duplicated Dsub mRNA IDs: "
            f"{sorted(duplicated_dsub_mrna)[:10]}"
        )

    cds_sequences, duplicate_cds_keys = read_fasta(
        config["cds_source"],
        config["cds_key_type"],
    )

    selected_cds_dir = base / "cds_longest_RBH"
    selected_cds_file = (
        selected_cds_dir
        / f"{comparison}.new_RBH.CDS.fasta"
    )
    selected_map_file = (
        selected_cds_dir
        / f"{comparison}.new_RBH.CDS.map.tsv"
    )
    selected_missing_file = (
        base / "checks"
        / f"{comparison}.new_RBH.CDS.missing.tsv"
    )

    outdir = (
        base
        / "combined_fasta"
        / comparison
    )

    if outdir.exists():
        shutil.rmtree(outdir)

    outdir.mkdir(parents=True, exist_ok=True)

    msa_dir = (
        base
        / "vcf2msa_longest_RBH"
        / comparison
    )

    if not msa_dir.is_dir():
        raise FileNotFoundError(msa_dir)

    rbh_by_dsub = {
        row["dsub_mrna"]: row
        for row in rbh_records
    }

    missing_cds = []
    selected_cds_count = 0

    with selected_cds_file.open("w") as fasta_out, \
         selected_map_file.open("w") as map_out:

        map_out.write(
            "comparison\toutgroup_protein\t"
            "outgroup_mRNA\tDsub_protein\t"
            "Dsub_mRNA\tCDS_length\n"
        )

        for record in rbh_records:
            outgroup_protein = record[
                "outgroup_protein"
            ]

            sequence = cds_sequences.get(
                outgroup_protein
            )

            if sequence is None:
                missing_cds.append(record)
                continue

            write_sequence(
                fasta_out,
                outgroup_protein,
                sequence,
            )

            map_out.write(
                "\t".join(
                    [
                        comparison,
                        outgroup_protein,
                        record["outgroup_mrna"],
                        record["dsub_protein"],
                        record["dsub_mrna"],
                        str(len(sequence)),
                    ]
                )
                + "\n"
            )

            selected_cds_count += 1

    with selected_missing_file.open("w") as out:
        out.write(
            "outgroup_protein\toutgroup_mRNA\t"
            "Dsub_protein\tDsub_mRNA\treason\n"
        )

        for record in missing_cds:
            out.write(
                "\t".join(
                    [
                        record["outgroup_protein"],
                        record["outgroup_mrna"],
                        record["dsub_protein"],
                        record["dsub_mrna"],
                        "no_CDS_sequence",
                    ]
                )
                + "\n"
            )

    combined_count = 0
    msa_count = 0
    msa_not_in_rbh = []
    msa_missing_cds = []
    invalid_combined = []

    for msa_file in sorted(msa_dir.glob("*.fasta")):
        msa_count += 1
        dsub_mrna = msa_file.stem

        record = rbh_by_dsub.get(dsub_mrna)

        if record is None:
            msa_not_in_rbh.append(dsub_mrna)
            continue

        outgroup_protein = record[
            "outgroup_protein"
        ]
        outgroup_sequence = cds_sequences.get(
            outgroup_protein
        )

        if outgroup_sequence is None:
            msa_missing_cds.append(dsub_mrna)
            continue

        dsub_records = read_fasta_records(msa_file)

        if not dsub_records:
            invalid_combined.append(
                (dsub_mrna, "empty_Dsub_MSA")
            )
            continue

        existing_headers = {
            header.split()[0]
            for header, sequence in dsub_records
        }

        if config["label"] in existing_headers:
            invalid_combined.append(
                (
                    dsub_mrna,
                    "outgroup_label_already_in_Dsub_MSA",
                )
            )
            continue

        output_file = (
            outdir
            / f"{dsub_mrna}.combined.fasta"
        )

        with output_file.open("w") as out:
            write_sequence(
                out,
                config["label"],
                outgroup_sequence,
            )

            for header, sequence in dsub_records:
                write_sequence(
                    out,
                    header,
                    sequence,
                )

        output_records = read_fasta_records(
            output_file
        )

        if len(output_records) != (
            len(dsub_records) + 1
        ):
            invalid_combined.append(
                (
                    dsub_mrna,
                    "unexpected_sequence_count",
                )
            )
            output_file.unlink(missing_ok=True)
            continue

        combined_count += 1

    check_file = (
        base
        / "checks"
        / f"{comparison}.combined_FASTA.check.tsv"
    )

    with check_file.open("w") as out:
        out.write("category\tidentifier_or_count\n")
        out.write(
            f"protein_RBH_pairs\t{len(protein_rows)}\n"
        )
        out.write(
            f"mRNA_RBH_pairs\t{len(mrna_rows)}\n"
        )
        out.write(
            f"selected_CDS\t{selected_cds_count}\n"
        )
        out.write(
            f"missing_CDS\t{len(missing_cds)}\n"
        )
        out.write(
            f"input_Dsub_MSA\t{msa_count}\n"
        )
        out.write(
            f"combined_FASTA\t{combined_count}\n"
        )
        out.write(
            f"MSA_not_in_RBH\t{len(msa_not_in_rbh)}\n"
        )
        out.write(
            f"MSA_missing_CDS\t{len(msa_missing_cds)}\n"
        )
        out.write(
            f"invalid_combined\t{len(invalid_combined)}\n"
        )
        out.write(
            "duplicate_CDS_keys_with_identical_sequence\t"
            f"{len(duplicate_cds_keys)}\n"
        )

    comparison_fatal = (
        len(missing_cds)
        + len(msa_not_in_rbh)
        + len(msa_missing_cds)
        + len(invalid_combined)
    )

    fatal_problems += comparison_fatal

    orco = "rna-XM_037873416.1"
    orco_combined = (
        outdir
        / f"{orco}.combined.fasta"
    )

    print(f"protein_RBH_pairs\t{len(protein_rows)}")
    print(f"mRNA_RBH_pairs\t{len(mrna_rows)}")
    print(f"selected_outgroup_CDS\t{selected_cds_count}")
    print(f"missing_outgroup_CDS\t{len(missing_cds)}")
    print(f"input_Dsub_MSA\t{msa_count}")
    print(f"combined_FASTA_written\t{combined_count}")
    print(f"MSA_not_in_RBH\t{len(msa_not_in_rbh)}")
    print(f"MSA_missing_CDS\t{len(msa_missing_cds)}")
    print(f"invalid_combined\t{len(invalid_combined)}")
    print(
        "Orco_combined\t"
        + (
            "present"
            if orco_combined.is_file()
            else "missing"
        )
    )
    print(f"output_directory\t{outdir}")

    summary_rows.append(
        (
            comparison,
            len(protein_rows),
            selected_cds_count,
            len(missing_cds),
            msa_count,
            combined_count,
            len(msa_not_in_rbh),
            len(msa_missing_cds),
            len(invalid_combined),
            orco_combined.is_file(),
        )
    )


summary_file = (
    base
    / "checks"
    / "outgroup_CDS_and_combined_FASTA.summary.tsv"
)

with summary_file.open("w") as out:
    out.write(
        "comparison\tRBH_pairs\tselected_CDS\t"
        "missing_CDS\tinput_Dsub_MSA\t"
        "combined_FASTA\tMSA_not_in_RBH\t"
        "MSA_missing_CDS\tinvalid_combined\t"
        "Orco_combined\n"
    )

    for row in summary_rows:
        out.write(
            "\t".join(map(str, row)) + "\n"
        )

print(f"\nsummary\t{summary_file}")

if fatal_problems:
    raise RuntimeError(
        f"Found {fatal_problems} combined-FASTA problems. "
        f"See checks/*.combined_FASTA.check.tsv and "
        f"checks/*.new_RBH.CDS.missing.tsv"
    )

print("Done.")
PY
